#include "WireCellSst/GeomDataSource.h"
//#include "WireCellSst/ToyuBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/ToyMatrixExclusive.h"
#include "WireCell2dToy/ToyMatrixKalman.h"
#include "WireCell2dToy/ToyMatrixIterate.h"
#include "WireCell2dToy/ToyMetric.h"

#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"

#include "WireCellData/GeomCluster.h"
//#include "WireCellNav/SliceDataSource.h"


#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/SimDataSource.h"
#include "WireCellNav/SliceDataSource.h"
#include "WireCellSst/Util.h"
#include "WireCellData/SimTruth.h"
#include "WireCell2dToy/ToyDepositor.h"
#include "WireCellNav/GenerativeFDS.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include <iostream>

using namespace WireCell;
using namespace std;




int main(int argc, char* argv[])
{
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root" << endl;
    return 1;
  }
  
  
  WireCellSst::GeomDataSource gds(argv[1]);
  std::vector<float> ex = gds.extent();
  cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;
  
  
  const char* root_file = argv[2];
  const char* tpath = "/Event/Sim";
  
  WireCell::FrameDataSource* fds = 0;
  fds = WireCellSst::make_fds(root_file);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  

  
 

 


  WireCell::ToyDepositor toydep(fds);
  const PointValueVector pvv = toydep.depositions(1);
  
  


  WireCell::GenerativeFDS gfds(toydep,gds,2400,5);
  gfds.jump(1);

  WireCellSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons

  

  const int N = 100000;
  Double_t x[N],y[N],z[N];
  Double_t xt[N],yt[N],zt[N];
  int ncount = 0;
  int ncount_t = 0;
  

  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
  WireCell2dToy::MergeToyTiling **mergetiling = new WireCell2dToy::MergeToyTiling*[2400];
  WireCell2dToy::TruthToyTiling **truthtiling = new WireCell2dToy::TruthToyTiling*[2400];
  WireCell2dToy::ToyMatrix **toymatrix = new WireCell2dToy::ToyMatrix*[2400];
  WireCell2dToy::ToyMatrixIterate **toymatrix_it = new WireCell2dToy::ToyMatrixIterate*[2400];
  
  WireCell2dToy::ToyMetric toymetric;


  //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  
  int ncount_mcell = 0;
  

  //  int i=454;{ // 46, 26
  // int i=329;{  // 18, 6,
  // int i=344;{ // 29 14
  int i=459;{ // 23, 8,   5e5 
    //for (int i=0;i!=sds.size();i++){
    // for (int i=365;i!=378;i++){
 
    sds.jump(i);
    WireCell::Slice slice = sds.get();
    if ( slice.group().size() >0){
      toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds);
      mergetiling[i] = new WireCell2dToy::MergeToyTiling(*toytiling[i],i);

      GeomCellSelection allcell = toytiling[i]->get_allcell();
      GeomWireSelection allwire = toytiling[i]->get_allwire();
      GeomCellSelection allmcell = mergetiling[i]->get_allcell();
      GeomWireSelection allmwire = mergetiling[i]->get_allwire();
      
      cout << i << " " << allmcell.size() << " " << allmwire.size() << endl;

      truthtiling[i] = new WireCell2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds);
      toymatrix[i] = new WireCell2dToy::ToyMatrix(*toytiling[i],*mergetiling[i]);
      if (toymatrix[i]->Get_Solve_Flag()==0)
	//	cout << i << endl;
	toymatrix_it[i] = new WireCell2dToy::ToyMatrixIterate(*toymatrix[i]);
      
      


      


      GeomCellSelection calmcell;
      for (int j=0;j!=allmcell.size();j++){
	MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
	double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
	double charge_err = toymatrix[i]->Get_Cell_Charge(mcell,2);
	if (charge + charge_err > 2000) calmcell.push_back(mcell);
      }
      

      // for (int j=0;j!=allcell.size();j++){
      // 	Point p = allcell[j]->center();
      // 	x[ncount] = i*0.32;
      // 	y[ncount] = p.y/units::cm;
      // 	z[ncount] = p.z/units::cm;
      // 	ncount ++;
      // }


      CellChargeMap ccmap = truthtiling[i]->ccmap();
      toymetric.Add(allmcell,*toymatrix[i],ccmap);
      Double_t charge_min = 10000;
      Double_t charge_max = 0;

      for (auto it = ccmap.begin();it!=ccmap.end(); it++){
      	Point p = it->first->center();
      	xt[ncount_t] = i*0.32;
      	yt[ncount_t] = p.y/units::cm;
      	zt[ncount_t] = p.z/units::cm;
      	ncount_t ++;

      	float charge = it->second;
      	if (charge > charge_max) charge_max = charge;
      	if (charge < charge_min) charge_min = charge;
       	// cout << it->second << endl;
      }

      
      // //loop through merged cell and compare with truth cells
      // for (int j=0;j!=allmcell.size();j++){
      // 	MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      // 	if (mcell->CheckContainTruthCell(ccmap)) 
      // 	  cout << "xin: " << toymatrix[i]->Get_mcindex(mcell) << endl;
      // 	//cout << mergetiling.wires(*allmcell[j]).size() << endl;
      // }
      
      // WireChargeMap wcmap = toytiling[i]->wcmap();



      // // Check Charge for Single Cell, construct an index map
      // WireIndexMap wimap;
      // CellIndexMap cimap;

      // int cindex = 0;
      // int windex = 0;
      // for (auto it = ccmap.begin(); it != ccmap.end(); it++){
      // 	const GeomCell *cell1 = it->first;
      // 	if (cimap.find(cell1) == cimap.end()){
      // 	  cimap[cell1] = cindex;
      // 	  cindex ++;

      // 	  //find the wires associated with this cell
      // 	  const GeomWireSelection wires = toytiling[i]->wires(*cell1);
      // 	  for (Int_t j=0;j!=wires.size();j++){
      // 	    const GeomWire *wire1 = wires[j];
      // 	    if (wimap.find(wire1) == wimap.end()){
      // 	      wimap[wire1] = windex;
      // 	      windex ++;
      // 	    }
      // 	  }
      // 	}
      // }

      // // cout << cindex << " " << ccmap.size() << " " << windex << " " << wimap.size() << endl;
      // TVectorD Vcell_charge(cindex);
      // TVectorD Vwire_tcharge(windex);
      // TVectorD Vwire_charge(windex);
      // TMatrixD A(windex,cindex);
      // for (auto it = ccmap.begin(); it != ccmap.end(); it++){
      // 	int index = cimap[it->first];
      // 	Vcell_charge[index] = it->second;
      // 	const GeomCell *cell1 = it->first;
      // 	const GeomWireSelection wires = toytiling[i]->wires(*cell1);
      // 	for (Int_t j=0;j!=wires.size();j++){
      // 	  const GeomWire *wire1 = wires[j];
      // 	  int index1 = wimap[wire1];
      // 	  A(index1,index) = 1;
      // 	}
      // }
      // for (auto it = wcmap.begin();it != wcmap.end(); it++){
      // 	int index = wimap[it->first];
      // 	Vwire_tcharge[index] = it->second;
      // }
      
      // //Vwire_tcharge.Print();
      // //A.Print();
      
      // Vwire_charge = A * Vcell_charge;
      // //Vwire_charge.Print();
      
      // double sum1 = 0, sum2 = 0, sum3 = 0;

      // for (int j=0;j!=windex;j++){
      // 	//cout << j << " " << Vwire_tcharge[j]<< " " << Vwire_charge[j] << " " << (Vwire_tcharge[j] - Vwire_charge[j]) << endl;
      // 	sum1 += Vwire_tcharge[j]; 
      // 	sum2 += Vwire_charge[j];
      // }
      // for (int j=0;j!=cindex;j++){
      // 	sum3 += Vcell_charge[j];
      // }

      // cout << "Single Cell Check: " << sum1/3. << " " << sum2/3. << " " << sum3 << endl;
      
      // // Vcell_charge.Print();


      // // Check Mapping of Merged Cell and Wire
      // WireIndexMap mwimap;
      // WireIndexMap swimap;
      // CellIndexMap mcimap;

      // int mcindex = 0;
      // int mwindex = 0;
      // int swindex = 0;

      // for (int j=0;j!=allmcell.size();j++){
      // 	//construct merged cell index
      // 	const MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      // 	if (mcimap.find(mcell) == mcimap.end()){
      // 	  mcimap[mcell] = mcindex;
      // 	  mcindex ++;

      // 	  const GeomWireSelection wires = mergetiling[i]->wires(*allmcell[j]);
      // 	  // cout << wires.size() << endl;
      // 	  for (int k=0;k!=wires.size();k++){
      // 	    //construct merged wire index
      // 	    const MergeGeomWire *mwire = (MergeGeomWire*)wires[k];
      // 	    if (mwimap.find(mwire) == mwimap.end()){
      // 	      mwimap[mwire] = mwindex;
      // 	      mwindex ++;

      // 	      //construct single wire index
      // 	      GeomWireSelection swires = mwire->get_allwire();
      // 	      for (int kk = 0; kk!=swires.size(); kk++){
      // 	       	const GeomWire* wire1 = swires[kk];
      // 	       	if (swimap.find(wire1) == swimap.end()){
      // 	       	  swimap[wire1] = swindex;
      // 	       	  swindex ++;
      // 	       	}
      // 	      }

      // 	    }
      // 	  }

      // 	}
      // }

      // // cout << mcindex << " " << mwindex << " " << swindex << " " << allmcell.size() << " " << allmwire.size() << " " << allwire.size() << endl;
      
      // TVectorD Vmwire_lcharge(mwindex);
      // TVectorD Vmwire_rcharge(mwindex);

      // TVectorD Vswire_charge(swindex);
      // TVectorD Vmcell_charge(mcindex);

      // TMatrixD MA(mwindex,mcindex);
      // TMatrixD MB(mwindex,swindex);

      // for (int j=0;j!=allwire.size();j++){
      // 	int index = swimap[allwire[j]];
      // 	float charge = wcmap[allwire[j]];
      // 	Vswire_charge[index] =charge;
      // }
      // // Vswire_charge.Print();
      // for (int j = 0;j!=allmcell.size();j++){
      // 	int index = mcimap[allmcell[j]];
      // 	float charge = 0;
      // 	for (int k=0; k!=((const MergeGeomCell*)allmcell[j])->get_allcell().size(); k++){
      // 	  charge += ccmap[((const MergeGeomCell*)allmcell[j])->get_allcell().at(k)];
      // 	}
      // 	Vmcell_charge[index] = charge;
      // }
      // // Vmcell_charge.Print();

      // for (int j=0;j!=allmwire.size();j++){
      // 	int index = mwimap[allmwire[j]];
      // 	//construct MA
      // 	for (int k=0; k!=mergetiling[i]->cells(*allmwire[j]).size();k++){
      // 	  int index1 = mcimap[mergetiling[i]->cells(*allmwire[j]).at(k)];
      // 	  MA(index,index1) = 1;
      // 	}

      // 	//construct MB
      // 	for (int k=0;k!=((MergeGeomWire*)allmwire[j])->get_allwire().size();k++){
      // 	  int index1 = swimap[((MergeGeomWire*)allmwire[j])->get_allwire().at(k)];
      // 	  MB(index,index1) = 1;
      // 	}

      // }
      // // MA.Print();
      // //MB.Print();

      // Vmwire_lcharge = MB * Vswire_charge;
      // Vmwire_rcharge = MA * Vmcell_charge;
      
      // sum1 = 0; sum2 = 0;
      // for (int j=0;j!= mwindex;j++){
      // 	sum1 += Vmwire_lcharge[j];
      // 	sum2 += Vmwire_rcharge[j];
      // 	//cout << Vmwire_lcharge[j] << " " << Vmwire_rcharge[j] << " " << Vmwire_lcharge[j] -Vmwire_rcharge[j] << endl;
      // }
      // cout << "Merged Cell Check: " << sum1/3. << " " << sum2/3. << endl;



      // for (auto it = wcmap.begin();it!=wcmap.end(); it++){
      // 	float charge = it->second;
      // 	if (charge > charge_max) charge_max = charge;
      // 	if (charge < charge_min) charge_min = charge;
      // }
    // int sum = 0;
    // for (int j=0;j!=allmcell.size();j++){
    //   sum += ((WireCell::MergeGeomCell*)allmcell[j])->get_allcell().size() ;
    // }
    // cout << allcell.size() << " " << allmcell.size() << " "  << sum << endl;

     // for (int j=0;j!=allmcell.size();j++){
     //   cout << j << " " << mergetiling.wires(*allmcell[j]).size() << endl;
     // }

    // for (int j=0;j!=allwire.size();j++){
    //   const GeomWire *wire = allwire[j];
    //   const GeomCellSelection targetcells =  mergetiling.cells(*wire);
    //   cout << j << " " << targetcells.size() << endl;
    // }

    

    // GeomCellSelection allcell = toytiling.get_allcell();
    // GeomWireSelection allwire = toytiling.get_allwire();
      
    //
    //}
  

  // cout << toytiling.wiremap[allwire.at(0)].size() << endl;
  // cout << toytiling.cellmap[allcell.at(0)].size() << endl;
  // //  

      
    // TApplication theApp("theApp",&argc,argv);
    // theApp.SetReturnFromRun(true);
    
    // TCanvas c1("ToyMC","ToyMC",800,600);
    // c1.Draw();
    
    // WireCell2dToy::ToyEventDisplay display(c1, gds);
    // display.charge_min = charge_min;
    // display.charge_max = charge_max;


    // gStyle->SetOptStat(0);
    
    // const Int_t NRGBs = 5;
    // const Int_t NCont = 255;
    // Int_t MyPalette[NCont];
    // Double_t stops[NRGBs] = {0.0, 0.34, 0.61, 0.84, 1.0};
    // Double_t red[NRGBs] = {0.0, 0.0, 0.87 ,1.0, 0.51};
    // Double_t green[NRGBs] = {0.0, 0.81, 1.0, 0.2 ,0.0};
    // Double_t blue[NRGBs] = {0.51, 1.0, 0.12, 0.0, 0.0};
    // Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    // gStyle->SetNumberContours(NCont);
    // for (int kk=0;kk!=NCont;kk++) MyPalette[kk] = FI+kk;
    // gStyle->SetPalette(NCont,MyPalette);

    

    // display.init(0,10.3698,-2.33/2.,2.33/2.);
    // //display.init(1.1,1.8,0.7,1.0);
    
    // display.draw_mc(1,WireCell::PointValueVector(),"colz");
    
    

    // // display.draw_slice(slice,"");
    // display.draw_cells(toytiling[i]->get_allcell(),"*same");
    // //display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",1); //0 is normal, 1 is only draw the ones containt the truth cell
    // display.draw_mergecells(calmcell,"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell
    // display.draw_truthcells(ccmap,"*same");
    
    // //display.draw_wires_charge(wcmap,"Fsame",FI);
    // // display.draw_cells_charge(toytiling.get_allcell(),"Fsame");
    // // display.draw_truthcells_charge(ccmap,"lFsame",FI);
    
    
    // theApp.Run();
    }
  }

 


 
  // TGraph2D *g = new TGraph2D(ncount,x,y,z);
  // TGraph2D *gt = new TGraph2D(ncount_t,xt,yt,zt);
  // TFile *file = new TFile("shower3D.root","RECREATE");
  // g->Write("shower3D");
  // gt->Write("shower3D_truth");

  toymetric.Print();

  return 0;
  
} // main()
