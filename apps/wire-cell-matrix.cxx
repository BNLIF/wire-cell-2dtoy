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
  //int i=459;{ // 23, 8,   5e5 
  //int i = 351;{
  //for (int i=0;i!=sds.size();i++){
  for (int i=348;i!=354;i++){
 
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
	toymatrix_it[i] = new WireCell2dToy::ToyMatrixIterate(*toymatrix[i]);
      
      cout << "chi2: " << toymatrix[i]->Get_Chi2() << endl;
      cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
      
      


      GeomCellSelection calmcell;
      for (int j=0;j!=allmcell.size();j++){
	MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
	double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
	double charge_err = toymatrix[i]->Get_Cell_Charge(mcell,2);
	
	//	cout << "Recon: " << j << " " << charge << " " << charge_err << endl;

	if (charge + charge_err > 2000) calmcell.push_back(mcell);
      }
      

      


      CellChargeMap ccmap = truthtiling[i]->ccmap();
      if (toymatrix[i]->Get_Solve_Flag()!=0)
	toymetric.Add(allmcell,*toymatrix[i],ccmap);

      toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());

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
      // 	  cout << "True: " << toymatrix[i]->Get_mcindex(mcell) << " " << mcell->GetTruthCharge()<< endl;
      // 	//cout << mergetiling.wires(*allmcell[j]).size() << endl;
      // }
      
      // WireChargeMap wcmap = toytiling[i]->wcmap();




      
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
    
    

    // // display.draw_slice(slice,""); // draw wire 
    // display.draw_cells(toytiling[i]->get_allcell(),"*same");
    // //display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",1); //0 is normal, 1 is only draw the ones containt the truth cell
    // display.draw_mergecells(calmcell,"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell
    // display.draw_truthcells(ccmap,"*same");
    
    // // display.draw_wires_charge(wcmap,"Fsame",FI);
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
