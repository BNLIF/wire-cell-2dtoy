#include "WireCellSst/GeomDataSource.h"
//#include "WireCellSst/ToyuBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCellData/MergeGeomCell.h"
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
  std::vector<double> ex = gds.extent();
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
  
  int max_events = 100;
  int eve_num = 18;

  
  WireCell::ToyDepositor toydep(fds);
  const PointValueVector pvv = toydep.depositions(eve_num);
  

  WireCell::GenerativeFDS gfds(toydep,gds,2400,max_events,2.0*1.6*units::millimeter);
  gfds.jump(eve_num);

  WireCellSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons

  

  const int N = 100000;
  Double_t x[N],y[N],z[N];
  Double_t xt[N],yt[N],zt[N];
  int ncount = 0;
  int ncount_t = 0;
  

  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
  WireCell2dToy::MergeToyTiling **mergetiling = new WireCell2dToy::MergeToyTiling*[2400];
  WireCell2dToy::TruthToyTiling **truthtiling = new WireCell2dToy::TruthToyTiling*[2400];
  
  //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  
  int ncount_mcell = 0;
  

  //int i=178;{
  int i=1172-800;{
    //int i=441;{
    // for (int i=0;i!=sds.size();i++){
    //for (int i=365;i!=378;i++){
 
    sds.jump(i);
    WireCell::Slice slice = sds.get();
    if ( slice.group().size() >0){
      cout << i << " " << slice.group().size() << endl;

      toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds);
      GeomCellSelection allcell = toytiling[i]->get_allcell();

      cout << allcell.size() << endl;
      
      mergetiling[i] = new WireCell2dToy::MergeToyTiling(*toytiling[i],i,2);
      truthtiling[i] = new WireCell2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds);
      
      

      // for (int j=0;j!=allcell.size();j++){
      // 	std::cout << toytiling[i]->wires(*allcell.at(j)).size() << std::endl;
      // }
      

      GeomCellSelection allmcell = mergetiling[i]->get_allcell();
      GeomWireSelection allmwire = mergetiling[i]->get_allwire();
      
      cout << allmcell.size() << endl;
     
      if (cluster_set.empty()){
  	// if cluster is empty, just insert all the mcell, each as a cluster
       	for (int j=0;j!=allmcell.size();j++){
  	  GeomCluster *cluster = new GeomCluster(*((MergeGeomCell*)allmcell[j]));
  	  cluster_set.insert(cluster);
  	}
      }else{
  	for (int j=0;j!=allmcell.size();j++){
  	  int flag = 0;
  	  int flag_save = 0;
  	  GeomCluster *cluster_save = 0;
	  
  	  cluster_delset.clear();

  	  // for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  	  //   if (i==318)
  	  //     cout << "b " << (*it)->get_allcell().size() << endl;
  	  // } 
	  

  	  // loop through merged cell
  	  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  	    //loop through clusters
	   
  	    flag += (*it)->AddCell(*((MergeGeomCell*)allmcell[j]));
  	    if (flag==1 && flag != flag_save){
  	      cluster_save = *it;
  	    }else if (flag>1 && flag != flag_save){
  	      cluster_save->MergeCluster(*(*it));
  	      cluster_delset.insert(*it);
  	    }
  	    flag_save = flag;
  	    // if (i==318)
  	    //   cout << "c " << flag << endl;
  	  }

  	  for (auto it = cluster_delset.begin();it!=cluster_delset.end();it++){
  	    cluster_set.erase(*it);
  	    delete (*it);
  	  }
	  
	  

  	  // if (i==318)
  	  //   cout << j << " " << flag << endl;
  	  if (flag==0){
  	    GeomCluster *cluster = new GeomCluster(*((MergeGeomCell*)allmcell[j]));
  	    cluster_set.insert(cluster);
  	  }

  	  // for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  	  //   if (i==318)
  	  //     cout << (*it)->get_allcell().size() << endl;
  	  // }
	  
  	}
      }
      

      for (int j=0;j!=allcell.size();j++){
  	Point p = allcell[j]->center();
  	x[ncount] = i*0.32;
  	y[ncount] = p.y/units::cm;
  	z[ncount] = p.z/units::cm;
  	ncount ++;
      }


      int ncount_mcell_cluster = 0;
      for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  	ncount_mcell_cluster += (*it)->get_allcell().size();
      }
      ncount_mcell += allmcell.size();
      
      
      int ncells_qx = 0;
      for (int j=0;j!=allmcell.size();j++){
	MergeGeomCell* mcell = (MergeGeomCell*)allmcell.at(j);
	ncells_qx += mcell->get_allcell().size();
      }

      cout << i << " " << allcell.size() << " " << ncells_qx << " " << allmcell.size() << " " << allmwire.size() << " " << cluster_set.size()  << endl;
      

      // for (int j=0;j!=allmwire.size();j++){
      // 	cout << mergetiling.cells(*allmwire[j]).size() << endl;
      // }


      CellChargeMap ccmap = truthtiling[i]->ccmap();

      Double_t charge_min = 10000;
      Double_t charge_max = 0;

      for (auto it = ccmap.begin();it!=ccmap.end(); it++){
  	Point p = it->first->center();
      	xt[ncount_t] = i*0.32;
      	yt[ncount_t] = p.y/units::cm;
      	zt[ncount_t] = p.z/units::cm;
      	ncount_t ++;

  	double charge = it->second;
  	if (charge > charge_max) charge_max = charge;
  	if (charge < charge_min) charge_min = charge;
       	// cout << it->second << endl;
      }


      //loop through merged cell and compare with truth cells
      for (int j=0;j!=allmcell.size();j++){
  	MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
  	mcell->CheckContainTruthCell(ccmap);
  	// 	cout << mergetiling.wires(*allmcell[j]).size() << endl;
      }


      
      // WireChargeMap wcmap = toytiling.wcmap();
      // for (auto it = wcmap.begin();it!=wcmap.end(); it++){
      // 	double charge = it->second;
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
    // cout << i << " " << allcell.size() << " " << allwire.size() << endl;
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
    
    

    // display.draw_slice(slice,"");
    // display.draw_cells(toytiling[i]->get_allcell(),"*same");
    // display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell
    // display.draw_truthcells(ccmap,"*same");
    
    // // display.draw_wires_charge(wcmap,"Fsame",FI);
    // // display.draw_cells_charge(toytiling.get_allcell(),"Fsame");
    // // display.draw_truthcells_charge(ccmap,"lFsame",FI);
    
    
    // theApp.Run();
    }
  }

  int ncount_mcell_cluster = 0;
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    ncount_mcell_cluster += (*it)->get_allcell().size();
  }


  cout << "Summary: " << ncount << " " << ncount_mcell << " " << ncount_mcell_cluster << endl;
  TGraph2D *g = new TGraph2D(ncount,x,y,z);
  TGraph2D *gt = new TGraph2D(ncount_t,xt,yt,zt);
  TFile *file = new TFile("shower3D.root","RECREATE");
  g->Write("shower3D");
  gt->Write("shower3D_truth");

  //save cluster
  int ncluster = 0;
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    ncount = 0;
    for (int i=0; i!=(*it)->get_allcell().size();i++){
      const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
      for (int j=0; j!=mcell->get_allcell().size();j++){
  	Point p = mcell->get_allcell().at(j)->center();
  	x[ncount] = mcell->GetTimeSlice()*0.32;
  	y[ncount] = p.y/units::cm;
  	z[ncount] = p.z/units::cm;
  	ncount ++;
      }
    }
    // cout << ncount << endl;
    TGraph2D *g1 = new TGraph2D(ncount,x,y,z);
    g1->Write(Form("cluster_%d",ncluster));
    ncluster ++;
  }

  file->Write();
  file->Close();

  return 0;
  
} // main()
