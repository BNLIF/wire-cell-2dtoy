#include "WCPSst/GeomDataSource.h"
//#include "WCPSst/ToyuBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/TruthToyTiling.h"
#include "WCP2dToy/SimpleBlobToyTiling.h"


#include "WCPData/MergeGeomCell.h"
#include "WCPData/GeomCluster.h"
//#include "WCPNav/SliceDataSource.h"
#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/ToyMatrixIterate.h"
#include "WCP2dToy/ToyMatrixMarkov.h"
#include "WCP2dToy/ToyMetric.h"
#include "WCP2dToy/BlobMetric.h"


#include "WCPNav/FrameDataSource.h"
#include "WCPNav/SimDataSource.h"
#include "WCPNav/SliceDataSource.h"
#include "WCPSst/Util.h"
#include "WCPData/SimTruth.h"
#include "WCP2dToy/ToyDepositor.h"
#include "WCPNav/GenerativeFDS.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TColor.h"
#include <iostream>

using namespace WCP;
using namespace std;

double rms(double a1, double a2, double a3){
  double ave = (a1 + a2 + a3)/3.;
  double rms1 = sqrt(pow(a1-ave,2) + pow(a2-ave,2));
  double rms2 = sqrt(pow(a2-ave,2) + pow(a3-ave,2));
  double rms3 = sqrt(pow(a1-ave,2) + pow(a3-ave,2));
  
  double rms;
  if (rms1 < rms2){
    rms = rms1;
  }else{
    rms = rms2;
  } 
  if (rms > rms3){
    rms = rms3;
  }
  
  //double rms  = sqrt((pow(a1-ave,2) + pow(a2-ave,2) + pow(a3-ave,2))/2.);
  return rms;
}


int main(int argc, char* argv[])
{
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root" << endl;
    return 1;
  }
  
  
  WCPSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;
  
  
  const char* root_file = argv[2];
  const char* tpath = "/Event/Sim";
  
  WCP::FrameDataSource* fds = 0;
  fds = WCPSst::make_fds(root_file);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  
  
  
  WCP::ToyDepositor toydep(fds);
  const PointValueVector pvv = toydep.depositions(1);
  
  
   
  WCP::GenerativeFDS gfds(toydep,gds,2400,5,2.0*1.6*units::millimeter);
  gfds.jump(1);

  
  WCPSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons

  

  const int N = 100000;
  Double_t x[N],y[N],z[N];
  Double_t xt[N],yt[N],zt[N];
  int ncount = 0;
  int ncount_t = 0;
  

  WCP2dToy::ToyTiling **toytiling = new WCP2dToy::ToyTiling*[2400];
  WCP2dToy::MergeToyTiling **mergetiling = new WCP2dToy::MergeToyTiling*[2400];
  WCP2dToy::TruthToyTiling **truthtiling = new WCP2dToy::TruthToyTiling*[2400];
  WCP2dToy::SimpleBlobToyTiling **blobtiling = new WCP2dToy::SimpleBlobToyTiling*[2400];
  
  WCP2dToy::ToyMatrix **toymatrix = new WCP2dToy::ToyMatrix*[2400];
  WCP2dToy::ToyMatrixIterate **toymatrix_it = new WCP2dToy::ToyMatrixIterate*[2400];
  WCP2dToy::ToyMatrixMarkov **toymatrix_markov = new WCP2dToy::ToyMatrixMarkov*[2400];
  
  //WCP2dToy::ToyMatrix **blobmatrix = new WCP2dToy::ToyMatrix*[2400];
  
  //WCP2dToy::ToyMatrixMarkov **blobmatrix_markov = new WCP2dToy::ToyMatrixMarkov*[2400];
  WCP2dToy::ToyMetric toymetric;
  WCP2dToy::BlobMetric blobmetric;
 
  //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  
  GeomCellSelection total_cells;
  GeomCellSelection total_recon_cells;
  GeomCellSelection total_corner_cells;
  GeomCellSelection total_blob_cells;
  CellChargeMap total_ccmap;
  
  CellChargeMap total_scrms;
  CellChargeMap total_scmap;


  Double_t charge_min = 10000;
  Double_t charge_max = 0;
    
  int ncount_mcell = 0;
  
  //simple cosmic
  int start_num =318;
  int end_num = 319;


  //nue cc 
  // int start_num =355;
  // int end_num = 357;
    
  //delta 
  // int start_num =678;
  // int end_num = 682;

  //complicated blob
  // int start_num = 454;
  // int end_num = 454;
  
  WCP::Slice slice;
  for (int i=start_num;i!=end_num+1;i++){
    sds.jump(i);
    slice = sds.get();
    
    
    toytiling[i] = new WCP2dToy::ToyTiling(slice,gds);
    mergetiling[i] = new WCP2dToy::MergeToyTiling(*toytiling[i],i);
    truthtiling[i] = new WCP2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds);
    GeomCellSelection allcell = toytiling[i]->get_allcell();
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    GeomWireSelection allmwire = mergetiling[i]->get_allwire();
    cout << i << " " << allcell.size() << " " << allmcell.size() << " " << allmwire.size()  << endl;
    
    for (int j=0;j!=allcell.size();j++){
      total_cells.push_back(allcell.at(j));
    }


    toymatrix[i] = new WCP2dToy::ToyMatrix(*toytiling[i],*mergetiling[i]);
    if (toymatrix[i]->Get_Solve_Flag()==0)
      //toymatrix_markov[i] = new WCP2dToy::ToyMatrixMarkov(toymatrix[i],&allmcell);    
      toymatrix_it[i] = new WCP2dToy::ToyMatrixIterate(*toymatrix[i]);
      
    cout << "chi2: " << toymatrix[i]->Get_Chi2() << endl;
    cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;

    toymatrix[i]->JudgeSimpleBlob(*toytiling[i],*mergetiling[i]);
    
    //save stuff for display ... 
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell =(MergeGeomCell*)allmcell.at(j);
      double charge =toymatrix[i]->Get_Cell_Charge(mcell);
      if (charge > 2000){
	for (int k=0;k!=mcell->get_allcell().size();k++){
	  total_recon_cells.push_back(mcell->get_allcell().at(k));
	}
      }
      if (mcell->IsSimpleBlob() == true){
	GeomCellSelection corners = mcell->get_cornercells();
	CellIndexMap indexmap = mcell->get_cornercells_index();
	for (int k=0;k!=corners.size();k++){
	  if (indexmap[corners.at(k)]>=2){
	    total_corner_cells.push_back(corners.at(k));
	  }
	}
     	//total_corner_cells.insert(total_corner_cells.end(),corners.begin(),corners.end());
      }
    }
    
    
    // blobmatrix[i] = new WCP2dToy::ToyMatrix(*toytiling[i],*blobtiling[i]);

    // GeomCellSelection ballmcell = blobtiling[i]->get_allcell();
    // GeomWireSelection ballmwire = blobtiling[i]->get_allwire();
    // cout << i << " " << ballmcell.size() << " " << ballmwire.size()  << endl;
    
    // if (blobmatrix[i]->Get_Solve_Flag()==0)
    //   blobmatrix_markov[i] = new WCP2dToy::ToyMatrixMarkov(blobmatrix[i],&ballmcell);
    
    // cout << "chi2: " << blobmatrix[i]->Get_Chi2() << endl;
    // cout << "NDF: " << blobmatrix[i]->Get_ndf() << endl;
    
    CellChargeMap ccmap = truthtiling[i]->ccmap();
    total_ccmap.insert(ccmap.begin(),ccmap.end());

    if (toymatrix[i]->Get_Solve_Flag()!=0)
	toymetric.Add(allmcell,*toymatrix[i],ccmap);
    toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
    
    
    
    
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
    
  }
  
  //use time information
  for (int i=start_num;i!=end_num+1;i++){
    cout << i << endl;
    if (toymatrix[i]->GetSimpleBlobReduction()){
      if (i==start_num){
	blobtiling[i] = new WCP2dToy::SimpleBlobToyTiling(*toytiling[i],*mergetiling[i],*toymatrix[i],*mergetiling[i+1],*toymatrix[i+1],*mergetiling[i+1],*toymatrix[i+1]);
      }else if (i==end_num){
	blobtiling[i] = new WCP2dToy::SimpleBlobToyTiling(*toytiling[i],*mergetiling[i],*toymatrix[i],*mergetiling[i-1],*toymatrix[i-1],*mergetiling[i-1],*toymatrix[i-1]);
      }else{
	blobtiling[i] = new WCP2dToy::SimpleBlobToyTiling(*toytiling[i],*mergetiling[i],*toymatrix[i],*mergetiling[i-1],*toymatrix[i-1],*mergetiling[i+1],*toymatrix[i+1]);
      }
      
      
      //save stuff
      GeomCellSelection blob_cells = blobtiling[i]->get_allcell();
      // cout << blob_cells.size() << endl;
      for (int k=0;k!=blob_cells.size();k++){
	total_blob_cells.push_back(blob_cells.at(k));
      }
      CellChargeMap ccmap = truthtiling[i]->ccmap();
      blobmetric.Add(*blobtiling[i],ccmap);
      blobmetric.Print();
    }
  }
  


  toymetric.Print();
  std::cout << std::endl;
  blobmetric.Print();
  

  //do display
  
  TApplication theApp("theApp",&argc,argv);
  theApp.SetReturnFromRun(true);
  
  TCanvas c1("ToyMC","ToyMC",800,600);
  //TCanvas c2("ToyMC","ToyMC",800,600);
  
  WCP2dToy::ToyEventDisplay display(c1, gds);
  //WCP2dToy::ToyEventDisplay display1(c2, gds);
  display.charge_min = charge_min;
  display.charge_max = charge_max;
  // display.charge_min = 0.;
  // display.charge_max = 1.;
  
  // display1.charge_min = charge_min;
  // display1.charge_max = charge_max;
  
  gStyle->SetOptStat(0);
  
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Int_t MyPalette[NCont];
  Double_t stops[NRGBs] = {0.0, 0.34, 0.61, 0.84, 1.0};
  Double_t red[NRGBs] = {0.0, 0.0, 0.87 ,1.0, 0.51};
  Double_t green[NRGBs] = {0.0, 0.81, 1.0, 0.2 ,0.0};
  Double_t blue[NRGBs] = {0.51, 1.0, 0.12, 0.0, 0.0};
  Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  for (int kk=0;kk!=NCont;kk++) MyPalette[kk] = FI+kk;
  gStyle->SetPalette(NCont,MyPalette);
  
  display.init(0,10.3698,-2.33/2.,2.33/2.);
  //display.init(8.61,8.66,0.0,0.06);
  // display1.init(1.4,1.65,0.7,1.0);
  
  
  display.draw_mc(1,WCP::PointValueVector(),"colz");
  
  display.draw_slice(slice,"");
  display.draw_cells(total_cells,"*same");
  //display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",1); //0 is normal, 1 is only draw the ones containt the truth cell
  
  display.draw_cells(total_recon_cells,"*same",1);
  //display.draw_truthcells_charge(total_ccmap,"*same",FI);
  //display.draw_truthcells_charge(total_scmap,"*same",FI);
  //display.draw_truthcells_charge(total_scrms,"*same",FI);
  display.draw_cells(total_blob_cells,"*same",6);
  // display.draw_cells(total_corner_cells,"*same",2);
  display.draw_truthcells(total_ccmap,"*same");
 


  //display.draw_reconcells(mergetiling[i]->get_allcell(),toymatrix[i],"*same",1);
  //display.draw_reconcells(blobtiling[i]->get_allcell(),blobmatrix[i],"*same",2);
  
  // display1.draw_mc(1,WCP::PointValueVector(),"colz");
  // display1.draw_truthcells_charge(total_ccmap,"*same",FI);
  

  theApp.Run();
  
  
  
  // TFile *file = new TFile("shower3D.root","RECREATE");
  // TGraph2D *g = new TGraph2D();
  // TGraph2D *gt = new TGraph2D();
  // TGraph2D *g_rec = new TGraph2D();
  // Double_t x_save, y_save, z_save;
  // Double_t charge_save;
  // Double_t chi2_save;
  // Double_t ndf_save;
  // Double_t ncell_save;
  // ncount = 0;
  // int ncount1 = 0;
  // ncount_t = 0;
  // //save results 
  // for (int i=start_num;i!=end_num+1;i++){
  //   //truth
  //   CellChargeMap ccmap = truthtiling[i]->ccmap();
  //   for (auto it = ccmap.begin();it!=ccmap.end(); it++){
  //     Point p = it->first->center();
  //     x_save = i*0.32;
  //     y_save = p.y/units::cm;
  //     z_save = p.z/units::cm;
  //     charge_save = it->second;
      
  //     gt->SetPoint(ncount_t,x_save,y_save,z_save);
  //     //t_true->Fill();
            
  //     ncount_t ++;
  //   }
    
  //   //recon 1
  //   GeomCellSelection allcell = toytiling[i]->get_allcell();
  //   for (int j=0;j!=allcell.size();j++){
      
  //     Point p = allcell[j]->center();
  //     x_save = i*0.32;
  //     y_save = p.y/units::cm;
  //     z_save = p.z/units::cm;
      

  //     g->SetPoint(ncount,x_save,y_save,z_save);
  //     //t_rec->Fill();

  //     ncount ++;
  //   }

  //   //recon 2 with charge
  //   GeomCellSelection allmcell = mergetiling[i]->get_allcell();
  //   for (int j=0;j!=allmcell.size();j++){
  //     MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
  //     double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
  //     if (charge>2000){
  // 	//truth
  // 	for (int k=0;k!=mcell->get_allcell().size();k++){
  // 	  Point p = mcell->get_allcell().at(k)->center();
  // 	  x_save = i*0.32;
  // 	  y_save = p.y/units::cm;
  // 	  z_save = p.z/units::cm;
  // 	  charge_save = charge/mcell->get_allcell().size();
  // 	  ncell_save = mcell->get_allcell().size();
  // 	  chi2_save = toymatrix[i]->Get_Chi2();
  // 	  ndf_save = toymatrix[i]->Get_ndf();

  // 	  g_rec->SetPoint(ncount1,x_save,y_save,z_save);
  // 	  //  t_rec_charge->Fill();
	  
  // 	  ncount1 ++;
  // 	}
  //     }
  //   }

  //   //save all results
  //   // file->Write(Form("toytiling_%d",i),toytiling[i]);
  //   // file->Write(Form("mergetiling_%d",i),mergetiling[i]);
  //   // file->Write(Form("truthtiling_%d",i),truthtiling[i]);
  //   // file->Write(Form("toymatrix_%d",i),toymatrix[i]);

  // }

  // g->Write("shower3D");
  // gt->Write("shower3D_truth");
  // g_rec->Write("shower3D_charge");
  // file->Write();
  // file->Close();
  
  
  return 0;
  
} // main()
