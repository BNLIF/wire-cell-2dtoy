#include "WireCellSst/GeomDataSource.h"
#include "WireCellSst/ToyuBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCell2dToy/SimpleBlobToyTiling.h"

#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/ToyMatrixExclusive.h"
#include "WireCell2dToy/ToyMatrixKalman.h"
#include "WireCell2dToy/ToyMatrixIterate.h"
#include "WireCell2dToy/ToyMatrixMarkov.h"
#include "WireCell2dToy/ToyMetric.h"
#include "WireCell2dToy/BlobMetric.h"

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
#include "WireCell2dToy/ToySignalSimu.h"
#include "WireCell2dToy/ToySignalSimuTrue.h"
#include "WireCell2dToy/ToySignalGaus.h"
#include "WireCell2dToy/ToySignalWien.h"

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
  

  int recon_threshold = 2000;
  int max_events = 5;
  int eve_num = 1;
  // WireCell::ToyDepositor toydep(fds);
  // const PointValueVector pvv = toydep.depositions(eve_num);
  
  // WireCell::GenerativeFDS gfds(toydep,gds,2400,max_events,2.0*1.6*units::millimeter);
  // gfds.jump(eve_num);

  // WireCellSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons

  WireCell::ToyDepositor toydep(fds);
  const PointValueVector pvv = toydep.depositions(eve_num);
  //WireCell::GenerativeFDS gfds(toydep,gds,9600,max_events,0.5*1.605723*units::millimeter); // 87 K at 0.5 kV/cm
  WireCell::GenerativeFDS gfds(toydep,gds,9600,max_events,0.5*1.60*units::millimeter); // 87 K at 0.5 kV/cm
  WireCell2dToy::ToySignalSimuFDS simu_fds(gfds,gds,9600,max_events,1.647,1.539+1.647,1); // time offset among different planes for the time electrons travel among different planes
  simu_fds.jump(eve_num);
  //simu_fds.Save();

  //  WireCell2dToy::ToySignalSimuTrueFDS st_fds(gfds,gds,9600/4,5); //truth
  WireCell::GenerativeFDS st_fds(toydep,gds,9600/4,max_events,0.5*1.60*units::millimeter); // 87 K at 0.5 kV/cm
  st_fds.jump(eve_num);
  // st_fds.Save();
  
  WireCell2dToy::ToySignalGausFDS gaus_fds(simu_fds,gds,9600/4,max_events,1.647,1.539+1.647); // gaussian smearing for charge estimation
  gaus_fds.jump(eve_num);
  //gaus_fds.Save();

   WireCell2dToy::ToySignalWienFDS wien_fds(simu_fds,gds,9600/4,max_events,1.647,1.539+1.647); // weiner smearing for hit identification
  wien_fds.jump(eve_num);
  //wien_fds.Save();
  
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  int nwire_u = wires_u.size();
  int nwire_v = wires_v.size();
  int nwire_w = wires_w.size();
  
  float threshold_u = 5.87819e+02 * 4.0;
  float threshold_v = 8.36644e+02 * 4.0;
  float threshold_w = 5.67974e+02 * 4.0;

  float threshold_ug = 755.96;
  float threshold_vg = 822.81;
  float threshold_wg = 510.84;
  
  // float threshold_u = 1000;
  // float threshold_v = 1000;
  // float threshold_w = 1000;
  

  WireCellSst::ToyuBooNESliceDataSource sds(wien_fds,gaus_fds,threshold_u, 
					    threshold_v, threshold_w, 
					    threshold_ug, 
					    threshold_vg, threshold_wg, 
					    nwire_u, 
					    nwire_v, nwire_w); 

  WireCellSst::ToyuBooNESliceDataSource sds_th(st_fds,st_fds,1, 
					    1, 1, 
					    threshold_ug, 
					    threshold_vg, threshold_wg, 
					    nwire_u, 
					    nwire_v, nwire_w); 

  // const int N = 100000;
  // Double_t x[N],y[N],z[N];
  // Double_t x1[N],y1[N],z1[N], charge_r1[N];
  // Double_t xt[N],yt[N],zt[N], charge_t[N];
  int ncount = 0;
  int ncount1 = 0;  
  int ncount2 = 0;

  int ncount_t = 0;
  

  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
  WireCell2dToy::MergeToyTiling **mergetiling = new WireCell2dToy::MergeToyTiling*[2400];
  WireCell2dToy::TruthToyTiling **truthtiling = new WireCell2dToy::TruthToyTiling*[2400];
  WireCell2dToy::SimpleBlobToyTiling **blobtiling = new WireCell2dToy::SimpleBlobToyTiling*[2400];

  WireCell2dToy::ToyMatrix **toymatrix = new WireCell2dToy::ToyMatrix*[2400];
  WireCell2dToy::ToyMatrixIterate **toymatrix_it = new WireCell2dToy::ToyMatrixIterate*[2400];
  WireCell2dToy::ToyMatrixMarkov **toymatrix_markov = new WireCell2dToy::ToyMatrixMarkov*[2400];
  
  //save truth ...
  WireCell2dToy::ToyTiling **toytiling_th = new WireCell2dToy::ToyTiling*[2400];
  WireCell2dToy::TruthToyTiling **truthtiling_th = new WireCell2dToy::TruthToyTiling*[2400];
 
  WireCell2dToy::ToyMetric toymetric;
  WireCell2dToy::BlobMetric blobmetric;

  //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  
  int ncount_mcell = 0;

  int start_num = 0 ;
  int end_num = sds.size()-1;

  // int start_num =1241;
  // int end_num = 1243;
  // int end_num = sds.size()-1;

  // int start_num = 400;
  // int end_num = 462;

  //  int i=454;{ // 46, 26
  // int i=329;{  // 18, 6,
  // int i=344;{ // 29 14
  //int i=459;{ // 23, 8,   5e5 
  //int i = 351;{
  //for (int i=0;i!=sds.size();i++){
  for (int i=start_num;i!=end_num+1;i++){
 
    sds.jump(i);
    sds_th.jump(i);
    WireCell::Slice slice = sds.get();
    WireCell::Slice slice_th = sds_th.get();
    //if ( slice.group().size() >0){
      
    toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
    

      mergetiling[i] = new WireCell2dToy::MergeToyTiling(*toytiling[i],i);

      GeomCellSelection allcell = toytiling[i]->get_allcell();
      GeomWireSelection allwire = toytiling[i]->get_allwire();
      GeomCellSelection allmcell = mergetiling[i]->get_allcell();
      GeomWireSelection allmwire = mergetiling[i]->get_allwire();
      
      cout << i << " " << allmcell.size() << " " << allmwire.size() << endl;

      truthtiling[i] = new WireCell2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds,800);
      toymatrix[i] = new WireCell2dToy::ToyMatrix(*toytiling[i],*mergetiling[i]);
      if (toymatrix[i]->Get_Solve_Flag()==0)
      	toymatrix_it[i] = new WireCell2dToy::ToyMatrixIterate(*toymatrix[i]);
      
      cout << "chi2: " << toymatrix[i]->Get_Chi2() << endl;
      cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
      
      toytiling_th[i] = new WireCell2dToy::ToyTiling(slice_th,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
      truthtiling_th[i] = new WireCell2dToy::TruthToyTiling(*toytiling_th[i],pvv,i,gds,800);

      // GeomCellSelection calmcell;
      // for (int j=0;j!=allmcell.size();j++){
      // 	MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      // 	double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
      // 	double charge_err = toymatrix[i]->Get_Cell_Charge(mcell,2);
	
      // 	//	cout << "Recon: " << j << " " << charge << " " << charge_err << endl;

      // 	if (charge > 2000) calmcell.push_back(mcell);
      // }
      
      //



      CellChargeMap ccmap = truthtiling[i]->ccmap();
      if (toymatrix[i]->Get_Solve_Flag()!=0)
	toymetric.Add(allmcell,*toymatrix[i],ccmap);

      toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());

      Double_t charge_min = 10000;
      Double_t charge_max = 0;

     
      
      
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
      // }
  }

  toymetric.Print();
  std::cout << "Starting MCMC" << std::endl;

  // //without  time information
  // for (int i=start_num;i!=end_num+1;i++){
  //   if (toymatrix[i]->Get_Solve_Flag()==0){
  //     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
  //     toymatrix_markov[i] = new WireCell2dToy::ToyMatrixMarkov(toymatrix[i],&allmcell);
  //     CellChargeMap ccmap = truthtiling[i]->ccmap();
  //     if (toymatrix[i]->Get_Solve_Flag()!=0)
  // 	toymetric.Add(allmcell,*toymatrix[i],ccmap);
  //     toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
  //     cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
  //     cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
  //   }
  // }
    



  //with time information
  if (start_num != end_num){
    int first_solve;
    for (int i=start_num; i!=end_num+1;i++){
      if (toymatrix[i]->Get_Solve_Flag()!=0){
  	first_solve = i;
  	break;
      }
    }
  

    for (int i=first_solve+1;i<=end_num-1;i++){
      if (toymatrix[i]->Get_Solve_Flag()==0){
  	GeomCellSelection allmcell = mergetiling[i]->get_allcell();
  	toymatrix_markov[i] = new WireCell2dToy::ToyMatrixMarkov(*toymatrix[i-1],*toymatrix[i],*toymatrix[i+1],*mergetiling[i-1],*mergetiling[i],*mergetiling[i+1],&allmcell);
  	CellChargeMap ccmap = truthtiling[i]->ccmap();
  	if (toymatrix[i]->Get_Solve_Flag()!=0)
  	  toymetric.Add(allmcell,*toymatrix[i],ccmap);
  	toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());

  	cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
  	cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;

  	//toymetric.Print();
      }
    }
     
    if (toymatrix[end_num]->Get_Solve_Flag()==0){
      GeomCellSelection allmcell = mergetiling[end_num]->get_allcell();
      toymatrix_markov[end_num] = new WireCell2dToy::ToyMatrixMarkov(*toymatrix[end_num-1],*toymatrix[end_num],*toymatrix[end_num-1],*mergetiling[end_num-1],*mergetiling[end_num],*mergetiling[end_num-1],&allmcell);

      
      CellChargeMap ccmap = truthtiling[end_num]->ccmap();
      if (toymatrix[end_num]->Get_Solve_Flag()!=0)
  	toymetric.Add(allmcell,*toymatrix[end_num],ccmap);
      toymetric.AddSolve(toymatrix[end_num]->Get_Solve_Flag());

      cout << "chi2: " << end_num << " " << toymatrix[end_num]->Get_Chi2() << endl;
      cout << "NDF: " << toymatrix[end_num]->Get_ndf() << endl;
    }
      
    // go to early ones 
     for (int i=first_solve-1;i>=start_num+1;i--){
      if (toymatrix[i]->Get_Solve_Flag()==0){
  	GeomCellSelection allmcell = mergetiling[i]->get_allcell();
  	toymatrix_markov[i] = new WireCell2dToy::ToyMatrixMarkov(*toymatrix[i-1],*toymatrix[i],*toymatrix[i+1],*mergetiling[i-1],*mergetiling[i],*mergetiling[i+1],&allmcell);
  	
  	CellChargeMap ccmap = truthtiling[i]->ccmap();
  	if (toymatrix[i]->Get_Solve_Flag()!=0)
  	  toymetric.Add(allmcell,*toymatrix[i],ccmap);
  	toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	
  	cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
  	cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
      }
    }
     
    if (toymatrix[start_num]->Get_Solve_Flag()==0){
      GeomCellSelection allmcell = mergetiling[start_num]->get_allcell();
      toymatrix_markov[start_num] = new WireCell2dToy::ToyMatrixMarkov(*toymatrix[start_num+1],*toymatrix[start_num],*toymatrix[start_num+1],*mergetiling[start_num+1],*mergetiling[start_num],*mergetiling[start_num+1],&allmcell);

      
      CellChargeMap ccmap = truthtiling[start_num]->ccmap();
      if (toymatrix[start_num]->Get_Solve_Flag()!=0)
  	toymetric.Add(allmcell,*toymatrix[start_num],ccmap);
      toymetric.AddSolve(toymatrix[start_num]->Get_Solve_Flag());

      cout << "chi2: " << start_num << " " << toymatrix[start_num]->Get_Chi2() << endl;
      cout << "NDF: " << toymatrix[start_num]->Get_ndf() << endl;
    }
  }


  


  //do blob thing ... 
  //use time information
  std::cout << "Reduce Blob" << std::endl; 
  for (int i=start_num;i!=end_num+1;i++){
    std::cout << "Check Blob " << i << std::endl;
    //std::cout << toymatrix[i]->GetSimpleBlobReduction() << std::endl;
    toymatrix[i]->JudgeSimpleBlob(*toytiling[i],*mergetiling[i]);
    //std::cout << toymatrix[i]->GetSimpleBlobReduction() << std::endl;
    if (toymatrix[i]->GetSimpleBlobReduction()){
      if (i==start_num){
  	blobtiling[i] = new WireCell2dToy::SimpleBlobToyTiling(*toytiling[i],*mergetiling[i],*toymatrix[i],*mergetiling[i+1],*toymatrix[i+1],*mergetiling[i+1],*toymatrix[i+1]);
      }else if (i==end_num){
  	blobtiling[i] = new WireCell2dToy::SimpleBlobToyTiling(*toytiling[i],*mergetiling[i],*toymatrix[i],*mergetiling[i-1],*toymatrix[i-1],*mergetiling[i-1],*toymatrix[i-1]);
      }else{
  	blobtiling[i] = new WireCell2dToy::SimpleBlobToyTiling(*toytiling[i],*mergetiling[i],*toymatrix[i],*mergetiling[i-1],*toymatrix[i-1],*mergetiling[i+1],*toymatrix[i+1]);
      }
      
      //save stuff
      CellChargeMap ccmap = truthtiling[i]->ccmap();
      blobmetric.Add(*blobtiling[i],ccmap);
      
      std::cout << "Check Blob " << i << std::endl;
      WireCell2dToy::BlobMetric tempblob;
      tempblob.Add(*blobtiling[i],ccmap);
      tempblob.Print();
    }
  }
  
  //do clustering ... 
   for (int i=start_num;i!=end_num+1;i++){
     GeomCellSelection pallmcell = mergetiling[i]->get_allcell();
     GeomCellSelection allmcell;
     for (int j=0;j!=pallmcell.size();j++){
       const GeomCell* mcell = pallmcell[j];
       if (toymatrix[i]->Get_Cell_Charge(mcell)> recon_threshold){
	 allmcell.push_back(mcell);
       }
     }
     
     
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
  	   
	 }
	 
	 for (auto it = cluster_delset.begin();it!=cluster_delset.end();it++){
	   cluster_set.erase(*it);
	   delete (*it);
	 }
	 
	 if (flag==0){
	   GeomCluster *cluster = new GeomCluster(*((MergeGeomCell*)allmcell[j]));
	   cluster_set.insert(cluster);
	 }
	 
       }
     }

     int ncount_mcell_cluster = 0;
      for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  	ncount_mcell_cluster += (*it)->get_allcell().size();
      }
      ncount_mcell += allmcell.size();
      cout << i << " " << allmcell.size()  << " " << cluster_set.size()  << endl;
   }

   int ncount_mcell_cluster = 0;
   for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
     ncount_mcell_cluster += (*it)->get_allcell().size();
   }
   cout << "Summary: " << ncount << " " << ncount_mcell << " " << ncount_mcell_cluster << endl;


  TFile *file = new TFile("shower3D_signal.root","RECREATE");
  TTree *t_true = new TTree("T_true","T_true");
  TTree *t_rec = new TTree("T_rec","T_rec");
  TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
  TTree *t_rec_charge_blob = new TTree("T_rec_charge_blob","T_rec_charge_blob");

  Double_t x_save, y_save, z_save;
  Double_t charge_save;
  Double_t ncharge_save;
  Double_t chi2_save;
  Double_t ndf_save;

  t_true->SetDirectory(file);
  t_true->Branch("x",&x_save,"x/D");
  t_true->Branch("y",&y_save,"y/D");
  t_true->Branch("z",&z_save,"z/D");
  t_true->Branch("q",&charge_save,"q/D");
  
  t_rec->SetDirectory(file);
  t_rec->Branch("x",&x_save,"x/D");
  t_rec->Branch("y",&y_save,"y/D");
  t_rec->Branch("z",&z_save,"z/D");
  
  t_rec_charge->SetDirectory(file);
  t_rec_charge->Branch("x",&x_save,"x/D");
  t_rec_charge->Branch("y",&y_save,"y/D");
  t_rec_charge->Branch("z",&z_save,"z/D");
  t_rec_charge->Branch("q",&charge_save,"q/D");
  t_rec_charge->Branch("nq",&ncharge_save,"nq/D");
  t_rec_charge->Branch("chi2",&chi2_save,"chi2/D");
  t_rec_charge->Branch("ndf",&ndf_save,"ndf/D");

  //blob stuff
  t_rec_charge_blob->SetDirectory(file);
  t_rec_charge_blob->Branch("x",&x_save,"x/D");
  t_rec_charge_blob->Branch("y",&y_save,"y/D");
  t_rec_charge_blob->Branch("z",&z_save,"z/D");
  t_rec_charge_blob->Branch("q",&charge_save,"q/D");
  t_rec_charge_blob->Branch("nq",&ncharge_save,"nq/D");
  
  TGraph2D *g = new TGraph2D();
  TGraph2D *gt = new TGraph2D();
  TGraph2D *g_rec = new TGraph2D();
  TGraph2D *g_rec_blob = new TGraph2D();

  //save results 
  for (int i=start_num;i!=end_num+1;i++){
    //truth
    CellChargeMap ccmap = truthtiling_th[i]->ccmap();
    for (auto it = ccmap.begin();it!=ccmap.end(); it++){
      Point p = it->first->center();
      x_save = i*0.32 - 256;
      y_save = p.y/units::cm;
      z_save = p.z/units::cm;
      charge_save = it->second;
      
      gt->SetPoint(ncount_t,x_save,y_save,z_save);
      t_true->Fill();
            
      ncount_t ++;
    }
    
    //recon 1
    GeomCellSelection allcell = toytiling[i]->get_allcell();
    for (int j=0;j!=allcell.size();j++){
      Point p = allcell[j]->center();
      x_save = i*0.32- 256;
      y_save = p.y/units::cm;
      z_save = p.z/units::cm;
      

      g->SetPoint(ncount,x_save,y_save,z_save);
      t_rec->Fill();

      ncount ++;
    }

    //recon 2 with charge
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
      if (charge> recon_threshold){
	//truth
	for (int k=0;k!=mcell->get_allcell().size();k++){
	  Point p = mcell->get_allcell().at(k)->center();
	  x_save = i*0.32- 256;
	  y_save = p.y/units::cm;
	  z_save = p.z/units::cm;
	  charge_save = charge/mcell->get_allcell().size();
	  ncharge_save = mcell->get_allcell().size();
	  chi2_save = toymatrix[i]->Get_Chi2();
	  ndf_save = toymatrix[i]->Get_ndf();

	  g_rec->SetPoint(ncount1,x_save,y_save,z_save);
	  t_rec_charge->Fill();
	  
	  ncount1 ++;
	}
      }
    }

    //recon 3 with charge and deblob
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
      if (charge> recon_threshold && !(mcell->IsSimpleBlob() && mcell->IsBlob())){
  	for (int k=0;k!=mcell->get_allcell().size();k++){
  	  Point p = mcell->get_allcell().at(k)->center();
  	  x_save = i*0.32-256;
  	  y_save = p.y/units::cm;
  	  z_save = p.z/units::cm;
  	  charge_save = charge/mcell->get_allcell().size();
  	  ncharge_save = mcell->get_allcell().size();
	  
  	  g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
  	  t_rec_charge_blob->Fill();
	  
  	  ncount2 ++;
  	}
      }
    }
     if (toymatrix[i]->GetSimpleBlobReduction()){
       for (int j=0;j!=blobtiling[i]->Get_Cells().size();j++){
  	 const GeomCell *cell = blobtiling[i]->Get_Cells().at(j);
  	 Point p = cell->center();
  	 x_save = i*0.32-256;
  	 y_save = p.y/units::cm;
  	 z_save = p.z/units::cm;
  	 charge_save = blobtiling[i]->Get_Cell_Charge(cell,1);
  	 ncharge_save = 1;
	 
  	 g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
  	 t_rec_charge_blob->Fill();
	  
  	 ncount2 ++;
       }
     }
    
    //save all results
    // file->Write(Form("toytiling_%d",i),toytiling[i]);
    // file->Write(Form("mergetiling_%d",i),mergetiling[i]);
    // file->Write(Form("truthtiling_%d",i),truthtiling[i]);
    // file->Write(Form("toymatrix_%d",i),toymatrix[i]);

  }
 
  

  g->Write("shower3D");
  gt->Write("shower3D_truth");
  g_rec->Write("shower3D_charge");
  g_rec_blob->Write("shower3D_charge_blob");
  
  const int N = 100000;
  Double_t x[N],y[N],z[N];
  //save cluster
  int ncluster = 0;
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    ncount = 0;
    for (int i=0; i!=(*it)->get_allcell().size();i++){
      const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
      for (int j=0; j!=mcell->get_allcell().size();j++){
  	Point p = mcell->get_allcell().at(j)->center();
  	x[ncount] = mcell->GetTimeSlice()*0.32- 256;
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

  // WireCell2dToy::ToyTiling* tt1;
  // WireCell2dToy::MergeToyTiling* tt2;
  // WireCell2dToy::TruthToyTiling* tt3;

  // WireCell2dToy::ToyMatrix* tt4;
  // WireCell2dToy::ToyMatrixIterate* tt5;
  // WireCell2dToy::ToyMatrixMarkov* tt6;
  // WireCell2dToy::SimpleBlobToyTiling *tt7;
  // int time_slice;
  
  // TTree* ttree = new TTree("T","T");
  // ttree->Branch("time_slice",&time_slice,"time_slice/I");
  // ttree->Branch("toytiling",&tt1);
  // //ttree->Branch("truthtiling",&tt3);
  // //ttree->Branch("mergetiling",&tt2);

  
  // //ttree->Branch("toymatrix",&tt4);
  // // ttree->Branch("toymatrix_ite",&tt5);
  // // ttree->Branch("toymatrix_markov",&tt6);
  // //ttree->Branch("blobtiling",&tt7);
  

  // ttree->SetDirectory(file);
  // for (int i=start_num;i!=end_num+1;i++){
  //   tt1 = toytiling[i];
  //   //tt3 = truthtiling[i];

  //   //tt2 = mergetiling[i];
  

  //   //tt4 = toymatrix[i];
  //   // tt5 = toymatrix_it[i];
  //   // tt6 = toymatrix_markov[i];
  //   //tt7 = blobtiling[i];
    
  //   time_slice = i;
  //   ttree->Fill();
  // }

  // // TTree *ttree1 = new TTree("TC","TC");
  // // GeomCluster *cluster;
  // // ttree1->Branch("cluster",&cluster);
  // // ttree1->SetDirectory(file);
  // // for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  // //   cluster = *it;
  // //   ttree1->Fill();
  // // }
  
 


  file->Write();
  file->Close();

  toymetric.Print();
  blobmetric.Print();

  return 0;
  
} // main()
