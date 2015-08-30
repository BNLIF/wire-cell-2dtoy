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
  if (argc < 4) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num" << endl;
    return 1;
  }
  
  
  WireCellSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;
  cout << "Pitch: " << gds.pitch(WirePlaneType_t(0)) 
       << " " << gds.pitch(WirePlaneType_t(1)) 
       << " " << gds.pitch(WirePlaneType_t(2))
       << endl;
  cout << "Angle: " << gds.angle(WirePlaneType_t(0)) 
       << " " << gds.angle(WirePlaneType_t(1)) 
       << " " << gds.angle(WirePlaneType_t(2))
       << endl;
  
  
  const char* root_file = argv[2];
  const char* tpath = "/Event/Sim";
  
  int recon_threshold = 2000;
  int max_events = 100;
  int eve_num = atoi(argv[3]);
  
  float unit_dis = 1.6;

  float toffset_1=1.647;
  float toffset_2=1.539+1.647;
  float toffset_3=0;

  int total_time_bin=9600;
  int frame_length = 3200;
  int nrebin = 4;


  float threshold_u = 5.87819e+02 * 4.0;
  float threshold_v = 8.36644e+02 * 4.0;
  float threshold_w = 5.67974e+02 * 4.0;

  float threshold_ug = 755.96;
  float threshold_vg = 822.81;
  float threshold_wg = 510.84;

  int time_offset = 0;

  TFile *tfile = TFile::Open(root_file);

  TTree* sst = dynamic_cast<TTree*>(tfile->Get(tpath));

  int run_no, subrun_no, event_no;
  sst->SetBranchAddress("eventNo",&event_no);
  sst->SetBranchAddress("runNo",&run_no);
  sst->SetBranchAddress("subRunNo",&subrun_no);
  sst->GetEntry(eve_num);
  
  cout << "Run No: " << run_no << " " << subrun_no << " " << eve_num << endl;


  WireCell::FrameDataSource* fds = 0;
  fds = WireCellSst::make_fds(*tfile);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  
  

  
  // WireCell::ToyDepositor toydep(fds);
  // const PointValueVector pvv = toydep.depositions(eve_num);
  
  // WireCell::GenerativeFDS gfds(toydep,gds,2400,max_events,2.0*unit_dis*units::millimeter);
  // gfds.jump(eve_num);

  // WireCellSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons

  WireCell::ToyDepositor toydep(fds,1,unit_dis);
  const PointValueVector& pvv = toydep.depositions(eve_num);
  const std::vector<int>& timeoffset = toydep.timeoffset();

  //cout << pvv.size() << endl;

  //WireCell::GenerativeFDS gfds(toydep,gds,total_time_bin,max_events,0.5*1.605723*units::millimeter); // 87 K at 0.5 kV/cm
  WireCell::GenerativeFDS gfds(toydep,gds,total_time_bin,max_events,0.5*unit_dis*units::millimeter); // 87 K at 0.5 kV/cm
  
  cout << "Put in Truth " << endl; 
  WireCell2dToy::ToySignalSimuTrueFDS st_fds(gfds,gds,total_time_bin/nrebin,max_events,0); //truth
  //WireCell::GenerativeFDS st_fds(toydep,gds,total_time_bin/nrebin,max_events,2.0*unit_dis*units::millimeter); // 87 K at 0.5 kV/cm
  st_fds.jump(eve_num);
  // st_fds.Save();
  
  cout << "Simulate Raw WaveForm " << endl; 
  WireCell2dToy::ToySignalSimuFDS simu_fds(gfds,gds,total_time_bin,max_events,toffset_1,toffset_2,1); // time offset among different planes for the time electrons travel among different planes
  simu_fds.jump(eve_num);
  //simu_fds.Save();

  cout << "Deconvolution with Gaussian filter" << endl;
  WireCell2dToy::ToySignalGausFDS gaus_fds(simu_fds,gds,total_time_bin/nrebin,max_events,toffset_1,toffset_2); // gaussian smearing for charge estimation
  gaus_fds.jump(eve_num);
  //gaus_fds.Save();

  cout << "Deconvolution with Wiener filter" << endl;
   WireCell2dToy::ToySignalWienFDS wien_fds(simu_fds,gds,total_time_bin/nrebin,max_events,toffset_1,toffset_2); // weiner smearing for hit identification
  wien_fds.jump(eve_num);
  //wien_fds.Save();
  
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  int nwire_u = wires_u.size();
  int nwire_v = wires_v.size();
  int nwire_w = wires_w.size();
  
 
  
  // float threshold_u = 1000;
  // float threshold_v = 1000;
  // float threshold_w = 1000;
  

  WireCellSst::ToyuBooNESliceDataSource sds(wien_fds,gaus_fds,threshold_u, 
  					    threshold_v, threshold_w, 
  					    threshold_ug, 
  					    threshold_vg, threshold_wg, 
  					    nwire_u, 
  					    nwire_v, nwire_w); 

  WireCellSst::ToyuBooNESliceDataSource sds_th(st_fds,st_fds,500, 
  					    500, 500, 
  					    threshold_ug, 
  					    threshold_vg, threshold_wg, 
  					    nwire_u, 
  					    nwire_v, nwire_w); 

 
  int ncount = 0;
  int ncount1 = 0;  
  int ncount2 = 0;

  int ncount_t = 0;
  

  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
  WireCell2dToy::MergeToyTiling **mergetiling = new WireCell2dToy::MergeToyTiling*[2400];
  WireCell2dToy::TruthToyTiling **truthtiling = new WireCell2dToy::TruthToyTiling*[2400];
  WireCell2dToy::SimpleBlobToyTiling **blobtiling = new WireCell2dToy::SimpleBlobToyTiling*[2400];

  WireCell2dToy::ToyMatrix **toymatrix = new WireCell2dToy::ToyMatrix*[2400];
    
  //save truth ...
  WireCell2dToy::ToyTiling **toytiling_th = new WireCell2dToy::ToyTiling*[2400];
  WireCell2dToy::TruthToyTiling **truthtiling_th = new WireCell2dToy::TruthToyTiling*[2400];
 
  WireCell2dToy::ToyMetric toymetric;
  WireCell2dToy::BlobMetric blobmetric;

  //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  
  int ncount_mcell = 0;

  delete fds;
  //tfile->Close();


  int start_num = 0 ;
  int end_num = sds.size()-1;

  cout << "Start the Reconstruction " << endl; 

  for (int i=start_num;i!=end_num+1;i++){
    
    sds.jump(i);
    sds_th.jump(i);
    WireCell::Slice slice = sds.get();
    WireCell::Slice slice_th = sds_th.get();
    cout << i << " " << slice.group().size() << " " << slice_th.group().size() << endl;
    
    toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
    GeomCellSelection allcell = toytiling[i]->get_allcell();
    GeomWireSelection allwire = toytiling[i]->get_allwire();
    cout << "Single Cell: " << i << " "  << allcell.size() << " " << allwire.size() << endl;
    

    mergetiling[i] = new WireCell2dToy::MergeToyTiling(*toytiling[i],i,3,1);
    
    
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    GeomWireSelection allmwire = mergetiling[i]->get_allwire();
    
    cout <<"Blob: " << i << " " << allmcell.size() << " " << allmwire.size() << endl;
    
    
    truthtiling[i] = new WireCell2dToy::TruthToyTiling(*toytiling[i],pvv,timeoffset,i,gds,unit_dis);
    //cout << "finish truth tiling " << endl; 
    toymatrix[i] = new WireCell2dToy::ToyMatrix(*toytiling[i],*mergetiling[i]);
    // cout << "start the iterate " << endl; 
    if (toymatrix[i]->Get_Solve_Flag()==0){
      WireCell2dToy::ToyMatrixIterate toymatrix_it(*toymatrix[i]);
    }

    cout << "chi2: " << toymatrix[i]->Get_Chi2() << endl;
    cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
    
    toytiling_th[i] = new WireCell2dToy::ToyTiling(slice_th,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
    truthtiling_th[i] = new WireCell2dToy::TruthToyTiling(*toytiling_th[i],pvv,timeoffset,i,gds,unit_dis);
    
    
    CellChargeMap ccmap = truthtiling[i]->ccmap();
    if (toymatrix[i]->Get_Solve_Flag()!=0)
      toymetric.Add(allmcell,*toymatrix[i],ccmap);
    
    toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
    
    Double_t charge_min = 10000;
    Double_t charge_max = 0;
    
  }

  toymetric.Print();
  std::cout << "Starting MCMC" << std::endl;

  // //without  time information
  // for (int i=start_num;i!=end_num+1;i++){
  //   if (toymatrix[i]->Get_Solve_Flag()==0){
  //     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
  //     WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i],&allmcell);
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
    int first_solve=-1;
    for (int i=start_num; i!=end_num+1;i++){
      if (toymatrix[i]->Get_Solve_Flag()!=0){
  	first_solve = i;
  	break;
      }
    }
    if (first_solve <0){
      for (int i=start_num;i!=end_num+1;i++){
  	if (toymatrix[i]->Get_Solve_Flag()==0){
  	  GeomCellSelection allmcell = mergetiling[i]->get_allcell();
  	  WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i],&allmcell);
  	  CellChargeMap ccmap = truthtiling[i]->ccmap();
  	  if (toymatrix[i]->Get_Solve_Flag()!=0)
  	    toymetric.Add(allmcell,*toymatrix[i],ccmap);
  	  toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
  	  cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
  	  cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
  	}
      }
    }else{
      for (int i=first_solve+1;i<=end_num-1;i++){
  	if (toymatrix[i]->Get_Solve_Flag()==0){
  	  GeomCellSelection allmcell = mergetiling[i]->get_allcell();
  	  WireCell2dToy::ToyMatrixMarkov toymatrix_markov(*toymatrix[i-1],*toymatrix[i],*toymatrix[i+1],*mergetiling[i-1],*mergetiling[i],*mergetiling[i+1],&allmcell);
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
  	WireCell2dToy::ToyMatrixMarkov toymatrix_markov(*toymatrix[end_num-1],*toymatrix[end_num],*toymatrix[end_num-1],*mergetiling[end_num-1],*mergetiling[end_num],*mergetiling[end_num-1],&allmcell);
	
	
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
  	  WireCell2dToy::ToyMatrixMarkov toymatrix_markov(*toymatrix[i-1],*toymatrix[i],*toymatrix[i+1],*mergetiling[i-1],*mergetiling[i],*mergetiling[i+1],&allmcell);
	  
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
  	WireCell2dToy::ToyMatrixMarkov toymatrix_markov(*toymatrix[start_num+1],*toymatrix[start_num],*toymatrix[start_num+1],*mergetiling[start_num+1],*mergetiling[start_num],*mergetiling[start_num+1],&allmcell);
	
	
  	CellChargeMap ccmap = truthtiling[start_num]->ccmap();
  	if (toymatrix[start_num]->Get_Solve_Flag()!=0)
  	  toymetric.Add(allmcell,*toymatrix[start_num],ccmap);
  	toymetric.AddSolve(toymatrix[start_num]->Get_Solve_Flag());
	
  	cout << "chi2: " << start_num << " " << toymatrix[start_num]->Get_Chi2() << endl;
  	cout << "NDF: " << toymatrix[start_num]->Get_ndf() << endl;
      }
    }
  }


  


  // //do blob thing ... 
  // //use time information
  // std::cout << "Reduce Blob" << std::endl; 
  // for (int i=start_num;i!=end_num+1;i++){
  //   std::cout << "Check Blob " << i << std::endl;
  //   //std::cout << toymatrix[i]->GetSimpleBlobReduction() << std::endl;
  //   if (!mergetiling[i]->GetRemerged()){
  //     toymatrix[i]->JudgeSimpleBlob(*toytiling[i],*mergetiling[i]);
  //   }

  //   //std::cout << toymatrix[i]->GetSimpleBlobReduction() << std::endl;
  //   if (toymatrix[i]->GetSimpleBlobReduction()){
  //     if (i==start_num){
  // 	blobtiling[i] = new WireCell2dToy::SimpleBlobToyTiling(*toytiling[i],*mergetiling[i],*toymatrix[i],*mergetiling[i+1],*toymatrix[i+1],*mergetiling[i+1],*toymatrix[i+1]);
  //     }else if (i==end_num){
  // 	blobtiling[i] = new WireCell2dToy::SimpleBlobToyTiling(*toytiling[i],*mergetiling[i],*toymatrix[i],*mergetiling[i-1],*toymatrix[i-1],*mergetiling[i-1],*toymatrix[i-1]);
  //     }else{
  // 	blobtiling[i] = new WireCell2dToy::SimpleBlobToyTiling(*toytiling[i],*mergetiling[i],*toymatrix[i],*mergetiling[i-1],*toymatrix[i-1],*mergetiling[i+1],*toymatrix[i+1]);
  //     }
  //   }
  //   if (toymatrix[i]->GetSimpleBlobReduction()){
  //     //save stuff
  //     CellChargeMap ccmap = truthtiling[i]->ccmap();
  //     blobmetric.Add(*blobtiling[i],ccmap);
      
  //     std::cout << "Check Blob " << i << std::endl;
  //     WireCell2dToy::BlobMetric tempblob;
  //     tempblob.Add(*blobtiling[i],ccmap);
  //     tempblob.Print();
  //   }
  // }
  
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


   TFile *file = new TFile(Form("shower3D_signal_%d.root",eve_num),"RECREATE");
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
      x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
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
      x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
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
  	  x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
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
    // if (toymatrix[i]->GetSimpleBlobReduction()){
    //   for (int j=0;j!=blobtiling[i]->Get_Cells().size();j++){
    // 	const GeomCell *cell = blobtiling[i]->Get_Cells().at(j);
    // 	Point p = cell->center();
    // 	x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
    // 	y_save = p.y/units::cm;
    // 	z_save = p.z/units::cm;
    // 	charge_save = blobtiling[i]->Get_Cell_Charge(cell,1);
    // 	ncharge_save = 1;
	
    // 	g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
    // 	t_rec_charge_blob->Fill();
	
    // 	ncount2 ++;
    //   }
    // }else{
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
      if (charge> recon_threshold ){
	for (int k=0;k!=mcell->get_allcell().size();k++){
	  Point p = mcell->get_allcell().at(k)->center();
	  x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
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
    // }
    
    
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
  
  // const int N = 100000;
  // Double_t x[N],y[N],z[N];
  Double_t x,y,z;
  //save cluster
  int ncluster = 0;
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    ncount = 0;
    TGraph2D *g1 = new TGraph2D();
    for (int i=0; i!=(*it)->get_allcell().size();i++){
      const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
      for (int j=0; j!=mcell->get_allcell().size();j++){
  	Point p = mcell->get_allcell().at(j)->center();
  	x = mcell->GetTimeSlice()*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
  	y = p.y/units::cm;
  	z = p.z/units::cm;
  	g1->SetPoint(ncount,x,y,z);
  	ncount ++;
      }
    }
    // cout << ncount << endl;
    g1->Write(Form("cluster_%d",ncluster));
    ncluster ++;
  }


  // save all the toy tiling stuff
  WireCell2dToy::ToyTiling* tt1 = 0;
  int time_slice;
  
  TTree* ttree = new TTree("T","T");
  ttree->Branch("time_slice",&time_slice,"time_slice/I");
  ttree->Branch("toytiling",&tt1);
  ttree->SetDirectory(file);
  for (int i=start_num;i!=end_num+1;i++){
    tt1 = toytiling[i];
    time_slice = i;
    ttree->Fill();
  }
  ttree->Write();

  TTree *ttree1 = new TTree("TC","TC");
  // To save cluster, we need to save
  // 1. time slice
  // 2. single cell
  // 3. charge
  // 4. cluster number
  const GeomCell* cell_save = 0;
  int cluster_num = -1;
  int mcell_id = -1;
  
  ttree1->Branch("time_slice",&time_slice,"time_slice/I"); // done
  ttree1->Branch("cell",&cell_save);
  ttree1->Branch("ncluster",&cluster_num,"cluster_num/I"); //done
  ttree1->Branch("mcell_id",&mcell_id,"mcell_id/I");
  ttree1->Branch("charge",&charge_save,"charge/D"); 

  double xx;
  ttree1->Branch("xx",&xx,"xx/D");    //don
  ttree1->Branch("x",&x,"x/D");    //done
  ttree1->Branch("y",&y,"y/D");
  ttree1->Branch("z",&z,"z/D");
  
   // save information to reconstruct the toytiling
  int u_index, v_index, w_index;
  double u_charge, v_charge, w_charge;
  double u_charge_err, v_charge_err, w_charge_err;
  
  ttree1->Branch("u_index",&u_index,"u_index/I");
  ttree1->Branch("v_index",&v_index,"v_index/I");
  ttree1->Branch("w_index",&w_index,"w_index/I");
  
  ttree1->Branch("u_charge",&u_charge,"u_charge/D");
  ttree1->Branch("v_charge",&v_charge,"v_charge/D");
  ttree1->Branch("w_charge",&w_charge,"w_charge/D");

  ttree1->Branch("u_charge_err",&u_charge_err,"u_charge_err/D");
  ttree1->Branch("v_charge_err",&v_charge_err,"v_charge_err/D");
  ttree1->Branch("w_charge_err",&w_charge_err,"w_charge_err/D");

  //end save 



  ttree1->SetDirectory(file);
  
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    cluster_num ++;
    //loop merged cell
    for (int i=0; i!=(*it)->get_allcell().size();i++){
      const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
      mcell_id ++;
      time_slice = mcell->GetTimeSlice();
      x = time_slice *nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
      xx = x;
      //loop single cell
      for (int j=0; j!=mcell->get_allcell().size();j++){
	cell_save = mcell->get_allcell().at(j);

	//fill the information needed for toytiling
	GeomWireSelection wires = toytiling[time_slice]->wires(*cell_save);
	
	//	if (i==0 && j==0) cout << "abc: " << time_slice << " " << toytiling[time_slice]->get_allcell().size() << " " << wires.size() << endl;

	for (int k=0;k!=wires.size();k++){
	  const GeomWire *wire = wires.at(k);
	  WirePlaneType_t plane = wire->plane();
	  if (plane==0){
	    u_index = wire->index();
	    u_charge = toytiling[time_slice]->wcmap()[wire];
	    u_charge_err = toytiling[time_slice]->wcemap()[wire];
	  }else if (plane==1){
	    v_index = wires.at(k)->index();
	    v_charge = toytiling[time_slice]->wcmap()[wire];
	    v_charge_err = toytiling[time_slice]->wcemap()[wire];
	  }else if (plane==2){
	    w_index = wire->index();
	    w_charge = toytiling[time_slice]->wcmap()[wire];
	    w_charge_err = toytiling[time_slice]->wcemap()[wire];
	  }
	}
	//end fill

	Point p = mcell->get_allcell().at(j)->center();
	charge_save = toymatrix[time_slice]->Get_Cell_Charge(mcell,1)/mcell->cross_section() * cell_save->cross_section();
	y = p.y/units::cm;
  	z = p.z/units::cm;
	ttree1->Fill();
	
      }
    }
  }
  ttree1->Write();
  
 
  TTree *Trun = new TTree("Trun","Trun");
  Trun->SetDirectory(file);
  Trun->Branch("eventNo",&event_no,"eventNo/I");
  Trun->Branch("runNo",&run_no,"runNo/I");
  Trun->Branch("subRunNo",&subrun_no,"runRunNo/I");

   Trun->Branch("unit_dis",&unit_dis,"unit_dis/F");
  Trun->Branch("toffset_uv",&toffset_1,"toffset_uv/F");
  Trun->Branch("toffset_uw",&toffset_2,"toffset_uw/F");
  Trun->Branch("toffset_u",&toffset_3,"toffset_u/F");
  Trun->Branch("total_time_bin",&total_time_bin,"total_time_bin/I");
  Trun->Branch("recon_threshold",&recon_threshold,"recon_threshold/I");
  Trun->Branch("frame_length",&frame_length,"frame_length/I");
  Trun->Branch("max_events",&max_events,"max_events/I");
  Trun->Branch("eve_num",&eve_num,"eve_num/I");
  Trun->Branch("nrebin",&nrebin,"nrebin/I");
  Trun->Branch("threshold_u",&threshold_u,"threshold_u/F");
  Trun->Branch("threshold_v",&threshold_v,"threshold_v/F");
  Trun->Branch("threshold_w",&threshold_w,"threshold_w/F");
  Trun->Branch("time_offset",&time_offset,"time_offset/I");

  Trun->Fill();

  file->Write();
  file->Close();

  toymetric.Print();
  blobmetric.Print();

  return 0;
  
} // main()
