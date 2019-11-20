#include "WCPSst/GeomDataSource.h"
#include "WCPSst/ToyuBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ClusterDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/TruthToyTiling.h"
#include "WCP2dToy/SimpleBlobToyTiling.h"

#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/ToyMatrixExclusive.h"
#include "WCP2dToy/ToyMatrixKalman.h"
#include "WCP2dToy/ToyMatrixIterate.h"
#include "WCP2dToy/ToyMatrixMarkov.h"
#include "WCP2dToy/ToyMetric.h"
#include "WCP2dToy/BlobMetric.h"
#include "WCPData/TPCParams.h"
#include "WCPData/Singleton.h"

#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"

#include "WCPData/GeomCluster.h"

#include "WCPNav/FrameDataSource.h"
#include "WCPNav/SimDataSource.h"
#include "WCPNav/SliceDataSource.h"
#include "WCPSst/Util.h"
#include "WCPData/SimTruth.h"
#include "WCP2dToy/ToyDepositor.h"
#include "WCPNav/GenerativeFDS.h"
#include "WCP2dToy/ToySignalSimu.h"
#include "WCP2dToy/ToySignalSimuTrue.h"
#include "WCP2dToy/ToySignalGaus.h"
#include "WCP2dToy/ToySignalWien.h"


#include "WCPData/SpaceCell.h"
#include "WCPData/MergeSpaceCell.h"
#include "WCP2dToy/ToyCrawler.h"

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

using namespace WCP;
using namespace std;

#define MAX_TRACKS 30000


int main(int argc, char* argv[])
{
  if (argc < 4) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num" << endl;
    return 1;
  }
  
  
  WCPSst::GeomDataSource gds(argv[1]);
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

  int save_image_outline_flag = 1; // prescale flag 
  
  
  int total_time_bin=9600;
  int frame_length = 3200;
  int nrebin = 4;
  
  TPCParams& mp = Singleton<TPCParams>::Instance();
  
  double pitch_u = gds.pitch(WirePlaneType_t(0));
  double pitch_v = gds.pitch(WirePlaneType_t(1));
  double pitch_w = gds.pitch(WirePlaneType_t(2));
  double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;

  mp.set_pitch_u(pitch_u);
  mp.set_pitch_v(pitch_v);
  mp.set_pitch_w(pitch_w);
  mp.set_ts_width(time_slice_width);
  
  std::cout << "Singleton: " << mp.get_pitch_u() << " " << mp.get_pitch_v() << " " << mp.get_pitch_w() << " " << mp.get_ts_width() << std::endl;


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
  
  //save MC truth ...
  int mc_Ntrack;  // number of tracks in MC
  int mc_id[MAX_TRACKS];  // track id; size == mc_Ntrack
  int mc_pdg[MAX_TRACKS];  // track particle pdg; size == mc_Ntrack
  int mc_process[MAX_TRACKS];  // track generation process code; size == mc_Ntrack
  int mc_mother[MAX_TRACKS];  // mother id of this track; size == mc_Ntrack
  float mc_startXYZT[MAX_TRACKS][4];  // start position of this track; size == mc_Ntrack
  float mc_endXYZT[MAX_TRACKS][4];  // end position of this track; size == mc_Ntrack
  float mc_startMomentum[MAX_TRACKS][4];  // start momentum of this track; size == mc_Ntrack
  float mc_endMomentum[MAX_TRACKS][4];  // end momentum of this track; size == mc_Ntrack
  std::vector<std::vector<int> > *mc_daughters =  new std::vector<std::vector<int> >;  // daughters id of this track; vector
  TObjArray* mc_trackPosition = new TObjArray();
  
  sst->SetBranchAddress("mc_Ntrack", &mc_Ntrack);  // number of tracks in MC
  sst->SetBranchAddress("mc_id", &mc_id);  // track id; size == mc_Ntrack
  sst->SetBranchAddress("mc_pdg", &mc_pdg);  // track particle pdg; size == mc_Ntrack
  sst->SetBranchAddress("mc_process", &mc_process);  // track generation process code; size == mc_Ntrack
  sst->SetBranchAddress("mc_mother", &mc_mother);  // mother id of this track; size == mc_Ntrack
  sst->SetBranchAddress("mc_daughters", &mc_daughters);  // daughters id of this track; vector
  sst->SetBranchAddress("mc_startXYZT", &mc_startXYZT);  // start position of this track; size == mc_Ntrack
  sst->SetBranchAddress("mc_endXYZT", &mc_endXYZT);  // start position of this track; size == mc_Ntrack
  sst->SetBranchAddress("mc_startMomentum", &mc_startMomentum);  // start momentum of this track; size == mc_Ntrack
  sst->SetBranchAddress("mc_endMomentum", &mc_endMomentum);  // start momentum of this track; size == mc_Ntrack
  sst->SetBranchAddress("mc_trackPosition",&mc_trackPosition);

  sst->GetEntry(eve_num);


  cout << "Run No: " << run_no << " " << subrun_no << " " << eve_num << endl;

  WCP::FrameDataSource* fds = 0;
  fds = WCPSst::make_fds(*tfile);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  
  TH1::AddDirectory(kFALSE);
  
  
 
  
  WCP::ToyDepositor *toydep = new WCP::ToyDepositor(fds,0,unit_dis,frame_length);
  const PointValueVector& pvv = toydep->depositions(eve_num);
  
  WCP::GenerativeFDS *gfds = new WCP::GenerativeFDS(*toydep,gds,total_time_bin,max_events,0.5*unit_dis*units::millimeter); // 87 K at 0.5 kV/cm
  gfds->jump(eve_num);

  
  int abc;
  
  tfile->Close("R");
  delete tfile;

  
  cout << "Put in Truth " << pvv.size() << endl; 
  WCP2dToy::ToySignalSimuTrueFDS *st_fds = new WCP2dToy::ToySignalSimuTrueFDS(*gfds,gds,total_time_bin/nrebin,max_events,0); //truth
  st_fds->jump(eve_num);
  
  cout << "Simulate Raw WaveForm " << endl; 
  WCP2dToy::ToySignalSimuFDS *simu_fds = new WCP2dToy::ToySignalSimuFDS(*gfds,gds,total_time_bin,max_events,toffset_1,toffset_2,1); // time offset among different planes for the time electrons travel among different planes
  simu_fds->jump(eve_num);
  
  cout << "Deconvolution with Gaussian filter" << endl;
  WCP2dToy::ToySignalGausFDS *gaus_fds = new WCP2dToy::ToySignalGausFDS(*simu_fds,gds,total_time_bin/nrebin,max_events,toffset_1,toffset_2); // gaussian smearing for charge estimation
  gaus_fds->jump(eve_num);
  
  cout << "Deconvolution with Wiener filter" << endl;
   WCP2dToy::ToySignalWienFDS *wien_fds = new WCP2dToy::ToySignalWienFDS(*simu_fds,gds,total_time_bin/nrebin,max_events,toffset_1,toffset_2); // weiner smearing for hit identification
  wien_fds->jump(eve_num);
  
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  int nwire_u = wires_u.size();
  int nwire_v = wires_v.size();
  int nwire_w = wires_w.size();
  
 
  
  // cin >> abc;

  WCPSst::ToyuBooNESliceDataSource *sds = new WCPSst::ToyuBooNESliceDataSource(*wien_fds,*gaus_fds,threshold_u, 
  					    threshold_v, threshold_w, 
  					    threshold_ug, 
  					    threshold_vg, threshold_wg, 
  					    nwire_u, 
  					    nwire_v, nwire_w); 

  WCPSst::ToyuBooNESliceDataSource *sds_th = new WCPSst::ToyuBooNESliceDataSource(*st_fds,*st_fds,500, 
  					    500, 500, 
  					    threshold_ug, 
  					    threshold_vg, threshold_wg, 
  					    nwire_u, 
  					    nwire_v, nwire_w); 

 
  int ncount = 0;
  int ncount1 = 0;  
  int ncount2 = 0;

  int ncount_t = 0;
  
  
  

  WCP2dToy::ToyTiling **toytiling = new WCP2dToy::ToyTiling*[2400];
  WCP2dToy::MergeToyTiling **mergetiling = new WCP2dToy::MergeToyTiling*[2400];
  WCP2dToy::TruthToyTiling **truthtiling = new WCP2dToy::TruthToyTiling*[2400];
  WCP2dToy::SimpleBlobToyTiling **blobtiling = new WCP2dToy::SimpleBlobToyTiling*[2400];
  WCP2dToy::ToyMatrix **toymatrix = new WCP2dToy::ToyMatrix*[2400];
  //save truth ...
  WCP2dToy::ToyTiling **toytiling_th = new WCP2dToy::ToyTiling*[2400];
  WCP2dToy::TruthToyTiling **truthtiling_th = new WCP2dToy::TruthToyTiling*[2400];
  WCP2dToy::ToyMetric toymetric;
  WCP2dToy::BlobMetric blobmetric;


   //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  
  int ncount_mcell = 0;
  
  int start_num = 0 ;
  int end_num = sds->size()-1;
  // int start_num = 800 ;
  // int end_num = 810;
  // GeomCellSelection total_cells;
  // GeomCellSelection total_edge_cells;

  cout << "Start the Tiling " << endl; 
  WCP::Slice slice;
  for (int i=start_num;i!=end_num+1;i++){
    
    sds->jump(i);
    sds_th->jump(i);
    slice = sds->get();
    WCP::Slice slice_th = sds_th->get();
    // cout << i << " " << slice.group().size() << " " << slice_th.group().size() << endl;
    
    toytiling[i] = new WCP2dToy::ToyTiling(slice,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
    GeomCellSelection allcell = toytiling[i]->get_allcell();
    
    GeomWireSelection allwire = toytiling[i]->get_allwire();
    cout << "Single Cell: " << i << " "  << allcell.size() << " " << allwire.size() << endl;
    //    mergetiling[i] = new WCP2dToy::MergeToyTiling(*toytiling[i],i,3,1);
    mergetiling[i] = new WCP2dToy::MergeToyTiling(*toytiling[i],i,3);
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    
    GeomWireSelection allmwire = mergetiling[i]->get_allwire();
    cout <<"Blob: " << i << " " << allmcell.size() << " " << allmwire.size() << endl;
    truthtiling[i] = new WCP2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds,frame_length/nrebin,unit_dis);
    
    
    toytiling_th[i] = new WCP2dToy::ToyTiling(slice_th,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
    truthtiling_th[i] = new WCP2dToy::TruthToyTiling(*toytiling_th[i],pvv,i,gds,frame_length/nrebin,unit_dis);
    //CellChargeMap ccmap = truthtiling[i]->ccmap();
        
    // for (int j=0;j!=allcell.size();j++){
    //   total_cells.push_back(allcell.at(j));
    // }
    // for (int j=0;j!=allmcell.size();j++){
    //   const MergeGeomCell *mcell = (const MergeGeomCell*)allmcell.at(j);
    //   GeomCellSelection edgecells = mcell->get_edgecells();
    //   for (int k=0;k!=edgecells.size();k++){
    // 	total_edge_cells.push_back(edgecells.at(k));
    //   }
    // }

  }
  delete sds;
  delete sds_th;
  
  delete simu_fds;
  delete gaus_fds;
  delete wien_fds;

  delete st_fds;
  delete gfds;
  delete toydep;
  delete fds;


  
    //do display
  
  // TApplication theApp("theApp",&argc,argv);
  // theApp.SetReturnFromRun(true);
  
  // TCanvas c1("ToyMC","ToyMC",800,600);
  // WCP2dToy::ToyEventDisplay display(c1, gds);
  
  // display.charge_min = 0.;
  // display.charge_max = 1.;
  
  
  
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
    
  // display.draw_mc(1,WCP::PointValueVector(),"colz");
  // display.draw_slice(slice,"");
  // display.draw_cells(total_cells,"*same");
  // display.draw_cells(total_edge_cells,"*same",2);
  // theApp.Run();



  


   for (int i=start_num;i!=end_num+1;i++){
     toymatrix[i] = new WCP2dToy::ToyMatrix(*toytiling[i],*mergetiling[i]);
    // cout << "start the iterate " << endl; 
    // if (toymatrix[i]->Get_Solve_Flag()==0){
    //   WCP2dToy::ToyMatrixIterate toymatrix_it(*toymatrix[i],2000,1e5);
    // }
    
    cout << i << " chi2: " << toymatrix[i]->Get_Chi2() <<
      " NDF: " << toymatrix[i]->Get_ndf() << endl;
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    CellChargeMap ccmap = truthtiling[i]->ccmap();
    if (toymatrix[i]->Get_Solve_Flag()!=0)
      toymetric.Add(allmcell,*toymatrix[i],ccmap);
    
    toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
   }


   toymetric.Print();

   double penalty = 5;
   std::cout << "Starting to use connectivitiy" << std::endl;
   std::list<int> solve_list;
   
   if (start_num != end_num){
     int first_solve=-1;
     for (int i=start_num; i!=end_num+1;i++){
       if (toymatrix[i]->Get_Solve_Flag()!=0){
   	 first_solve = i;
   	 break;
       }
     }
     if (first_solve == -1) first_solve = start_num;

    
     for (int i=first_solve+1;i<=end_num-1;i++){
       if (toymatrix[i]->Get_Solve_Flag()==0 && toymatrix[i]->Get_mcindex()>0){ 
  	 if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
  	   WCP2dToy::ToyMatrixIterate toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty);
	   
  	   GeomCellSelection allmcell = mergetiling[i]->get_allcell();
  	   CellChargeMap ccmap = truthtiling[i]->ccmap();
  	   if (toymatrix[i]->Get_Solve_Flag()!=0)
  	     toymetric.Add(allmcell,*toymatrix[i],ccmap);
  	   toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	   
  	   cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
  	   cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
  	 }else{
  	   solve_list.push_back(i); 
  	 }
       }
     }
     
     for (int i=first_solve-1;i>=start_num+1;i--){
       if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
  	 if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
  	   WCP2dToy::ToyMatrixIterate toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty);
	 
  	   GeomCellSelection allmcell = mergetiling[i]->get_allcell();
  	   CellChargeMap ccmap = truthtiling[i]->ccmap();
  	   if (toymatrix[i]->Get_Solve_Flag()!=0)
  	     toymetric.Add(allmcell,*toymatrix[i],ccmap);
  	   toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	 
  	   cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
  	   cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
  	 }else{
  	   solve_list.push_back(i);
  	 }
       }
     }
   }
   
   // start second round ...
   // std::cout << "Connectivitiy rest " << solve_list.size() << std::endl;
   int prev_count = 0;
   while (solve_list.size() >0){
     int curr_count = solve_list.size();
     std::cout << "Connectivitiy rest " << solve_list.size() << std::endl;

     if (curr_count == prev_count){
       int i = solve_list.front(); // pick the first element ... 
       if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
  	 WCP2dToy::ToyMatrixIterate toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty);
  	 GeomCellSelection allmcell = mergetiling[i]->get_allcell();
  	 CellChargeMap ccmap = truthtiling[i]->ccmap();
  	 if (toymatrix[i]->Get_Solve_Flag()!=0)
  	   toymetric.Add(allmcell,*toymatrix[i],ccmap);
  	 toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
  	 cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
  	 cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
  	 solve_list.erase(solve_list.begin());
       }
     }else{
       for (auto it = solve_list.begin(); it!= solve_list.end(); it++){
  	 int i = *it;
  	 if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
  	   if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
  	     WCP2dToy::ToyMatrixIterate toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty);
  	     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
  	     CellChargeMap ccmap = truthtiling[i]->ccmap();
  	     if (toymatrix[i]->Get_Solve_Flag()!=0)
  	       toymetric.Add(allmcell,*toymatrix[i],ccmap);
  	     toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
  	     cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
  	     cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
  	     it = solve_list.erase(it);
  	   }
  	 }
       }
     }
     
     prev_count = curr_count;
   }


   
   // by the end do the final two
   if (toymatrix[end_num]->Get_Solve_Flag()==0&& toymatrix[end_num]->Get_mcindex()>0){
     WCP2dToy::ToyMatrixIterate toymatrix_it(toymatrix[end_num-1],toymatrix[end_num],0,mergetiling[end_num-1],mergetiling[end_num],0,2000,1e5,penalty);
     GeomCellSelection allmcell = mergetiling[end_num]->get_allcell();
     CellChargeMap ccmap = truthtiling[end_num]->ccmap();
     if (toymatrix[end_num]->Get_Solve_Flag()!=0)
       toymetric.Add(allmcell,*toymatrix[end_num],ccmap);
     toymetric.AddSolve(toymatrix[end_num]->Get_Solve_Flag());
     
     cout << "chi2: " << end_num << " " << toymatrix[end_num]->Get_Chi2() << endl;
     cout << "NDF: " << toymatrix[end_num]->Get_ndf() << endl;
   }
   
   if (toymatrix[start_num]->Get_Solve_Flag()==0){
     WCP2dToy::ToyMatrixIterate toymatrix_it(0,toymatrix[start_num],toymatrix[start_num+1],0,mergetiling[start_num],mergetiling[start_num+1],2000,1e5,penalty);
     
     GeomCellSelection allmcell = mergetiling[start_num]->get_allcell();
     CellChargeMap ccmap = truthtiling[start_num]->ccmap();
     if (toymatrix[start_num]->Get_Solve_Flag()!=0&& toymatrix[start_num]->Get_mcindex()>0)
       toymetric.Add(allmcell,*toymatrix[start_num],ccmap);
     toymetric.AddSolve(toymatrix[start_num]->Get_Solve_Flag());
     
     cout << "chi2: " << start_num << " " << toymatrix[start_num]->Get_Chi2() << endl;
     cout << "NDF: " << toymatrix[start_num]->Get_ndf() << endl;
   }
   
   toymetric.Print();
   std::cout << "Starting MCMC" << std::endl;
   solve_list.clear();
   
   //without  time information
   // for (int i=start_num;i!=end_num+1;i++){
   //   if (toymatrix[i]->Get_Solve_Flag()==0){
   //     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   //     WCP2dToy::ToyMatrixMarkov toymatrix_markov(*toymatrix[i],*mergetiling[i],&allmcell,Good_MCells.at(i-start_num));
   //     //WCP2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i],&allmcell);
   //     CellChargeMap ccmap = truthtiling[i]->ccmap();
   //     if (toymatrix[i]->Get_Solve_Flag()!=0)
   // 	toymetric.Add(allmcell,*toymatrix[i],ccmap);
   //     toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
   //     cout << " chi2: " << i << " " << toymatrix[i]->Get_Chi2() 
   // 	   << " NDF: " << toymatrix[i]->Get_ndf() << endl;
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
     if (first_solve == -1) first_solve = start_num;
     
     
     for (int i=first_solve+1;i<=end_num-1;i++){
       if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
   	 if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
   	   GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   	   WCP2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell,1500,2000,penalty);
   	   CellChargeMap ccmap = truthtiling[i]->ccmap();
   	   if (toymatrix[i]->Get_Solve_Flag()!=0)
   	     toymetric.Add(allmcell,*toymatrix[i],ccmap);
   	   toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	   
   	   cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
   	   cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
   	 }else{
   	   solve_list.push_back(i);
   	 }
   	 //toymetric.Print();
       }
     }
     
     // go to early ones 
     for (int i=first_solve-1;i>=start_num+1;i--){
       if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
   	 if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
   	   GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   	   WCP2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell,1500,2000,penalty);
	   
   	   CellChargeMap ccmap = truthtiling[i]->ccmap();
   	   if (toymatrix[i]->Get_Solve_Flag()!=0)
   	     toymetric.Add(allmcell,*toymatrix[i],ccmap);
   	   toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	   
   	   cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
   	   cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
   	 }else{
   	   solve_list.push_back(i);
   	 }
       }
     }
      
     // do the while ... 
     int prev_count = 0;
   while (solve_list.size() >0){
     int curr_count = solve_list.size();
     std::cout << "Connectivitiy rest " << solve_list.size() << std::endl;

     if (curr_count == prev_count){
       int i = solve_list.front(); // pick the first element ... 
       if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
   	 GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   	 WCP2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell,1500,2000,penalty);
	 
   	 CellChargeMap ccmap = truthtiling[i]->ccmap();
   	 if (toymatrix[i]->Get_Solve_Flag()!=0)
   	   toymetric.Add(allmcell,*toymatrix[i],ccmap);
   	 toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
   	 cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
   	 cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
   	 solve_list.erase(solve_list.begin());
       }
     }else{
       for (auto it = solve_list.begin(); it!= solve_list.end(); it++){
   	 int i = *it;
   	 if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
   	   if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
   	     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   	     WCP2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell,1500,2000,penalty);
	     
   	     CellChargeMap ccmap = truthtiling[i]->ccmap();
   	     if (toymatrix[i]->Get_Solve_Flag()!=0)
   	       toymetric.Add(allmcell,*toymatrix[i],ccmap);
   	     toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
   	     cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
   	     cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
   	     it = solve_list.erase(it);
   	   }
   	 }
       }
     }
     
     prev_count = curr_count;
   }


     //deal with the start/end ones ... 
     if (toymatrix[end_num]->Get_Solve_Flag()==0&& toymatrix[end_num]->Get_mcindex()>0){
       GeomCellSelection allmcell = mergetiling[end_num]->get_allcell();
       WCP2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[end_num-1],toymatrix[end_num],0,mergetiling[end_num-1],mergetiling[end_num],0,&allmcell,1500,2000,penalty);
       
       
       CellChargeMap ccmap = truthtiling[end_num]->ccmap();
       if (toymatrix[end_num]->Get_Solve_Flag()!=0)
   	 toymetric.Add(allmcell,*toymatrix[end_num],ccmap);
       toymetric.AddSolve(toymatrix[end_num]->Get_Solve_Flag());
       
       cout << "chi2: " << end_num << " " << toymatrix[end_num]->Get_Chi2() << endl;
       cout << "NDF: " << toymatrix[end_num]->Get_ndf() << endl;
     }
     
     if (toymatrix[start_num]->Get_Solve_Flag()==0&& toymatrix[start_num]->Get_mcindex()>0){
       GeomCellSelection allmcell = mergetiling[start_num]->get_allcell();
       WCP2dToy::ToyMatrixMarkov toymatrix_markov(0,toymatrix[start_num],toymatrix[start_num+1],0,mergetiling[start_num],mergetiling[start_num+1],&allmcell,1500,2000,penalty);
       
       
       CellChargeMap ccmap = truthtiling[start_num]->ccmap();
       if (toymatrix[start_num]->Get_Solve_Flag()!=0)
   	 toymetric.Add(allmcell,*toymatrix[start_num],ccmap);
       toymetric.AddSolve(toymatrix[start_num]->Get_Solve_Flag());
       
       cout << "chi2: " << start_num << " " << toymatrix[start_num]->Get_Chi2() << endl;
       cout << "NDF: " << toymatrix[start_num]->Get_ndf() << endl;
     }
   }

   // form a map to illustrate connectivities 
   GeomCellCellsMap cell_prev_map;
   GeomCellCellsMap cell_next_map;
   
   for (int i=start_num;i!=end_num;i++){
     GeomCellSelection curr_mcell = mergetiling[i]->get_allcell();
     GeomCellSelection next_mcell = mergetiling[i+1]->get_allcell();
     for (int j=0;j!=curr_mcell.size();j++){
       const MergeGeomCell *curr_cell = (MergeGeomCell*)curr_mcell.at(j);
       for (int k=0;k!=next_mcell.size();k++){
	 const MergeGeomCell *next_cell = (MergeGeomCell*)next_mcell.at(k);
	 if (curr_cell->Overlap(*next_cell)){

	   if (cell_next_map.find(curr_cell)==cell_next_map.end()){
	     GeomCellSelection cells;
	     cell_next_map[curr_cell] = cells;
	   }
	   cell_next_map[curr_cell].push_back(next_cell);
	   
	   if (cell_prev_map.find(next_cell)==cell_prev_map.end()){
	     GeomCellSelection cells;
	     cell_prev_map[next_cell] = cells;
	   }
	   cell_prev_map[next_cell].push_back(curr_cell);
	   
	 }
       }
     }
   }
   // 
   
    // save good cluster_cells;
   GeomCellSelection good_cluster_cells;
   for (int i=start_num;i!=end_num+1;i++){
     GeomCellSelection curr_mcell = mergetiling[i]->get_allcell();
     for (int j=0;j!=curr_mcell.size();j++){
       const MergeGeomCell *curr_cell = (MergeGeomCell*)curr_mcell.at(j);
       int flag_good_cluster_cell = 0;
     
       if (cell_next_map.find(curr_cell) != cell_next_map.end()){
	 for (int k=0;k!=cell_next_map[curr_cell].size();k++){
	   if (toymatrix[i+1]->Get_Cell_Charge(cell_next_map[curr_cell].at(k))> recon_threshold || toymatrix[i+1]->Get_Solve_Flag()==0){
	     flag_good_cluster_cell = 1;
	     break;
	   }
	 }
       }
       if (flag_good_cluster_cell == 0){
	 if (cell_prev_map.find(curr_cell) != cell_prev_map.end()){
	   for (int k=0;k!=cell_prev_map[curr_cell].size();k++){
	     if (toymatrix[i-1]->Get_Cell_Charge(cell_prev_map[curr_cell].at(k))> recon_threshold || toymatrix[i-1]->Get_Solve_Flag()==0){
	       flag_good_cluster_cell = 1;
	       break;
	     }
	   }
	 }
       }

       if (flag_good_cluster_cell == 1)
	 good_cluster_cells.push_back(curr_cell);
     }
   }

 

  // reset existing cluster and redo it ...
  std::cout << "Start Clustering after solving equations " << std::endl;
  
  
  ncount_mcell = 0;

     
  //do clustering ... 
   for (int i=start_num;i!=end_num+1;i++){
     GeomCellSelection pallmcell = mergetiling[i]->get_allcell();
     GeomCellSelection allmcell;
     for (int j=0;j!=pallmcell.size();j++){
       const GeomCell* mcell = pallmcell[j];
       
       // if (toymatrix[i]->Get_Cell_Charge(mcell)> recon_threshold || toymatrix[i]->Get_Solve_Flag()==0){
       // 	 allmcell.push_back(mcell);
       // }
       
        int flag_save_cell = 0;

      if (toymatrix[i]->Get_Solve_Flag()==0){
      	flag_save_cell = 1;
      }else{
      	if (toymatrix[i]->Get_Cell_Charge(mcell)> recon_threshold){
      	  flag_save_cell = 1;
      	}else{
      	  if (i == start_num || i == end_num + 1) continue;
	  
      	  // there are good cells from the prev and next
      	  flag_save_cell = 0;
      	  //	  std::cout << "Xin: " << cell_next_map[mcell].size() <<  " " << cell_prev_map[mcell].size() << std::endl;

      	  for (int k=0;k!=cell_next_map[mcell].size();k++){
	    if (find(good_cluster_cells.begin(),good_cluster_cells.end(),cell_next_map[mcell].at(k))!=good_cluster_cells.end()){
	      //if (toymatrix[i+1]->Get_Cell_Charge(cell_next_map[mcell].at(k))> recon_threshold){
      	      flag_save_cell = 1;
      	      break;
      	    }
      	  }
      	  if (flag_save_cell==1){
      	    flag_save_cell = 0;
      	    for (int k=0;k!=cell_prev_map[mcell].size();k++){
	      if (find(good_cluster_cells.begin(),good_cluster_cells.end(),cell_prev_map[mcell].at(k))!=good_cluster_cells.end()){
		//if (toymatrix[i-1]->Get_Cell_Charge(cell_prev_map[mcell].at(k))> recon_threshold){
      		flag_save_cell = 1;
      		break;
      	      }
      	    }
      	  }
      	}
      }
      
      
      if (flag_save_cell == 1)
      	allmcell.push_back(mcell);
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


  
  // start crawler
   cout << "Start Crawling after solving equations " << endl;
   std::vector<WCP2dToy::ToyCrawler*> crawlers;
   MergeSpaceCellSelection all_msc_cells;
   SpaceCellSelection all_sc_cells;
   int ncluster = 0;
   for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
     
     MergeSpaceCellSelection mscells;
     for (int i=0; i!=(*it)->get_allcell().size();i++){
       const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
       MergeSpaceCell *mscell = new MergeSpaceCell();
       all_msc_cells.push_back(mscell);
       mscell->set_mcell(mcell);
       for (int j=0;j!=mcell->get_allcell().size();j++){
	 const GeomCell *cell = mcell->get_allcell().at(j);
	 SpaceCell *space_cell = new SpaceCell(ncluster,*cell,(mcell->GetTimeSlice()*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.)*units::cm,1,unit_dis/10.*nrebin/2.*units::cm);
	 all_sc_cells.push_back(space_cell);
	 mscell->AddSpaceCell(space_cell);
       }
       mscells.push_back(mscell);
     }
     WCP2dToy::ToyCrawler* toycrawler = new WCP2dToy::ToyCrawler(mscells,2);
     crawlers.push_back(toycrawler);
     
     // std::cout << ncluster << " " << toycrawler->Get_mcells_map().size() << " " << toycrawler->Get_allCT().size() << " " << toycrawler->Get_allMCT().size()  << std::endl;
     
     // if (toycrawler->Get_mcells_map().size()>200){
     //   TApplication theApp("theApp",&argc,argv);
     //   theApp.SetReturnFromRun(true);
     
     //   TCanvas c1("ToyMC","ToyMC",800,600);
     //   c1.Draw();
     
     //   WCP2dToy::ClusterDisplay display(c1);
     //   //display.DrawCluster(cells);
     //   display.DrawCluster(mscells);
     
     
     //   display.DrawCrawler(*toycrawler,"psame",1);
     
     //   theApp.Run();
     // }
     
     
     ncluster ++;
   }
  
  
  
  //check # of clusters 
  // sum = 0 ;
  // for (int i=0;i!=crawlers.size();i++){
  //   for (auto it = crawlers.at(i)->Get_mcells_map().begin(); it!= crawlers.at(i)->Get_mcells_map().end();it++){
  //     MergeSpaceCell *mcell1 = it->first;
  //     sum += mcell1->Get_all_spacecell().size();
  //   }
  // }
  
  // std::cout << "Check: " << crawlers.size() << " "  << sum << std::endl;


  //save files
  TFile *file = new TFile(Form("shower3D_cluster_%d.root",eve_num),"RECREATE");
  TTree *t_true = new TTree("T_true","T_true");
  TTree *t_rec = new TTree("T_rec","T_rec");
  TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
  TTree *t_rec_charge_blob = new TTree("T_rec_charge_blob","T_rec_charge_blob");

  Double_t x_save, y_save, z_save;
  Double_t charge_save;
  Double_t ncharge_save;
  Double_t chi2_save;
  Double_t ndf_save;
  Double_t type_save;
  
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

  if (save_image_outline_flag==1){
    t_rec->Branch("type",&type_save,"type/D");
    t_rec_charge->Branch("type",&type_save,"type/D");
    t_rec_charge_blob->Branch("type",&type_save,"type/D");
  }
  
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
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];

      if (save_image_outline_flag==1){
	Point p = mcell->get_allcell().at(0)->center();
	x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
	y_save = p.y/units::cm;
	z_save = p.z/units::cm;
	g->SetPoint(ncount,x_save,y_save,z_save);
	ncount ++;
	type_save = 1; //center
	t_rec->Fill();
	
	for (int k=0;k!=mcell->get_edgecells().size();k++){
	  Point p = mcell->get_edgecells().at(k)->center();
	  y_save = p.y/units::cm;
    	  z_save = p.z/units::cm;
	  g->SetPoint(ncount,x_save,y_save,z_save);
	  ncount ++;
	  type_save = 2; //boundary ...
	  t_rec->Fill();
	}
	
	for (int k=0;k!=mcell->get_allcell().size();k++){
    	  Point p = mcell->get_allcell().at(k)->center();
	  x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
	  y_save = p.y/units::cm;
    	  z_save = p.z/units::cm;
	  type_save = 0; //original
	  t_rec->Fill();
	}


      }else{
	for (int k=0;k!=mcell->get_allcell().size();k++){
    	  Point p = mcell->get_allcell().at(k)->center();
	  x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
	  y_save = p.y/units::cm;
    	  z_save = p.z/units::cm;
	  g->SetPoint(ncount,x_save,y_save,z_save);
	  ncount ++;
	  t_rec->Fill();
	 
	}
      }
    }

    // GeomCellSelection allcell = toytiling[i]->get_allcell();
    // for (int j=0;j!=allcell.size();j++){
    //   Point p = allcell[j]->center();
    //   x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
    //   y_save = p.y/units::cm;
    //   z_save = p.z/units::cm;
      

    //   g->SetPoint(ncount,x_save,y_save,z_save);
    //   t_rec->Fill();

    //   ncount ++;
    // }
    
    //recon 2 with charge
    //GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
      if (charge> recon_threshold || toymatrix[i]->Get_Solve_Flag()==0){

    	if (toymatrix[i]->Get_Solve_Flag()==0)
	  charge = toytiling[i]->get_ave_charge();

	if (save_image_outline_flag==1){
	  Point p = mcell->get_allcell().at(0)->center();
	  x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
	  y_save = p.y/units::cm;
	  z_save = p.z/units::cm;
	  charge_save = charge;
	  ncharge_save = mcell->get_allcell().size();
	  chi2_save = toymatrix[i]->Get_Chi2();
	  ndf_save = toymatrix[i]->Get_ndf();
	  g_rec->SetPoint(ncount1,x_save,y_save,z_save);
	  ncount1 ++;
	  type_save = 1; //center
	  t_rec_charge->Fill();

	  charge_save = 0;
	  ncharge_save = 0;
	  chi2_save= 0;
	  ndf_save = 0;

	  for (int k=0;k!=mcell->get_edgecells().size();k++){
	    Point p = mcell->get_edgecells().at(k)->center();
	    y_save = p.y/units::cm;
	    z_save = p.z/units::cm;
	    g_rec->SetPoint(ncount1,x_save,y_save,z_save);
	    ncount1 ++;
	    type_save = 2; //boundary ...
	    t_rec_charge->Fill();
	  }

	  for (int k=0;k!=mcell->get_allcell().size();k++){
	    Point p = mcell->get_allcell().at(k)->center();
	    x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
	    y_save = p.y/units::cm;
	    z_save = p.z/units::cm;
	    charge_save = charge/mcell->get_allcell().size();
	    ncharge_save = mcell->get_allcell().size();
	    chi2_save = toymatrix[i]->Get_Chi2();
	    ndf_save = toymatrix[i]->Get_ndf();
	    type_save = 0; // original
	    t_rec_charge->Fill();
	  }
	}else{
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
    }
    
    // for (int j=0;j!=allmcell.size();j++){
    //   MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
    //   double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
    //   if (charge> recon_threshold || toymatrix[i]->Get_Solve_Flag()==0){
    // 	if (toymatrix[i]->Get_Solve_Flag()==0)
    // 	  charge = toytiling[i]->get_ave_charge();

    // 	for (int k=0;k!=mcell->get_allcell().size();k++){
    // 	  Point p = mcell->get_allcell().at(k)->center();
    // 	  x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
    // 	  y_save = p.y/units::cm;
    // 	  z_save = p.z/units::cm;
    // 	  charge_save = charge/mcell->get_allcell().size();
    // 	  ncharge_save = mcell->get_allcell().size();
	  
    // 	  g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
    // 	  t_rec_charge_blob->Fill();
	  
    // 	  ncount2 ++;
    // 	}
    //   }
    // }


  }
 
  

  
  
  // // const int N = 100000;
  // // Double_t x[N],y[N],z[N];
  // Double_t x,y,z;
  // //save cluster
  // ncluster = 0;
  // for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  //   ncount = 0;
  //   TGraph2D *g1 = new TGraph2D();
  //   for (int i=0; i!=(*it)->get_allcell().size();i++){
  //     const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
  //     for (int j=0; j!=mcell->get_allcell().size();j++){
  // 	Point p = mcell->get_allcell().at(j)->center();
  // 	x = mcell->GetTimeSlice()*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
  // 	y = p.y/units::cm;
  // 	z = p.z/units::cm;
  // 	g1->SetPoint(ncount,x,y,z);
  // 	ncount ++;
  //     }
  //   }
  //   // cout << ncount << endl;
  //   g1->Write(Form("cluster_%d",ncluster));
  //   ncluster ++;
  // }


  // save all the toy tiling stuff
  WCP2dToy::ToyTiling* tt1 = 0;
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
  double xx,yy,zz;
  ttree1->Branch("xx",&xx,"xx/D");    //don
  ttree1->Branch("yy",&yy,"yy/D");    //don
  ttree1->Branch("zz",&zz,"zz/D");    //don

  // ttree1->Branch("x",&x_save,"x/D");    //done
  //ttree1->Branch("y",&y_save,"y/D");
  //ttree1->Branch("z",&z_save,"z/D");
  

   // save information to reconstruct the toytiling
  int u_index, v_index, w_index;
  double u_charge, v_charge, w_charge;
  double u_charge_err, v_charge_err, w_charge_err;
  
  int apa_no=0, cryostat_no=0;
  int face = 0;
  ttree1->Branch("face",&face,"face/I");
  ttree1->Branch("apa_no",&apa_no,"apa_no/I");
  ttree1->Branch("cryostat_no",&cryostat_no,"cryostat_no/I");

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
      x_save = time_slice *nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
      xx = x_save;
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
	//hack for now
	charge_save = toymatrix[time_slice]->Get_Cell_Charge(mcell,1)/mcell->cross_section() * cell_save->cross_section();
	//charge_save = 1;
	
	//std::cout << time_slice << " " << x_save << std::endl;
	yy = p.y/units::cm;
  	zz = p.z/units::cm;
	ttree1->Fill();
	
	if (save_image_outline_flag==0){
	  //save the g_rec_blob tree ... 
	  x_save = xx;
	  y_save = yy;
	  z_save = zz;
	  ncharge_save = mcell->get_allcell().size();
	  g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
	  ncount2++;
	  t_rec_charge_blob->Fill();
	}else if (save_image_outline_flag==1){
	  x_save = xx;
	  y_save = yy;
	  z_save = zz;
	  type_save = 0;
	  ncharge_save = mcell->get_allcell().size();
	  t_rec_charge_blob->Fill();
	}
      }

      if (save_image_outline_flag==1){
	Point p = mcell->get_allcell().at(0)->center();
	x_save = time_slice *nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
	y_save = mcell->center().y/units::cm;
	z_save = mcell->center().z/units::cm;
	charge_save = toymatrix[time_slice]->Get_Cell_Charge(mcell,1);
	ncharge_save = mcell->get_allcell().size();
	g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
	ncount2++;
	type_save = 1;
	t_rec_charge_blob->Fill();
	
	charge_save = 0;
	ncharge_save = 0;
	
	for (int k=0;k!=mcell->get_edgecells().size();k++){
	  Point p = mcell->get_edgecells().at(k)->center();
	  y_save = p.y/units::cm;
	  z_save = p.z/units::cm;
	  g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
	  ncount2 ++;
	  type_save = 2; //boundary ...
	  t_rec_charge_blob->Fill();
	}
      }

    }
  }
  ttree1->Write();
  g->Write("shower3D");
  gt->Write("shower3D_truth");
  g_rec->Write("shower3D_charge");
  g_rec_blob->Write("shower3D_charge_blob");
 
  TTree *Trun = new TTree("Trun","Trun");
  Trun->SetDirectory(file);
  
  int detector = 0; // MicroBooNE
  Trun->Branch("detector",&detector,"detector/I");

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

  pitch_u = pitch_u/units::cm;
  pitch_v = pitch_v/units::cm;
  pitch_w = pitch_w/units::cm;
  Trun->Branch("pitch_u",&pitch_u,"pitch_u/D");
  Trun->Branch("pitch_v",&pitch_v,"pitch_v/D");
  Trun->Branch("pitch_w",&pitch_w,"pitch_w/D");

  Trun->Fill();

  TTree *TMC = new TTree("TMC","TMC");
  TMC->SetDirectory(file);
 
  TMC->Branch("mc_Ntrack", &mc_Ntrack);  // number of tracks in MC
  TMC->Branch("mc_id", &mc_id, "mc_id[mc_Ntrack]/I");  // track id; size == mc_Ntrack
  TMC->Branch("mc_pdg", &mc_pdg, "mc_id[mc_Ntrack]/I");  // track particle pdg; size == mc_Ntrack
  TMC->Branch("mc_process", &mc_process, "mc_process[mc_Ntrack]/I");  // track generation process code; size == mc_Ntrack
  TMC->Branch("mc_mother", &mc_mother, "mc_mother[mc_Ntrack]/I");  // mother id of this track; size == mc_Ntrack
  TMC->Branch("mc_daughters", mc_daughters);  // daughters id of this track; vector
  TMC->Branch("mc_startXYZT", &mc_startXYZT, "mc_startXYZT[mc_Ntrack][4]/F");  // start position of this track; size == mc_Ntrack
  TMC->Branch("mc_endXYZT", &mc_endXYZT, "mc_endXYZT[mc_Ntrack][4]/F");  // start position of this track; size == mc_Ntrack
  TMC->Branch("mc_startMomentum", &mc_startMomentum, "mc_startMomentum[mc_Ntrack][4]/F");  // start momentum of this track; size == mc_Ntrack
  TMC->Branch("mc_endMomentum", &mc_endMomentum, "mc_endMomentum[mc_Ntrack][4]/F");  // start momentum of this track; size == mc_Ntrack
  TMC->Branch("mc_trackPosition",mc_trackPosition);


  TMC->Fill();



  file->Write();
  file->Close();

  toymetric.Print();
  blobmetric.Print();

  // cin >> abc;
  
  return 0;
  
} // main()
