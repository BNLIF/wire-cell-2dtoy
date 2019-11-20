#include "WCPSst/GeomDataSource.h"
#include "WCPSst/DatauBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCPSst/uBooNESliceDataSource.h"

#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/BadTiling.h"
#include "WCP2dToy/LowmemTiling.h"
#include "WCP2dToy/uBooNE_L1SP.h"
#include "WCP2dToy/WCPHolder.h"

#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/TruthToyTiling.h"
#include "WCP2dToy/SimpleBlobToyTiling.h"

#include "WCP2dToy/ChargeSolving.h"
#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/ToyMatrixExclusive.h"
#include "WCP2dToy/ToyMatrixKalman.h"
#include "WCP2dToy/ToyMatrixIterate.h"
#include "WCP2dToy/ToyMatrixIterate_SingleWire.h"
#include "WCP2dToy/ToyMatrixIterate_Only.h"

#include "WCP2dToy/ToyMatrixMarkov.h"
#include "WCP2dToy/ToyMetric.h"
#include "WCP2dToy/BlobMetric.h"
#include "WCPData/TPCParams.h"
#include "WCPData/Singleton.h"

#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"

#include "WCPData/Slim3DCluster.h"
//#include "WCPNav/SliceDataSource.h"


#include "WCPNav/FrameDataSource.h"
#include "WCPNav/SimDataSource.h"
#include "WCPNav/SliceDataSource.h"
#include "WCPSst/Util.h"
#include "WCPData/SimTruth.h"
#include "WCP2dToy/ToyDepositor.h"
#include "WCPNav/GenerativeFDS.h"
#include "WCP2dToy/ToySignalSimu.h"
#include "WCP2dToy/ToySignalSimuTrue.h"
#include "WCP2dToy/DataSignalGaus_ROI.h"
#include "WCP2dToy/DataSignalWien_ROI.h"

#include "WCP2dToy/uBooNE_Data_2D_Deconvolution.h"
#include "WCP2dToy/uBooNE_Data_ROI.h"
#include "WCP2dToy/uBooNE_Data_After_ROI.h"
#include "WCP2dToy/uBooNE_Data_After_ROI_gaus.h"
#include "WCP2dToy/pd_Data_FDS.h"
#include "WCP2dToy/uBooNE_Data_Error.h"
#include "WCP2dToy/ExecMon.h"

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

bool GeomWireSelectionCompare(GeomWireSelection a, GeomWireSelection b) {
  if (a.size() > b.size()){
    return true;
  }else if (a.size() < b.size()){
    return false;
  }else{
    for (int i=0;i!=a.size();i++){
      if (a.at(i)->ident() > b.at(i)->ident()){
        return true;
      }else if (a.at(i)->ident() < b.at(i)->ident()){
        return false;
      }
    }
    return false;
  }
  return false;
}


int main(int argc, char* argv[])
{
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root -t[0,1] -s[0,1,2]" << endl;
    return 1;
  }

  TH1::AddDirectory(kFALSE);

  int two_plane = 0;
  int save_file = 0;
  int nt_off1 = 0;
  int nt_off2 = 0;
  int solve_charge = 1;
  for(Int_t i = 1; i != argc; i++){
     switch(argv[i][1]){
     case 't':
       two_plane = atoi(&argv[i][2]); 
       break;
     case 's':
       save_file = atoi(&argv[i][2]); 
       break;
     case 'd':
       solve_charge = atoi(&argv[i][2]); 
       break;
     case 'a':
       nt_off1 = atoi(&argv[i][2]);
       break;
     case 'b':
       nt_off2 = atoi(&argv[i][2]);
       break;
     }
  }
  
  if (two_plane)
    cout << "Enable Two Plane Reconstruction " << endl; 
  
  ExecMon em("starting");
  cerr << em("load geometry") << endl;

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
      
  //Note: the above one is still on the high end, 
  float unit_dis = 1.101; // match 256 cm

  //**** time offset for 58kV ****// 
  float toffset_1=0.0; //(nt_off1 * 0.2 - 1.0 );  // time offset between u/v 
  float toffset_2=0.0; //(nt_off2 * 0.2 - 1.0); // time offset between u/w
  float toffset_3=0.0;
  
  int save_image_outline_flag = 0; // prescale flag 
  
  int total_time_bin = 9592;
  int recon_threshold = 2000;
  int frame_length = 3200;
  int max_events = 100;
  int eve_num;
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
  

  int time_offset = 4; // Now the time offset is taken care int he signal processing, so we just need the overall offset ... 
  

  const char* root_file = argv[2];  
  int run_no, subrun_no, event_no;
  
  cerr << em("load data") << endl;

  // load Trun
  TFile *file1 = new TFile(root_file);
  TTree *Trun = (TTree*)file1->Get("Trun");
  
  Trun->SetBranchAddress("eventNo",&event_no);
  Trun->SetBranchAddress("runNo",&run_no);
  Trun->SetBranchAddress("subRunNo",&subrun_no);
  Trun->SetBranchAddress("unit_dis",&unit_dis);
  Trun->SetBranchAddress("toffset_uv",&toffset_1);
  Trun->SetBranchAddress("toffset_uw",&toffset_2);
  Trun->SetBranchAddress("toffset_u",&toffset_3);
  Trun->SetBranchAddress("total_time_bin",&total_time_bin);
  // Trun->SetBranchAddress("recon_threshold",&recon_threshold);
  Trun->SetBranchAddress("frame_length",&frame_length);
  //Trun->SetBranchAddress("max_events",&max_events);
  Trun->SetBranchAddress("eve_num",&eve_num);
  Trun->SetBranchAddress("nrebin",&nrebin);
  //Trun->SetBranchAddress("threshold_u",&threshold_u);
  //Trun->SetBranchAddress("threshold_v",&threshold_v);
  //Trun->SetBranchAddress("threshold_w",&threshold_w);
  //Trun->SetBranchAddress("threshold_ug",&threshold_ug);
  //Trun->SetBranchAddress("threshold_vg",&threshold_vg);
  //Trun->SetBranchAddress("threshold_wg",&threshold_wg);
  Trun->SetBranchAddress("time_offset",&time_offset);
  
  Trun->GetEntry(0);

  //  TTree *T_op = (TTree*)file1->Get("T_op");
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  int nwire_u = wires_u.size();
  int nwire_v = wires_v.size();
  int nwire_w = wires_w.size();


  cout << "Run No: " << run_no << " " << subrun_no << " " << event_no << endl;
  ChirpMap uplane_map;
  ChirpMap vplane_map;
  ChirpMap wplane_map;
  TTree *T_chirp = (TTree*)file1->Get("T_bad");
  Int_t chid=0, plane=0;
  Int_t start_time=0, end_time=0;
  T_chirp->SetBranchAddress("chid",&chid);
  T_chirp->SetBranchAddress("plane",&plane);
  T_chirp->SetBranchAddress("start_time",&start_time);
  T_chirp->SetBranchAddress("end_time",&end_time);
  for (Int_t i=0;i!=T_chirp->GetEntries();i++){
    T_chirp->GetEntry(i);
    
    std::pair<int,int> abc(start_time,end_time);
    if (plane == 0){
      uplane_map[chid] = abc;
    }else if (plane == 1){
      vplane_map[chid-nwire_u] = abc;
    }else if (plane == 2){
      wplane_map[chid-nwire_u-nwire_v] = abc;
    }
  }
  
  std::vector<float> uplane_rms;
  std::vector<float> vplane_rms;
  std::vector<float> wplane_rms;
  // Note, there is a mismatch here
  // These RMS values are for single ticks
  // Later they are applied to the rebinned data
  // Probably OK as the TPC signal processing already took care the fake hits ... 
  TH1F *hu_threshold = (TH1F*)file1->Get("hu_threshold");
  TH1F *hv_threshold = (TH1F*)file1->Get("hv_threshold");
  TH1F *hw_threshold = (TH1F*)file1->Get("hw_threshold");
  for (int i=0;i!=hu_threshold->GetNbinsX();i++){
    uplane_rms.push_back(hu_threshold->GetBinContent(i+1));
  }
  for (int i=0;i!=hv_threshold->GetNbinsX();i++){
    vplane_rms.push_back(hv_threshold->GetBinContent(i+1));
  }
  for (int i=0;i!=hw_threshold->GetNbinsX();i++){
    wplane_rms.push_back(hw_threshold->GetBinContent(i+1));
  }
  
  TH2F *hu_decon = (TH2F*)file1->Get("hu_decon");
  TH2F *hv_decon = (TH2F*)file1->Get("hv_decon");
  TH2F *hw_decon = (TH2F*)file1->Get("hw_decon");

  TH2F *hu_raw = (TH2F*)file1->Get("hu_raw");
  TH2F *hv_raw = (TH2F*)file1->Get("hv_raw");
  TH2F *hw_raw = (TH2F*)file1->Get("hw_raw");
  
  WCP2dToy::pdDataFDS roi_fds(gds,hu_decon,hv_decon,hw_decon,eve_num);
  roi_fds.jump(eve_num);

  TH2F *hu_decon_g = (TH2F*)file1->Get("hu_decon");
  TH2F *hv_decon_g = (TH2F*)file1->Get("hv_decon");
  TH2F *hw_decon_g = (TH2F*)file1->Get("hw_decon");

  //The last input is the time offset
  // if it is + 2.8, it will shift L1SP decon result to later time by 2.8 us
  // if it is -2.8, it will shift L1SP decon result to early time by 2.8 us
  // This is to take into account potential time difference between
  // WCT standard SP and WCP L1SP (based on WCP standard SP) ...
  WCP2dToy::uBooNE_L1SP l1sp(hv_raw,hv_decon,hv_decon_g,nrebin,0);
  //  WCP2dToy::uBooNE_L1SP l1sp(hv_raw,hv_decon,hv_decon_g,nrebin,-2.8);
  
  WCP2dToy::pdDataFDS roi_gaus_fds(gds,hu_decon_g,hv_decon_g,hw_decon_g,eve_num);
  roi_gaus_fds.jump(eve_num);

  WCP2dToy::uBooNEDataError error_fds(gds,hu_decon_g, hv_decon_g, hw_decon_g, eve_num, nrebin);
  error_fds.jump(eve_num);
  
  


  WCPSst::uBooNESliceDataSource sds(roi_fds,roi_gaus_fds,error_fds,
					 threshold_u, threshold_v, threshold_w,
					 nwire_u, nwire_v, nwire_w,
					 &uplane_rms, &vplane_rms, &wplane_rms); 
  
 
  
  cerr << em("begin tiling") << endl;

  
  int ncount = 0;
  int ncount1 = 0;
  int ncount2 = 0;

  int ncount_t = 0;
  

 
  WCP2dToy::LowmemTiling **lowmemtiling = new WCP2dToy::LowmemTiling*[2400];
  WCP2dToy::ChargeSolving **chargesolver = new WCP2dToy::ChargeSolving*[2400];
  
  WCP2dToy::WCPHolder WCholder;

  //add in cluster
  Slim3DClusterSet cluster_set, cluster_delset;
  
  int ncount_mcell = 0;

 
  int start_num = 0 ;
  int end_num = sds.size()-1;

  // start_num = 1925;
  // end_num = 1925;

  // start_num = 1665;
  // end_num = 1665;
  
  // start_num = 6600/4.;
  // end_num = 7500/4.;
  
  // start_num = 100;
  // end_num = 105;
  // end_num = 150;
  
  TFile *file = new TFile(Form("result_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
  
  TTree *Twc = new TTree("Twc","Twc");
  Twc->SetDirectory(file);
  Int_t n_cells;
  Int_t n_good_wires;
  Int_t n_bad_wires;
  Int_t time_slice;
  Int_t n_single_cells;
  Int_t ndirect_solved;
  Int_t nL1_solved;
  Int_t total_ndf;
  Double_t total_chi2;

  Int_t L1_ndf[1000];
  Int_t direct_ndf[1000];
  Double_t L1_chi2_base[1000];
  Double_t L1_chi2_penalty[1000];
  Double_t direct_chi2[1000];
  
  Twc->Branch("time_slice",&time_slice,"time_slice/I");
  Twc->Branch("n_cells",&n_cells,"n_cells/I");
  Twc->Branch("n_good_wires",&n_good_wires,"n_good_wires/I");  
  Twc->Branch("n_bad_wires",&n_bad_wires,"n_bad_wires/I");
  Twc->Branch("n_single_cells",&n_single_cells,"n_single_cells/I");
  
  Twc->Branch("ndirect_solved",&ndirect_solved,"ndirect_solved/I");
  Twc->Branch("nL1_solved",&nL1_solved,"nL1_solved/I");
  Twc->Branch("total_ndf",&total_ndf,"total_ndf/I");
  Twc->Branch("total_chi2",&total_chi2,"total_chi2/D");

  Twc->Branch("L1_ndf",L1_ndf,"L1_ndf[nL1_solved]/I");
  Twc->Branch("L1_chi2_base",L1_chi2_base,"L1_chi2_base[nL1_solved]/D");
  Twc->Branch("L1_chi2_penalty",L1_chi2_penalty,"L1_chi2_penalty[nL1_solved]/D");
  Twc->Branch("direct_ndf",direct_ndf,"direct_ndf[ndirect_solved]/I");
  Twc->Branch("direct_chi2",direct_chi2,"direct_chi2[ndirect_solved]/D");
  
  //test 
  // uplane_map.begin()->second.second=5000;

  for (int i=start_num;i!=end_num+1;i++){
    if (i%200==0)
      std::cout << "Tiling: " << i << std::endl;

    sds.jump(i);
    WCP::Slice slice = sds.get();
    WCP::Slice slice_err = sds.get_error();

    // WCP::Channel::Group group = slice.group();
    // WCP::Channel::Group group_err = slice_err.group();
    //double sum = 0;
    // for (int i=0;i!=group.size();i++){
    //   std::cout << group.at(i).first << " " << group.at(i).second << " " << group_err.at(i).first << " " << group_err.at(i).second << std::endl;
    // }
    
    lowmemtiling[i] = new WCP2dToy::LowmemTiling(i,nrebin,gds,WCholder);
    if (i==start_num){
      lowmemtiling[i]->init_bad_cells(uplane_map,vplane_map,wplane_map);
    }else{
      lowmemtiling[i]->check_bad_cells(lowmemtiling[i-1],uplane_map,vplane_map,wplane_map);
    }
    lowmemtiling[i]->init_good_cells(slice,slice_err,uplane_rms,vplane_rms,wplane_rms);

    // std::cout << i << " " << lowmemtiling[i]->get_three_good_wire_cells().size() << " " << lowmemtiling[i]->get_two_good_wire_cells().size() << " " << lowmemtiling[i]->get_all_good_wires().size() << " " << lowmemtiling[i]->get_all_bad_wires().size() << " " << lowmemtiling[i]->get_cell_wires_map().size() << " " << lowmemtiling[i]->get_wire_cells_map().size() << " " << lowmemtiling[i]->get_wire_charge_map().size() << " " << lowmemtiling[i]->get_wire_charge_error_map().size() << " " << lowmemtiling[i]->get_two_bad_wire_cells().size() << " " << lowmemtiling[i]->get_wire_type_map().size() << " " << lowmemtiling[i]->get_wire_pwire_map().size() << " " << lowmemtiling[i]->get_pwire_wires_map().size() << std::endl;
    // lowmemtiling[i]->reset_cells();
    
    // if (i==start_num){
    //   lowmemtiling[i]->init_bad_cells(uplane_map,vplane_map,wplane_map);
    // }else{
    //   lowmemtiling[i]->check_bad_cells(lowmemtiling[i-1],uplane_map,vplane_map,wplane_map);
    // }
    // lowmemtiling[i]->init_good_cells(slice,slice_err,uplane_rms,vplane_rms,wplane_rms);
    
    // std::cout << i << " " << lowmemtiling[i]->get_three_good_wire_cells().size() << " " << lowmemtiling[i]->get_two_good_wire_cells().size() << " " << lowmemtiling[i]->get_all_good_wires().size() << " " << lowmemtiling[i]->get_all_bad_wires().size() << " " << lowmemtiling[i]->get_cell_wires_map().size() << " " << lowmemtiling[i]->get_wire_cells_map().size() << " " << lowmemtiling[i]->get_wire_charge_map().size() << " " << lowmemtiling[i]->get_wire_charge_error_map().size() << " " << lowmemtiling[i]->get_two_bad_wire_cells().size() << " " << lowmemtiling[i]->get_wire_type_map().size() << " " << lowmemtiling[i]->get_wire_pwire_map().size() << " " << lowmemtiling[i]->get_pwire_wires_map().size()<< std::endl;
    
    GeomWireSelection wires = lowmemtiling[i]->find_L1SP_wires();
    l1sp.AddWires(i,wires);
    
  }
   cerr << em("finish tiling") << endl;
   
  l1sp.AddWireTime_Raw();
  cerr << em("finish raw signal examination") << endl;
  
  l1sp.Form_rois(6);
  cerr << em("finish L1SP") << endl;
  
  
  roi_fds.refresh(hu_decon,hv_decon,hw_decon,eve_num);
  roi_gaus_fds.refresh(hu_decon_g,hv_decon_g,hw_decon_g,eve_num);
  error_fds.refresh(hu_decon_g, hv_decon_g, hw_decon_g, eve_num);
  
  std::set<int> time_slice_set = l1sp.get_time_slice_set();

  //std::cout << time_slice_set.size() << std::endl;
  for (auto it = time_slice_set.begin(); it!= time_slice_set.end(); it++){
    int time_slice = *it;
    if (time_slice >= start_num && time_slice <=end_num){

      
      sds.jump(time_slice);
      WCP::Slice slice = sds.get();
      WCP::Slice slice_err = sds.get_error();

      std::cout << time_slice << " " << slice.group().size() << " " << slice_err.group().size() << std::endl;
      // if (time_slice == 2070){
      // 	WCP::Channel::Group group = slice.group();
      // 	for (int k=0;k!=group.size();k++){
      // 	  std::cout << k << " " << group.at(k).first << " " << group.at(k).second << std::endl;
      // 	}
      // }
      //  std::cout << time_slice << " " << lowmemtiling[time_slice]->get_three_good_wire_cells().size() << " " << lowmemtiling[time_slice]->get_two_good_wire_cells().size() << " " << lowmemtiling[time_slice]->get_all_good_wires().size() << " " << lowmemtiling[time_slice]->get_all_bad_wires().size() << " " << lowmemtiling[time_slice]->get_cell_wires_map().size() << " " << lowmemtiling[time_slice]->get_wire_cells_map().size() << " " << lowmemtiling[time_slice]->get_wire_charge_map().size() << " " << lowmemtiling[time_slice]->get_two_bad_wire_cells().size() << std::endl;
       lowmemtiling[time_slice]->reset_cells();
       if (time_slice==start_num){
	 lowmemtiling[time_slice]->init_bad_cells(uplane_map,vplane_map,wplane_map);
       }else{
	 lowmemtiling[time_slice]->check_bad_cells(lowmemtiling[time_slice-1],uplane_map,vplane_map,wplane_map);
       }
       lowmemtiling[time_slice]->init_good_cells(slice,slice_err,uplane_rms,vplane_rms,wplane_rms);
      // std::cout << time_slice << " " << lowmemtiling[time_slice]->get_three_good_wire_cells().size() << " " << lowmemtiling[time_slice]->get_two_good_wire_cells().size() << " " << lowmemtiling[time_slice]->get_all_good_wires().size() << " " << lowmemtiling[time_slice]->get_all_bad_wires().size() << " " << lowmemtiling[time_slice]->get_cell_wires_map().size() << " " << lowmemtiling[time_slice]->get_wire_cells_map().size() << " " << lowmemtiling[time_slice]->get_wire_charge_map().size() << " " << lowmemtiling[time_slice]->get_two_bad_wire_cells().size() << std::endl;
      // std::cout << time_slice << " " << slice.group().size() << " ";
      // slice = sds.get();
      // slice_err = sds.get_error();
      // std::cout << slice.group().size() << std::endl;
      
    }
  }
  // for (int i=start_num;i!=end_num+1;i++){
  //    sds.jump(i);
  //    WCP::Slice slice = sds.get();
  //    WCP::Slice slice_err = sds.get_error();
     
  //    lowmemtiling[i]->init_good_cells(slice,slice_err,uplane_rms,vplane_rms,wplane_rms);
  // }
  
  cerr << em("finish tiling") << endl;

  TFile *file2 = new TFile("temp_l1sp.root","RECREATE");
  hu_raw->SetDirectory(file2);
  hv_raw->SetDirectory(file2);
  hw_raw->SetDirectory(file2);
  hu_decon->SetDirectory(file2);
  hv_decon->SetDirectory(file2);
  hw_decon->SetDirectory(file2);
  hu_decon_g->SetDirectory(file2);
  hv_decon_g->SetDirectory(file2);
  hw_decon_g->SetDirectory(file2);
  hu_threshold->SetDirectory(file2);
  hv_threshold->SetDirectory(file2);
  hw_threshold->SetDirectory(file2);
  file2->Write();
  file2->Close();
  
  

  
} // main()
