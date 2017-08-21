#include "WireCellSst/GeomDataSource.h"
#include "WireCellSst/DatauBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCellSst/uBooNESliceDataSource.h"

#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/BadTiling.h"
#include "WireCell2dToy/LowmemTiling.h"
#include "WireCell2dToy/WireCellHolder.h"

#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCell2dToy/SimpleBlobToyTiling.h"

#include "WireCell2dToy/ChargeSolving.h"
#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/ToyMatrixExclusive.h"
#include "WireCell2dToy/ToyMatrixKalman.h"
#include "WireCell2dToy/ToyMatrixIterate.h"
#include "WireCell2dToy/ToyMatrixIterate_SingleWire.h"
#include "WireCell2dToy/ToyMatrixIterate_Only.h"


#include "WireCell2dToy/ToyMatrixMarkov.h"
#include "WireCell2dToy/ToyMetric.h"
#include "WireCell2dToy/BlobMetric.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"

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
#include "WireCell2dToy/DataSignalGaus_ROI.h"
#include "WireCell2dToy/DataSignalWien_ROI.h"

#include "WireCell2dToy/uBooNE_Data_2D_Deconvolution.h"
#include "WireCell2dToy/uBooNE_Data_ROI.h"
#include "WireCell2dToy/uBooNE_Data_After_ROI.h"
#include "WireCell2dToy/uBooNE_Data_After_ROI_gaus.h"
#include "WireCell2dToy/pd_Data_FDS.h"
#include "WireCell2dToy/uBooNE_Data_Error.h"
#include "WireCell2dToy/ExecMon.h"

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
  Trun->SetBranchAddress("recon_threshold",&recon_threshold);
  Trun->SetBranchAddress("frame_length",&frame_length);
  Trun->SetBranchAddress("max_events",&max_events);
  Trun->SetBranchAddress("eve_num",&eve_num);
  Trun->SetBranchAddress("nrebin",&nrebin);
  Trun->SetBranchAddress("threshold_u",&threshold_u);
  Trun->SetBranchAddress("threshold_v",&threshold_v);
  Trun->SetBranchAddress("threshold_w",&threshold_w);
  Trun->SetBranchAddress("threshold_ug",&threshold_ug);
  Trun->SetBranchAddress("threshold_vg",&threshold_vg);
  Trun->SetBranchAddress("threshold_wg",&threshold_wg);
  Trun->SetBranchAddress("time_offset",&time_offset);
  
  Trun->GetEntry(0);

  TTree *T_op = (TTree*)file1->Get("T_op");
  
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
  TTree *T_chirp = (TTree*)file1->Get("T_chirp");
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
  
  TH2I *hu_decon = (TH2I*)file1->Get("hu_decon");
  TH2I *hv_decon = (TH2I*)file1->Get("hv_decon");
  TH2I *hw_decon = (TH2I*)file1->Get("hw_decon");
  
  WireCell2dToy::pdDataFDS roi_fds(gds,hu_decon,hv_decon,hw_decon,eve_num);
  roi_fds.jump(eve_num);

  TH2I *hu_decon_g = (TH2I*)file1->Get("hu_decon_g");
  TH2I *hv_decon_g = (TH2I*)file1->Get("hv_decon_g");
  TH2I *hw_decon_g = (TH2I*)file1->Get("hw_decon_g");
  
  WireCell2dToy::pdDataFDS roi_gaus_fds(gds,hu_decon_g,hv_decon_g,hw_decon_g,eve_num);
  roi_gaus_fds.jump(eve_num);

  WireCell2dToy::uBooNEDataError error_fds(gds,hu_decon_g, hv_decon_g, hw_decon_g, eve_num, nrebin);
  error_fds.jump(eve_num);
  
  
  // WireCellSst::ToyuBooNESliceDataSource sds(roi_fds,roi_gaus_fds,threshold_u, 
  // 					    threshold_v, threshold_w, 
  // 					    threshold_ug, 
  // 					    threshold_vg, threshold_wg, 
  // 					    nwire_u, 
  // 					    nwire_v, nwire_w,
  // 					    &uplane_rms, &vplane_rms, &wplane_rms); 

  WireCellSst::uBooNESliceDataSource sds(roi_fds,roi_gaus_fds,error_fds,
					 threshold_u, threshold_v, threshold_w,
					 nwire_u, nwire_v, nwire_w,
					 &uplane_rms, &vplane_rms, &wplane_rms); 
  
  // sds.jump(100);
  // full_sds.jump(100);
  // WireCell::Slice slice = sds.get();
  // WireCell::Slice slice1 = full_sds.get();
  // WireCell::Slice slice2 = full_sds.get_error();

  // WireCell::Channel::Group group = slice.group();
  // WireCell::Channel::Group group1 = slice1.group();
  // WireCell::Channel::Group group2 = slice2.group();
  // for (int i=0;i!=group.size();i++){
  //   std::cout << group.at(i).second << " " << group1.at(i).second << " " << group2.at(i).second << std::endl;
  // }
  
  cerr << em("begin tiling") << endl;

  
  int ncount = 0;
  int ncount1 = 0;
  int ncount2 = 0;

  int ncount_t = 0;
  

  // WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
  // WireCell2dToy::BadTiling **badtiling = new WireCell2dToy::BadTiling*[2400];
  // WireCell2dToy::MergeToyTiling **mergetiling = new WireCell2dToy::MergeToyTiling*[2400];
  // WireCell2dToy::ToyMatrix **toymatrix = new WireCell2dToy::ToyMatrix*[2400];
  WireCell2dToy::LowmemTiling **lowmemtiling = new WireCell2dToy::LowmemTiling*[2400];
  WireCell2dToy::ChargeSolving **chargesolver = new WireCell2dToy::ChargeSolving*[2400];
  
  WireCell2dToy::WireCellHolder WCholder;

  //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  
  int ncount_mcell = 0;

 
  int start_num = 0 ;
  int end_num = sds.size()-1;

  // start_num = 1925;
  // end_num = 1925;

  // start_num = 650;
  // end_num = 650;

  // end_num = 500;
  
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
  Twc->Branch("time_slice",&time_slice,"time_slice/I");
  Twc->Branch("n_cells",&n_cells,"n_cells/I");
  Twc->Branch("n_good_wires",&n_good_wires,"n_good_wires/I");  
  Twc->Branch("n_bad_wires",&n_bad_wires,"n_bad_wires/I");
  Twc->Branch("n_single_cells",&n_single_cells,"n_single_cells/I");
  Twc->Branch("ndirect_solved",&ndirect_solved,"ndirect_solved/I");
  Twc->Branch("nL1_solved",&nL1_solved,"nL1_solved/I");
  //test 
  // uplane_map.begin()->second.second=5000;

  for (int i=start_num;i!=end_num+1;i++){
    if (i%50==0)
      std::cout << i << std::endl;

    sds.jump(i);
    WireCell::Slice slice = sds.get();
    WireCell::Slice slice_err = sds.get_error();
    
    lowmemtiling[i] = new WireCell2dToy::LowmemTiling(i,nrebin,gds,WCholder);
    if (i==start_num){
      lowmemtiling[i]->init_bad_cells(uplane_map,vplane_map,wplane_map);
    }else{
      lowmemtiling[i]->check_bad_cells(lowmemtiling[i-1],uplane_map,vplane_map,wplane_map);
    }
    lowmemtiling[i]->init_good_cells(slice,slice_err,uplane_rms,vplane_rms,wplane_rms);
    
    //std::cout << lowmemtiling[i]->get_cell_wires_map().size() << " " << lowmemtiling[i]->get_wire_cells_map().size() << std::endl;

    // refine the merge cells 
    //lowmemtiling[i]->DivideWires(3,0);
    lowmemtiling[i]->MergeWires();

    // create individual cells ...
    //GeomCellSelection single_cells = lowmemtiling[i]->create_single_cells();
    time_slice = i;
    n_cells = lowmemtiling[i]->get_cell_wires_map().size();//get_all_cell_centers().size();
    n_good_wires = lowmemtiling[i]->get_all_good_wires().size();
    n_bad_wires = lowmemtiling[i]->get_all_bad_wires().size();
    //n_single_cells = single_cells.size();
    

    // L1 solving
    chargesolver[i] = new WireCell2dToy::ChargeSolving(gds, *lowmemtiling[i]);
    ndirect_solved = chargesolver[i]->get_ndirect_solved();
    nL1_solved = chargesolver[i]->get_nL1_solved();

    Twc->Fill();
    
    // std::cout << lowmemtiling[i]->get_cell_wires_map().size()<< " " << lowmemtiling[i]->get_all_good_wires().size() << " " << lowmemtiling[i]->get_all_bad_wires().size() << " " << lowmemtiling[i]->get_all_cell_centers().size() << " " << lowmemtiling[i]->get_wire_cells_map().size() << " " << single_cells.size() << std::endl;
    
    //std::cout << lowmemtiling[i]->get_cell_wires_map().size() << " " << lowmemtiling[i]->get_wire_cells_map().size() << std::endl;

    
    
    // //    std::cout << lowmemtiling[i]->get_three_good_wire_cells().size() << std::endl;
    // std::vector<GeomWireSelection> vec1_wires;
    // std::vector<GeomWireSelection> vec2_wires;
    // // GeomWireSelection dwires;
    // for (int j=0;j!=lowmemtiling[i]->get_three_good_wire_cells().size();j++){
    //   //std::cout << "N: " << ((SlimMergeGeomCell*)lowmemtiling[i]->get_three_good_wire_cells().at(j))->get_uwires().size() << " "
    //   //<< ((SlimMergeGeomCell*)lowmemtiling[i]->get_three_good_wire_cells().at(j))->get_vwires().size() << " " 
    //   //<< ((SlimMergeGeomCell*)lowmemtiling[i]->get_three_good_wire_cells().at(j))->get_wwires().size() << " " << std::endl;
    //   GeomWireSelection dwires;
    //   //   if (j==0){
    //   for (int k=0;k!=((SlimMergeGeomCell*)lowmemtiling[i]->get_three_good_wire_cells().at(j))->get_uwires().size();k++){
    // 	dwires.push_back(((SlimMergeGeomCell*)lowmemtiling[i]->get_three_good_wire_cells().at(j))->get_uwires().at(k));
    //   }
    //   for (int k=0;k!=((SlimMergeGeomCell*)lowmemtiling[i]->get_three_good_wire_cells().at(j))->get_vwires().size();k++){
    // 	dwires.push_back(((SlimMergeGeomCell*)lowmemtiling[i]->get_three_good_wire_cells().at(j))->get_vwires().at(k));
    //   }
    //   for (int k=0;k!=((SlimMergeGeomCell*)lowmemtiling[i]->get_three_good_wire_cells().at(j))->get_wwires().size();k++){
    // 	dwires.push_back(((SlimMergeGeomCell*)lowmemtiling[i]->get_three_good_wire_cells().at(j))->get_wwires().at(k));
    //   }
    //   // std::cout << dwires.size() << std::endl;
    //   // for (int j=0;j!=dwires.size();j++){
    //   // 	std::cout << dwires.at(j)->plane() << " " << dwires.at(j)->index() << std::endl;
    //   // }
    //   sort_by_ident(dwires);
    //   //}
      
    //   vec1_wires.push_back(dwires);
    // 	//      dwires.insert(dwires.begin(),((SlimMergeGeomCell*)lowmemtiling[i]->get_three_good_wire_cells().at(j))->get_uwires().begin(),((SlimMergeGeomCell*)lowmemtiling[i]->get_three_good_wire_cells().at(j))->get_uwires().end());
    // }


    // toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds,0.15,0.2,0.1,threshold_ug,threshold_vg, threshold_wg, &uplane_rms, &vplane_rms, &wplane_rms);
    
    // // toytiling[i]->twoplane_tiling(i,nrebin,gds,uplane_rms,vplane_rms,wplane_rms, uplane_map, vplane_map, wplane_map);


    // // // GeomCellSelection allcell = toytiling[i]->get_allcell();
    // // // GeomWireSelection allwire = toytiling[i]->get_allwire();
    // // // cout << i << " " << allcell.size() << " " << allwire.size() << endl;

    // mergetiling[i] = new WireCell2dToy::MergeToyTiling(*toytiling[i],i,3);
    
    // if (lowmemtiling[i]->get_cell_wires_map().size()!= mergetiling[i]->get_allcell().size() ||
    // 	lowmemtiling[i]->get_wire_cells_map().size()!= mergetiling[i]->get_wire_all().size()){

    //   std::cout << i << " " << lowmemtiling[i]->get_cell_wires_map().size() << " " << lowmemtiling[i]->get_wire_cells_map().size()<< " " << mergetiling[i]->get_allcell().size() << " " << mergetiling[i]->get_wire_all().size() << std::endl;
    // }

    // for (int j=0;j!=mergetiling[i]->get_allcell().size();j++){
    //   // std::cout << "O: " << ((MergeGeomCell*)mergetiling[i]->get_allcell().at(j))->get_uwires().size() << " " 
    //   // 		<< ((MergeGeomCell*)mergetiling[i]->get_allcell().at(j))->get_vwires().size() << " " 
    //   // 		<< ((MergeGeomCell*)mergetiling[i]->get_allcell().at(j))->get_wwires().size() << " " << std::endl;
    //   GeomWireSelection dwires;
    //   for (int k = 0; k!= ((MergeGeomCell*)mergetiling[i]->get_allcell().at(j))->get_uwires().size(); k++){
    // 	dwires.push_back( ((MergeGeomCell*)mergetiling[i]->get_allcell().at(j))->get_uwires().at(k));
    //   }
    //   for (int k = 0; k!= ((MergeGeomCell*)mergetiling[i]->get_allcell().at(j))->get_vwires().size(); k++){
    // 	dwires.push_back( ((MergeGeomCell*)mergetiling[i]->get_allcell().at(j))->get_vwires().at(k));
    //   }
    //   for (int k = 0; k!= ((MergeGeomCell*)mergetiling[i]->get_allcell().at(j))->get_wwires().size(); k++){
    // 	dwires.push_back( ((MergeGeomCell*)mergetiling[i]->get_allcell().at(j))->get_wwires().at(k));
    //   }
    //   sort_by_ident(dwires);
    //   // std::cout << dwires.size() << std::endl;
    //   vec2_wires.push_back(dwires);
    // }
    
    // bool flag_test = false;
    // if (vec1_wires.size()!=vec2_wires.size())
    //   flag_test = true;
    
    // if (!flag_test){
    //   sort(vec1_wires.begin(),vec1_wires.end(),GeomWireSelectionCompare);
    //   sort(vec2_wires.begin(),vec2_wires.end(),GeomWireSelectionCompare);
    //   for (int j=0;j!=vec1_wires.size();j++){
    // 	//std::cout << j << " " << vec1_wires.at(j).size() << " " << vec2_wires.at(j).size() << " " << vec1_wires.at(j).at(0)->index() << " " << vec2_wires.at(j).at(0)->index() << std::endl;
    // 	if (vec1_wires.at(j).size()!=vec2_wires.at(j).size()){
    // 	  flag_test = true;
    // 	  //break;
    // 	}
	
    //   }
    // }
    // if(flag_test){
    //   std::cout << "Bad " << i << std::endl;
    // }
    

    // if (i==0){
    //   badtiling[i] = new WireCell2dToy::BadTiling(i,nrebin,uplane_map,vplane_map,wplane_map,gds,0,1); // 2 plane bad tiling
    //   // badtiling[i] = new WireCell2dToy::BadTiling(i,nrebin,uplane_map,vplane_map,wplane_map,gds,1,1); // 1 plane bad tiling
    // }

    //badtiling[i] = new WireCell2dToy::BadTiling(i,nrebin,uplane_map,vplane_map,wplane_map,gds);

    // // if (toymatrix[i]->Get_Solve_Flag()!=0){
    // //   toymatrix[i]->Update_pred();
    // //   toymatrix[i]->Print();
    // // }

    // //draw ... 
    // TApplication theApp("theApp",&argc,argv);
    // theApp.SetReturnFromRun(true);
    
    // TCanvas c1("ToyMC","ToyMC",800,600);
    // c1.Draw();
    
    // WireCell2dToy::ToyEventDisplay display(c1, gds);
    // display.charge_min = 0;
    // display.charge_max = 5e4;


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
    // display.draw_mc(1,WireCell::PointValueVector(),"colz");
    // display.draw_slice(slice,""); // draw wire 
    // // display.draw_wires(vec1_wires.at(64),"same"); // draw wire 
    // // // display.draw_bad_region(uplane_map,i,nrebin,0,"same");
    // // // display.draw_bad_region(vplane_map,i,nrebin,1,"same");
    // // // display.draw_bad_region(wplane_map,i,nrebin,2,"same");
    // // display.draw_bad_cell(badtiling[i]->get_cell_all());
    // display.draw_cells(single_cells,"*same");
    // //display.draw_points(lowmemtiling[i]->get_all_cell_centers(),"*");
    // //display.draw_merged_wires(lowmemtiling[i]->get_all_good_wires(),"same",2);
    // //display.draw_merged_wires(lowmemtiling[i]->get_all_bad_wires(),"same",1);
    
    // //display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell
    
    // // display.draw_wires_charge(toytiling[i]->wcmap(),"Fsame",FI);
    // // display.draw_cells_charge(toytiling[i]->get_allcell(),"Fsame");
    //  theApp.Run();
  }
  
  cerr << em("finish tiling") << endl;


  TGraph2D *g = new TGraph2D();
  TGraph2D *g_rec = new TGraph2D();
  TTree *t_rec_simple = new TTree("T_rec","T_rec");
  TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
  Double_t x_save, y_save, z_save;
  Double_t charge_save;
  Double_t ncharge_save;
  Double_t chi2_save;
  Double_t ndf_save;
  
  t_rec_simple->SetDirectory(file);
  t_rec_simple->Branch("x",&x_save,"x/D");
  t_rec_simple->Branch("y",&y_save,"y/D");
  t_rec_simple->Branch("z",&z_save,"z/D");
  
  t_rec_charge->SetDirectory(file);
  t_rec_charge->Branch("x",&x_save,"x/D");
  t_rec_charge->Branch("y",&y_save,"y/D");
  t_rec_charge->Branch("z",&z_save,"z/D");
  t_rec_charge->Branch("q",&charge_save,"q/D");
  t_rec_charge->Branch("nq",&ncharge_save,"nq/D");
  t_rec_charge->Branch("chi2",&chi2_save,"chi2/D");
  t_rec_charge->Branch("ndf",&ndf_save,"ndf/D");
  
  
  for (int i=start_num; i!=end_num+1;i++){
    GeomCellMap cell_wires_map = lowmemtiling[i]->get_cell_wires_map();
    chi2_save = chargesolver[i]->get_chi2();
    ndf_save = chargesolver[i]->get_ndf();
    
    for (auto it = cell_wires_map.begin(); it!= cell_wires_map.end(); it++){
      MergeGeomCell *mcell = (MergeGeomCell*) it->first;
      GeomCellSelection temp_cells = lowmemtiling[i]->create_single_cells((SlimMergeGeomCell*)it->first);
      for (auto it1 = temp_cells.begin(); it1!=temp_cells.end(); it1++){
	Point p = (*it1)->center();
	x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
	y_save = p.y/units::cm;
	z_save = p.z/units::cm;
	g->SetPoint(ncount,x_save, y_save, z_save);
	t_rec_simple->Fill();
	ncount ++;
      }
      
      // fill the charge ... 
      if (chargesolver[i]->get_mcell_charge(mcell)>300){
      	charge_save = chargesolver[i]->get_mcell_charge(mcell) / (temp_cells.size() * 1.0) ;
      	ncharge_save = temp_cells.size();
	
      	for (auto it1 = temp_cells.begin(); it1!=temp_cells.end(); it1++){
      	  Point p = (*it1)->center();
      	  x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
      	  y_save = p.y/units::cm;
      	  z_save = p.z/units::cm;
      	  g_rec->SetPoint(ncount1,x_save, y_save, z_save);
	  
      	  t_rec_charge->Fill();
	  
      	  ncount1 ++;
      	}
      }
    }
    // 
  }
  g->Write("g");
  g_rec->Write("g_rec");
  
  TTree *t_bad = new TTree("T_bad","T_bad");
  t_bad->SetDirectory(file);
  Int_t bad_npoints;
  Double_t bad_y[100],bad_z[100];
  t_bad->Branch("bad_npoints",&bad_npoints,"bad_npoints/I");
  t_bad->Branch("bad_y",bad_y,"bad_y[bad_npoints]/D");
  t_bad->Branch("bad_z",bad_z,"bad_z[bad_npoints]/D");

  for (int i=0; i!=lowmemtiling[start_num]->get_two_bad_wire_cells().size();i++){
    const SlimMergeGeomCell *cell = (SlimMergeGeomCell*)lowmemtiling[start_num]->get_two_bad_wire_cells().at(i);
    PointVector ps = cell->boundary();
    bad_npoints = ps.size();
    for (int j=0;j!=bad_npoints;j++){
      bad_y[j] = ps.at(j).y/units::cm;
      bad_z[j] = ps.at(j).z/units::cm;
    }
    t_bad->Fill();
  }
  
  if (T_op!=0){
    T_op->CloneTree()->Write();
  }
  Trun->CloneTree()->Write();
  file->Write();
  file->Close();
  
  return 0;
  
} // main()
