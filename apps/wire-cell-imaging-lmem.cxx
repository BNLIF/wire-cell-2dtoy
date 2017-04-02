#include "WireCellSst/GeomDataSource.h"
#include "WireCellSst/DatauBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/BadTiling.h"
#include "WireCell2dToy/LowmemTiling.h"

#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCell2dToy/SimpleBlobToyTiling.h"

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
#include "WireCell2dToy/pd_Data_FDS.h"
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
  
  int total_time_bin=9594;
  int recon_threshold = 2000;
  int frame_length = 3200;
  int max_events = 100;
  int eve_num;
  int nrebin = 6;

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

  WireCellSst::ToyuBooNESliceDataSource sds(roi_fds,roi_fds,threshold_u, 
  					    threshold_v, threshold_w, 
  					    threshold_ug, 
  					    threshold_vg, threshold_wg, 
  					    nwire_u, 
  					    nwire_v, nwire_w,
  					    &uplane_rms, &vplane_rms, &wplane_rms); 
    
  
  cerr << em("begin tiling") << endl;

  
  int ncount = 0;
  int ncount1 = 0;
  int ncount2 = 0;

  int ncount_t = 0;
  

  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
  WireCell2dToy::BadTiling **badtiling = new WireCell2dToy::BadTiling*[2400];
  WireCell2dToy::MergeToyTiling **mergetiling = new WireCell2dToy::MergeToyTiling*[2400];
  WireCell2dToy::ToyMatrix **toymatrix = new WireCell2dToy::ToyMatrix*[2400];
  WireCell2dToy::LowmemTiling **lowmemtiling = new WireCell2dToy::LowmemTiling*[2400];


  //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  
  int ncount_mcell = 0;

 
  int start_num = 0 ;
  int end_num = sds.size()-1;

  start_num = 555;
  end_num = 555;
  

  //test 
  uplane_map.begin()->second.second=5000;

  for (int i=start_num;i!=end_num+1;i++){
 
    sds.jump(i);
    WireCell::Slice slice = sds.get();
    
    lowmemtiling[i] = new WireCell2dToy::LowmemTiling(i,nrebin,slice,gds,uplane_rms,vplane_rms,wplane_rms);
    if (i==start_num){
      lowmemtiling[i]->init_bad_cells(uplane_map,vplane_map,wplane_map);
    }else{
      lowmemtiling[i]->check_bad_cells(lowmemtiling[i-1],uplane_map,vplane_map,wplane_map);
    }


    //toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds,0.15,0.2,0.1,threshold_ug,threshold_vg, threshold_wg, &uplane_rms, &vplane_rms, &wplane_rms);

    //    if (two_plane)
    //  toytiling[i]->twoplane_tiling(i,nrebin,gds,uplane_rms,vplane_rms,wplane_rms, uplane_map, vplane_map, wplane_map);


    //    GeomCellSelection allcell = toytiling[i]->get_allcell();
    //GeomWireSelection allwire = toytiling[i]->get_allwire();

    //    cout << i << " " << allcell.size() << " " << allwire.size() << endl;

    //    mergetiling[i] = new WireCell2dToy::MergeToyTiling(*toytiling[i],i,3);
    
   
    
    //    if (i==0){
    // badtiling[i] = new WireCell2dToy::BadTiling(i,nrebin,uplane_map,vplane_map,wplane_map,gds,0,1); // 2 plane bad tiling
      // badtiling[i] = new WireCell2dToy::BadTiling(i,nrebin,uplane_map,vplane_map,wplane_map,gds,1,1); // 1 plane bad tiling
    //}

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
    // // // display.draw_bad_region(uplane_map,i,nrebin,0,"same");
    // // // display.draw_bad_region(vplane_map,i,nrebin,1,"same");
    // // // display.draw_bad_region(wplane_map,i,nrebin,2,"same");
    // // display.draw_bad_cell(badtiling[i]->get_cell_all());
  
    // display.draw_cells(toytiling[i]->get_allcell(),"*same");
    // display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell
    
    // // display.draw_wires_charge(toytiling[i]->wcmap(),"Fsame",FI);
    // // display.draw_cells_charge(toytiling[i]->get_allcell(),"Fsame");
    //  theApp.Run();
  }
  
  cerr << em("finish tiling") << endl;


  TFile *file = new TFile(Form("result_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
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
