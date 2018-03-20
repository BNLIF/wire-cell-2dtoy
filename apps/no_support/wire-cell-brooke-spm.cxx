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
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"

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

//#define MAX_TRACKS 30000


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
  float unit_dis = 1.6;  // 
  float toffset_1=1.647;
  float toffset_2=1.539+1.647;
  //float toffset_1 = 0;
  //float toffset_2 = 0;
  float toffset_3=0;

  int save_image_outline_flag = 1; // prescale flag 
  

  int total_time_bin=9600;
  int frame_length = 3200;
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
  std::vector<std::vector<int> > *mc_daughters= new std::vector<std::vector<int> >;  // daughters id of this track; vector
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


  WireCell::FrameDataSource* fds = 0;
  fds = WireCellSst::make_fds(*tfile);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  
  ChirpMap uplane_map;
  ChirpMap vplane_map;
  ChirpMap wplane_map;


  // WireCell::ToyDepositor toydep(fds);
  // const PointValueVector pvv = toydep.depositions(eve_num);
  
  // WireCell::GenerativeFDS gfds(toydep,gds,2400,max_events,2.0*unit_dis*units::millimeter);
  // gfds.jump(eve_num);


  
 

  WireCell::ToyDepositor toydep(fds,0,unit_dis,frame_length);
  const PointValueVector& pvv = toydep.depositions(eve_num);
  
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
  
 
   TFile *file = new TFile(Form("2D_display_%d_%d_%d.root",run_no,subrun_no,eve_num),"RECREATE");


  TH2I *hu_orig = new TH2I("hu_orig","hu_orig",nwire_u,-0.5,nwire_u-0.5,total_time_bin/nrebin,0,total_time_bin);
  TH2I *hv_orig = new TH2I("hv_orig","hv_orig",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin/nrebin,0,total_time_bin);
  TH2I *hw_orig = new TH2I("hw_orig","hw_orig",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin/nrebin,0,total_time_bin);

  TH1I *hu_baseline = new TH1I("hu_baseline","hu_threshold",nwire_u,-0.5,-0.5+nwire_u);
  TH1I *hv_baseline = new TH1I("hv_baseline","hv_threshold",nwire_v,-0.5+nwire_u,-0.5+nwire_u+nwire_v);
  TH1I *hw_baseline = new TH1I("hw_baseline","hw_threshold",nwire_w,-0.5+nwire_u+nwire_v,-0.5+nwire_u+nwire_v+nwire_w);

  
  TH1I *hu_threshold = new TH1I("hu_threshold","hu_basline",nwire_u,-0.5,-0.5+nwire_u);
  TH1I *hv_threshold = new TH1I("hv_threshold","hv_basline",nwire_v,-0.5+nwire_u,-0.5+nwire_u+nwire_v);
  TH1I *hw_threshold = new TH1I("hw_threshold","hw_basline",nwire_w,-0.5+nwire_u+nwire_v,-0.5+nwire_u+nwire_v+nwire_w);

  TH2F *hu_raw = new TH2F("hu_raw","hu_raw",nwire_u,-0.5,nwire_u-0.5,total_time_bin,0,total_time_bin);
  TH2F *hv_raw = new TH2F("hv_raw","hv_raw",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin,0,total_time_bin);
  TH2F *hw_raw = new TH2F("hw_raw","hw_raw",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin,0,total_time_bin);

  TH2F *hu_decon = new TH2F("hu_decon","hu_decon",nwire_u,-0.5,nwire_u-0.5,total_time_bin/nrebin,0,total_time_bin);
  TH2F *hv_decon = new TH2F("hv_decon","hv_decon",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin/nrebin,0,total_time_bin);
  TH2F *hw_decon = new TH2F("hw_decon","hw_decon",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin/nrebin,0,total_time_bin);
  

  TH2F *htemp;
  TH2F *htemp1;
  
  const Frame& frame = simu_fds.get();
  size_t ntraces = frame.traces.size();
  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    WirePlaneType_t plane = gds.by_channel(chid).at(0)->plane();
    if (plane == WirePlaneType_t(0)){
      htemp = hu_raw;
      
    }else if (plane == WirePlaneType_t(1)){
      htemp = hv_raw;
      
      chid -= nwire_u;
    }else if (plane == WirePlaneType_t(2)){
      htemp = hw_raw;
     
      chid -= nwire_u + nwire_v;
    }
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp->SetBinContent(chid+1,tt,trace.charge.at(i));
    }
  }
  
  const Frame& frame1 = gaus_fds.get();
  ntraces = frame1.traces.size();
  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    WirePlaneType_t plane = gds.by_channel(chid).at(0)->plane();
    if (plane == WirePlaneType_t(0)){
      htemp1 = hu_decon;
    }else if (plane == WirePlaneType_t(1)){
      htemp1 = hv_decon;
      chid -= nwire_u;
    }else if (plane == WirePlaneType_t(2)){
      htemp1 = hw_decon;
      chid -= nwire_u + nwire_v;
    }
     for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp1->SetBinContent(chid+1,tt,trace.charge.at(i));
    }
  }

  const Frame& frame2 = st_fds.get();
  size_t ntraces2 = frame2.traces.size();
  TH2I *htemp2;
  for (size_t ind=0; ind<ntraces2; ++ind) {
    const Trace& trace = frame2.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    WirePlaneType_t plane = gds.by_channel(chid).at(0)->plane();
    if (plane == WirePlaneType_t(0)){
      htemp2 = hu_orig;
      
    }else if (plane == WirePlaneType_t(1)){
      htemp2 = hv_orig;
      
      chid -= nwire_u;
    }else if (plane == WirePlaneType_t(2)){
      htemp2 = hw_orig;
      chid -= nwire_u + nwire_v;
    }
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp2->SetBinContent(chid+1,tt,trace.charge.at(i)/500.);
    }
  }

  
  // save original data ... 
  //const char* tpath = "/Event/Sim";
  //TFile tfile(root_file,"read");
  //TTree* tree = dynamic_cast<TTree*>(tfile.Get(tpath));
  // tree->SetBranchStatus("*",0);
  // std::vector<int> *channelid = new std::vector<int>;
  // TClonesArray* esignal = new TClonesArray;
  
  // tree->SetBranchStatus("raw_channelId",1);
  // tree->SetBranchAddress("raw_channelId", &channelid);
  // tree->SetBranchStatus("raw_wf",1);
  // tree->SetBranchAddress("raw_wf", &esignal);
  // tree->GetEntry(eve_num);
  
  

  std::vector<float>& uplane_rms = wien_fds.get_uplane_rms();
  std::vector<float>& vplane_rms = wien_fds.get_vplane_rms();
  std::vector<float>& wplane_rms = wien_fds.get_wplane_rms();
  for (Int_t i=0;i!=uplane_rms.size();i++){
    hu_threshold->SetBinContent(i+1,uplane_rms.at(i)*3.6*nrebin );
  }
  for (Int_t i=0;i!=vplane_rms.size();i++){
    hv_threshold->SetBinContent(i+1,vplane_rms.at(i)*3.6*nrebin );
  }
  for (Int_t i=0;i!=wplane_rms.size();i++){
    hw_threshold->SetBinContent(i+1,wplane_rms.at(i)*3.6*nrebin );
  }

  // finish saving

  TTree *Trun = new TTree("Trun","Trun");
  Trun->SetDirectory(file);

  int detector = 0; // MicroBooNE
  Trun->Branch("detector",&detector,"detector/I");

  Trun->Branch("eventNo",&event_no,"eventNo/I");
  Trun->Branch("runNo",&run_no,"runNo/I");
  Trun->Branch("subRunNo",&subrun_no,"subRunNo/I");
  
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

  TTree *T_lf = new TTree("T_lf","T_lf");
  Int_t channel;
  T_lf->SetDirectory(file);
  T_lf->Branch("channel",&channel,"channel/I");
  std::set<int> lf_noisy_channels;// = data_fds.get_lf_noisy_channels();
  for (auto it = lf_noisy_channels.begin(); it!= lf_noisy_channels.end(); it++){
    channel = *it;
    T_lf->Fill();
  }


  TTree *T_bad = new TTree("T_bad","T_bad");
  Int_t chid, plane;
  Int_t start_time,end_time;
  T_bad->Branch("chid",&chid,"chid/I");
  T_bad->Branch("plane",&plane,"plane/I");
  T_bad->Branch("start_time",&start_time,"start_time/I");
  T_bad->Branch("end_time",&end_time,"end_time/I");
  T_bad->SetDirectory(file);
  for (auto it = uplane_map.begin(); it!=uplane_map.end();it++){
    chid = it->first;
    plane = 0;
    start_time = it->second.first;
    end_time = it->second.second;
    T_bad->Fill();
  }
  for (auto it = vplane_map.begin(); it!=vplane_map.end();it++){
    chid = it->first + nwire_u;
    plane = 1;
    start_time = it->second.first;
    end_time = it->second.second;
    T_bad->Fill();
  }
  for (auto it = wplane_map.begin(); it!=wplane_map.end();it++){
    chid = it->first + nwire_u + nwire_v;
    plane = 2;
    start_time = it->second.first;
    end_time = it->second.second;
    T_bad->Fill();
  }
  



  file->Write();
  file->Close();
} // main()
