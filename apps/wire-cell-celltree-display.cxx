#include "WireCellSst/GeomDataSource.h"
#include "WireCellSst/DatauBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/BadTiling.h"

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

#include "WireCell2dToy/uBooNE_Data_2D_Deconvolution.h"
#include "WireCell2dToy/uBooNE_Data_ROI.h"
#include "WireCell2dToy/uBooNE_Data_After_ROI.h"
#include "WireCell2dToy/uBooNE_Data_After_ROI_gaus.h"

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
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num " << endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);
  
  WireCellSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  
  float unit_dis = 1.119;  // 70 KV @ 273 V/cm
  //final offset after time scan (70kV)
  float toffset_1=0.0; //-0.787;
  float toffset_2=0.0;//-0.603;
  float toffset_3=0.0;
  int total_time_bin=9592;
  int frame_length = 3200;
  int nrebin = 4; // 6 is default
  int eve_num  = atoi(argv[3]);
  int time_offset = -92.;
  
  const char* root_file = argv[2];
  int run_no, subrun_no, event_no;
  const char* tpath = "/Event/Sim";
  TFile tfile(root_file,"read");
  TTree* tree = dynamic_cast<TTree*>(tfile.Get(tpath));
  tree->SetBranchAddress("eventNo" , &event_no);
  tree->SetBranchAddress("runNo"   , &run_no);
  tree->SetBranchAddress("subRunNo", &subrun_no);

  Int_t nf_nChannel;
  std::vector<int> *nf_channelid = new std::vector<int>;
  TClonesArray* nf_esignal = new TClonesArray;
  tree->SetBranchAddress("nf_nChannel",&nf_nChannel);
  tree->SetBranchAddress("nf_channelId", &nf_channelid);
  tree->SetBranchAddress("nf_wf", &nf_esignal);

  Short_t nf_shift, nf_scale;
  tree->SetBranchAddress("nf_shift",&nf_shift);
  tree->SetBranchAddress("nf_scale",&nf_scale);

  Int_t calibWiener_nChannel;
  std::vector<int> *calibWiener_channelid = new std::vector<int>;
  TClonesArray* calibWiener_esignal = new TClonesArray;
  tree->SetBranchAddress("calibWiener_nChannel",&calibWiener_nChannel);
  tree->SetBranchAddress("calibWiener_channelId", &calibWiener_channelid);
  tree->SetBranchAddress("calibWiener_wf", &calibWiener_esignal);


  
  std::vector<double> *channelThreshold = new std::vector<double>;
  tree->SetBranchAddress("channelThreshold",&channelThreshold);

  std::vector<int> *badChannel = new std::vector<int>;
  std::vector<int> *badBegin = new std::vector<int>;
  std::vector<int> *badEnd = new std::vector<int>;
  tree->SetBranchAddress("badChannel",&badChannel);
  tree->SetBranchAddress("badBegin",&badBegin);
  tree->SetBranchAddress("badEnd",&badEnd);
  
  tree->GetEntry(eve_num);
 
  cout << "Run No: " << run_no << " " << subrun_no << " " << event_no << endl;

  // std::cout << nf_nChannel << std::endl;
  
  TFile *file = new TFile(Form("nsp_2D_display_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");

  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  Int_t nwire_u = wires_u.size();
  Int_t nwire_v = wires_v.size();
  Int_t nwire_w = wires_w.size();

  TH1I *hu_threshold = new TH1I("hu_threshold","hu_threshold",nwire_u,-0.5,-0.5+nwire_u);
  TH1I *hv_threshold = new TH1I("hv_threshold","hv_threshold",nwire_v,-0.5+nwire_u,-0.5+nwire_u+nwire_v);
  TH1I *hw_threshold = new TH1I("hw_threshold","hw_threshold",nwire_w,-0.5+nwire_u+nwire_v,-0.5+nwire_u+nwire_v+nwire_w);

  hu_threshold->SetDirectory(file);
  hv_threshold->SetDirectory(file);
  hw_threshold->SetDirectory(file);
  for (size_t i=0; i!= channelThreshold->size();i++){
    if (i<nwire_u){
      hu_threshold->SetBinContent(i+1,int(channelThreshold->at(i)*nrebin));
    }else if (i< nwire_u+nwire_v){
      hv_threshold->SetBinContent(i+1-nwire_u,int(channelThreshold->at(i)*nrebin));
    }else{
      hw_threshold->SetBinContent(i+1-nwire_u-nwire_v,int(channelThreshold->at(i)*nrebin));
    }
  }
  //std::cout << channelThreshold->size() << std::endl;


  

  //std::cout << signal->GetBinContent(1) << " " << signal->GetNbinsX() << std::endl;
  TH2F *hu_raw = new TH2F("hu_raw","hu_raw",nwire_u,-0.5,nwire_u-0.5,total_time_bin,0,total_time_bin);
  TH2F *hv_raw = new TH2F("hv_raw","hv_raw",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin,0,total_time_bin);
  TH2F *hw_raw = new TH2F("hw_raw","hw_raw",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin,0,total_time_bin);

  hu_raw->SetDirectory(file);
  hv_raw->SetDirectory(file);
  hw_raw->SetDirectory(file);
  
  for (size_t i=0; i!=nf_channelid->size(); i++){
    int chid = nf_channelid->at(i);
    TH1S* signal = dynamic_cast<TH1S*>(nf_esignal->At(i));
    if (chid < nwire_u){
      for (int j=0;j!=total_time_bin;j++){
	hu_raw->SetBinContent(chid+1,j+1, (signal->GetBinContent(j+1)+nf_shift)*1.0/nf_scale);
      }
    }else if (chid < nwire_v + nwire_u){
      for (int j=0;j!=total_time_bin;j++){
	hv_raw->SetBinContent(chid+1-nwire_u,j+1, (signal->GetBinContent(j+1)+nf_shift)*1.0/nf_scale);
      }
    }else{
      for (int j=0;j!=total_time_bin;j++){
	hw_raw->SetBinContent(chid+1-nwire_u-nwire_v,j+1, (signal->GetBinContent(j+1)+nf_shift)*1.0/nf_scale);
      }
    }
  
  }
 
  
  TH2F *hu_decon = new TH2F("hu_decon","hu_decon",nwire_u,-0.5,nwire_u-0.5,total_time_bin/nrebin,0,total_time_bin);
  TH2F *hv_decon = new TH2F("hv_decon","hv_decon",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin/nrebin,0,total_time_bin);
  TH2F *hw_decon = new TH2F("hw_decon","hw_decon",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin/nrebin,0,total_time_bin);
  hu_decon->SetDirectory(file);
  hv_decon->SetDirectory(file);
  hw_decon->SetDirectory(file);

  {
    TH1F* signal = dynamic_cast<TH1F*>(calibWiener_esignal->At(0));
    std::cout << signal->GetNbinsX() << " " << signal->GetBinContent(2398) << std::endl;
  }
  
  for (size_t i=0; i!=calibWiener_channelid->size(); i++){
    int chid = calibWiener_channelid->at(i);
    TH1F* signal = dynamic_cast<TH1F*>(calibWiener_esignal->At(i));
    if (chid < nwire_u){
      for (int j=0;j!=total_time_bin/nrebin;j++){
    	hu_decon->SetBinContent(chid+1,j+1, signal->GetBinContent(j+1));
      }
    }else if (chid < nwire_v + nwire_u){
      for (int j=0;j!=total_time_bin/nrebin;j++){
    	hv_decon->SetBinContent(chid+1-nwire_u,j+1, signal->GetBinContent(j+1));
      }
    }else{
      for (int j=0;j!=total_time_bin/nrebin;j++){
    	hw_decon->SetBinContent(chid+1-nwire_u-nwire_v,j+1, signal->GetBinContent(j+1));
      }
    }
  }
  
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
  //Trun->Branch("recon_threshold",&recon_threshold,"recon_threshold/I");
  Trun->Branch("frame_length",&frame_length,"frame_length/I");
  //Trun->Branch("max_events",&max_events,"max_events/I");
  Trun->Branch("eve_num",&eve_num,"eve_num/I");
  Trun->Branch("nrebin",&nrebin,"nrebin/I");
  //Trun->Branch("threshold_u",&threshold_u,"threshold_u/F");
  //Trun->Branch("threshold_v",&threshold_v,"threshold_v/F");
  //Trun->Branch("threshold_w",&threshold_w,"threshold_w/F");
  Trun->Branch("time_offset",&time_offset,"time_offset/I");
  Trun->Fill();

  TTree *T_lf = new TTree("T_lf","T_lf");
  Int_t channel;
  T_lf->SetDirectory(file);
  T_lf->Branch("channel",&channel,"channel/I");
  

  TTree *T_bad = new TTree("T_bad","T_bad");
  Int_t chid, plane;
  Int_t start_time,end_time;
  T_bad->Branch("chid",&chid,"chid/I");
  T_bad->Branch("plane",&plane,"plane/I");
  T_bad->Branch("start_time",&start_time,"start_time/I");
  T_bad->Branch("end_time",&end_time,"end_time/I");
  T_bad->SetDirectory(file);

  for (size_t i=0; i!=badChannel->size();i++){
    chid = badChannel->at(i);
    if (chid < nwire_u){
      plane = 0;
    }else if (chid < nwire_u + nwire_v){
      plane = 1;
    }else{
      plane = 2;
    }
    start_time = badBegin->at(i);
    end_time = badEnd->at(i);
    T_bad->Fill();
  }
  


  file->Write();
  file->Close();
  
  
}
