#include "WCPSst/GeomDataSource.h"
#include "WCPSst/DatauBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/BadTiling.h"

#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/TruthToyTiling.h"
#include "WCP2dToy/SimpleBlobToyTiling.h"

#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/ToyMatrixExclusive.h"
#include "WCP2dToy/ToyMatrixKalman.h"
#include "WCP2dToy/ToyMatrixIterate.h"
#include "WCP2dToy/ToyMatrixIterate_SingleWire.h"
#include "WCP2dToy/ToyMatrixIterate_Only.h"


#include "WCP2dToy/ToyMatrixMarkov.h"
#include "WCP2dToy/ToyMetric.h"
#include "WCP2dToy/BlobMetric.h"

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

//#include "WCP2dToy/DataSignalGaus.h"
//#include "WCP2dToy/DataSignalWien_ROI.h"

#include "WCP2dToy/uBooNE_Data_2D_Deconvolution.h"
#include "WCP2dToy/uBooNE_Data_ROI.h"
#include "WCP2dToy/uBooNE_Data_After_ROI.h"

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


int main(int argc, char* argv[])
{
  if (argc < 4) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt decon.root raw.root" << endl;
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

  TString filename = argv[2];
  TFile *file = new TFile(filename);
  TTree *Trun = (TTree*)file->Get("Trun");
  TTree *T_bad = (TTree*)file->Get("T_bad");
  TTree *T_lf = (TTree*)file->Get("T_lf");

  TH2I *hu_decon = (TH2I*)file->Get("hu_decon");
  TH2I *hv_decon = (TH2I*)file->Get("hv_decon");
  TH2I *hw_decon = (TH2I*)file->Get("hw_decon");
  
  filename = argv[3];
  TFile *file2 = new TFile(filename);
  TH2F *hu_raw = (TH2F*)file2->Get("hu_raw");
  TH2F *hv_raw = (TH2F*)file2->Get("hv_raw");
  TH2F *hw_raw = (TH2F*)file2->Get("hw_raw");
 
  WCPSst::DatauBooNEFrameDataSource data_fds(hu_raw,hv_raw,hw_raw,T_bad,T_lf,Trun,gds);

  WCP2dToy::uBooNEData2DDeconvolutionFDS wien_fds(hu_decon,hv_decon,hw_decon,T_bad, gds);
  ChirpMap& uplane_map = wien_fds.get_u_cmap();
  ChirpMap& vplane_map = wien_fds.get_v_cmap();
  ChirpMap& wplane_map = wien_fds.get_w_cmap();
  std::set<int>& lf_noisy_channels = data_fds.get_lf_noisy_channels();

  int run_no;
  int subrun_no;
  int event_no;
 
  Trun->SetBranchAddress("eventNo",&event_no);
  Trun->SetBranchAddress("runNo",&run_no);
  Trun->SetBranchAddress("subRunNo",&subrun_no);
  Trun->GetEntry(0);


  //std::cout << uplane_map.size() << " " << vplane_map.size() << " " << wplane_map.size() <<  " " << uplane_map[880].first << " " << uplane_map[880].second << std::endl;
  int rebin = 4;
  WCP2dToy::uBooNEDataROI uboone_rois(data_fds,wien_fds,gds,uplane_map,vplane_map,wplane_map,lf_noisy_channels);
  WCP2dToy::uBooNEDataAfterROI roi_fds(wien_fds,gds,uboone_rois,rebin);
  roi_fds.jump(0);
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));
  

  Int_t nwire_u = wires_u.size();
  Int_t nwire_v = wires_v.size();
  Int_t nwire_w = wires_w.size();

  int  total_time_bin = hu_decon->GetNbinsY();

  TFile *file1 = new TFile(Form("nsp3_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");

  TH2F *hu_roi = new TH2F("hu_roi","hu_roi",nwire_u,-0.5,nwire_u-0.5,total_time_bin/rebin,0,total_time_bin);
  TH2F *hv_roi = new TH2F("hv_roi","hv_roi",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin/rebin,0,total_time_bin);
  TH2F *hw_roi = new TH2F("hw_roi","hw_roi",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin/rebin,0,total_time_bin);
  TH2F *htemp1;
  
  const Frame& frame1 = roi_fds.get();
  int ntraces = frame1.traces.size();
  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    WirePlaneType_t plane = gds.by_channel(chid).at(0)->plane();
    if (plane == WirePlaneType_t(0)){
      htemp1 = hu_roi;
    }else if (plane == WirePlaneType_t(1)){
      htemp1 = hv_roi;
      chid -= nwire_u;
    }else if (plane == WirePlaneType_t(2)){
      htemp1 = hw_roi;
      chid -= nwire_u + nwire_v;
    }
     for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp1->SetBinContent(chid+1,tt,trace.charge.at(i));
    }
  }
  
  file1->Write();
  file1->Close();
  


  // std::vector<std::pair<int,int>>& rois = uboone_rois.get_self_rois(500);
  // std::cout << rois.size() << std::endl;
  // for (int i=0;i!=rois.size();i++){
  //   std::cout << i << " S " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  // }
  // rois = uboone_rois.get_others_rois(500);
  // std::cout << rois.size() << std::endl;
  // for (int i=0;i!=rois.size();i++){
  //   std::cout << i << " O " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  // }
  // rois = uboone_rois.get_combined_rois(500);
  // std::cout << rois.size() << std::endl;
  // for (int i=0;i!=rois.size();i++){
  //   std::cout << i << " C " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  // }

  // rois = uboone_rois.get_self_rois(3500);
  // std::cout << rois.size() << std::endl;
  // for (int i=0;i!=rois.size();i++){
  //   std::cout << i << " S " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  // }
  // rois = uboone_rois.get_others_rois(3500);
  // std::cout << rois.size() << std::endl;
  // for (int i=0;i!=rois.size();i++){
  //   std::cout << i << " O " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  // }
  // rois = uboone_rois.get_combined_rois(3500);
  // std::cout << rois.size() << std::endl;
  // for (int i=0;i!=rois.size();i++){
  //   std::cout << i << " C " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  // }


  // rois = uboone_rois.get_combined_rois(7000);
  // std::cout << rois.size() << std::endl;
  // for (int i=0;i!=rois.size();i++){
  //   std::cout << i << " C " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  // }


  
  
}
