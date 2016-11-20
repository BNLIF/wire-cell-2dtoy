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


#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/SimDataSource.h"
#include "WireCellNav/SliceDataSource.h"
#include "WireCellSst/Util.h"
#include "WireCellData/SimTruth.h"
#include "WireCell2dToy/ToyDepositor.h"
#include "WireCellNav/GenerativeFDS.h"
#include "WireCell2dToy/ToySignalSimu.h"
#include "WireCell2dToy/ToySignalSimuTrue.h"

//#include "WireCell2dToy/DataSignalGaus.h"
//#include "WireCell2dToy/DataSignalWien_ROI.h"

#include "WireCell2dToy/uBooNE_Data_2D_Deconvolution.h"

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
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root " << endl;
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

  TString filename = argv[2];
  TFile *file = new TFile(filename);
  TTree *Trun = (TTree*)file->Get("Trun");
  TTree *T_bad = (TTree*)file->Get("T_bad");
  TTree *T_lf = (TTree*)file->Get("T_lf");

  TTree *Trun1 = Trun->CloneTree();
  TTree *T_bad1 = T_bad->CloneTree();
  TTree *T_lf1 = T_lf->CloneTree();

  TH2F *hu_raw = (TH2F*)file->Get("hu_raw");
  TH2F *hv_raw = (TH2F*)file->Get("hv_raw");
  TH2F *hw_raw = (TH2F*)file->Get("hw_raw");

  

 
  

  
  
  WireCellSst::DatauBooNEFrameDataSource data_fds(hu_raw,hv_raw,hw_raw,T_bad,T_lf,Trun,gds);


  //data_fds.jump(0);

  int run_no = data_fds.get_run_no();
  int subrun_no = data_fds.get_subrun_no();
  int event_no = data_fds.get_event_no();
  cout << "Run No: " << run_no << " " << subrun_no << " " << event_no << endl;
  
 
  


  ChirpMap& uplane_map = data_fds.get_u_cmap();
  ChirpMap& vplane_map = data_fds.get_v_cmap();
  ChirpMap& wplane_map = data_fds.get_w_cmap();
  cout << "Bad Channels: " << uplane_map.size() << " " << vplane_map.size() << " " << wplane_map.size() << endl;

  // const Frame& frame1 = data_fds.get();
  // size_t ntraces = frame1.traces.size();
  // cout << ntraces << endl;

  
  cout << "Deconvolution with Wiener filter" << endl; 
  WireCell2dToy::uBooNEData2DDeconvolutionFDS wien_fds(data_fds,gds,uplane_map, vplane_map, wplane_map,100,0.,0.,0.);
  wien_fds.jump(0);
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  Int_t nwire_u = wires_u.size();
  Int_t nwire_v = wires_v.size();
  Int_t nwire_w = wires_w.size();


  
  
  
  TFile *file1 = new TFile(Form("nsp2_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
  Trun1->SetDirectory(file1);
  T_bad1->SetDirectory(file1);
  T_lf1->SetDirectory(file1);

  Trun1->Write();
  T_bad1->Write();
  T_lf1->Write();

  int total_time_bin = wien_fds.Get_Bins_Per_Frame();

  
  TH2I *hu_decon = new TH2I("hu_decon","hu_decon",nwire_u,-0.5,nwire_u-0.5,total_time_bin,0,total_time_bin);
  TH2I *hv_decon = new TH2I("hv_decon","hv_decon",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin,0,total_time_bin);
  TH2I *hw_decon = new TH2I("hw_decon","hw_decon",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin,0,total_time_bin);

  //save histograms ... 
  TH2I *htemp1;
  const Frame& frame1 = wien_fds.get();
  int ntraces = frame1.traces.size();
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

    
  file1->Write();
  file1->Close();


}
