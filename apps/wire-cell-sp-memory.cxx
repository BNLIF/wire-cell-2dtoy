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
  if (argc < 4) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num -t[0,1] -s[0,1,2]" << endl;
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
  cerr << em("...done") << endl;
  
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
  
  
 
  float unit_dis = 1.101; // match 256 cm

   float toffset_1=0.0; //(nt_off1 * 0.2 - 1.0 );  // time offset between u/v 
  float toffset_2=0.0; //(nt_off2 * 0.2 - 1.0); // time offset between u/w
  float toffset_3=0.0;
  

  
  int save_image_outline_flag = 0; // prescale flag 
  

  int total_time_bin=9594;
  int recon_threshold = 2000;
  int frame_length = 3200;
  int max_events = 100;
  int eve_num  = atoi(argv[3]);
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
 
  
  cerr << em("load data frame and noise filtering") << endl;
  WireCellSst::DatauBooNEFrameDataSource *data_fds = new WireCellSst::DatauBooNEFrameDataSource(root_file,gds,total_time_bin);
  
  cerr << em("prejump") << endl;
  delete data_fds;
  
  cerr << em("... done") << endl;


  // if (save_file != 2){
  //   data_fds.jump(eve_num);
  //   if (save_file == 1)
  //     data_fds.Save();
  // }
  //cerr << em("... done") << endl;
  
  // run_no = data_fds.get_run_no();
  // subrun_no = data_fds.get_subrun_no();
  // event_no = data_fds.get_event_no();
  
  // cout << "Run No: " << run_no << " " << subrun_no << " " << event_no << endl;

  // ChirpMap& uplane_map = data_fds.get_u_cmap();
  // ChirpMap& vplane_map = data_fds.get_v_cmap();
  // ChirpMap& wplane_map = data_fds.get_w_cmap();
  
  // std::set<int>& lf_noisy_channels = data_fds.get_lf_noisy_channels();
  


  // cout << "Deconvolution with Wiener filter" << endl; 

  // cerr << em("do TPC signal processing") << endl;
  // WireCell2dToy::uBooNEData2DDeconvolutionFDS wien_fds(data_fds,gds,uplane_map, vplane_map, wplane_map,100,toffset_1,toffset_2,toffset_3);
  // wien_fds.jump(eve_num);
  // WireCell2dToy::uBooNEDataROI uboone_rois(data_fds,wien_fds,gds,uplane_map,vplane_map,wplane_map,lf_noisy_channels);
  // WireCell2dToy::uBooNEDataAfterROI roi_fds(wien_fds,gds,uboone_rois,nrebin);
  // roi_fds.jump(eve_num);
  // cerr << em("... done") << endl;
  
  // cerr << em("Clear data inside data_fds") << endl;
  // data_fds.Clear();
  // cerr << em("... done") << endl;


  // std::vector<float>& uplane_rms = uboone_rois.get_uplane_rms();
  // std::vector<float>& vplane_rms = uboone_rois.get_vplane_rms();
  // std::vector<float>& wplane_rms = uboone_rois.get_wplane_rms();


  // GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  // GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  // GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  // int nwire_u = wires_u.size();
  // int nwire_v = wires_v.size();
  // int nwire_w = wires_w.size();
  
  
  // cerr << em("creat time slice ") << endl;
  // WireCellSst::ToyuBooNESliceDataSource sds(roi_fds,roi_fds,threshold_u, 
  // 					    threshold_v, threshold_w, 
  // 					    threshold_ug, 
  // 					    threshold_vg, threshold_wg, 
  // 					    nwire_u, 
  // 					    nwire_v, nwire_w,
  // 					    &uplane_rms, &vplane_rms, &wplane_rms); 
    
  // cerr << em("... done") << endl;
  


  

  return 0;
  
} // main()
