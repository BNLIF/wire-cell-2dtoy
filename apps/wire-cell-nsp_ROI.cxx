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
#include "WireCell2dToy/uBooNE_Data_ROI.h"

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

  TH2I *hu_decon = (TH2I*)file->Get("hu_decon");
  TH2I *hv_decon = (TH2I*)file->Get("hv_decon");
  TH2I *hw_decon = (TH2I*)file->Get("hw_decon");
  
 
  WireCell2dToy::uBooNEData2DDeconvolutionFDS wien_fds(hu_decon,hv_decon,hw_decon,T_bad, gds);

  ChirpMap& uplane_map = wien_fds.get_u_cmap();
  ChirpMap& vplane_map = wien_fds.get_v_cmap();
  ChirpMap& wplane_map = wien_fds.get_w_cmap();
  
  // std::cout << uplane_map.size() << " " << vplane_map.size() << " " << wplane_map.size() <<  " " << wplane_map[7500].first << " " << wplane_map[7500].second << std::endl;

  WireCell2dToy::uBooNEDataROI uboone_rois(wien_fds,gds,uplane_map,vplane_map,wplane_map);

  std::vector<std::pair<int,int>>& rois = uboone_rois.get_self_rois(500);
  std::cout << rois.size() << std::endl;
  for (int i=0;i!=rois.size();i++){
    std::cout << i << " S " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  }
  rois = uboone_rois.get_others_rois(500);
  std::cout << rois.size() << std::endl;
  for (int i=0;i!=rois.size();i++){
    std::cout << i << " O " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  }


  rois = uboone_rois.get_self_rois(3500);
  std::cout << rois.size() << std::endl;
  for (int i=0;i!=rois.size();i++){
    std::cout << i << " S " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  }
  rois = uboone_rois.get_others_rois(3500);
  std::cout << rois.size() << std::endl;
  for (int i=0;i!=rois.size();i++){
    std::cout << i << " O " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  }

  rois = uboone_rois.get_self_rois(7000);
  std::cout << rois.size() << std::endl;
  for (int i=0;i!=rois.size();i++){
    std::cout << i << " S " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  }


  
  
}
