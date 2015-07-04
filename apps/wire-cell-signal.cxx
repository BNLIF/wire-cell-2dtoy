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
#include "WireCell2dToy/ToySignalGaus.h"
#include "WireCell2dToy/ToySignalWien.h"

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
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root" << endl;
    return 1;
  }
  
  
  WireCellSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;
  
  
  const char* root_file = argv[2];
  const char* tpath = "/Event/Sim";
  
  WireCell::FrameDataSource* fds = 0;
  fds = WireCellSst::make_fds(root_file);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  

  // TFile tfile(root_file,"read");
  // TTree* sst = dynamic_cast<TTree*>(tfile.Get(tpath));
  // WireCellSst::ToyuBooNEFrameDataSource fds_data(*sst,gds);
  // fds_data.jump(1);
  // fds_data.Save();
  // WireCell2dToy::ToySignalPreFDS pre_fds(fds_data,gds,9600/4,5);
  // pre_fds.jump(1);
  // pre_fds.Save();
  


  WireCell::ToyDepositor toydep(fds);
  const PointValueVector pvv = toydep.depositions(1);
  WireCell::GenerativeFDS gfds(toydep,gds,9600,5,0.5*1.605723*units::millimeter); // 87 K at 0.5 kV/cm
  WireCell2dToy::ToySignalSimuFDS simu_fds(gfds,gds,9600,5);
  simu_fds.jump(1);
  simu_fds.Save();
  
  WireCell2dToy::ToySignalGausFDS gaus_fds(simu_fds,gds,9600/4,5);
  gaus_fds.jump(1);
  gaus_fds.Save();
  
  WireCell2dToy::ToySignalWienFDS wien_fds(simu_fds,gds,9600/4,5);
  wien_fds.jump(1);
  wien_fds.Save();
  
  

  // WireCellSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons
  // int start_num =180;
  // int end_num = 181;
  // for (int i=start_num;i!=end_num+1;i++){
  //   sds.jump(i);
  //   WireCell::Slice slice = sds.get();
  // }

 
  return 0;
  
} // main()
