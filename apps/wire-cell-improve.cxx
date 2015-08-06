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
  
  TFile *tfile = TFile::Open(root_file);

  WireCell::FrameDataSource* fds = 0;
  fds = WireCellSst::make_fds(*tfile);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  
  TH1::AddDirectory(kFALSE);
  

  int recon_threshold = 2000;
  int max_events = 100;
  int eve_num = atoi(argv[3]);
  
  WireCell::ToyDepositor *toydep = new WireCell::ToyDepositor(fds);
  const PointValueVector& pvv = toydep->depositions(eve_num);
  
  WireCell::GenerativeFDS *gfds = new WireCell::GenerativeFDS(*toydep,gds,9600,max_events,0.5*1.60*units::millimeter); // 87 K at 0.5 kV/cm
  gfds->jump(eve_num);

  
  int abc;
  
  tfile->Close("R");
  delete tfile;

  
  cout << "Put in Truth " << endl; 
  WireCell2dToy::ToySignalSimuTrueFDS *st_fds = new WireCell2dToy::ToySignalSimuTrueFDS(*gfds,gds,9600/4,max_events,0); //truth
  st_fds->jump(eve_num);
  
  cout << "Simulate Raw WaveForm " << endl; 
  WireCell2dToy::ToySignalSimuFDS *simu_fds = new WireCell2dToy::ToySignalSimuFDS(*gfds,gds,9600,max_events,1.647,1.539+1.647,1); // time offset among different planes for the time electrons travel among different planes
  simu_fds->jump(eve_num);
  
  cout << "Deconvolution with Gaussian filter" << endl;
  WireCell2dToy::ToySignalGausFDS *gaus_fds = new WireCell2dToy::ToySignalGausFDS(*simu_fds,gds,9600/4,max_events,1.647,1.539+1.647); // gaussian smearing for charge estimation
  gaus_fds->jump(eve_num);
  
  cout << "Deconvolution with Wiener filter" << endl;
   WireCell2dToy::ToySignalWienFDS *wien_fds = new WireCell2dToy::ToySignalWienFDS(*simu_fds,gds,9600/4,max_events,1.647,1.539+1.647); // weiner smearing for hit identification
  wien_fds->jump(eve_num);
  
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  int nwire_u = wires_u.size();
  int nwire_v = wires_v.size();
  int nwire_w = wires_w.size();
  
  float threshold_u = 5.87819e+02 * 4.0;
  float threshold_v = 8.36644e+02 * 4.0;
  float threshold_w = 5.67974e+02 * 4.0;

  float threshold_ug = 755.96;
  float threshold_vg = 822.81;
  float threshold_wg = 510.84;
  
  // cin >> abc;

  WireCellSst::ToyuBooNESliceDataSource *sds = new WireCellSst::ToyuBooNESliceDataSource(*wien_fds,*gaus_fds,threshold_u, 
  					    threshold_v, threshold_w, 
  					    threshold_ug, 
  					    threshold_vg, threshold_wg, 
  					    nwire_u, 
  					    nwire_v, nwire_w); 

  WireCellSst::ToyuBooNESliceDataSource *sds_th = new WireCellSst::ToyuBooNESliceDataSource(*st_fds,*st_fds,500, 
  					    500, 500, 
  					    threshold_ug, 
  					    threshold_vg, threshold_wg, 
  					    nwire_u, 
  					    nwire_v, nwire_w); 

 
  int ncount = 0;
  int ncount1 = 0;  
  int ncount2 = 0;

  int ncount_t = 0;
  
  
  

  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
  WireCell2dToy::MergeToyTiling **mergetiling = new WireCell2dToy::MergeToyTiling*[2400];
  WireCell2dToy::TruthToyTiling **truthtiling = new WireCell2dToy::TruthToyTiling*[2400];
  WireCell2dToy::SimpleBlobToyTiling **blobtiling = new WireCell2dToy::SimpleBlobToyTiling*[2400];
  WireCell2dToy::ToyMatrix **toymatrix = new WireCell2dToy::ToyMatrix*[2400];
  //save truth ...
  WireCell2dToy::ToyTiling **toytiling_th = new WireCell2dToy::ToyTiling*[2400];
  WireCell2dToy::TruthToyTiling **truthtiling_th = new WireCell2dToy::TruthToyTiling*[2400];
  WireCell2dToy::ToyMetric toymetric;
  WireCell2dToy::BlobMetric blobmetric;

  //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  int ncount_mcell = 0;

  
  int start_num = 0 ;
  int end_num = sds->size()-1;

  cout << "Start the Reconstruction " << endl; 

  for (int i=start_num;i!=end_num+1;i++){
    
    sds->jump(i);
    sds_th->jump(i);
    WireCell::Slice slice = sds->get();
    WireCell::Slice slice_th = sds_th->get();
    cout << i << " " << slice.group().size() << " " << slice_th.group().size() << endl;
    
    toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
    GeomCellSelection allcell = toytiling[i]->get_allcell();
    GeomWireSelection allwire = toytiling[i]->get_allwire();
    cout << "Single Cell: " << i << " "  << allcell.size() << " " << allwire.size() << endl;
    mergetiling[i] = new WireCell2dToy::MergeToyTiling(*toytiling[i],i,3,1);
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    GeomWireSelection allmwire = mergetiling[i]->get_allwire();
    cout <<"Blob: " << i << " " << allmcell.size() << " " << allmwire.size() << endl;
    truthtiling[i] = new WireCell2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds,800);
    
    
    toytiling_th[i] = new WireCell2dToy::ToyTiling(slice_th,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
    truthtiling_th[i] = new WireCell2dToy::TruthToyTiling(*toytiling_th[i],pvv,i,gds,800);
    CellChargeMap ccmap = truthtiling[i]->ccmap();
        
  }

  
  delete sds;
  delete sds_th;
  
  delete simu_fds;
  delete gaus_fds;
  delete wien_fds;

  delete st_fds;
  delete gfds;
  delete toydep;
  delete fds;
  
  
  

  cin >> abc;
  
  return 0;
  
} // main()
