#include "WCPSst/GeomDataSource.h"
#include "WCPSst/ToyuBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/TruthToyTiling.h"
#include "WCP2dToy/SimpleBlobToyTiling.h"

#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/ToyMatrixExclusive.h"
#include "WCP2dToy/ToyMatrixKalman.h"
#include "WCP2dToy/ToyMatrixIterate.h"
#include "WCP2dToy/ToyMatrixMarkov.h"
#include "WCP2dToy/ToyMetric.h"
#include "WCP2dToy/BlobMetric.h"

#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"

#include "WCPData/GeomCluster.h"
//#include "WCPNav/SliceDataSource.h"


#include "WCPNav/FrameDataSource.h"
#include "WCPNav/SimDataSource.h"
#include "WCPNav/SliceDataSource.h"
#include "WCPSst/Util.h"
#include "WCPData/SimTruth.h"
#include "WCP2dToy/ToyDepositor.h"
#include "WCPNav/GenerativeFDS.h"
#include "WCP2dToy/ToySignalSimu.h"
#include "WCP2dToy/ToySignalSimuTrue.h"
#include "WCP2dToy/ToySignalGaus.h"
#include "WCP2dToy/ToySignalWien.h"

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
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root" << endl;
    return 1;
  }
  
  
  WCPSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;
  
  
  const char* root_file = argv[2];
  const char* tpath = "/Event/Sim";
  
  WCP::FrameDataSource* fds = 0;
  fds = WCPSst::make_fds(root_file);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  

  // TFile tfile(root_file,"read");
  // TTree* sst = dynamic_cast<TTree*>(tfile.Get(tpath));
  // WCPSst::ToyuBooNEFrameDataSource fds_data(*sst,gds);
  // fds_data.jump(1);
  // fds_data.Save();
  // WCP2dToy::ToySignalPreFDS pre_fds(fds_data,gds,9600/4,5);
  // pre_fds.jump(1);
  // pre_fds.Save();
  


  WCP::ToyDepositor toydep(fds);
  const PointValueVector pvv = toydep.depositions(1);

  // float sum_charge = 0;
  // for (int i=0;i!=pvv.size();i++){
  //   sum_charge += pvv.at(i).second;
  // }
  // cout << sum_charge << endl;


  //WCP::GenerativeFDS gfds(toydep,gds,9600,5,0.5*1.605723*units::millimeter); // 87 K at 0.5 kV/cm
  WCP::GenerativeFDS gfds(toydep,gds,9600,5,0.5*1.60*units::millimeter); // 87 K at 0.5 kV/cm
  WCP2dToy::ToySignalSimuFDS simu_fds(gfds,gds,9600,5,1.647,1.539+1.647,1); // time offset among different planes for the time electrons travel among different planes
  simu_fds.jump(1);
  //simu_fds.Save();

  WCP2dToy::ToySignalSimuTrueFDS st_fds(gfds,gds,9600/4,5); //truth
  st_fds.jump(1);
  //st_fds.Save();
  
  WCP2dToy::ToySignalGausFDS gaus_fds(simu_fds,gds,9600/4,5,1.647,1.539+1.647); // gaussian smearing for charge estimation
  gaus_fds.jump(1);
  //gaus_fds.Save();
  
  WCP2dToy::ToySignalWienFDS wien_fds(simu_fds,gds,9600/4,5,1.647,1.539+1.647); // weiner smearing for hit identification
  wien_fds.jump(1);
  //wien_fds.Save();
  
  
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

  // float threshold_u = 1000;
  // float threshold_v = 1000;
  // float threshold_w = 1000;
  

  WCPSst::ToyuBooNESliceDataSource sds(wien_fds,gaus_fds,threshold_u, 
					    threshold_v, threshold_w, 
					    threshold_ug, 
					    threshold_vg, threshold_wg, 
					    nwire_u, 
					    nwire_v, nwire_w); 
  
  WCPSst::ToyuBooNESliceDataSource sds_th(st_fds,st_fds,1, 
					    1, 1, 
					    threshold_ug, 
					    threshold_vg, threshold_wg, 
					    nwire_u, 
					    nwire_v, nwire_w); 


  WCP2dToy::ToyTiling **toytiling = new WCP2dToy::ToyTiling*[2400];
  WCP2dToy::MergeToyTiling **mergetiling = new WCP2dToy::MergeToyTiling*[2400];
  WCP2dToy::TruthToyTiling **truthtiling = new WCP2dToy::TruthToyTiling*[2400];

  GeomCellSelection total_cells;
  GeomCellSelection total_recon_cells;
  GeomCellSelection total_corner_cells;
  GeomCellSelection total_blob_cells;
  CellChargeMap total_ccmap;


   
  int start_num =184 + 800;
  int end_num = 186 + 800;
  for (int i=start_num;i!=end_num+1;i++){
    sds.jump(i);
    sds_th.jump(i);
    WCP::Slice slice = sds.get();
    WCP::Slice slice_th = sds_th.get();
    
    toytiling[i] = new WCP2dToy::ToyTiling(slice,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
    mergetiling[i] = new WCP2dToy::MergeToyTiling(*toytiling[i],i);
    
    GeomCellSelection allcell = toytiling[i]->get_allcell();
    GeomWireSelection allwire = toytiling[i]->get_allwire();
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    GeomWireSelection allmwire = mergetiling[i]->get_allwire();
    cout << i << " " << allmcell.size() << " " << allmwire.size() << " " << slice_th.group().size() << endl;
    truthtiling[i] = new WCP2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds,800);
    
    // for (int j=0;j!=allmwire.size();j++){
    //   cout << mergetiling[i]->cells(*allmwire.at(j)).size() << endl;
    // }

        
    for (int j=0;j!=allcell.size();j++){
      total_cells.push_back(allcell.at(j));
    }
    CellChargeMap ccmap = truthtiling[i]->ccmap();
    total_ccmap.insert(ccmap.begin(),ccmap.end());
  }

  Double_t charge_min = 10000;
  Double_t charge_max = 0;
  
  TApplication theApp("theApp",&argc,argv);
  theApp.SetReturnFromRun(true);
  
  TCanvas c1("ToyMC","ToyMC",800,600);
  c1.Draw();
  
  WCP2dToy::ToyEventDisplay display(c1, gds);
  display.charge_min = charge_min;
  display.charge_max = charge_max;
  
  
  gStyle->SetOptStat(0);
  
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Int_t MyPalette[NCont];
  Double_t stops[NRGBs] = {0.0, 0.34, 0.61, 0.84, 1.0};
  Double_t red[NRGBs] = {0.0, 0.0, 0.87 ,1.0, 0.51};
  Double_t green[NRGBs] = {0.0, 0.81, 1.0, 0.2 ,0.0};
  Double_t blue[NRGBs] = {0.51, 1.0, 0.12, 0.0, 0.0};
  Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  for (int kk=0;kk!=NCont;kk++) MyPalette[kk] = FI+kk;
  gStyle->SetPalette(NCont,MyPalette);
  
  
  
  display.init(0,10.3698,-2.33/2.,2.33/2.);
  //display.init(1.1,1.8,0.7,1.0);
  display.draw_mc(1,WCP::PointValueVector(),"colz");
  
  
  
  // display.draw_slice(slice,""); // draw wire 
  display.draw_cells(total_cells,"*same");
  //display.draw_mergecells(,"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell
  //display.draw_mergecells(calmcell,"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell
  display.draw_truthcells(total_ccmap,"*same");
  
  // display.draw_wires_charge(wcmap,"Fsame",FI);
  // display.draw_cells_charge(toytiling.get_allcell(),"Fsame");
  // display.draw_truthcells_charge(ccmap,"lFsame",FI);
  
  
  theApp.Run();


 
  return 0;
  
} // main()
