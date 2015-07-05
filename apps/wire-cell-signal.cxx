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
  //WireCell::GenerativeFDS gfds(toydep,gds,9600,5,0.5*1.605723*units::millimeter); // 87 K at 0.5 kV/cm
  WireCell::GenerativeFDS gfds(toydep,gds,9600,5,0.5*1.60*units::millimeter); // 87 K at 0.5 kV/cm
  WireCell2dToy::ToySignalSimuFDS simu_fds(gfds,gds,9600,5,1.647,1.539+1.647); // time offset among different planes for the time electrons travel among different planes
  simu_fds.jump(1);
  //simu_fds.Save();

  // WireCell2dToy::ToySignalSimuTrueFDS st_fds(gfds,gds,9600,5); //truth
  // st_fds.jump(1);
  // st_fds.Save();
  
  WireCell2dToy::ToySignalGausFDS gaus_fds(simu_fds,gds,9600/4,5,1.647,1.539+1.647); // gaussian smearing for charge estimation
  gaus_fds.jump(1);
  gaus_fds.Save();
  
  WireCell2dToy::ToySignalWienFDS wien_fds(simu_fds,gds,9600/4,5,1.647,1.539+1.647); // weiner smearing for hit identification
  wien_fds.jump(1);
  wien_fds.Save();
  
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  int nwire_u = wires_u.size();
  int nwire_v = wires_v.size();
  int nwire_w = wires_w.size();
  
  float threshold_u = 5.87819e+02 * 3.5;
  float threshold_v = 8.36644e+02 * 3.5;
  float threshold_w = 5.67974e+02 * 3.5;

  float threshold_ug = 410.543*2.5;
  float threshold_vg = 631.936*2.5;
  float threshold_wg = 315.031*2.5;

  // float threshold_u = 1000;
  // float threshold_v = 1000;
  // float threshold_w = 1000;
  

  WireCellSst::ToyuBooNESliceDataSource sds(wien_fds,gaus_fds,threshold_u, 
					    threshold_v, threshold_w, 
					    threshold_ug, 
					    threshold_vg, threshold_wg, 
					    nwire_u, 
					    nwire_v, nwire_w); 
  
  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
  WireCell2dToy::MergeToyTiling **mergetiling = new WireCell2dToy::MergeToyTiling*[2400];
  WireCell2dToy::TruthToyTiling **truthtiling = new WireCell2dToy::TruthToyTiling*[2400];
   
  int start_num =459 + 800;
  int end_num = 459 + 800;
  for (int i=start_num;i!=end_num+1;i++){
    sds.jump(i);
    WireCell::Slice slice = sds.get();
    toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds);
    mergetiling[i] = new WireCell2dToy::MergeToyTiling(*toytiling[i],i);
    
    GeomCellSelection allcell = toytiling[i]->get_allcell();
    GeomWireSelection allwire = toytiling[i]->get_allwire();
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    GeomWireSelection allmwire = mergetiling[i]->get_allwire();
    cout << i << " " << allmcell.size() << " " << allmwire.size() << endl;
    truthtiling[i] = new WireCell2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds,800);
    
    // for (int j=0;j!=allmwire.size();j++){
    //   cout << mergetiling[i]->cells(*allmwire.at(j)).size() << endl;
    // }

    CellChargeMap ccmap = truthtiling[i]->ccmap();
    
    Double_t charge_min = 10000;
    Double_t charge_max = 0;

    TApplication theApp("theApp",&argc,argv);
    theApp.SetReturnFromRun(true);
    
    TCanvas c1("ToyMC","ToyMC",800,600);
    c1.Draw();
    
    WireCell2dToy::ToyEventDisplay display(c1, gds);
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
    display.draw_mc(1,WireCell::PointValueVector(),"colz");
    
    

    display.draw_slice(slice,""); // draw wire 
    display.draw_cells(toytiling[i]->get_allcell(),"*same");
    display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell
    //display.draw_mergecells(calmcell,"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell
    display.draw_truthcells(ccmap,"*same");
    
    // display.draw_wires_charge(wcmap,"Fsame",FI);
    // display.draw_cells_charge(toytiling.get_allcell(),"Fsame");
    // display.draw_truthcells_charge(ccmap,"lFsame",FI);
    
    
    theApp.Run();
  }

 
  return 0;
  
} // main()
