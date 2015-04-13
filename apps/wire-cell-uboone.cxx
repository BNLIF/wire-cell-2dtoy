#include "WireCellSst/GeomDataSource.h"
#include "WireCellSst/ToyuBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
//#include "WireCellNav/SliceDataSource.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
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
  std::vector<float> ex = gds.extent();
  cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;


  const char* root_file = argv[2];
  const char* tpath = "/Event/Sim";
  
  TFile tfile(root_file,"read");
  TTree* sst = dynamic_cast<TTree*>(tfile.Get(tpath));
  WireCellSst::ToyuBooNEFrameDataSource fds(*sst);
  std::cerr << "Got " << fds.size() 
	    << " frames from " << tpath 
	    << " in " << root_file << std::endl;
  
  fds.jump(1);
  WireCell::Frame frame = fds.get();
  
  WireCellSst::ToyuBooNESliceDataSource sds(fds,1);
  
  int i=1137;
  //int i=1146;
  // for (int i=0;i!=sds.size();i++){
    sds.jump(i);
    WireCell::Slice slice = sds.get();
    if ( slice.group().size() >0){
      WireCell::ToyTiling toytiling(slice,gds);
      
      GeomCellSelection allcell = toytiling.get_allcell();
      GeomWireSelection allwire = toytiling.get_allwire();
      cout << i << " " << allcell.size() << " " << allwire.size() << endl;
      //
    //}
  

  // cout << toytiling.wiremap[allwire.at(0)].size() << endl;
  // cout << toytiling.cellmap[allcell.at(0)].size() << endl;
  // //  

  
  TApplication theApp("theApp",&argc,argv);
  theApp.SetReturnFromRun(true);

  TCanvas c1("ToyMC","ToyMC",800,600);
  c1.Draw();

  WireCell2dToy::ToyEventDisplay display(c1, gds);
  
  gStyle->SetOptStat(0);

  //display.init(0,10.3698,-2.33/2.,2.33/2.);
  display.init(0.6,1.0,0.0,0.3);
  //display.init(0.6,0.7,0.07,0.12);
  //display.init();
  display.draw_mc(1,WireCell::PointValueVector(),"");
  //display.draw_mc(1,fds.mctruth,"");
  //display.draw_mc(2,fds.mctruth,"TEXT");
  
  
  display.draw_slice(slice,"");
 
  display.draw_cells(toytiling.get_allcell(),"*same");
  //display.draw_mc(3,fds.mctruth,"*same");

  theApp.Run();
    }

  return 0;
  
} // main()
