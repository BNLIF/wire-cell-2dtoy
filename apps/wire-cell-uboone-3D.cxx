#include "WireCellSst/GeomDataSource.h"
//#include "WireCellSst/ToyuBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/GeomCluster.h"
//#include "WireCellNav/SliceDataSource.h"

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/SimDataSource.h"
#include "WireCellNav/SliceDataSource.h"
#include "WireCellSst/Util.h"
#include "WireCellData/SimTruth.h"
#include "WireCell2dToy/ToyDepositor.h"
#include "WireCellNav/GenerativeFDS.h"


#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph2D.h"
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
  

  
  WireCell::ToyDepositor toydep(fds);
  const PointValueVector pvv = toydep.depositions(1);
  

  WireCell::GenerativeFDS gfds(toydep,gds,2400,5);
  gfds.jump(1);

  WireCellSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons
  
  const int N = 100000;
  Double_t x[N],y[N],z[N];
  
  int ncount = 0;
  //int i=1140;{
  
  for (int i=0; i!=sds.size();i++){
    sds.jump(i);
    WireCell::Slice slice = sds.get();
    if ( slice.group().size() >0){
      WireCell2dToy::ToyTiling toytiling(slice,gds);
      //WireCell2dToy::MergeToyTiling mergetiling(toytiling);
      
      GeomCellSelection allcell = toytiling.get_allcell();
      
      for (int j=0;j!=allcell.size();j++){
	Point p = allcell[j]->center();
	x[ncount] = i*0.32;
	y[ncount] = p.y/units::cm;
	z[ncount] = p.z/units::cm;
	ncount ++;
      }
      
      // GeomCellSelection allmcell = mergetiling.get_allcell();
      // GeomWireSelection allwire = mergetiling.get_allwire();

      // cout << i << " " << allmcell.size() << " " << allwire.size() << endl;

      // int sum = 0;
      // for (int j=0;j!=allmcell.size();j++){
      //   sum += ((WireCell::MergeGeomCell*)allmcell[j])->get_allcell().size() ;
      // }
      // cout << allcell.size() << " " << allmcell.size() << " "  << sum << endl;
      
      
    }
  }


  TApplication theApp("theApp",&argc,argv);
  theApp.SetReturnFromRun(true);
  
  TCanvas c1("ToyMC","ToyMC",800,600);
  c1.Draw();
  gStyle->SetOptStat(0);

  TGraph2D *g = new TGraph2D(ncount,x,y,z);
  g->Draw("p");

  TFile *file = new TFile("shower3D.root","RECREATE");
  g->Write("shower3D");
  file->Write();
  file->Close();

  //theApp.Run();
  return 0;
  
} // main()
