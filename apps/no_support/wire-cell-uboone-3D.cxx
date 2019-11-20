#include "WCPSst/GeomDataSource.h"
//#include "WCPSst/ToyuBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/TruthToyTiling.h"
#include "WCPData/MergeGeomCell.h"
#include "WCPData/GeomCluster.h"
//#include "WCPNav/SliceDataSource.h"

#include "WCPNav/FrameDataSource.h"
#include "WCPNav/SimDataSource.h"
#include "WCPNav/SliceDataSource.h"
#include "WCPSst/Util.h"
#include "WCPData/SimTruth.h"
#include "WCP2dToy/ToyDepositor.h"
#include "WCPNav/GenerativeFDS.h"


#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TFile.h"
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
  

  
  WCP::ToyDepositor toydep(fds);
  const PointValueVector pvv = toydep.depositions(1);
  

  WCP::GenerativeFDS gfds(toydep,gds,2400,5,2.0*1.6*units::millimeter);
  gfds.jump(1);

  WCPSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons
  
  const int N = 100000;
  Double_t x[N],y[N],z[N];
  
  int ncount = 0;
  //int i=1140;{
  
  for (int i=0; i!=sds.size();i++){
    sds.jump(i);
    WCP::Slice slice = sds.get();
    if ( slice.group().size() >0){
      WCP2dToy::ToyTiling toytiling(slice,gds);
      //WCP2dToy::MergeToyTiling mergetiling(toytiling);
      
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
      //   sum += ((WCP::MergeGeomCell*)allmcell[j])->get_allcell().size() ;
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
