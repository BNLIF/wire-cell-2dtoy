#include "WireCellSst/GeomDataSource.h"
//#include "WireCellSst/ToyuBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCellData/MergeGeomCell.h"
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
  
  WireCell::FrameDataSource* fds = 0;
  fds = WireCellSst::make_fds(root_file);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  

  
 

  // WireCell::SimDataSource* sim = dynamic_cast<WireCell::SimDataSource*>(fds);
  // if (!sim) {
  //   cerr << "ERROR: the FDS is not also an SimDS " << endl;
  //   return 2;
  // }
  // fds->jump(1);
  // WireCell::SimTruthSelection sts = sim->truth();
  // cerr << "Got " << sts.size() << " true hits" << endl;
  // for (size_t itruth = 0; itruth < sts.size(); ++itruth) {
  //   const WireCell::SimTruth* st = sts[itruth];
  //   cerr << "Hit: "
  // 	 << " @ (" << st->x() << " " << st->y() << " " << st->z() << ")"
  // 	 << " q=" << st->charge()
  // 	 << " tdc=" << st->tdc()
  // 	 << endl;
  // }
  // cout << units::cm << endl;


  WireCell::ToyDepositor toydep(fds);
  const PointValueVector pvv = toydep.depositions(1);
  
  // for (int itruth = 0; itruth < pvv.size(); ++itruth){
  //   cout << pvv[itruth].first.x << " " << pvv[itruth].first.y << " " << pvv[itruth].first.z << " " << pvv[itruth].second << endl;
  // }


  WireCell::GenerativeFDS gfds(toydep,gds,800,5);
  gfds.jump(1);

  WireCellSst::ToyuBooNESliceDataSource sds(gfds,2000); //set threshold at 2000 electrons





  const int N = 100000;
  Double_t x[N],y[N],z[N];
  int ncount = 0;
  
  int i=454;
  //for (int i=0;i!=sds.size();i++){
    sds.jump(i);
    WireCell::Slice slice = sds.get();
    if ( slice.group().size() >0){
      WireCell2dToy::ToyTiling toytiling(slice,gds);
      WireCell2dToy::MergeToyTiling mergetiling(toytiling);
      
      GeomCellSelection allcell = toytiling.get_allcell();
      GeomCellSelection allmcell = mergetiling.get_allcell();
      GeomWireSelection allwire = mergetiling.get_allwire();

      for (int j=0;j!=allcell.size();j++){
	Point p = allcell[j]->center();
	x[ncount] = i*0.32*units::cm;
	y[ncount] = p.y/units::cm;
	z[ncount] = p.z/units::cm;
	ncount ++;
      }

      // cout << i << " " << allcell.size() << endl;

    // int sum = 0;
    // for (int j=0;j!=allmcell.size();j++){
    //   sum += ((WireCell::MergeGeomCell*)allmcell[j])->get_allcell().size() ;
    // }
    // cout << allcell.size() << " " << allmcell.size() << " "  << sum << endl;

     // for (int j=0;j!=allmcell.size();j++){
     //   cout << j << " " << mergetiling.wires(*allmcell[j]).size() << endl;
     // }

    // for (int j=0;j!=allwire.size();j++){
    //   const GeomWire *wire = allwire[j];
    //   const GeomCellSelection targetcells =  mergetiling.cells(*wire);
    //   cout << j << " " << targetcells.size() << endl;
    // }

    

    // GeomCellSelection allcell = toytiling.get_allcell();
    // GeomWireSelection allwire = toytiling.get_allwire();
    // cout << i << " " << allcell.size() << " " << allwire.size() << endl;
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
    
    display.init(0,10.3698,-2.33/2.,2.33/2.);
    //display.init(0.6,1.0,0.0,0.3);
    //display.init(0.6,0.7,0.07,0.12);
    //display.init();
    display.draw_mc(1,WireCell::PointValueVector(),"");
    //display.draw_mc(1,fds.mctruth,"");
    //display.draw_mc(2,fds.mctruth,"TEXT");
    
    
    display.draw_slice(slice,"");
    
    display.draw_cells(toytiling.get_allcell(),"*same");
    display.draw_mergecells(mergetiling.get_allcell(),"*same");
    //display.draw_mc(3,fds.mctruth,"*same");
    
    theApp.Run();
    //}
  }

    TGraph2D *g = new TGraph2D(ncount,x,y,z);
    TFile *file = new TFile("shower3D.root","RECREATE");
  g->Write("shower3D");
  file->Write();
  file->Close();

  return 0;
  
} // main()
