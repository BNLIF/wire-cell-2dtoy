#include "WCPSst/GeomDataSource.h"
#include "WCPSst/ToyuBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/MergeToyTiling.h"

#include "WCPData/MergeGeomCell.h"
//#include "WCPNav/SliceDataSource.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
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
  
  TFile tfile(root_file,"read");
  TTree* sst = dynamic_cast<TTree*>(tfile.Get(tpath));
  WCPSst::ToyuBooNEFrameDataSource fds(*sst,gds);
  std::cerr << "Got " << fds.size() 
	    << " frames from " << tpath 
	    << " in " << root_file << std::endl;
  
  fds.jump(1);
  WCP::Frame frame = fds.get();
  
  WCPSst::ToyuBooNESliceDataSource sds(fds,1);
  
  // int i=1129;{
  int i=331;{
  //for (int i=0;i!=sds.size();i++){
  //for (int i=1143;i!=1145;i++){
  sds.jump(i);
  WCP::Slice slice = sds.get();
  if ( slice.group().size() >0){
    WCP2dToy::ToyTiling toytiling(slice,gds);
    GeomCellSelection allcell = toytiling.get_allcell();

   
    

    WCP2dToy::MergeToyTiling mergetiling(toytiling);    
    GeomCellSelection allmcell = mergetiling.get_allcell();
    GeomWireSelection allwire = mergetiling.get_allwire();
    //cout << i << endl;
    // int sum = 0;
    // for (int j=0;j!=allmcell.size();j++){
    //   GeomCellSelection cells = ((WCP::MergeGeomCell*)allmcell[j])->get_allcell();
    //   for (int k=0;k!=cells.size();k++){
    // 	GeomWireSelection wires = toytiling.wires(*cells[k]);
    // 	if( wires[0]->ident()==0 || wires[1]->ident()==0 || wires[2]->ident()==0){
    // 	  cout << i << " Wrong!!" << endl;
    // 	}
    //   }
    //   //   sum += ((WCP::MergeGeomCell*)allmcell[j])->get_allcell().size() ;
    // }

    //  //debug the toytiling itself. 
    // for (int j=0;j!=allcell.size();j++){
    //   GeomWireSelection wires = toytiling.wires(*allcell[j]);
    //   if( wires[0]->ident()==0 || wires[1]->ident()==0 || wires[2]->ident()==0){
    // 	cout << i << " Wrong!!" << endl;
    //   }  
    // }



    if (allcell.size()>0){
      cout << i << " " << allmcell.size() << " "  << allwire.size() << endl;
    }

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
    
    WCP2dToy::ToyEventDisplay display(c1, gds);
    
    gStyle->SetOptStat(0);
    
    display.init(0,10.3698,-2.33/2.,2.33/2.);
    //display.init(0.6,1.0,0.0,0.3);
    //display.init(0.6,0.7,0.07,0.12);
    //display.init();
    display.draw_mc(1,WCP::PointValueVector(),"");
    //display.draw_mc(1,fds.mctruth,"");
    //display.draw_mc(2,fds.mctruth,"TEXT");
    
    
    display.draw_slice(slice,"");
    
    display.draw_cells(toytiling.get_allcell(),"*same");
    display.draw_mergecells(mergetiling.get_allcell(),"*same");
    //display.draw_mc(3,fds.mctruth,"*same");
    
    theApp.Run();
  }
  }

  return 0;
  
} // main()
