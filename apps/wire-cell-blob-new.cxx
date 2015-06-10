#include "WireCellSst/GeomDataSource.h"
//#include "WireCellSst/ToyuBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCell2dToy/BlobToyTiling.h"
#include "WireCell2dToy/CaveToyTiling.h"


#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/GeomCluster.h"
//#include "WireCellNav/SliceDataSource.h"
#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/ToyMatrixIterate.h"
#include "WireCell2dToy/ToyMatrixMarkov.h"

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
#include "TColor.h"
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
  

  
  WireCell::ToyDepositor toydep(fds);
  const PointValueVector pvv = toydep.depositions(1);
  

  WireCell::GenerativeFDS gfds(toydep,gds,2400,5);
  gfds.jump(1);

  WireCellSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons

  

  const int N = 100000;
  Double_t x[N],y[N],z[N];
  Double_t xt[N],yt[N],zt[N];
  int ncount = 0;
  int ncount_t = 0;
  

  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
  WireCell2dToy::MergeToyTiling **mergetiling = new WireCell2dToy::MergeToyTiling*[2400];
  WireCell2dToy::TruthToyTiling **truthtiling = new WireCell2dToy::TruthToyTiling*[2400];
  WireCell2dToy::CaveToyTiling **cavetiling = new WireCell2dToy::CaveToyTiling*[2400];
  
  WireCell2dToy::ToyMatrix **toymatrix = new WireCell2dToy::ToyMatrix*[2400];
  WireCell2dToy::ToyMatrixMarkov **toymatrix_markov = new WireCell2dToy::ToyMatrixMarkov*[2400];
  // WireCell2dToy::ToyMatrix **blobmatrix = new WireCell2dToy::ToyMatrix*[2400];
  
  // WireCell2dToy::ToyMatrixMarkov **blobmatrix_markov = new WireCell2dToy::ToyMatrixMarkov*[2400];
 
  //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  
  int ncount_mcell = 0;
  
  int i=454;{
    //for (int i=0;i!=sds.size();i++){
    
    sds.jump(i);
    WireCell::Slice slice = sds.get();
          
    toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds);
    mergetiling[i] = new WireCell2dToy::MergeToyTiling(*toytiling[i],i);
    truthtiling[i] = new WireCell2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds);
    GeomCellSelection allcell = toytiling[i]->get_allcell();
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    GeomWireSelection allmwire = mergetiling[i]->get_allwire();
    cout << i << " " << allcell.size() << " " << allmcell.size() << " " << allmwire.size()  << endl;
    
    toymatrix[i] = new WireCell2dToy::ToyMatrix(*toytiling[i],*mergetiling[i]);
    
    if (toymatrix[i]->Get_Solve_Flag()==0){
      toymatrix_markov[i] = new WireCell2dToy::ToyMatrixMarkov(toymatrix[i],&allmcell);
    }
    
    cout << "chi2: " << toymatrix[i]->Get_Chi2() << endl;
    cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;

    
    // for now put this part here
    cavetiling[i] = new WireCell2dToy::CaveToyTiling(toytiling[i],*mergetiling[i],*toymatrix[i]);
    
    //blobtiling[i] = new WireCell2dToy::BlobToyTiling(*toytiling[i],*mergetiling[i],*toymatrix[i],i,5);
    //blobmatrix[i] = new WireCell2dToy::ToyMatrix(*toytiling[i],*blobtiling[i]);

    // GeomCellSelection ballmcell = blobtiling[i]->get_allcell();
    // GeomWireSelection ballmwire = blobtiling[i]->get_allwire();
    // cout << i << " " << ballmcell.size() << " " << ballmwire.size()  << endl;
    
    // if (blobmatrix[i]->Get_Solve_Flag()==0)
    //   blobmatrix_markov[i] = new WireCell2dToy::ToyMatrixMarkov(blobmatrix[i],&ballmcell);
    
    // cout << "chi2: " << blobmatrix[i]->Get_Chi2() << endl;
    // cout << "NDF: " << blobmatrix[i]->Get_ndf() << endl;
    
    CellChargeMap ccmap = truthtiling[i]->ccmap();
    
    Double_t charge_min = 10000;
    Double_t charge_max = 0;
    
    for (auto it = ccmap.begin();it!=ccmap.end(); it++){
      Point p = it->first->center();
      xt[ncount_t] = i*0.32;
      yt[ncount_t] = p.y/units::cm;
      zt[ncount_t] = p.z/units::cm;
      ncount_t ++;
      
      float charge = it->second;
      if (charge > charge_max) charge_max = charge;
      if (charge < charge_min) charge_min = charge;
      // cout << it->second << endl;
    }
    
    
    //loop through merged cell and compare with truth cells
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      mcell->CheckContainTruthCell(ccmap);
      // 	cout << mergetiling.wires(*allmcell[j]).size() << endl;
    }
    
    
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
      //display.init(1.4,1.7,0.5,1.2);
      
      display.draw_mc(1,WireCell::PointValueVector(),"colz");
      
      display.draw_slice(slice,"");
      //display.draw_cells(toytiling[i]->get_allcell(),"*same");
      //display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",1); //0 is normal, 1 is only draw the ones containt the truth cell
      //display.draw_reconcells(blobtiling[i]->get_allcell(),blobmatrix[i],"*same");
      display.draw_reconcells(mergetiling[i]->get_allcell(),toymatrix[i],"*same");
      display.draw_truthcells(ccmap,"*same");
      
      // display.draw_wires_charge(wcmap,"Fsame",FI);
      // display.draw_cells_charge(toytiling.get_allcell(),"Fsame");
      // display.draw_truthcells_charge(ccmap,"lFsame",FI);
      
      theApp.Run();
    
  }

  

  return 0;
  
} // main()
