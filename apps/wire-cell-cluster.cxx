#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/ClusterDisplay.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"
#include "WireCellData/SpaceCell.h"


#include "TApplication.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TTree.h"
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
    cerr << "usage: wire-cell-cluster /path/to/shower_3D.root cluster#" << endl;
    return 1;
  }
   
  TString filename = argv[1];
  int ncluster = atoi(argv[2]);


  TFile *file = new TFile(filename);
  TTree *T = (TTree*)file->Get("T");
  TTree *TC = (TTree*)file->Get("TC");
  
  //cout << T->GetEntries() << " " << TC->GetEntries() << endl;
  
  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
  int time_slice;
  WireCell2dToy::ToyTiling* tt = 0;
  
  T->SetBranchAddress("time_slice",&time_slice);
  T->SetBranchAddress("toytiling",&tt);

  T->GetEntry(855);
  // cout << tt->get_allwire().size() << " " << tt->get_allcell().size() << endl;

  const GeomCell *cell = tt->get_allcell().at(5);
  const GeomWire *wire = tt->get_allwire().at(3);
  //cout << cell->cross_section() << " " << cell->center().y << endl;

  // GeomCellMap cellmap = tt->cmap();
  // GeomWireMap wiremap = tt->wmap();
  // WireChargeMap wirechargemap = tt->wcmap();

  // cout << wirechargemap[wire] << endl;
  //GeomWireSelection wires = cellmap[cell];
  //cout << wires.size() << " " << endl;
  
  double charge, x,y,z;
  int cluster_num;
  TC->SetBranchAddress("time_slice",&time_slice);
  TC->SetBranchAddress("charge",&charge);
  TC->SetBranchAddress("x",&x);
  TC->SetBranchAddress("y",&y);
  TC->SetBranchAddress("z",&z);
  TC->SetBranchAddress("ncluster",&cluster_num);
  TC->SetBranchAddress("cell",&cell);
  
  SpaceCellSelection cells;

  //int flag = 0;

  for (int i=0;i!=TC->GetEntries();i++){
    TC->GetEntry(i);

    // if (flag == 0 && cluster_num == ncluster){
    //   flag = 1;
    // }
    // if (flag==1 && cluster_num != ncluster){
    //   break;
    // }

    if (cluster_num == ncluster){
      SpaceCell *space_cell = new SpaceCell(cluster_num,*cell,x*units::cm,charge,0.32*units::cm);
      cells.push_back(space_cell);
    }
  }

  TApplication theApp("theApp",&argc,argv);
  theApp.SetReturnFromRun(true);
  
  TCanvas c1("ToyMC","ToyMC",800,600);
  c1.Draw();
  
  WireCell2dToy::ClusterDisplay display(c1);
  display.DrawCluster(cells);

  Point p;
  p.x = cells.at(0)->x();
  p.y = cells.at(0)->y();
  p.z = cells.at(0)->z();
  display.DrawHough(cells,p,-1,10*units::m);
  

  theApp.Run();
  //std::cout << cells.size() << std::endl;
  //successfully read the TC tree 
  // TC->GetEntry(0);
  //cout << x << " " << y << " " << z << " " << charge << " " << time_slice << " " << cluster_num << " " << cell->cross_section() << " " << cell->center().y << endl;
    
    
}
