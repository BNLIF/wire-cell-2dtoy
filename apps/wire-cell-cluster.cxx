#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/ClusterDisplay.h"
#include "WireCell2dToy/ToyCrawler.h"
#include "WireCell2dToy/ToyTracking.h"

#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"
#include "WireCellData/SpaceCell.h"
#include "WireCellData/MergeSpaceCell.h"



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

  
  TString filename = argv[2];
  int ncluster = atoi(argv[3]);


  TFile *file = new TFile(filename);
  TTree *T = (TTree*)file->Get("T");
  TTree *TC = (TTree*)file->Get("TC");
  TTree *Trun = (TTree*)file->Get("Trun");
  
  float unit_dis;
  int nrebin;
  int total_time_bin;
  Trun->SetBranchAddress("nrebin",&nrebin);
  Trun->SetBranchAddress("unit_dis",&unit_dis);
  Trun->SetBranchAddress("total_time_bin",&total_time_bin);
  Trun->GetEntry(0);
  
  // unit_dis = 1.6;
  // total_time_bin = 9600;
  // nrebin = 4;

  
  //cout << nrebin << " " << unit_dis << " " << total_time_bin << endl;
  const int ntime = total_time_bin/nrebin;
  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[ntime];
  for (int i=0;i!=ntime;i++){
    toytiling[i] = new WireCell2dToy::ToyTiling();
  }


  int time_slice;
  // WireCell2dToy::ToyTiling* tt = 0;
  // T->SetBranchAddress("time_slice",&time_slice);
  // T->SetBranchAddress("toytiling",&tt);

  

  //T->GetEntry(855);
  // cout << tt->get_allwire().size() << " " << tt->get_allcell().size() << endl;

  const GeomCell *cell = 0;//tt->get_allcell().at(5);
  const GeomWire *wire = 0;//tt->get_allwire().at(3);
  //cout << cell->cross_section() << " " << cell->center().y << endl;

  // GeomCellMap cellmap = tt->cmap();
  // GeomWireMap wiremap = tt->wmap();
  // WireChargeMap wirechargemap = tt->wcmap();

  // cout << wirechargemap[wire] << endl;
  //GeomWireSelection wires = cellmap[cell];
  //cout << wires.size() << " " << endl;
  
  double charge, x,y,z;
  int cluster_num;
  int mcell_id;
  TC->SetBranchAddress("time_slice",&time_slice);
  TC->SetBranchAddress("charge",&charge);
  TC->SetBranchAddress("xx",&x);
  TC->SetBranchAddress("yy",&y);
  TC->SetBranchAddress("zz",&z);
  TC->SetBranchAddress("ncluster",&cluster_num);
  TC->SetBranchAddress("mcell_id",&mcell_id);
  TC->SetBranchAddress("cell",&cell);
  
  int u_index, v_index, w_index;
  double u_charge, v_charge, w_charge;
  double u_charge_err, v_charge_err, w_charge_err;

  TC->SetBranchAddress("u_index",&u_index);
  TC->SetBranchAddress("v_index",&v_index);
  TC->SetBranchAddress("w_index",&w_index);
  TC->SetBranchAddress("u_charge",&u_charge);
  TC->SetBranchAddress("v_charge",&v_charge);
  TC->SetBranchAddress("w_charge",&w_charge);
  TC->SetBranchAddress("u_charge_err",&u_charge_err);
  TC->SetBranchAddress("v_charge_err",&v_charge_err);
  TC->SetBranchAddress("w_charge_err",&w_charge_err);


  

  
  int prev_mcell_id = -1;
  MergeSpaceCellSelection mcells; // save all the cells
  
  int flag = 0;
  SpaceCellSelection cells;


  MergeSpaceCell *mcell;
  
  for (int i=0;i!=TC->GetEntries();i++){
    TC->GetEntry(i);
   
    //cluster starting at 0
    if (cluster_num == ncluster){
      
      if (flag == 0){
	mcell = new MergeSpaceCell();
	flag = 1;
      }else if (flag==1 && mcell_id!=prev_mcell_id){
	mcells.push_back(mcell);
	mcell = new MergeSpaceCell();
      }

      //cout << x << " " << y << " " << z << " " << charge << endl;
      
      GeomCell *cell1 = new GeomCell(cell);
      toytiling[time_slice]->AddCell(gds,cell1,u_index,v_index,w_index,u_charge,v_charge,w_charge,u_charge_err,v_charge_err,w_charge_err);
      

      SpaceCell *space_cell = new SpaceCell(cluster_num,*cell1,x*units::cm,charge,unit_dis/10.*nrebin/2.*units::cm);
      mcell->AddSpaceCell(space_cell);
      cells.push_back(space_cell);
      

      prev_mcell_id = mcell_id;
    }

  }
  mcells.push_back(mcell);
  

  // for (int i=0;i!=ntime;i++){
  //   GeomCellSelection allcell = toytiling[i]->get_allcell();
  //   GeomWireSelection allwire = toytiling[i]->get_allwire();
  //   cout << "Single Cell: " << i << " "  << allcell.size() << " " << allwire.size() << endl;
  // }


  //cout << mcells.size() << endl;

  // do the Toy Crawler
  std::cout << "Crawling " << std::endl;
  WireCell2dToy::ToyCrawler toycrawler(mcells);

  std::cout << "Tracking " << std::endl;
  WireCell2dToy::ToyTracking toytracking(toycrawler);
  


  std::cout << "Drawing " << std::endl; 
  TApplication theApp("theApp",&argc,argv);
  theApp.SetReturnFromRun(true);
  
  TCanvas c1("ToyMC","ToyMC",800,600);
  c1.Draw();
  
  WireCell2dToy::ClusterDisplay display(c1);
  // display.DrawCluster(cells);
  //display.DrawCluster(mcells);
  display.DrawCluster(mcells,toytracking);
  // display.DrawCrawler(toycrawler,"psame",1);

  
  WCVertexSelection& vertices = toytracking.get_vertices();
  display.DrawVertex(vertices,"psame");
  
  WCTrackSelection& bad_tracks = toytracking.get_bad_tracks();
  display.DrawTracks(bad_tracks,"same",2);

  WCTrackSelection& short_tracks = toytracking.get_short_tracks();
  display.DrawTracks(short_tracks,"psame",4);

  // Point p;
  // p.x = cells.at(0)->x();
  // p.y = cells.at(0)->y();
  // p.z = cells.at(0)->z();
  // display.DrawHough(cells,p,-1,10*units::m);
  

  theApp.Run();
  //std::cout << cells.size() << std::endl;
  //successfully read the TC tree 
  // TC->GetEntry(0);
  //cout << x << " " << y << " " << z << " " << charge << " " << time_slice << " " << cluster_num << " " << cell->cross_section() << " " << cell->center().y << endl;
    
    
}
