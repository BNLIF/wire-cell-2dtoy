#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/ClusterDisplay.h"
#include "WireCell2dToy/ToyCrawler.h"

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
   if (argc < 2) {
    cerr << "usage: wire-cell-allcluster /path/to/shower_3D.root " << endl;
    return 1;
  }
   
  TString filename = argv[1];
  

  TFile *file = new TFile(filename);
  TTree *T = (TTree*)file->Get("T");
  TTree *TC = (TTree*)file->Get("TC");
  
  //cout << T->GetEntries() << " " << TC->GetEntries() << endl;
  
  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
  int time_slice;
  WireCell2dToy::ToyTiling* tt = 0;
  
  T->SetBranchAddress("time_slice",&time_slice);
  T->SetBranchAddress("toytiling",&tt);

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
  TC->SetBranchAddress("x",&x);
  TC->SetBranchAddress("y",&y);
  TC->SetBranchAddress("z",&z);
  TC->SetBranchAddress("ncluster",&cluster_num);
  TC->SetBranchAddress("mcell_id",&mcell_id);
  TC->SetBranchAddress("cell",&cell);
  
  //save all the crawlers ... 
  std::vector<WireCell2dToy::ToyCrawler*> crawlers;
  

  int prev_mcell_id = -1;
  int prev_cluster_num = -1;
  
  MergeSpaceCellSelection mcells; // save all the cells
  int flag = 0;
  MergeSpaceCell *mcell;
  SpaceCellSelection cells;
  
  for (int i=0;i!=TC->GetEntries();i++){
    TC->GetEntry(i);

    if (cluster_num != prev_cluster_num){
      if (prev_cluster_num!=-1){
	mcells.push_back(mcell);  
	WireCell2dToy::ToyCrawler* toycrawler = new WireCell2dToy::ToyCrawler(mcells);
	toycrawler->FormGraph();

	// if (cluster_num==3){
	//   TApplication theApp("theApp",&argc,argv);
	//   theApp.SetReturnFromRun(true);	  
	//   TCanvas c1("ToyMC","ToyMC",800,600);
	//   c1.Draw();
	//   WireCell2dToy::ClusterDisplay display(c1);
	//   display.DrawCluster(cells);
	//   display.DrawCluster(mcells);
	  
	
	//   display.DrawCrawler(*toycrawler,"psame",1);
	  
	//   theApp.Run();
	// }
	


	crawlers.push_back(toycrawler);
	mcells.clear();
	cells.clear();
	flag = 0;
      }
    }


    if (flag == 0){
      mcell = new MergeSpaceCell();
      flag = 1;
    }else if (flag==1 && (mcell_id!=prev_mcell_id)){
      mcells.push_back(mcell);
      mcell = new MergeSpaceCell();
    }
    GeomCell *cell1 = new GeomCell(cell);
    SpaceCell *space_cell = new SpaceCell(cluster_num,*cell1,x*units::cm,charge,0.32*units::cm);
    mcell->AddSpaceCell(space_cell);
    cells.push_back(space_cell);
    
    prev_cluster_num = cluster_num;
    prev_mcell_id = mcell_id;
   
  }

  mcells.push_back(mcell);  
  WireCell2dToy::ToyCrawler* toycrawler = new WireCell2dToy::ToyCrawler(mcells);
  toycrawler->FormGraph();
  crawlers.push_back(toycrawler);

  //check # of clusters 
  int sum = 0 ;
  for (int i=0;i!=crawlers.size();i++){
    for (auto it = crawlers.at(i)->Get_mcells_map().begin(); it!= crawlers.at(i)->Get_mcells_map().end();it++){
      MergeSpaceCell *mcell1 = it->first;
      sum += mcell1->Get_all_spacecell().size();
      //    sum += .size();
    }
  }
  
  std::cout << "Check: " << crawlers.size() << " " << TC->GetEntries() << " " << sum << std::endl;


  //start the prepare the important merge cell vectors
  // int start_num = 0 ;
  // int end_num = 2399;
  // std::vector<GeomCellSelection> Good_MCells;
  
  // for (int i=start_num;i!=end_num+1;i++){
  // GeomCellSelection cells;
  //Good_MCells.push_back(cells);
  // }

  MergeSpaceCellSelection ms_cells;
  
  
  for (int i=0;i!=crawlers.size();i++){
    WireCell2dToy::ToyCrawler *toycrawler = crawlers.at(i);
    for (int j=0;j!=toycrawler->Get_allMCT().size();j++){
      MergeClusterTrack *mct = toycrawler->Get_allMCT().at(j);
      int ntime = mct->Get_TimeLength();
      if (ntime >=5){
	// do something
	
	for (int k=0;k!=ntime;k++){
	  MergeSpaceCellSelection cells = mct->Get_MSCS(k);
	  if (cells.size()==1){
	    ms_cells.push_back(cells.at(0));
	  }else if (cells.size()>1){
	    MergeSpaceCell *cell = cells.at(0);
	    for (int i1 = 1; i1!=cells.size();i1++){
	      if (cell->Get_all_spacecell().size() < cells.at(i1)->Get_all_spacecell().size()){
		cell = cells.at(i1);
	      }
	    }
	    ms_cells.push_back(cell);
	  }
	}
	

      }
    }
  }
  
  // plot it
  
  TApplication theApp("theApp",&argc,argv);
  theApp.SetReturnFromRun(true);
  
  TCanvas c1("ToyMC","ToyMC",800,600);
  c1.Draw();
  
  WireCell2dToy::ClusterDisplay display(c1);
  display.DrawCluster(ms_cells);

  theApp.Run();
}
