#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/ClusterDisplay.h"
#include "WireCell2dToy/ToyCrawler.h"
#include "WireCell2dToy/ToyTracking.h"
#include "WireCell2dToy/ToyWalking.h"

#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"

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

  TPCParams& mp = Singleton<TPCParams>::Instance();
  
  double pitch_u = gds.pitch(WirePlaneType_t(0));
  double pitch_v = gds.pitch(WirePlaneType_t(1));
  double pitch_w = gds.pitch(WirePlaneType_t(2));
  double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;

  mp.set_pitch_u(pitch_u);
  mp.set_pitch_v(pitch_v);
  mp.set_pitch_w(pitch_w);
  mp.set_ts_width(time_slice_width);
  
  std::cout << "Singleton: " << mp.get_pitch_u() << " " << mp.get_pitch_v() << " " << mp.get_pitch_w() << " " << mp.get_ts_width() << std::endl;

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


  MergeSpaceCell *mcell=0;
  
  for (int i=0;i!=TC->GetEntries();i++){
    //std::cout << i << std::endl;
    TC->GetEntry(i);
   
    //std::cout << i << std::endl;

    //cluster starting at 0
    if (cluster_num == ncluster){
      
      if (flag == 0){
	mcell = new MergeSpaceCell();
	flag = 1;
      }else if (flag==1 && mcell_id!=prev_mcell_id){
	if (mcell->Get_all_spacecell().size()>0)
	  mcells.push_back(mcell);
	mcell = new MergeSpaceCell();
      }

      // cout << x << " " << y << " " << z << " " << charge << endl;
      
      GeomCell *cell1 = new GeomCell(cell);
      toytiling[time_slice]->AddCell(gds,cell1,u_index,v_index,w_index,u_charge,v_charge,w_charge,u_charge_err,v_charge_err,w_charge_err);
      

      SpaceCell *space_cell = new SpaceCell(cluster_num,*cell1,x*units::cm,charge,unit_dis/10.*nrebin/2.*units::cm);
      mcell->AddSpaceCell(space_cell);
      cells.push_back(space_cell);
      

      prev_mcell_id = mcell_id;
    }

  }
  if (mcell!=0){
    if (mcell->Get_all_spacecell().size()>0)
      mcells.push_back(mcell);
  }

  // for (int i=0;i!=ntime;i++){
  //   GeomCellSelection allcell = toytiling[i]->get_allcell();
  //   GeomWireSelection allwire = toytiling[i]->get_allwire();
  //   cout << "Single Cell: " << i << " "  << allcell.size() << " " << allwire.size() << endl;
  // }


  //cout << mcells.size() << endl;

  // do the Toy Crawler
  std::cout << "Crawling " << std::endl;
  WireCell2dToy::ToyCrawler toycrawler(mcells);
  MergeSpaceCellMap& mcells_map = toycrawler.Get_mcells_map();
  // find the center
  Point center(0,0,0);
  int sum =0;
  for (int i=0;i!=mcells.size();i++){
    center.x +=mcells.at(i)->Get_Center().x *mcells.at(i)->Get_all_spacecell().size();
    center.y +=mcells.at(i)->Get_Center().y *mcells.at(i)->Get_all_spacecell().size();
    center.z +=mcells.at(i)->Get_Center().z *mcells.at(i)->Get_all_spacecell().size();
    sum += mcells.at(i)->Get_all_spacecell().size();
  }
  center.x/= sum;
  center.y/= sum;
  center.z/= sum;
  
  float max_dis = 0;
  MergeSpaceCell *mcell1;
  // find the furthest one
  for (int i=0;i!=mcells.size();i++){
    float dis = sqrt(pow(mcells.at(i)->Get_Center().x-center.x,2) +
		     pow(mcells.at(i)->Get_Center().y-center.y,2) +
		     pow(mcells.at(i)->Get_Center().z-center.z,2));
    if (dis > max_dis){
      max_dis = dis;
      mcell1 = mcells.at(i);
    }
  }
  // find the furthest one to this one
  max_dis = 0;
  MergeSpaceCell *mcell2;
  for (int i=0;i!=mcells.size();i++){
    float dis = sqrt(pow(mcells.at(i)->Get_Center().x-mcell1->Get_Center().x,2) +
		     pow(mcells.at(i)->Get_Center().y-mcell1->Get_Center().y,2) +
		     pow(mcells.at(i)->Get_Center().z-mcell1->Get_Center().z,2));
    if (dis > max_dis){
      max_dis = dis;
      mcell2 = mcells.at(i);
    }
  }
  std::cout << mcell1->Get_Center().x/units::cm << " " << mcell1->Get_Center().y/units::cm << " " << mcell1->Get_Center().z/units::cm << std::endl;
  std::cout << mcell2->Get_Center().x/units::cm << " " << mcell2->Get_Center().y/units::cm << " " << mcell2->Get_Center().z/units::cm << std::endl;

  Point start_point = mcell1->Get_Center();
  Point end_point = mcell2->Get_Center();
  
  // do toywalking test ... 
  WireCell2dToy::ToyWalking walking(mcell1,start_point,mcell2,end_point,mcells_map,5000);
  MergeSpaceCellSelection mcells1 = walking.get_cells();
  std::cout << "xin: " << mcells.size() << " " << mcells1.size() << " " << walking.get_counter() << " " << walking.get_global_counter() << std::endl;

  Point sp_point(163.508*units::cm, -16.2064*units::cm, 363.983*units::cm);
  float min_dis = 1e9;
  MergeSpaceCell *mcell3;
  for (int i=0;i!=mcells.size();i++){
    float dis = sqrt(pow(mcells.at(i)->Get_Center().x-sp_point.x,2) +
		     pow(mcells.at(i)->Get_Center().y-sp_point.y,2) +
		     pow(mcells.at(i)->Get_Center().z-sp_point.z,2));
    if (dis < min_dis){
      min_dis = dis;
      mcell3 = mcells.at(i);
    }
  }
  std::cout << mcell3->Get_Center().x/units::cm << " " << mcell3->Get_Center().y/units::cm << " " << mcell3->Get_Center().z/units::cm << " " << mcell3->Get_all_spacecell().size() << std::endl;
  
  for (int i=0;i!=mcells_map[mcell3].size();i++){
    std::cout << mcells_map[mcell3].at(i)->Get_Center().x/units::cm << " " << mcells_map[mcell3].at(i)->Get_Center().y/units::cm << " " << mcells_map[mcell3].at(i)->Get_Center().z/units::cm << " " << std::endl;
  }
  
  std::cout << std::endl;
  MergeSpaceCell *mcell4 = 0;
  float max_ncells = 0;
  for (int i=0;i!=mcells.size();i++){
    if (fabs(mcells.at(i)->Get_Center().x - 163.84*units::cm)<1*units::mm){
      if (max_ncells < mcells.at(i)->Get_all_spacecell().size()){
	max_ncells = mcells.at(i)->Get_all_spacecell().size();
	mcell4 = mcells.at(i);
      }
      std::cout << mcells.at(i)->Get_Center().x/units::cm << " " << mcells.at(i)->Get_Center().y/units::cm << " " <<  mcells.at(i)->Get_Center().z/units::cm << " " << mcells.at(i)->Get_all_spacecell().size() <<  std::endl;
    } 
  }
  std::cout << mcell4->Get_Center().x/units::cm << " " << mcell4->Get_Center().y/units::cm << " " << mcell4->Get_Center().z/units::cm << " " << mcell4->Get_all_spacecell().size() << std::endl;

  auto it1 = find(mcells_map[mcell4].begin(),mcells_map[mcell4].end(),mcell3);
  auto it2 = find(mcells_map[mcell3].begin(),mcells_map[mcell3].end(),mcell4);
  std::cout << it1 - mcells_map[mcell4].end() << " " << it2 - mcells_map[mcell3].end() << " " << mcell3->Overlap(*mcell4) << " " << mcell4->Overlap(*mcell3) << " " << fabs(mcell4->Get_Center().x - mcell3->Get_Center().x) <<  std::endl;

  //WireCell2dToy::ToyCrawler toycrawler(mcells,1,2); //cosmic tune?

 //  // test
 //  std::cout << "Tracking " << std::endl;
 //  WireCell2dToy::ToyTracking toytracking(toycrawler);
 //  MergeSpaceCellMap& mcells_map = toycrawler.Get_mcells_map();
 //  //WireCell2dToy::ToyTracking toytracking(toycrawler,1); //cosmic tune?
 //  toytracking.IterateMergeTracks(mcells_map);
  
 // std:cout << mcells.size() << " " << mcells_map.size() << std::endl;
 //  // for (int i=0;i!=mcells.size();i++){
 //  //   if (mcells_map.find(mcells.at(i)) == mcells_map.end()){
 //  //     std::cout << i << std::endl;
 //  //   }else{
 //  //     std::cout << i << " " << mcells_map[mcells.at(i)].size() << std::endl;
 //  //   }
 //  // }
 //  // for (auto it = mcells_map.begin();it!=mcells_map.end();it++){
 //  //   auto it1 = find(mcells.begin(),mcells.end(),it->first);
 //  //   if (it1 == mcells.end()){
 //  //     std::cout << it->first->Get_all_spacecell().size() << " " 
 //  //   		<< it->first->Get_Center().x/units::cm << " " 
 //  //   		<< it->first->Get_Center().y/units::cm << " " 
 //  //   		<< it->first->Get_Center().z/units::cm << " " 
 //  // 		<< it->second.size() << std::endl;
 //  //   }
 //  // }


 //  std::cout << "Good Tracks:     " << toytracking.get_good_tracks().size() <<std::endl;
 //  std::cout << "Good Vertices:        " << toytracking.get_good_vertices().size() << std::endl;
 //  std::cout << "Bad Tracks:      " << toytracking.get_bad_tracks().size() << std::endl;
 //  std::cout << "Parallel Tracks: " << toytracking.get_parallel_tracks().size() << std::endl;
 //  std::cout << "Showers:         " << toytracking.get_showers().size() << std::endl;


 //  std::cout << "Drawing " << std::endl; 
 //  TApplication theApp("theApp",&argc,argv);
 //  theApp.SetReturnFromRun(true);
  
 //  TCanvas c1("ToyMC","ToyMC",800,600);
 //  c1.Draw();
  
 //  WireCell2dToy::ClusterDisplay display(c1);
 //  // display.DrawCluster(cells);
 //  display.DrawCluster(mcells);
 //  //display.DrawCluster(mcells,toytracking);
 //  //display.DrawCrawler(toycrawler,"psame",1);

 //  WCVertexSelection& vertices = toytracking.get_good_vertices();
 //  //WCVertexSelection& vertices = toytracking.get_vertices();
 //  display.DrawVertex(vertices,"psame");
  


 //  WCTrackSelection& bad_tracks = toytracking.get_bad_tracks();
 //  //display.DrawTracks(bad_tracks,"same",2);

 //  WCTrackSelection& short_tracks = toytracking.get_short_tracks();
 //  //display.DrawTracks(short_tracks,"psame",4);

 //  WCShowerSelection& showers =toytracking.get_showers();
 //  if (showers.size() > 0)
 //    display.DrawShower(showers.at(0),"psame",8);
 //  // Point p;
 //  // p.x = cells.at(0)->x();
 //  // p.y = cells.at(0)->y();
 //  // p.z = cells.at(0)->z();
 //  // display.DrawHough(cells,p,-1,10*units::m);
  

 //  theApp.Run();
 //  //std::cout << cells.size() << std::endl;
 //  //successfully read the TC tree 
 //  // TC->GetEntry(0);
 //  //cout << x << " " << y << " " << z << " " << charge << " " << time_slice << " " << cluster_num << " " << cell->cross_section() << " " << cell->center().y << endl;
    
    
}
