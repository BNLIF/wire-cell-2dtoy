#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/ClusterDisplay.h"
#include "WireCell2dToy/ToyCrawler.h"
#include "WireCell2dToy/ToyTracking.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"

#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"
#include "WireCellData/SpaceCell.h"
#include "WireCellData/MergeSpaceCell.h"

#include "WireCellSst/MCTruth.h"

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

  if (argc < 5){
    cerr << "usage: wire-cell-dune-det-op /path/to/celltree.root cluster_num -o[0,1](do_rotation) -p[0,1](is_3mm) " << endl;
    return 1;
  }

  bool rotate_90deg=false;
  bool is_3mm=false;
  bool random_vertices=false;
  int seed=0;
  for (Int_t i = 1; i != argc; i++){
     switch(argv[i][1]){
     case 'o':
       rotate_90deg = atoi(&argv[i][2]); 
       break;
     case 'p':
       is_3mm = atoi(&argv[i][2]); 
       break;
     case 'r':
       random_vertices = atoi(&argv[i][2]);
       break;
     case 's':
       seed = atoi(&argv[i][2]);
       break;
     }
  }


  if (rotate_90deg) cout<<"Beam is perpendicular to wire planes. ";
  else cout<<"Beam is parallel to wire planes (default). ";
  if (is_3mm) cout<<"Wire pitch is 3 mm."<<endl;
  else cout<<"Wire pitch is 5 mm (default)."<<endl;
  
  double pitchU, pitchV, pitchW;
  if (!is_3mm) {
    pitchU = 0.4667*units::cm;
    pitchV = 0.4667*units::cm;
    pitchW = 0.479*units::cm;
  } else {
    pitchU = 0.3*units::cm;
    pitchV = 0.3*units::cm;
    pitchW = 0.3*units::cm;
  }
  
  //build GDS ... 
  DetectorGDS gds;
  gds.set_ncryos(1);
  gds.set_napas(0,4);
  Vector center0(-0*units::cm, -300.05*units::cm, 115.318875*units::cm);
  Vector halves0(3.995355*units::cm, 299.5*units::cm, 115.318875*units::cm);
  gds.set_apa(0, 0, 35.71*units::deg, 35.71*units::deg, pitchU, pitchV, pitchW, center0, halves0);
  Vector center1(-0*units::cm, 300.05*units::cm, 115.318875*units::cm);
  Vector halves1(3.995355*units::cm, 299.5*units::cm, 115.318875*units::cm);
  gds.set_apa(0, 1, 35.71*units::deg, 35.71*units::deg, pitchU, pitchV, pitchW, center1, halves1);
  Vector center2(-0*units::cm, -300.05*units::cm, 347.709*units::cm);
  Vector halves2(3.995355*units::cm, 299.5*units::cm, 115.318875*units::cm);
  gds.set_apa(0, 2, 35.71*units::deg, 35.71*units::deg, pitchU, pitchV, pitchW, center2, halves2);
  Vector center3(-0*units::cm, 300.05*units::cm, 347.709*units::cm);
  Vector halves3(3.995355*units::cm, 299.5*units::cm, 115.318875*units::cm);
  gds.set_apa(0, 3, 35.71*units::deg, 35.71*units::deg, pitchU, pitchV, pitchW, center3, halves3);
  gds.buildGDS();

  
  TString filename = argv[1];
  int ncluster = atoi(argv[2]);


  TFile *file = new TFile(filename);
  TTree *T = (TTree*)file->Get("T");
  TTree *TC = (TTree*)file->Get("TC");
  TTree *Trun = (TTree*)file->Get("Trun");
  TTree *TMC = (TTree*)file->Get("TMC");
  
  TGraph2D *shower3D = (TGraph2D*)file->Get("shower3D");

  float unit_dis;
  int nrebin;
  int total_time_bin;
  Trun->SetBranchAddress("nrebin",&nrebin);
  Trun->SetBranchAddress("unit_dis",&unit_dis);
  Trun->SetBranchAddress("total_time_bin",&total_time_bin);
  Trun->GetEntry(0);
  
  TPCParams& mp = Singleton<TPCParams>::Instance();
  
  double pitch_u = pitchU;
  double pitch_v = pitchV;
  double pitch_w = pitchW;
  double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;

  mp.set_pitch_u(pitch_u);
  mp.set_pitch_v(pitch_v);
  mp.set_pitch_w(pitch_w);
  mp.set_ts_width(time_slice_width);
  
  std::cout << "Singleton: " << mp.get_pitch_u() << " " << mp.get_pitch_v() << " " << mp.get_pitch_w() << " " << mp.get_ts_width() << std::endl;

  const int ntime = total_time_bin/nrebin;
  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[ntime];
  for (int i=0;i!=ntime;i++){
    toytiling[i] = new WireCell2dToy::ToyTiling();
  }

  int time_slice;

  const GeomCell *cell = 0;//tt->get_allcell().at(5);
  const GeomWire *wire = 0;//tt->get_allwire().at(3);
  
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

  int apa_no=0, cryostat_no=0;
  int face = 0;
  TC->SetBranchAddress("face",&face);
  TC->SetBranchAddress("apa_no",&apa_no);
  TC->SetBranchAddress("cryostat_no",&cryostat_no);
    
  
  int prev_mcell_id = -1;
  int prev_cluster_num = -1;

  MergeSpaceCellSelection mcells; // save all the cells
  MergeSpaceCellSelection mcells_all; // save all the cells
  int flag = 0;
  MergeSpaceCell *mcell=0;
  SpaceCellSelection cells;
  
  for (int i=0;i!=TC->GetEntries();i++){
    TC->GetEntry(i);
    
    if (cluster_num != prev_cluster_num){
      if (prev_cluster_num!=-1){
	if (mcell->Get_all_spacecell().size()>0){
	  mcells.push_back(mcell);  
	  mcells_all.push_back(mcell);
	}
	// this is a cluster
	mcells.clear();
	cells.clear();
	flag = 0;
      }
    }
    
    if (flag == 0){
      mcell = new MergeSpaceCell();
      mcell->set_id(mcell_id);
      flag = 1;
    }else if (flag==1 && (mcell_id!=prev_mcell_id)){
      mcells.push_back(mcell);
      mcells_all.push_back(mcell);
      mcell = new MergeSpaceCell();
      mcell->set_id(mcell_id);
    }
    GeomCell *cell1 = new GeomCell(cell);

    toytiling[time_slice]->AddCell(gds,cryostat_no,apa_no,cell1,u_index,v_index,w_index,u_charge,v_charge,w_charge,u_charge_err,v_charge_err,w_charge_err);

    SpaceCell *space_cell = new SpaceCell(cluster_num,*cell1,x*units::cm,charge,0.32*units::cm);
    mcell->AddSpaceCell(space_cell);
    cells.push_back(space_cell);
    
    prev_cluster_num = cluster_num;
    prev_mcell_id = mcell_id;

  }

  // last cluster ... 
  if (mcell!=0){
    if (mcell->Get_all_spacecell().size()>0){
      mcells.push_back(mcell);  
      mcells_all.push_back(mcell);
    }
  }
   
  
  // for (int i=0;i!=ntime;i++){
  //   GeomCellSelection allcell = toytiling[i]->get_allcell();
  //   GeomWireSelection allwire = toytiling[i]->get_allwire();
  //   cout << "Single Cell: " << i << " "  << allcell.size() << " " << allwire.size() << endl;
  // }
  
  
  TGraph2D *grec = new TGraph2D();
  int num = 0;
  for (int i=0;i!=mcells_all.size();i++){
    // for (int j=0;j!=mcells_all.at(i)->Get_all_spacecell().size();j++){
    //   grec->SetPoint(num,mcells_all.at(i)->Get_all_spacecell().at(j)->x()/units::cm,
    // 		     mcells_all.at(i)->Get_all_spacecell().at(j)->y()/units::cm,
    // 		     mcells_all.at(i)->Get_all_spacecell().at(j)->z()/units::cm
    // 		     );
    //   num ++;
    // }
    grec->SetPoint(i,mcells_all.at(i)->Get_Center().x/units::cm,
    		    mcells_all.at(i)->Get_Center().y/units::cm,
    		    mcells_all.at(i)->Get_Center().z/units::cm);
    
		    
  }

  std::cout << TC->GetEntries() << " " << num << std::endl;

  
  // deal with MC truth ... 
  WireCellSst::MCTruth *mctruth = new WireCellSst::MCTruth(TMC);
  //mctruth->GetEntry(0);

  // Need the (original) primary vertex to be in the fiducial volume ...
  Point neutrino_vertex = mctruth->find_neutrino_vertex();
  bool contained_flag = false;
  if (fabs(neutrino_vertex[0]) < 360-10 && fabs(neutrino_vertex[1]) < 600-10 && neutrino_vertex[2]>10 && neutrino_vertex[2] < 360-20){
    contained_flag = true;
  }
  std::cout << "neutrino Vertex: " << neutrino_vertex[0] << " " << neutrino_vertex[1] << " " << neutrino_vertex[2] << " " << contained_flag << std::endl;

  // Now need to find the primary electron trajectory ...
  Point primary_vertex = mctruth->find_primary_vertex();
  std::cout << "Primary Vertex: " << primary_vertex[0] << " " << primary_vertex[1] << " " << primary_vertex[2] << std::endl;
  
  // use the vertex and the main point
  // use the direction of the electron
  MCParticle* electron = mctruth->find_primary_electron();
  if (electron != 0){
    std::cout << electron->startMomentum[0] << " " << electron->startMomentum[1] << " " << electron->startMomentum[2] << std::endl;
  }

  // Need to find all the seconary gammas' (energy deposition) starting point and trajectory... 
  // position, the energy of secondary electron etc
  // direction is just gamma's direction?
  MCParticleSelection photons = mctruth->find_primary_photons();
  std::cout << photons.size() << std::endl;
  
  std::cout << mctruth->find_neutrino_true_energy() << " " << mctruth->find_neutrino_visible_energy() << std::endl;
  
  // save wire and cell directly? (in addition to the point)
  // gap can be defined as a line between the position and vertex location 

  // Given a position, need to find the corresponding merged cells' location (through merged wire)... 
  // can also cut distance betweem wires ... 
  
  // Given a merged cell, need to find all the wires and then energy ... 
  // dE/dx through projection? 
  
  
  TApplication theApp("theApp",&argc,argv);
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  shower3D->Draw("p");
  grec->SetMarkerColor(2);
  grec->Draw("psame");
  grec->SetMarkerStyle(21);
  grec->SetMarkerSize(0.5);
  theApp.Run();


  // //cout << mcells.size() << endl;

  // // do the Toy Crawler
  // std::cout << "Crawling " << std::endl;
  // WireCell2dToy::ToyCrawler toycrawler(mcells);
  // //WireCell2dToy::ToyCrawler toycrawler(mcells,1,2); //cosmic tune?

  // // test
  // std::cout << "Tracking " << std::endl;
  // WireCell2dToy::ToyTracking toytracking(toycrawler);
  // MergeSpaceCellMap& mcells_map = toycrawler.Get_mcells_map();
  // //WireCell2dToy::ToyTracking toytracking(toycrawler,1); //cosmic tune?
  // toytracking.IterateMergeTracks(mcells_map);
  
  // //std:cout << mcells.size() << " " << mcells_map.size() << std::endl;
  // // for (int i=0;i!=mcells.size();i++){
  // //   if (mcells_map.find(mcells.at(i)) == mcells_map.end()){
  // //     std::cout << i << std::endl;
  // //   }else{
  // //     std::cout << i << " " << mcells_map[mcells.at(i)].size() << std::endl;
  // //   }
  // // }
  // // for (auto it = mcells_map.begin();it!=mcells_map.end();it++){
  // //   auto it1 = find(mcells.begin(),mcells.end(),it->first);
  // //   if (it1 == mcells.end()){
  // //     std::cout << it->first->Get_all_spacecell().size() << " " 
  // //   		<< it->first->Get_Center().x/units::cm << " " 
  // //   		<< it->first->Get_Center().y/units::cm << " " 
  // //   		<< it->first->Get_Center().z/units::cm << " " 
  // // 		<< it->second.size() << std::endl;
  // //   }
  // // }


  // std::cout << "Good Tracks:     " << toytracking.get_good_tracks().size() <<std::endl;
  // std::cout << "Good Vertices:        " << toytracking.get_good_vertices().size() << std::endl;
  // std::cout << "Bad Tracks:      " << toytracking.get_bad_tracks().size() << std::endl;
  // std::cout << "Parallel Tracks: " << toytracking.get_parallel_tracks().size() << std::endl;
  // std::cout << "Showers:         " << toytracking.get_showers().size() << std::endl;


  // std::cout << "Drawing " << std::endl; 
  // TApplication theApp("theApp",&argc,argv);
  // theApp.SetReturnFromRun(true);
  
  // TCanvas c1("ToyMC","ToyMC",800,600);
  // c1.Draw();
  
  // WireCell2dToy::ClusterDisplay display(c1);
  // shower3D_charge->Draw("p");

  // // display.DrawCluster(cells);
  // // display.DrawCluster(mcells);
  // //display.DrawCluster(mcells,toytracking);
  // //display.DrawCrawler(toycrawler,"psame",1);

  // WCVertexSelection& vertices = toytracking.get_good_vertices();
  // //WCVertexSelection& vertices = toytracking.get_vertices();
  // display.DrawVertex(vertices,"psame");
  
  // WCTrackSelection& bad_tracks = toytracking.get_bad_tracks();
  // //display.DrawTracks(bad_tracks,"same",2);

  // WCTrackSelection& short_tracks = toytracking.get_short_tracks();
  // //display.DrawTracks(short_tracks,"psame",4);

  // WCShowerSelection& showers =toytracking.get_showers();
  // if (showers.size() > 0)
  //   display.DrawShower(showers.at(0),"psame",8);
  // // Point p;
  // // p.x = cells.at(0)->x();
  // // p.y = cells.at(0)->y();
  // // p.z = cells.at(0)->z();
  // // display.DrawHough(cells,p,-1,10*units::m);
  

  // theApp.Run();
  // //std::cout << cells.size() << std::endl;
  // //successfully read the TC tree 
  // // TC->GetEntry(0);
  // //cout << x << " " << y << " " << z << " " << charge << " " << time_slice << " " << cluster_num << " " << cell->cross_section() << " " << cell->center().y << endl;
    
    
}
