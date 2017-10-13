#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/ClusterDisplay.h"
#include "WireCell2dToy/ToyCrawler.h"
#include "WireCell2dToy/ToyTracking.h"
#include "WireCell2dToy/ToyCosmic.h"
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
   if (argc < 4) {
    cerr << "usage: wire-cell-allcluster /path/to/ChannelWireGeometry.txt /path/to/shower_3D.root eve_num" << endl;
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
  
  int eve_no = atoi(argv[3]);

  TFile *file = new TFile(filename);
  //TTree *T = (TTree*)file->Get("T");
  TTree *TC = (TTree*)file->Get("TC");
  TTree *Trun = (TTree*)file->Get("Trun");
  
  float unit_dis;
  int nrebin;
  int total_time_bin;
  Trun->SetBranchAddress("nrebin",&nrebin);
  Trun->SetBranchAddress("unit_dis",&unit_dis);
  Trun->SetBranchAddress("total_time_bin",&total_time_bin);
  int event_no, run_no, subrun_no;
  Trun->SetBranchAddress("eventNo",&event_no);
  Trun->SetBranchAddress("runNo",&run_no);
  Trun->SetBranchAddress("subRunNo",&subrun_no);
  Trun->GetEntry(0);

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

  //cout << T->GetEntries() << " " << TC->GetEntries() << endl;
  const int ntime = total_time_bin/nrebin;
  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[ntime];
  for (int i=0;i!=ntime;i++){
    toytiling[i] = new WireCell2dToy::ToyTiling();
  }
  int time_slice;
  //  WireCell2dToy::ToyTiling* tt = 0;
  
  // T->SetBranchAddress("time_slice",&time_slice);
  // T->SetBranchAddress("toytiling",&tt);

  //T->GetEntry(855);
  // cout << tt->get_allwire().size() << " " << tt->get_allcell().size() << endl;

  const GeomCell *cell = 0;//tt->get_allcell().at(5);
  //const GeomWire *wire = 0;//tt->get_allwire().at(3);
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



  //save all the crawlers ... 
  std::vector<WireCell2dToy::ToyCrawler*> crawlers;
  std::vector<WireCell2dToy::ToyTracking*> trackings;

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
	WireCell2dToy::ToyCrawler* toycrawler = new WireCell2dToy::ToyCrawler(mcells,1,2); // cosmic tune
	//toycrawler->FormGraph();

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
      mcell->set_id(mcell_id);
      flag = 1;
    }else if (flag==1 && (mcell_id!=prev_mcell_id)){
      mcells.push_back(mcell);
      mcell = new MergeSpaceCell();
      mcell->set_id(mcell_id);
    }
    GeomCell *cell1 = new GeomCell(cell);

    // if (time_slice == 1041){
      
    //   cout << "Single Cell: " << i << " "  << toytiling[time_slice]->get_allcell().size() << " " << toytiling[time_slice]->get_allwire().size() << endl;

      toytiling[time_slice]->AddCell(gds,cell1,u_index,v_index,w_index,u_charge,v_charge,w_charge,u_charge_err,v_charge_err,w_charge_err);
    
    //   cout << u_index << " " << v_index << " " << w_index << " " << u_charge << " " << v_charge << " " << w_charge << " " <<u_charge_err << " " << v_charge_err << " " << w_charge_err << endl;
      
    //   cout << "Single Cell1: " << i << " "  << toytiling[time_slice]->get_allcell().size() << " " << toytiling[time_slice]->get_allwire().size() << endl;
    
    // }
    

    SpaceCell *space_cell = new SpaceCell(cluster_num,*cell1,x*units::cm,charge,0.32*units::cm);
    mcell->AddSpaceCell(space_cell);
    cells.push_back(space_cell);
    
    prev_cluster_num = cluster_num;
    prev_mcell_id = mcell_id;
   
  }

  
  // for (int i=0;i!=ntime;i++){
  //   GeomCellSelection allcell = toytiling[i]->get_allcell();
  //   GeomWireSelection allwire = toytiling[i]->get_allwire();
  //   cout << "Single Cell: " << i << " "  << allcell.size() << " " << allwire.size() << endl;
  // }


  mcells.push_back(mcell);  
  WireCell2dToy::ToyCrawler* toycrawler = new WireCell2dToy::ToyCrawler(mcells);
  toycrawler->FormGraph();
  crawlers.push_back(toycrawler);






  //check # of clusters 
  int sum = 0 ;
  
  for (int i=0;i!=crawlers.size();i++){
    WireCell2dToy::ToyTracking* toytracking = new WireCell2dToy::ToyTracking(*crawlers.at(i),1); // cosmic tune
    trackings.push_back(toytracking);

    std::cout << "Cluster          " << i << std::endl;
    std::cout << "Good Tracks:     " << toytracking->get_good_tracks().size() <<std::endl;
    std::cout << "Vertices:        " << toytracking->get_vertices().size() << std::endl;
    std::cout << "Bad Tracks:      " << toytracking->get_bad_tracks().size() << std::endl;
    std::cout << "Parallel Tracks: " << toytracking->get_parallel_tracks().size() << std::endl;
    std::cout << "Showers:         " << toytracking->get_showers().size() << std::endl;

    for (auto it = crawlers.at(i)->Get_mcells_map().begin(); it!= crawlers.at(i)->Get_mcells_map().end();it++){
      MergeSpaceCell *mcell1 = it->first;
      sum += mcell1->Get_all_spacecell().size();
      //    sum += .size();
    }
  }
  
  std::cout << "Check: " << crawlers.size() << " " << TC->GetEntries() << " " << sum << std::endl;

  
  WireCell2dToy::ToyCosmic toycosmic(trackings);
  
  // //Check tracking ... 
  // for (int i=0;i!=trackings.size();i++){
  // }



  TFile *file1 = new TFile(Form("./rootfiles/cosmic_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
  TTree *T1 = new TTree("T_goodtrack","T_goodtrack");
  TTree *T2 = new TTree("T_vertex","T_vertex");
  TTree *T3 = new TTree("T_badtrack","T_badtrack");
  TTree *T4 = new TTree("T_shorttrack","T_shortrack");
  TTree *T5 = new TTree("T_paratrack","T_paratrack");
  TTree *T6 = new TTree("T_shower","T_shower");

  TTree *T7 = new TTree("T_cosmic","T_cosmic");
  TTree *T8 = new TTree("T_neutrino","T_neutrino");

  T1->SetDirectory(file1);
  T2->SetDirectory(file1);
  T3->SetDirectory(file1);
  T4->SetDirectory(file1);
  T5->SetDirectory(file1);
  T6->SetDirectory(file1);
  T7->SetDirectory(file1);
  T8->SetDirectory(file1);
  
  

  Int_t ntracks;
  Int_t nshowers;
  Int_t npoints;
  Double_t xx[10000],yy[10000],zz[10000];
  Double_t theta[10000],phi[10000];
  Double_t energy[10000],dedx[10000];
  Int_t msc_id[10000];
  Int_t trackid;
  Int_t vtrack_id[100];
  Int_t showerid;
  Int_t vshower_id[100];
  
  T1->Branch("npoints",&npoints,"npoints/I");
  T1->Branch("trackid",&trackid,"trackid/I");
  T1->Branch("x",xx,"x[npoints]/D");
  T1->Branch("y",yy,"y[npoints]/D");
  T1->Branch("z",zz,"z[npoints]/D");
  T1->Branch("theta",theta,"theta[npoints]/D");
  T1->Branch("phi",phi,"phi[npoints]/D");
  T1->Branch("energy",energy,"energy[npoints]/D");
  T1->Branch("dedx",dedx,"dedx[npoints]/D");
  T1->Branch("msc_id",msc_id,"msc_id[npoints]/I");


  T2->Branch("x",xx,"x/D");
  T2->Branch("y",yy,"y/D");
  T2->Branch("z",zz,"z/D");
  T2->Branch("ntracks",&ntracks,"ntracks/I");
  T2->Branch("vtrack_id",vtrack_id,"vtrack_id[ntracks]/I");
  T2->Branch("nshowers",&nshowers,"nshowers/I");
  T2->Branch("vshower_id",vshower_id,"vshower_id[nshowers]/I");
  

  T3->Branch("npoints",&npoints,"npoints/I");
  T3->Branch("trackid",&trackid,"trackid/I");
  T3->Branch("x",xx,"x[npoints]/D");
  T3->Branch("y",yy,"y[npoints]/D");
  T3->Branch("z",zz,"z[npoints]/D");
  T3->Branch("theta",theta,"theta[npoints]/D");
  T3->Branch("phi",phi,"phi[npoints]/D");
  T3->Branch("energy",energy,"energy[npoints]/D");
  T3->Branch("dedx",dedx,"dedx[npoints]/D");
  T3->Branch("msc_id",msc_id,"msc_id[npoints]/I");

  T4->Branch("npoints",&npoints,"npoints/I");
  T4->Branch("trackid",&trackid,"trackid/I");
  T4->Branch("x",xx,"x[npoints]/D");
  T4->Branch("y",yy,"y[npoints]/D");
  T4->Branch("z",zz,"z[npoints]/D");
  T4->Branch("theta",theta,"theta[npoints]/D");
  T4->Branch("phi",phi,"phi[npoints]/D");
  T4->Branch("energy",energy,"energy[npoints]/D");
  T4->Branch("dedx",dedx,"dedx[npoints]/D");
  T4->Branch("msc_id",msc_id,"msc_id[npoints]/I");

  T5->Branch("npoints",&npoints,"npoints/I");
  T5->Branch("trackid",&trackid,"trackid/I");
  T5->Branch("x",xx,"x[npoints]/D");
  T5->Branch("y",yy,"y[npoints]/D");
  T5->Branch("z",zz,"z[npoints]/D");
  T5->Branch("theta",theta,"theta[npoints]/D");
  T5->Branch("phi",phi,"phi[npoints]/D");
  T5->Branch("energy",energy,"energy[npoints]/D");
  T5->Branch("dedx",dedx,"dedx[npoints]/D");
  T5->Branch("msc_id",msc_id,"msc_id[npoints]/I");
  
  T6->Branch("showerid",&showerid,"showerid/I");
  T6->Branch("npoints",&npoints,"npoints/I");
  T6->Branch("energy",energy,"energy[npoints]/D");
  T6->Branch("msc_id",msc_id,"msc_id[npoints]/I");
  T6->Branch("vertex_x",xx,"vertex_x/D");
  T6->Branch("vertex_y",yy,"vertex_y/D");
  T6->Branch("vertex_z",zz,"vertex_z/D");
  
  int cosmic_flag;
  T7->Branch("cosmic_flag",&cosmic_flag,"cosmic_flag/I");
  T7->Branch("trackid",&trackid,"trackid/I");
  T7->Branch("npoints",&npoints,"npoints/I");
  T7->Branch("x",xx,"x[npoints]/D");
  T7->Branch("y",yy,"y[npoints]/D");
  T7->Branch("z",zz,"z[npoints]/D");
  
  double length;
  T8->Branch("npoints",&npoints,"npoints/I");
  T8->Branch("trackid",&trackid,"trackid/I");
  T8->Branch("length",&length,"length/D");
  T8->Branch("x",xx,"x[npoints]/D");
  T8->Branch("y",yy,"y[npoints]/D");
  T8->Branch("z",zz,"z[npoints]/D");
  T8->Branch("theta",theta,"theta[npoints]/D");
  T8->Branch("phi",phi,"phi[npoints]/D");
  T8->Branch("energy",energy,"energy[npoints]/D");
  T8->Branch("dedx",dedx,"dedx[npoints]/D");
  T8->Branch("msc_id",msc_id,"msc_id[npoints]/I");
  


  WCTrackSelection all_tracks;
  WCTrackSelection good_tracks;
  WCTrackSelection bad_tracks;
  WCTrackSelection short_tracks;
  WCTrackSelection parallel_tracks;
  WCVertexSelection vertices;
  WCShowerSelection showers;
  
  
  TGraph2D *g = new TGraph2D();
  int ncount = 0;

  for (int i=0;i!=trackings.size();i++){
    for (int j=0;j!=trackings.at(i)->get_good_tracks().size();j++){
      good_tracks.push_back(trackings.at(i)->get_good_tracks().at(j));
      all_tracks.push_back(trackings.at(i)->get_good_tracks().at(j));
      // std::cout <<"Track: " << i << " " << j << " " << good_tracks.size() - 1<< std::endl;
    }
    for (int j=0;j!=trackings.at(i)->get_bad_tracks().size();j++){
      bad_tracks.push_back(trackings.at(i)->get_bad_tracks().at(j));
      all_tracks.push_back(trackings.at(i)->get_bad_tracks().at(j));
    }
    for (int j=0;j!=trackings.at(i)->get_short_tracks().size();j++){
      short_tracks.push_back(trackings.at(i)->get_short_tracks().at(j));
      all_tracks.push_back(trackings.at(i)->get_short_tracks().at(j));
    }
    for (int j=0;j!=trackings.at(i)->get_parallel_tracks().size();j++){
      parallel_tracks.push_back(trackings.at(i)->get_parallel_tracks().at(j));
      all_tracks.push_back(trackings.at(i)->get_parallel_tracks().at(j));
    }

    for (int j=0;j!=trackings.at(i)->get_vertices().size();j++){
      vertices.push_back(trackings.at(i)->get_vertices().at(j));
      // std::cout << i << " " << j << " " << trackings.at(i)->get_vertices().at(j)->Center().x/units::cm << " " <<
      // 	trackings.at(i)->get_vertices().at(j)->Center().y/units::cm << " " <<
      // 	trackings.at(i)->get_vertices().at(j)->Center().z/units::cm << " " << 
      // 	std::endl;
    }
    for (int j=0;j!=trackings.at(i)->get_showers().size();j++){
      showers.push_back(trackings.at(i)->get_showers().at(j));
    }
  }
  
  
  //fill T1
  for (int i = 0; i!=good_tracks.size();i++){
    WCTrack *track = good_tracks.at(i);
    
    npoints = track->get_centerVP().size();
    trackid = find(all_tracks.begin(),all_tracks.end(),track) - all_tracks.begin();
    for (int j=0;j!=npoints;j++){
      xx[j] = track->get_centerVP().at(j).x/units::cm;
      yy[j] = track->get_centerVP().at(j).y/units::cm;
      zz[j] = track->get_centerVP().at(j).z/units::cm;
      theta[j] = track->get_centerVP_theta().at(j);
      phi[j] = track->get_centerVP_phi().at(j);
      energy[j] = track->get_centerVP_energy().at(j);
      dedx[j] = track->get_centerVP_dedx().at(j);
      msc_id[j] = track->get_centerVP_cells().at(j)->get_id();
      g->SetPoint(ncount,xx[j],yy[j],zz[j]);
      ncount ++;
    }
    T1->Fill();
  }

  //fill T2  //Vertex
  for (int i=0;i!=vertices.size();i++){
    ntracks = 0;
    WCVertex *vertex = vertices.at(i);
    xx[0] = vertex->Center().x/units::cm;
    yy[0] = vertex->Center().y/units::cm;
    zz[0] = vertex->Center().z/units::cm;
    for (int j=0;j!=vertex->get_ntracks();j++){
      WCTrack *track = vertex->get_tracks().at(j);
      auto it = find(all_tracks.begin(),all_tracks.end(),track);
      if (it != all_tracks.end()){
	vtrack_id[ntracks] = it-all_tracks.begin();
	ntracks ++;
      }
    }
    nshowers = 0;
    for (int j=0;j!=showers.size();j++){
      if (vertex == showers.at(j)->get_vertex()){
	vshower_id[nshowers] = j;
	nshowers++;
      }
    }

    T2->Fill();
  }

  //fill T3
  for (int i = 0; i!=bad_tracks.size();i++){
    WCTrack *track = bad_tracks.at(i);
    
    npoints = track->get_centerVP().size();
    trackid = find(all_tracks.begin(),all_tracks.end(),track) - all_tracks.begin();
    for (int j=0;j!=npoints;j++){
      xx[j] = track->get_centerVP().at(j).x/units::cm;
      yy[j] = track->get_centerVP().at(j).y/units::cm;
      zz[j] = track->get_centerVP().at(j).z/units::cm;
      theta[j] = track->get_centerVP_theta().at(j);
      phi[j] = track->get_centerVP_phi().at(j);
      energy[j] = track->get_centerVP_energy().at(j);
      dedx[j] = track->get_centerVP_dedx().at(j);
      msc_id[j] = track->get_centerVP_cells().at(j)->get_id();
      g->SetPoint(ncount,xx[j],yy[j],zz[j]);
      ncount ++;
    }
    T3->Fill();
  }
  
  //fill T4
  for (int i = 0; i!=short_tracks.size();i++){
    WCTrack *track = short_tracks.at(i);
    npoints = track->get_all_cells().size();
    trackid = find(all_tracks.begin(),all_tracks.end(),track) - all_tracks.begin();
    for (int j=0;j!=npoints;j++){
      xx[j] = track->get_all_cells().at(j)->Get_Center().x/units::cm;
      yy[j] = track->get_all_cells().at(j)->Get_Center().y/units::cm;
      zz[j] = track->get_all_cells().at(j)->Get_Center().z/units::cm;
      theta[j] = 0;
      phi[j] = 0;
      energy[j] = track->get_all_cells().at(j)->Get_Charge();
      dedx[j] = 0;
      msc_id[j] = track->get_all_cells().at(j)->get_id();
      g->SetPoint(ncount,xx[j],yy[j],zz[j]);
      ncount ++;
    }
    T4->Fill();
  }


  //fill T5
  for (int i = 0; i!=parallel_tracks.size();i++){
    WCTrack *track = parallel_tracks.at(i);
    
    npoints = track->get_centerVP().size();
    trackid = find(all_tracks.begin(),all_tracks.end(),track) - all_tracks.begin();
    for (int j=0;j!=npoints;j++){
      xx[j] = track->get_centerVP().at(j).x/units::cm;
      yy[j] = track->get_centerVP().at(j).y/units::cm;
      zz[j] = track->get_centerVP().at(j).z/units::cm;
      theta[j] = track->get_centerVP_theta().at(j);
      phi[j] = track->get_centerVP_phi().at(j);
      energy[j] = track->get_centerVP_energy().at(j);
      dedx[j] = track->get_centerVP_dedx().at(j);
      msc_id[j] = track->get_centerVP_cells().at(j)->get_id();
      g->SetPoint(ncount,xx[j],yy[j],zz[j]);
      ncount ++;
    }
    T5->Fill();
  }

  //fill T6
  for (int i=0;i!=showers.size();i++){
    WCShower *shower = showers.at(i);
    npoints = shower->get_all_cells().size();
    showerid = find(showers.begin(),showers.end(),shower) - showers.begin();
    xx[0] = shower->get_vertex()->Center().x/units::cm;
    yy[0] = shower->get_vertex()->Center().y/units::cm;
    zz[0] = shower->get_vertex()->Center().z/units::cm;
    for (int j=0;j!=npoints;j++){
      g->SetPoint(ncount,
		  shower->get_all_cells().at(j)->Get_Center().x/units::cm,
		  shower->get_all_cells().at(j)->Get_Center().y/units::cm,
		  shower->get_all_cells().at(j)->Get_Center().z/units::cm);
      ncount ++;
      energy[j] = shower->get_all_cells().at(j)->Get_Charge();
      msc_id[j] = shower->get_all_cells().at(j)->get_id();
    }
    T6->Fill();
  }
  
  WireCell2dToy::WCCosmicSelection& cosmics = toycosmic.get_cosmics();
  for (int i=0;i!=cosmics.size();i++){
    trackid = i;
    npoints = 0;
    WireCell2dToy::WCCosmic *cosmic = cosmics.at(i);
    if (cosmic->IsCosmic()){
      cosmic_flag = 1;
    }else{
      cosmic_flag = 0;
    }
    for (int j=0;j!=cosmic->get_points().size();j++){
      xx[npoints] = cosmic->get_points().at(j).x/units::cm;
      yy[npoints] = cosmic->get_points().at(j).y/units::cm;
      zz[npoints] = cosmic->get_points().at(j).z/units::cm;
      npoints ++;

    }
    T7->Fill();
  }


  //fill T8
  WCTrackSelection neutrino_tracks;
  for (int i=0;i!=toycosmic.get_neutrinos().size();i++){
    for(int j=0;j!=toycosmic.get_neutrinos().at(i)->get_good_tracks().size();j++){
      neutrino_tracks.push_back(toycosmic.get_neutrinos().at(i)->get_good_tracks().at(j));
    }
  }

  for (int i = 0; i!=neutrino_tracks.size();i++){
    WCTrack *track = neutrino_tracks.at(i);
    
    npoints = track->get_centerVP().size();
    trackid = find(all_tracks.begin(),all_tracks.end(),track) - all_tracks.begin();
    length = track-> get_range()/units::cm;
    for (int j=0;j!=npoints;j++){
      xx[j] = track->get_centerVP().at(j).x/units::cm;
      yy[j] = track->get_centerVP().at(j).y/units::cm;
      zz[j] = track->get_centerVP().at(j).z/units::cm;
      theta[j] = track->get_centerVP_theta().at(j);
      phi[j] = track->get_centerVP_phi().at(j);
      energy[j] = track->get_centerVP_energy().at(j);
      dedx[j] = track->get_centerVP_dedx().at(j);
      msc_id[j] = track->get_centerVP_cells().at(j)->get_id();
      g->SetPoint(ncount,xx[j],yy[j],zz[j]);
      ncount ++;
    }
    T8->Fill();
  }


  // for (int i=0;i!=cosmics.size();i++){
  //   trackid = i;
  //   npoints = 0;
  //   WireCell2dToy::WCCosmic *cosmic = cosmics.at(i);
  //   for (int j=0;j!=cosmic->get_mcells().size();j++){
  //     xx[npoints] = cosmic->get_mcells().at(j)->Get_Center().x/units::cm;
  //     yy[npoints] = cosmic->get_mcells().at(j)->Get_Center().y/units::cm;
  //     zz[npoints] = cosmic->get_mcells().at(j)->Get_Center().z/units::cm;
  //     npoints ++;
  //   }
  //   T7->Fill();
  // }
  // std::vector<WireCell2dToy::ToyTrackingSelection>& cosmics = toycosmic.get_raw_candidates();
  // //int sum1 = 0;
  // for (int i=0;i!=cosmics.size();i++){
  //   trackid = i;
  //   npoints = 0;
  //   for (int j=0;j!=cosmics.at(i).size();j++){
  //     WireCell2dToy::ToyTracking *tracking = cosmics.at(i).at(j);
  //     // sum1 ++;
  //     WCTrackSelection tracking_tracks = tracking->get_good_tracks();
  //     for (int k=0;k!=tracking_tracks.size();k++){
  // 	WCTrack *track = tracking_tracks.at(k);
  // 	for (int k1=0;k1!= track->get_centerVP_cells().size();k1++){
  // 	  xx[npoints] = track->get_centerVP_cells().at(k1)->Get_Center().x/units::cm;
  // 	  yy[npoints] = track->get_centerVP_cells().at(k1)->Get_Center().y/units::cm;
  // 	  zz[npoints] = track->get_centerVP_cells().at(k1)->Get_Center().z/units::cm;
  // 	  npoints ++;
  // 	}
  //     }
  //   }
  //   T7->Fill();
  // }

  

  
  //std::cout << sum1 << std::endl;

  g->Write("shower3D");

  TC->CloneTree()->Write();
  Trun->CloneTree()->Write();

  file1->Write();
  file1->Close();
  
  



  //start the prepare the important merge cell vectors
  // int start_num = 0 ;
  // int end_num = 2399;
  // std::vector<GeomCellSelection> Good_MCells;
  
  // for (int i=start_num;i!=end_num+1;i++){
  // GeomCellSelection cells;
  //Good_MCells.push_back(cells);
  // }

  // MergeSpaceCellSelection ms_cells;
  
  
  // for (int i=0;i!=crawlers.size();i++){
  //   WireCell2dToy::ToyCrawler *toycrawler = crawlers.at(i);
  //   for (int j=0;j!=toycrawler->Get_allMCT().size();j++){
  //     MergeClusterTrack *mct = toycrawler->Get_allMCT().at(j);
  //     int ntime = mct->Get_TimeLength();
  //     if (ntime >=5){
  // 	// do something
	
  // 	for (int k=0;k!=ntime;k++){
  // 	  MergeSpaceCellSelection cells = mct->Get_MSCS(k);
  // 	  if (cells.size()==1){
  // 	    ms_cells.push_back(cells.at(0));
  // 	  }else if (cells.size()>1){
  // 	    MergeSpaceCell *cell = cells.at(0);
  // 	    for (int i1 = 1; i1!=cells.size();i1++){
  // 	      if (cell->Get_all_spacecell().size() < cells.at(i1)->Get_all_spacecell().size()){
  // 		cell = cells.at(i1);
  // 	      }
  // 	    }
  // 	    ms_cells.push_back(cell);
  // 	  }
  // 	}
	

  //     }
  //   }
  // }
  
  // // plot it
  
  // TApplication theApp("theApp",&argc,argv);
  // theApp.SetReturnFromRun(true);
  
  // TCanvas c1("ToyMC","ToyMC",800,600);
  // c1.Draw();
  
  // WireCell2dToy::ClusterDisplay display(c1);
  // display.DrawCluster(ms_cells);

  // theApp.Run();
}
