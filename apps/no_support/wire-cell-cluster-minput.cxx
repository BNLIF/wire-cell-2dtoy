#include "WireCellSst/GeomDataSource.h"
#include "WireCellData/PR3DCluster.h"
#include "WireCellData/SlimMergeGeomCell.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"
#include "WireCell2dToy/ExecMon.h"
#include "WireCell2dToy/CalcPoints.h"
#include "WireCell2dToy/ToyClustering.h"


#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace WireCell;
using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 4) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root entry_num (input: wire-cell-imaging-lmem-celltree -s2 mode output)" << endl;
    return 1;
  }

  ExecMon em("starting");
  //cerr << em("load geometry") << endl;

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



  // test geometry ...
  const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),0);
  const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),0);
  const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),0);
  double first_u_dis = gds.wire_dist(*uwire) ;
  double first_v_dis = gds.wire_dist(*vwire) ;
  double first_w_dis = gds.wire_dist(*wwire) ; 
  
  
  TString filename = argv[2];
  TFile *file = new TFile(filename);
  TTree *Trun = (TTree*)file->Get("Trun");

  int run_no, subrun_no, event_no;
  int time_offset;
  int nrebin;
  int frame_length;
  int eve_num;
  int entry_num = atoi(argv[3]); 
  float unit_dis;
  
  Trun->SetBranchAddress("eventNo",&event_no);
  Trun->SetBranchAddress("runNo",&run_no);
  Trun->SetBranchAddress("subRunNo",&subrun_no);
  Trun->SetBranchAddress("unit_dis",&unit_dis);
  Trun->SetBranchAddress("frame_length",&frame_length);
  Trun->SetBranchAddress("eve_num",&eve_num);
  Trun->SetBranchAddress("nrebin",&nrebin);
  Trun->SetBranchAddress("time_offset",&time_offset);
  Trun->GetEntry(entry_num);

  // define singleton ... 
  TPCParams& mp = Singleton<TPCParams>::Instance();
  
  double pitch_u = gds.pitch(WirePlaneType_t(0));
  double pitch_v = gds.pitch(WirePlaneType_t(1));
  double pitch_w = gds.pitch(WirePlaneType_t(2));
  double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;

  double angle_u = gds.angle(WirePlaneType_t(0));
  double angle_v = gds.angle(WirePlaneType_t(1));
  double angle_w = gds.angle(WirePlaneType_t(2));

  //std::cout << angle_u << " " << angle_v << " " << angle_w << std::endl;
  
  mp.set_pitch_u(pitch_u);
  mp.set_pitch_v(pitch_v);
  mp.set_pitch_w(pitch_w);
  mp.set_angle_u(angle_u);
  mp.set_angle_v(angle_v);
  mp.set_angle_w(angle_w);
  mp.set_ts_width(time_slice_width);
  mp.set_first_u_dis(first_u_dis);
  mp.set_first_v_dis(first_v_dis);
  mp.set_first_w_dis(first_w_dis);
  


  std::map<int,std::pair<double,double>> dead_u_index;
  std::map<int,std::pair<double,double>> dead_v_index;
  std::map<int,std::pair<double,double>> dead_w_index;
  
  // load mcell
  
  TTree *TC = (TTree*)file->Get("TC");
  std::vector<int> *cluster_id_vec = new std::vector<int>;
  std::vector<int> *time_slice_vec = new std::vector<int>;
  std::vector<double> *q_vec = new std::vector<double>;
  std::vector<double> *uq_vec = new std::vector<double>;
  std::vector<double> *vq_vec = new std::vector<double>;
  std::vector<double> *wq_vec = new std::vector<double>;
  std::vector<double> *udq_vec = new std::vector<double>;
  std::vector<double> *vdq_vec = new std::vector<double>;
  std::vector<double> *wdq_vec = new std::vector<double>;

  std::vector<int> *nwire_u_vec = new  std::vector<int>;
  std::vector<int> *nwire_v_vec = new  std::vector<int>;
  std::vector<int> *nwire_w_vec = new  std::vector<int>;
  std::vector<int> *flag_u_vec = new  std::vector<int>;
  std::vector<int> *flag_v_vec = new  std::vector<int>;
  std::vector<int> *flag_w_vec = new  std::vector<int>;

  std::vector<std::vector<int>> *wire_index_u_vec = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *wire_index_v_vec = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *wire_index_w_vec = new std::vector<std::vector<int>>;
  std::vector<std::vector<double>> *wire_charge_u_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_v_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_w_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_err_u_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_err_v_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_err_w_vec = new std::vector<std::vector<double>>;


  TC->SetBranchAddress("cluster_id",&cluster_id_vec);
  TC->SetBranchAddress("time_slice",&time_slice_vec);
  TC->SetBranchAddress("q",&q_vec);
  TC->SetBranchAddress("uq",&uq_vec);
  TC->SetBranchAddress("vq",&vq_vec);
  TC->SetBranchAddress("wq",&wq_vec);
  TC->SetBranchAddress("udq",&udq_vec);
  TC->SetBranchAddress("vdq",&vdq_vec);
  TC->SetBranchAddress("wdq",&wdq_vec);
  TC->SetBranchAddress("nwire_u",&nwire_u_vec);
  TC->SetBranchAddress("nwire_v",&nwire_v_vec);
  TC->SetBranchAddress("nwire_w",&nwire_w_vec);
  TC->SetBranchAddress("flag_u",&flag_u_vec);
  TC->SetBranchAddress("flag_v",&flag_v_vec);
  TC->SetBranchAddress("flag_w",&flag_w_vec);
  TC->SetBranchAddress("wire_index_u",&wire_index_u_vec);
  TC->SetBranchAddress("wire_index_v",&wire_index_v_vec);
  TC->SetBranchAddress("wire_index_w",&wire_index_w_vec);
  TC->SetBranchAddress("wire_charge_u",&wire_charge_u_vec);
  TC->SetBranchAddress("wire_charge_v",&wire_charge_v_vec);
  TC->SetBranchAddress("wire_charge_w",&wire_charge_w_vec);
  TC->SetBranchAddress("wire_charge_err_u",&wire_charge_err_u_vec);
  TC->SetBranchAddress("wire_charge_err_v",&wire_charge_err_v_vec);
  TC->SetBranchAddress("wire_charge_err_w",&wire_charge_err_w_vec);


  //load mcell
  
  TTree *TDC = (TTree*)file->Get("TDC");
  std::vector<int> *ntime_slice_vec = new std::vector<int>;
  std::vector<std::vector<int>> *time_slices_vec = new std::vector<std::vector<int>>;

  TDC->SetBranchAddress("cluster_id",&cluster_id_vec);
  TDC->SetBranchAddress("ntime_slice",&ntime_slice_vec);
  TDC->SetBranchAddress("time_slice",&time_slices_vec);
  
  TDC->SetBranchAddress("nwire_u",&nwire_u_vec);
  TDC->SetBranchAddress("nwire_v",&nwire_v_vec);
  TDC->SetBranchAddress("nwire_w",&nwire_w_vec);
  TDC->SetBranchAddress("flag_u",&flag_u_vec);
  TDC->SetBranchAddress("flag_v",&flag_v_vec);
  TDC->SetBranchAddress("flag_w",&flag_w_vec);
  TDC->SetBranchAddress("wire_index_u",&wire_index_u_vec);
  TDC->SetBranchAddress("wire_index_v",&wire_index_v_vec);
  TDC->SetBranchAddress("wire_index_w",&wire_index_w_vec);

  // load cells ... 
  GeomCellSelection mcells;
  PR3DClusterSelection live_clusters;
  PR3DClusterSelection dead_clusters;
  PR3DCluster *cluster;
  int prev_cluster_id=-1;
  int ident = 0;
  TC->GetEntry(entry_num);
  for (int i=0;i!=cluster_id_vec->size();i++){
    int cluster_id = cluster_id_vec->at(i);
    SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
    int time_slice = time_slice_vec->at(i);
    std::vector<int> wire_index_u = wire_index_u_vec->at(i);
    std::vector<int> wire_index_v = wire_index_v_vec->at(i);
    std::vector<int> wire_index_w = wire_index_w_vec->at(i);
    std::vector<double> wire_charge_u = wire_charge_u_vec->at(i);
    std::vector<double> wire_charge_v = wire_charge_v_vec->at(i);
    std::vector<double> wire_charge_w = wire_charge_w_vec->at(i);
    std::vector<double> wire_charge_err_u = wire_charge_err_u_vec->at(i);
    std::vector<double> wire_charge_err_v = wire_charge_err_v_vec->at(i);
    std::vector<double> wire_charge_err_w = wire_charge_err_w_vec->at(i);

    mcell->SetTimeSlice(time_slice);

    mcell->set_uq(uq_vec->at(i));
    mcell->set_vq(vq_vec->at(i));
    mcell->set_wq(wq_vec->at(i));

    mcell->set_udq(udq_vec->at(i));
    mcell->set_vdq(vdq_vec->at(i));
    mcell->set_wdq(wdq_vec->at(i));

    mcell->set_q(q_vec->at(i));

    double temp_x = (time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm;

    
    if (flag_u_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(0));
      for (int j=0;j!=nwire_u_vec->at(i);j++){
	if (dead_u_index.find(wire_index_u.at(j))==dead_u_index.end()){
	  dead_u_index[wire_index_u.at(j)] = std::make_pair(temp_x-0.1*units::cm,temp_x+0.1*units::cm);
	}else{
	  if (temp_x-0.1*units::cm < dead_u_index[wire_index_u.at(j)].first){
	    dead_u_index[wire_index_u.at(j)].first = temp_x - 0.1*units::cm;
	  }else if (temp_x+0.1*units::cm > dead_u_index[wire_index_u.at(j)].second){
	    dead_u_index[wire_index_u.at(j)].second = temp_x + 0.1*units::cm;
	  }
	}
	
      }
    }
    if (flag_v_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(1));
      for (int j=0;j!=nwire_v_vec->at(i);j++){
	if (dead_v_index.find(wire_index_v.at(j))==dead_v_index.end()){
	  dead_v_index[wire_index_v.at(j)] = std::make_pair(temp_x-0.1*units::cm,temp_x+0.1*units::cm);
	}else{
	  if (temp_x-0.1*units::cm < dead_v_index[wire_index_v.at(j)].first){
	    dead_v_index[wire_index_v.at(j)].first = temp_x-0.1*units::cm;
	  }else if (temp_x+0.1*units::cm > dead_v_index[wire_index_v.at(j)].second){
	    dead_v_index[wire_index_v.at(j)].second = temp_x + 0.1*units::cm;
	  }
	}
      }
    }
    if (flag_w_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(2));
      for (int j=0;j!=nwire_w_vec->at(i);j++){
	if (dead_w_index.find(wire_index_w.at(j))==dead_w_index.end()){
	  dead_w_index[wire_index_w.at(j)] = std::make_pair(temp_x-0.1*units::cm,temp_x+0.1*units::cm);
	}else{
	  if (temp_x-0.1*units::cm < dead_w_index[wire_index_w.at(j)].first){
	    dead_w_index[wire_index_w.at(j)].first = temp_x-0.1*units::cm;
	  }else if (temp_x+0.1*units::cm > dead_w_index[wire_index_w.at(j)].second){
	    dead_w_index[wire_index_w.at(j)].second = temp_x + 0.1*units::cm;
	  }
	}
      }
    }
    for (int j=0;j!=nwire_u_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.at(j));
      mcell->AddWire(wire,WirePlaneType_t(0),wire_charge_u.at(j),wire_charge_err_u.at(j));
    }
    for (int j=0;j!=nwire_v_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.at(j));
      mcell->AddWire(wire,WirePlaneType_t(1),wire_charge_v.at(j),wire_charge_err_v.at(j));
    }
    for (int j=0;j!=nwire_w_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.at(j));
      mcell->AddWire(wire,WirePlaneType_t(2),wire_charge_w.at(j),wire_charge_err_w.at(j));
    }
    mcells.push_back(mcell);

    if (cluster_id != prev_cluster_id){
      cluster = new PR3DCluster(cluster_id);
      live_clusters.push_back(cluster);
    }
    cluster->AddCell(mcell,time_slice);

    prev_cluster_id = cluster_id;
    ident++;
  }
  //  std::cout << live_clusters.size() << std::endl;

  prev_cluster_id = -1;
  // TDC
   cluster_id_vec->clear();
   wire_index_u_vec->clear();
   wire_index_v_vec->clear();
   wire_index_w_vec->clear();
   nwire_u_vec->clear();
   nwire_v_vec->clear();
   nwire_w_vec->clear();
   flag_u_vec->clear();
   flag_v_vec->clear();
   flag_w_vec->clear();

   TDC->GetEntry(entry_num);
   for (int i=0;i!=cluster_id_vec->size();i++){
    int cluster_id = cluster_id_vec->at(i);
    SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
    std::vector<int> time_slices = time_slices_vec->at(i);
    std::vector<int> wire_index_u = wire_index_u_vec->at(i);
    std::vector<int> wire_index_v = wire_index_v_vec->at(i);
    std::vector<int> wire_index_w = wire_index_w_vec->at(i);

    mcell->SetTimeSlice(time_slices.at(0));

    if (flag_u_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(0));
      // for (int i=0;i!=nwire_u;i++){
      // 	dead_u_index.insert(wire_index_u[i]);
      // }
    }
    if (flag_v_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(1));
      // for (int i=0;i!=nwire_v;i++){
      // 	dead_v_index.insert(wire_index_v[i]);
      // }
    }
    if (flag_w_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(2));
      // for (int i=0;i!=nwire_w;i++){
      // 	dead_w_index.insert(wire_index_w[i]);
      // }
    }
    for (int j=0;j!=nwire_u_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.at(j));
      mcell->AddWire(wire,WirePlaneType_t(0));
    }
    for (int j=0;j!=nwire_v_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.at(j));
      mcell->AddWire(wire,WirePlaneType_t(1));
    }
    for (int j=0;j!=nwire_w_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.at(j));
      mcell->AddWire(wire,WirePlaneType_t(2));
    }
    mcells.push_back(mcell);

    if (cluster_id!=prev_cluster_id){
      cluster = new PR3DCluster(cluster_id);
      dead_clusters.push_back(cluster);
    }
    for (int j=0;j!=ntime_slice_vec->at(i);j++){
      cluster->AddCell(mcell,time_slices.at(j));
    }
    
    prev_cluster_id=cluster_id;
    ident++;
  }

   // std::cout << live_clusters.size() << std::endl;
   // for (size_t i=0;i!=live_clusters.size();i++){
   //   std::cout << live_clusters.at(i)->get_cluster_id() << " " 
   // 	       << live_clusters.at(i)->get_num_mcells() << " "
   // 	       << live_clusters.at(i)->get_num_time_slices() << std::endl;
   // }
   // std::cout << dead_clusters.size() << std::endl;
   for (size_t i=0;i!=dead_clusters.size();i++){
     dead_clusters.at(i)->Remove_duplicated_mcells();
     // std::cout << dead_clusters.at(i)->get_cluster_id() << " " 
     // 	       << dead_clusters.at(i)->get_num_mcells() << " "
     // 	       << dead_clusters.at(i)->get_num_time_slices() << std::endl;
   }

   // for (auto it = dead_u_index.begin(); it!=dead_u_index.end(); it++){
   //   std::cout << it->first << " " << it->second.first/units::cm << " " << it->second.second/units::cm << std::endl;
   // }
  
   
   //   cerr << em("load clusters from file") << endl;

  

   // Start to add X, Y, Z points
   // form boundaries for the bad cells ... 
   for (size_t j = 0; j!= dead_clusters.size(); j++){
     WireCell2dToy::calc_boundary_points_dead(gds,dead_clusters.at(j));
   }
   // form sampling points for the normal cells ...
   for (size_t i=0; i!=live_clusters.size();i++){
     WireCell2dToy::calc_sampling_points(gds,live_clusters.at(i),nrebin, frame_length, unit_dis);
     // live_clusters.at(i)->Calc_PCA();
   }
   //cerr << em("Add X, Y, Z points") << std::endl;

   // create global point cloud and mcell to cluster map ...
   // ToyPointCloud *global_point_cloud =  new ToyPointCloud();
   //std::map<SlimMergeGeomCell*,PR3DCluster*> mcell_cluster_map;
   DynamicToyPointCloud global_point_cloud(angle_u,angle_v,angle_w);
   for (size_t i=0;i!=live_clusters.size();i++){
     live_clusters.at(i)->Create_point_cloud();
     global_point_cloud.AddPoints(live_clusters.at(i),0);
     // live_clusters.at(i)->Update_mcell_cluster_map(mcell_cluster_map);
   }
   //   global_point_cloud->build_kdtree_index();
   cerr << em("Build local point clouds") << std::endl;
   

   WireCell2dToy::Clustering_jump_gap_cosmics(live_clusters, dead_clusters, dead_u_index, dead_v_index, dead_w_index, global_point_cloud);
   cerr << em("Clustering to jump gap in cosmics") << std::endl;


   
   // //for (size_t i=0;i!=live_clusters.size();i++){
   //   // live_clusters.at(i)->Create_point_cloud();
   //   // std::cout << i << " "<< live_clusters.at(i)->get_point_cloud()->get_num_points() << std::endl;
   // //}
   
   for (size_t i=0;i!=live_clusters.size();i++){
     //    std::cout << live_clusters.at(i)->get_mcells().size() << " " << live_clusters.at(i)->get_num_time_slices() << std::endl;
     live_clusters.at(i)->Create_graph();
     std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = live_clusters.at(i)->get_highest_lowest_wcps();
     live_clusters.at(i)->dijkstra_shortest_paths(wcps.first);
     live_clusters.at(i)->cal_shortest_path(wcps.second);
     live_clusters.at(i)->fine_tracking();
   }
   
   // cerr << em("Trajectory fit in all clusters") << std::endl;
   
   // Point p1(337.346*units::cm,87.0524*units::cm,697.899*units::cm);
   // for (int i=0;i!=live_clusters.size();i++){
   //   if (live_clusters.at(i)->get_cluster_id()==10){
   //     live_clusters.at(i)->Create_point_cloud();
   //     Point p = live_clusters.at(i)->calc_ave_pos(p1,30*units::cm);
   //     std::cout << "R: " << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << " " << sqrt(pow(p1.x-p.x,2)+pow(p1.y-p.y,2)+pow(p1.z-p.z,2))/units::cm << std::endl;
   //   }
   // }
   
   
   
   
   TFile *file1 = new TFile(Form("pr_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");

   TTree *t_bad = new TTree("T_bad","T_bad");
   t_bad->SetDirectory(file1);
   Int_t bad_npoints;
   Int_t ncluster;
   Double_t bad_y[100],bad_z[100];
   t_bad->Branch("cluster_id",&ncluster,"cluster_id/I");
   t_bad->Branch("bad_npoints",&bad_npoints,"bad_npoints/I");
   t_bad->Branch("bad_y",bad_y,"bad_y[bad_npoints]/D");
   t_bad->Branch("bad_z",bad_z,"bad_z[bad_npoints]/D");
   
   for (size_t j = 0; j!= dead_clusters.size(); j++){
     SMGCSelection& mcells = dead_clusters.at(j)->get_mcells();
     ncluster = dead_clusters.at(j)->get_cluster_id();
     for (size_t i=0;i!=mcells.size();i++){
       PointVector ps = mcells.at(i)->boundary();
       bad_npoints = ps.size();
       for (int k=0;k!=bad_npoints;k++){
	 bad_y[k] = ps.at(k).y/units::cm;
	 bad_z[k] = ps.at(k).z/units::cm;
       }
       t_bad->Fill();
     }
   }
   
   TTree *T_cluster ;
   Double_t x,y,z,q,nq;
   
   T_cluster = new TTree("T_cluster","T_cluster");
   T_cluster->Branch("cluster_id",&ncluster,"cluster_id/I");
   T_cluster->Branch("x",&x,"x/D");
   T_cluster->Branch("y",&y,"y/D");
   T_cluster->Branch("z",&z,"z/D");
   T_cluster->Branch("q",&q,"q/D");
   T_cluster->Branch("nq",&nq,"nq/D");
	 
   T_cluster->SetDirectory(file1);

   TTree *T_rec = new TTree("T_rec","T_rec");
   T_rec->Branch("x",&x,"x/D");
   T_rec->Branch("y",&y,"y/D");
   T_rec->Branch("z",&z,"z/D");
   T_rec->SetDirectory(file1);

   Double_t charge_save=1, ncharge_save=1, chi2_save=1, ndf_save=1;
   TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
   t_rec_charge->SetDirectory(file1);
   t_rec_charge->Branch("x",&x,"x/D");
   t_rec_charge->Branch("y",&y,"y/D");
   t_rec_charge->Branch("z",&z,"z/D");
   t_rec_charge->Branch("q",&charge_save,"q/D");
   t_rec_charge->Branch("nq",&ncharge_save,"nq/D");
   t_rec_charge->Branch("chi2",&chi2_save,"chi2/D");
   t_rec_charge->Branch("ndf",&ndf_save,"ndf/D");
   
   // test ... 
   // ncluster = 0;
   // x=  0;
   // y=0;
   // z=0;
   // T_cluster->Fill();
   
   for (size_t j = 0; j!= live_clusters.size(); j++){
     SMGCSelection& mcells = live_clusters.at(j)->get_mcells();
     ncluster = live_clusters.at(j)->get_cluster_id();
     for (size_t i=0;i!=mcells.size();i++){
       PointVector ps = mcells.at(i)->get_sampling_points();
       int time_slice = mcells.at(i)->GetTimeSlice();
       if (ps.size()==0){
	 std::cout << "zero sampling points!" << std::endl;
       }else{
	 q = mcells.at(i)->get_q() / ps.size();
	 nq = ps.size();
	 for (int k=0;k!=ps.size();k++){
	   x = (ps.at(k).x)/units::cm ;
	   y = ps.at(k).y/units::cm;
	   z = ps.at(k).z/units::cm;
	   T_cluster->Fill();
	 }
       }

     }

     // save wcps
     std::list<WCPointCloud<double>::WCPoint>& wcps_list = live_clusters.at(j)->get_path_wcps();
     //ncluster = -1 * ncluster-100;
     for (auto it = wcps_list.begin(); it!=wcps_list.end(); it++){
       x = (*it).x/units::cm;
       y = (*it).y/units::cm;
       z = (*it).z/units::cm;
       T_rec->Fill();
     }

     PointVector& pts = live_clusters.at(j)->get_fine_tracking_path();
     for (size_t i=0; i!=pts.size(); i++){
       x = pts.at(i).x/units::cm;
       y = pts.at(i).y/units::cm;
       z = pts.at(i).z/units::cm;
       t_rec_charge->Fill();
     }

     

     // // save mcells
     // std::list<SlimMergeGeomCell*>& mcells_list = live_clusters.at(j)->get_path_mcells();
     // ncluster = -1 * ncluster-100;
     // for (auto it = mcells_list.begin(); it!=mcells_list.end(); it++){
     //   Point p = (*it)->center();
     //   x = p.x/units::cm;
     //   y = p.y/units::cm;
     //   z = p.z/units::cm;
     //   T_cluster->Fill();
     // }
     
     
     // if (live_clusters.at(j)->get_num_mcells()>30){
     //   // add PCA axis point
     //   Vector center = live_clusters.at(j)->get_center();
     //   Vector dir = live_clusters.at(j)->get_PCA_axis(0);
     //   for (int i=-200;i!=200;i++){
     // 	 x = (center.x + dir.x *(i*units::cm) )/units::cm;
     // 	 y = (center.y + dir.y *(i*units::cm) )/units::cm;
     // 	 z = (center.z + dir.z *(i*units::cm) )/units::cm;
     // 	 T_cluster->Fill();
     //   }
     // }
    }

   // ncluster = 0;
   // for (auto it = dead_live_cluster_mapping.begin(); it!= dead_live_cluster_mapping.end(); it++){
   //   std::vector<PR3DCluster*> clusters = (*it).second;
   //   if (clusters.size()>1){
   //     //std::cout << clusters.size() << std::endl;
   //     for (auto it1 = clusters.begin(); it1!=clusters.end(); it1++){
   // 	 PR3DCluster* cluster = (*it1);
   // 	 ncluster = cluster->get_cluster_id();
   // 	 SMGCSelection& mcells = cluster->get_mcells();
   // 	 for (size_t i=0;i!=mcells.size();i++){
   // 	   PointVector ps = mcells.at(i)->get_sampling_points();
   // 	   for (int k=0;k!=ps.size();k++){
   // 	     x = ps.at(k).x/units::cm;//time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
   // 	     y = ps.at(k).y/units::cm;
   // 	     z = ps.at(k).z/units::cm;
   // 	     T_cluster->Fill();
   // 	   }
   // 	 }
   //     }
   //   }
   //   // ncluster++;
   // }
   // for (auto it = dead_live_mcells_mapping.begin(); it!= dead_live_mcells_mapping.end(); it++){
   //   std::vector<std::vector<SlimMergeGeomCell*>> mcellss = (*it).second;
   //   // std::cout << mcellss.size() << std::endl;
   //   if (mcellss.size()>1){
   //     for (auto it1 = mcellss.begin(); it1!=mcellss.end(); it1++){
   // 	 std::vector<SlimMergeGeomCell*> mcells = (*it1);
   // 	 for (size_t i=0;i!=mcells.size();i++){
   // 	   PointVector ps = mcells.at(i)->get_sampling_points();
   // 	   for (int k=0;k!=ps.size();k++){
   // 	     x = ps.at(k).x/units::cm;//time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
   // 	     y = ps.at(k).y/units::cm;
   // 	     z = ps.at(k).z/units::cm;
   // 	     T_cluster->Fill();
   // 	   }
   // 	 }
   //     }
   //   }
   //   ncluster++;
   // }
   
   Trun->CloneTree()->Write();
   file1->Write();
   file1->Close();
}
