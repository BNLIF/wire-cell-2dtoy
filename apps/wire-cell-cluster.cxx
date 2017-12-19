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
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root " << endl;
    return 1;
  }

  ExecMon em("starting");
  cerr << em("load geometry") << endl;

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
  float unit_dis;
  
  Trun->SetBranchAddress("eventNo",&event_no);
  Trun->SetBranchAddress("runNo",&run_no);
  Trun->SetBranchAddress("subRunNo",&subrun_no);
  Trun->SetBranchAddress("unit_dis",&unit_dis);
  Trun->SetBranchAddress("frame_length",&frame_length);
  Trun->SetBranchAddress("eve_num",&eve_num);
  Trun->SetBranchAddress("nrebin",&nrebin);
  Trun->SetBranchAddress("time_offset",&time_offset);
  Trun->GetEntry(0);

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
  
  
  // load mcell
  
  TTree *TC = (TTree*)file->Get("TC");
  Int_t cluster_id;
  Int_t time_slice;
  Double_t q, uq, vq, wq, udq, vdq, wdq;
  Int_t nwire_u=0, flag_u; //number of wires, dead?
  Int_t nwire_v=0, flag_v;
  Int_t nwire_w=0, flag_w;
  Int_t wire_index_u[2400];
  Int_t wire_index_v[2400];
  Int_t wire_index_w[3256];
  Double_t wire_charge_u[2400];
  Double_t wire_charge_v[2400];
  Double_t wire_charge_w[2400];
  Double_t wire_charge_err_u[2400];
  Double_t wire_charge_err_v[2400];
  Double_t wire_charge_err_w[2400];

  TC->SetBranchAddress("cluster_id",&cluster_id);
  TC->SetBranchAddress("time_slice",&time_slice);
  TC->SetBranchAddress("q",&q);
  TC->SetBranchAddress("uq",&uq);
  TC->SetBranchAddress("vq",&vq);
  TC->SetBranchAddress("wq",&wq);
  TC->SetBranchAddress("udq",&udq);
  TC->SetBranchAddress("vdq",&vdq);
  TC->SetBranchAddress("wdq",&wdq);
  TC->SetBranchAddress("nwire_u",&nwire_u);
  TC->SetBranchAddress("nwire_v",&nwire_v);
  TC->SetBranchAddress("nwire_w",&nwire_w);
  TC->SetBranchAddress("flag_u",&flag_u);
  TC->SetBranchAddress("flag_v",&flag_v);
  TC->SetBranchAddress("flag_w",&flag_w);
  TC->SetBranchAddress("wire_index_u",wire_index_u);
  TC->SetBranchAddress("wire_index_v",wire_index_v);
  TC->SetBranchAddress("wire_index_w",wire_index_w);
  TC->SetBranchAddress("wire_charge_u",wire_charge_u);
  TC->SetBranchAddress("wire_charge_v",wire_charge_v);
  TC->SetBranchAddress("wire_charge_w",wire_charge_w);
  TC->SetBranchAddress("wire_charge_err_u",wire_charge_err_u);
  TC->SetBranchAddress("wire_charge_err_v",wire_charge_err_v);
  TC->SetBranchAddress("wire_charge_err_w",wire_charge_err_w);


  //load mcell
  
  TTree *TDC = (TTree*)file->Get("TDC");
  int ntime_slice = 0,time_slices[2400];
  TDC->SetBranchAddress("cluster_id",&cluster_id);
  TDC->SetBranchAddress("ntime_slice",&ntime_slice);
  TDC->SetBranchAddress("time_slice",time_slices);
  
  TDC->SetBranchAddress("nwire_u",&nwire_u);
  TDC->SetBranchAddress("nwire_v",&nwire_v);
  TDC->SetBranchAddress("nwire_w",&nwire_w);
  TDC->SetBranchAddress("flag_u",&flag_u);
  TDC->SetBranchAddress("flag_v",&flag_v);
  TDC->SetBranchAddress("flag_w",&flag_w);
  TDC->SetBranchAddress("wire_index_u",wire_index_u);
  TDC->SetBranchAddress("wire_index_v",wire_index_v);
  TDC->SetBranchAddress("wire_index_w",wire_index_w);

  // load cells ... 
  GeomCellSelection mcells;
  PR3DClusterSelection live_clusters;
  PR3DClusterSelection dead_clusters;
  PR3DCluster *cluster;
  int prev_cluster_id=-1;
  int ident = 0;
  for (int i=0;i!=TC->GetEntries();i++){
    TC->GetEntry(i);
    SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
    mcell->SetTimeSlice(time_slice);

    mcell->set_uq(uq);
    mcell->set_vq(vq);
    mcell->set_wq(wq);

    mcell->set_udq(udq);
    mcell->set_vdq(vdq);
    mcell->set_wdq(wdq);

    mcell->set_q(q);
    if (flag_u==0){
      mcell->add_bad_planes(WirePlaneType_t(0));
    }
    if (flag_v==0){
      mcell->add_bad_planes(WirePlaneType_t(1));
    }
    if (flag_w==0){
      mcell->add_bad_planes(WirePlaneType_t(2));
    }
    for (int i=0;i!=nwire_u;i++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(0),wire_index_u[i]);
      mcell->AddWire(wire,WirePlaneType_t(0),wire_charge_u[i],wire_charge_err_u[i]);
    }
    for (int i=0;i!=nwire_v;i++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(1),wire_index_v[i]);
      mcell->AddWire(wire,WirePlaneType_t(1),wire_charge_v[i],wire_charge_err_v[i]);
    }
    for (int i=0;i!=nwire_w;i++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(2),wire_index_w[i]);
      mcell->AddWire(wire,WirePlaneType_t(2),wire_charge_w[i],wire_charge_err_w[i]);
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
   for (int i=0;i!=TDC->GetEntries();i++){
    TDC->GetEntry(i);

    SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
    mcell->SetTimeSlice(time_slices[0]);

    if (flag_u==0){
      mcell->add_bad_planes(WirePlaneType_t(0));
    }
    if (flag_v==0){
      mcell->add_bad_planes(WirePlaneType_t(1));
    }
    if (flag_w==0){
      mcell->add_bad_planes(WirePlaneType_t(2));
    }
    for (int i=0;i!=nwire_u;i++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(0),wire_index_u[i]);
      mcell->AddWire(wire,WirePlaneType_t(0));
    }
    for (int i=0;i!=nwire_v;i++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(1),wire_index_v[i]);
      mcell->AddWire(wire,WirePlaneType_t(1));
    }
    for (int i=0;i!=nwire_w;i++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(2),wire_index_w[i]);
      mcell->AddWire(wire,WirePlaneType_t(2));
    }
    mcells.push_back(mcell);

    if (cluster_id!=prev_cluster_id){
      cluster = new PR3DCluster(cluster_id);
      dead_clusters.push_back(cluster);
    }
    for (int i=0;i!=ntime_slice;i++){
      cluster->AddCell(mcell,time_slices[i]);
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
  
   
   cerr << em("load clusters from file") << endl;

  

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
   cerr << em("Add X, Y, Z points") << std::endl;

   // create global point cloud and mcell to cluster map ...
   ToyPointCloud *global_point_cloud =  new ToyPointCloud();
   std::map<SlimMergeGeomCell*,PR3DCluster*> mcell_cluster_map;
   for (size_t i=0;i!=live_clusters.size();i++){
     live_clusters.at(i)->Create_point_cloud(global_point_cloud);
     live_clusters.at(i)->Update_mcell_cluster_map(mcell_cluster_map);
   }
   global_point_cloud->build_kdtree_index();

   // std::cout << mcell_cluster_map.size() << " " << global_point_cloud->get_num_points() << std::endl;
   
   cerr << em("Build global and local point clouds") << std::endl;
   
   WireCell2dToy::Clustering_live_dead(live_clusters, dead_clusters);
   cerr << em("Clustering live and dead clusters") << std::endl;

   WireCell2dToy::Clustering_jump_gap_cosmics(live_clusters);
   cerr << em("Clustering to jump gap in cosmics") << std::endl;


   
   // //for (size_t i=0;i!=live_clusters.size();i++){
   //   // live_clusters.at(i)->Create_point_cloud();
   //   // std::cout << i << " "<< live_clusters.at(i)->get_point_cloud()->get_num_points() << std::endl;
   // //}
   
   // for (size_t i=0;i!=live_clusters.size();i++){
   //   //    std::cout << live_clusters.at(i)->get_mcells().size() << " " << live_clusters.at(i)->get_num_time_slices() << std::endl;
   //   live_clusters.at(i)->Create_graph();
   //   std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = live_clusters.at(i)->get_highest_lowest_wcps();
   //   live_clusters.at(i)->dijkstra_shortest_paths(wcps.first);
   //   live_clusters.at(i)->cal_shortest_path(wcps.second);
   //   live_clusters.at(i)->fine_tracking();
   // }
   
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
   Double_t x,y,z;
   
   T_cluster = new TTree("T_cluster","T_cluster");
   T_cluster->Branch("cluster_id",&ncluster,"cluster_id/I");
   T_cluster->Branch("x",&x,"x/D");
   T_cluster->Branch("y",&y,"y/D");
   T_cluster->Branch("z",&z,"z/D");
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
       if (ps.size()==0) std::cout << "zero sampling points!" << std::endl;
       for (int k=0;k!=ps.size();k++){
    	 x = ps.at(k).x/units::cm;//time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
    	 y = ps.at(k).y/units::cm;
    	 z = ps.at(k).z/units::cm;
    	 T_cluster->Fill();
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
