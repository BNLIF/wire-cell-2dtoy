#include "WireCell2dToy/ToyClustering.h"

#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"
#include "WireCellData/Line.h"

using namespace WireCell;
using namespace std;

#include "ToyClustering_dead_live.h"
#include "ToyClustering_reg.h"
#include "ToyClustering_para_prol.h"
#include "ToyClustering_close.h"
#include "ToyClustering_extend.h"
#include "ToyClustering_isolated.h"

#include "WireCell2dToy/ExecMon.h"

double WireCell2dToy::cal_proj_angle_diff(TVector3& dir1, TVector3& dir2, double plane_angle){
  TVector3 temp_dir1;
  TVector3 temp_dir2;

  temp_dir1.SetXYZ(dir1.X(),0,-sin(plane_angle)*dir1.Y()+cos(plane_angle)*dir1.Z());
  temp_dir2.SetXYZ(dir2.X(),0,-sin(plane_angle)*dir2.Y()+cos(plane_angle)*dir2.Z());
  
  return temp_dir1.Angle(temp_dir2);
}

bool WireCell2dToy::is_angle_consistent(TVector3& dir1, TVector3& dir2, bool same_direction, double angle_cut, double uplane_angle, double vplane_angle, double wplane_angle, int num_cut){
  double angle_u = WireCell2dToy::cal_proj_angle_diff(dir1,dir2,uplane_angle);
  double angle_v = WireCell2dToy::cal_proj_angle_diff(dir1,dir2,vplane_angle);
  double angle_w = WireCell2dToy::cal_proj_angle_diff(dir1,dir2,wplane_angle);
  int num = 0;
  //input is degrees ...
  angle_cut *= 3.1415926/180.;
  
  if (same_direction){
    if (angle_u <= angle_cut) num++;
    if (angle_v <= angle_cut) num++;
    if (angle_w <= angle_cut) num++;
  }else{
    if ((3.1415926-angle_u) <= angle_cut) num++;
    if ((3.1415926-angle_v) <= angle_cut) num++;
    if ((3.1415926-angle_w) <= angle_cut) num++;
  }
  
  if (num>=num_cut ) return true;
  return false;
}


double WireCell2dToy::Find_Closeset_Points(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2,double length_1, double length_2, double length_cut, SlimMergeGeomCell *mcell1_save, SlimMergeGeomCell *mcell2_save, Point& p1_save, Point &p2_save){
  double dis_save = 1e9;

  // pick the first point from the small one
  SlimMergeGeomCell *prev_mcell1 = 0;
  SlimMergeGeomCell *prev_mcell2 = 0;
  SlimMergeGeomCell *mcell1 = 0; 
  Point p1;//
  SlimMergeGeomCell *mcell2=0;
  Point p2;
  double dis;

  if (length_1 < length_2){
    mcell1 = *(cluster1->get_time_cells_set_map().begin()->second.begin());
    p1 = mcell1->center();

    while(mcell1!=prev_mcell1 || mcell2!=prev_mcell2){
      prev_mcell1 = mcell1;
      prev_mcell2 = mcell2;
      
      // find the closest point and merged cell in cluster2
      std::pair<SlimMergeGeomCell*,Point> temp_results = cluster2->get_closest_point_mcell(p1);
      p2 = temp_results.second;
      mcell2 = temp_results.first;
      // find the closest point and merged cell in cluster1
      temp_results = cluster1->get_closest_point_mcell(p2);
      p1 = temp_results.second;
      mcell1 = temp_results.first;
    }
    dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));

    if (dis < dis_save){
      dis_save = dis;
      mcell1_save = mcell1;
      mcell2_save = mcell2;
      p1_save = p1;
      p2_save = p2;
    }

    prev_mcell1 = 0;
    prev_mcell2 = 0;

    mcell1 = *(cluster1->get_time_cells_set_map().rbegin()->second.begin());
    p1 = mcell1->center();

    while(mcell1!=prev_mcell1 || mcell2!=prev_mcell2){
      prev_mcell1 = mcell1;
      prev_mcell2 = mcell2;
      
      // find the closest point and merged cell in cluster2
      std::pair<SlimMergeGeomCell*,Point> temp_results = cluster2->get_closest_point_mcell(p1);
      p2 = temp_results.second;
      mcell2 = temp_results.first;
      // find the closest point and merged cell in cluster1
      temp_results = cluster1->get_closest_point_mcell(p2);
      p1 = temp_results.second;
      mcell1 = temp_results.first;
    }
    dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));

    if (dis < dis_save){
      dis_save = dis;
      mcell1_save = mcell1;
      mcell2_save = mcell2;
      p1_save = p1;
      p2_save = p2;
    }
    
  }else{

    mcell2 = *(cluster2->get_time_cells_set_map().begin()->second.begin());
    p2 = mcell2->center();

    while(mcell1!=prev_mcell1 || mcell2!=prev_mcell2){
      prev_mcell1 = mcell1;
      prev_mcell2 = mcell2;
      
      // find the closest point and merged cell in cluster2
      std::pair<SlimMergeGeomCell*,Point> temp_results = cluster1->get_closest_point_mcell(p2);
      p1 = temp_results.second;
      mcell1 = temp_results.first;
      // find the closest point and merged cell in cluster1
      temp_results = cluster2->get_closest_point_mcell(p1);
      p2 = temp_results.second;
      mcell2 = temp_results.first;
    }
    dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));

    if (dis < dis_save){
      dis_save = dis;
      mcell1_save = mcell1;
      mcell2_save = mcell2;
      p1_save = p1;
      p2_save = p2;
    }

    
    prev_mcell1 = 0;
    prev_mcell2 = 0;
    

    mcell2 = *(cluster2->get_time_cells_set_map().rbegin()->second.begin());
    p2 = mcell2->center();

    while(mcell1!=prev_mcell1 || mcell2!=prev_mcell2){
      prev_mcell1 = mcell1;
      prev_mcell2 = mcell2;
      
      // find the closest point and merged cell in cluster2
      std::pair<SlimMergeGeomCell*,Point> temp_results = cluster1->get_closest_point_mcell(p2);
      p1 = temp_results.second;
      mcell1 = temp_results.first;
      // find the closest point and merged cell in cluster1
      temp_results = cluster2->get_closest_point_mcell(p1);
      p2 = temp_results.second;
      mcell2 = temp_results.first;
    }
    dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));

    if (dis < dis_save){
      dis_save = dis;
      mcell1_save = mcell1;
      mcell2_save = mcell2;
      p1_save = p1;
      p2_save = p2;
    }
    
  }
  
  
  return dis_save;
}


void WireCell2dToy::Clustering_jump_gap_cosmics(WireCell::PR3DClusterSelection& live_clusters, WireCell::PR3DClusterSelection& dead_clusters){

  
  //  ExecMon em("starting");


  // include some parallel or prolonged, no need to do track fitting
  std::map<PR3DCluster*,double> cluster_length_map;

  std::set<PR3DCluster*> cluster_connected_dead;
  
  // calculate the length ...
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  //std::vector<double> cluster_length_vec;
  //estimated length

  for (size_t i=0; i!= live_clusters.size(); i++){
    PR3DCluster* cluster_1 = live_clusters.at(i);
    std::vector<int> range_v1 = cluster_1->get_uvwt_range();
    double length_1 = sqrt(2./3. * (pow(pitch_u*range_v1.at(0),2) + pow(pitch_v*range_v1.at(1),2) + pow(pitch_w*range_v1.at(2),2)) + pow(time_slice_width*range_v1.at(3),2));
    //cluster_length_vec.push_back(length_1);
    cluster_length_map[cluster_1] = length_1;
  }
  
  
  //cluster live dead ...
  Clustering_live_dead(live_clusters, dead_clusters, cluster_length_map, cluster_connected_dead);
  // try to do the large ones immediate ... 
  Clustering_extend(live_clusters, cluster_length_map,cluster_connected_dead,4,60*units::cm,0,15*units::cm);
  // cerr << em("live dead") << endl;

  // std::cout << cluster_connected_dead.size() << std::endl;
  
  // first round clustering
  Clustering_regular(live_clusters, cluster_length_map,cluster_connected_dead,60*units::cm,false);

  //  cerr << em("1st regular") << endl;
  Clustering_regular(live_clusters, cluster_length_map,cluster_connected_dead,30*units::cm,true); // do extension

  // cerr << em("2nd regular") << endl;

  
  //dedicated one dealing with parallel and prolonged track
  Clustering_parallel_prolong(live_clusters, cluster_length_map,cluster_connected_dead,35*units::cm);

  // cerr << em("parallel prolong") << endl;
  
  //clustering close distance ones ... 
  Clustering_close(live_clusters, cluster_length_map,cluster_connected_dead,1.2*units::cm);

  // cerr << em("close") << endl;

  // std::cout << cluster_connected_dead.size() << std::endl;

  
  //std::cout << "Num. of clusters: " << live_clusters.size() << std::endl;

  int num_try =3;
  // for very busy events do less ... 
  if (live_clusters.size() > 1100 ) num_try = 1;

  
  
  for (int i=0;i!= num_try ;i++){
    //extend the track ...
    // deal with prolong case
    Clustering_extend(live_clusters, cluster_length_map,cluster_connected_dead,1,150*units::cm,0);
    //  cerr << em("extend prolong") << endl;
    // deal with parallel case 
    Clustering_extend(live_clusters, cluster_length_map,cluster_connected_dead,2,30*units::cm,0);
    // cerr << em("extend parallel") << endl;
    // extension regular case
    Clustering_extend(live_clusters, cluster_length_map,cluster_connected_dead,3,15*units::cm,0);
    //  cerr << em("extend regular") << endl;
    // extension ones connected to dead region ...
    Clustering_extend(live_clusters, cluster_length_map,cluster_connected_dead,4,60*units::cm,i);
    
    //  cerr << em("extend dead") << endl;
  }

  
  
  //std::cout << cluster_connected_dead.size() << std::endl;
  
  
}







