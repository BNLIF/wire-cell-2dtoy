#ifndef WCP2dToy_TOYCLUSTERING_H
#define WCP2dToy_TOYCLUSTERING_H

#include "WCPData/PR3DCluster.h"
#include "WCPData/DynamicToyPointCloud.h"
#include "WCPData/ToyCTPointCloud.h"
//#include "WCPSst/GeomDataSource.h"

namespace WCP2dToy{
  double Find_Closeset_Points(WCP::PR3DCluster *cluster1, WCP::PR3DCluster *cluster2,double length_1, double length_2, double length_cut, WCP::SlimMergeGeomCell *mcell1, WCP::SlimMergeGeomCell *mcell2, WCP::Point& p1, WCP::Point &p2);

  
  double cal_proj_angle_diff(TVector3& dir1, TVector3& dir2, double plane_angle);
  bool is_angle_consistent(TVector3& dir1, TVector3& dir2, bool same_direction, double angle_cut, double angle_u, double angle_v, double angle_w, int num_cut = 2);
  
  
  std::map<WCP::PR3DCluster*,std::vector<std::pair<WCP::PR3DCluster*,double>>> Clustering_jump_gap_cosmics(WCP::PR3DClusterSelection& live_clusters, WCP::PR3DClusterSelection& dead_clusters, std::map<int,std::pair<double,double>>& dead_u_index, std::map<int,std::pair<double,double>>& dead_v_index, std::map<int,std::pair<double,double>>& dead_w_index, WCP::DynamicToyPointCloud& global_point_cloud, WCP::ToyCTPointCloud& ct_point_cloud, bool flag_neutrino = true);

  std::map<WCP::PR3DCluster*,std::vector<std::pair<WCP::PR3DCluster*,double>>> Clustering_jump_gap_cosmics(WCP::PR3DClusterSelection& live_clusters, WCP::PR3DClusterSelection& dead_clusters, std::map<int,std::pair<double,double>>& dead_u_index, std::map<int,std::pair<double,double>>& dead_v_index, std::map<int,std::pair<double,double>>& dead_w_index, WCP::DynamicToyPointCloud& global_point_cloud, bool flag_neutrino = true);

  
  // clustering live clusters associated with the same dead cluster
  void Clustering_live_dead(WCP::PR3DClusterSelection& live_clusters, WCP::PR3DClusterSelection& dead_clusters, std::map<WCP::PR3DCluster*,double>& cluster_length_map, std::set<WCP::PR3DCluster*>& cluster_connected_dead);
  // bool IsCrossing(WCP::Point& p, double theta, double phi, WCP::SlimMergeGeomCell *mcell, WCPSst::GeomDataSource& gds);

  
  

  
  

  void Clustering_regular(WCP::PR3DClusterSelection& live_clusters,  std::map<WCP::PR3DCluster*,double>& cluster_length_map, std::set<WCP::PR3DCluster*>& cluster_connected_dead, double length_cut = 45*units::cm, bool flag_enable_extend = true);
  bool Clustering_1st_round(WCP::PR3DCluster *cluster1, WCP::PR3DCluster *cluster2, double length_1, double length_2, double length_cut = 45*units::cm, bool flag_enable_extend = true);
 
  void Clustering_parallel_prolong(WCP::PR3DClusterSelection& live_clusters, std::map<WCP::PR3DCluster*,double>& cluster_length_map, std::set<WCP::PR3DCluster*>& cluster_connected_dead, double length_cut = 35*units::cm);
  bool Clustering_2nd_round(WCP::PR3DCluster *cluster1, WCP::PR3DCluster *cluster2, double length_1, double length_2, double length_cut = 35*units::cm);

  void Clustering_close(WCP::PR3DClusterSelection& live_clusters, std::map<WCP::PR3DCluster*,double>& cluster_length_map, std::set<WCP::PR3DCluster*>& cluster_connected_dead, double length_cut = 1*units::cm);
  bool Clustering_3rd_round(WCP::PR3DCluster *cluster1, WCP::PR3DCluster *cluster2, double length_1, double length_2, double length_cut = 1*units::cm);

  void Clustering_extend(WCP::PR3DClusterSelection& live_clusters, std::map<WCP::PR3DCluster*,double>& cluster_length_map, std::set<WCP::PR3DCluster*>& cluster_connected_dead, int flag, double length_cut = 150*units::cm, int num_try = 0, double length_2_cut = 3*units::cm, int num_dead_try=3);
  bool Clustering_4th_prol(WCP::PR3DCluster *cluster1, WCP::PR3DCluster *cluster2, double length_2, WCP::Point& earliest_p, TVector3& dir_earlp, double length_cut);
  bool Clustering_4th_para(WCP::PR3DCluster *cluster1, WCP::PR3DCluster *cluster2, double length_1, double length_2, WCP::Point& earliest_p, TVector3& dir_earlp, double length_cut);
  bool Clustering_4th_reg(WCP::PR3DCluster *cluster1, WCP::PR3DCluster *cluster2, double length_1, double length_2, WCP::Point p1, double length_cut);
  bool Clustering_4th_dead(WCP::PR3DCluster *cluster1, WCP::PR3DCluster *cluster2, double length_1, double length_2, double length_cut, int num_dead_try=3);
  
  
  //void Clustering_prolong(WCP::PR3DClusterSelection& live_clusters, std::map<WCP::PR3DCluster*,double>& cluster_length_map, double length_cut = 45*units::cm);

  
  void Clustering_separate(WCP::PR3DClusterSelection& live_clusters, std::map<WCP::PR3DCluster*,double>& cluster_length_map, std::map<int,std::pair<double,double>>& dead_u_index, std::map<int,std::pair<double,double>>& dead_v_index, std::map<int,std::pair<double,double>>& dead_w_index);
  void Clustering_separate(WCP::PR3DClusterSelection& live_clusters, std::map<WCP::PR3DCluster*,double>& cluster_length_map, std::map<int,std::pair<double,double>>& dead_u_index, std::map<int,std::pair<double,double>>& dead_v_index, std::map<int,std::pair<double,double>>& dead_w_index, WCP::ToyCTPointCloud& ct_point_cloud);

  
  bool JudgeSeparateDec_1(WCP::PR3DCluster* cluster, TVector3& drift_dir, double length, double time_slice_length);

  bool JudgeSeparateDec_2(WCP::PR3DCluster* cluster, TVector3& drift_dir ,std::vector<WCP::WCPointCloud<double>::WCPoint>& boundary_points, std::vector<WCP::WCPointCloud<double>::WCPoint>& independent_points, double cluster_length = 100*units::cm);
  
  std::vector<WCP::PR3DCluster*> Separate_1(WCP::PR3DCluster* cluster,std::vector<WCP::WCPointCloud<double>::WCPoint>& boundary_points, std::vector<WCP::WCPointCloud<double>::WCPoint>& independent_points, std::map<int,std::pair<double,double>>& dead_u_index, std::map<int,std::pair<double,double>>& dead_v_index, std::map<int,std::pair<double,double>>& dead_w_index, double length=0);
  std::vector<WCP::PR3DCluster*> Separate_1(WCP::ToyCTPointCloud& ct_point_cloud, WCP::PR3DCluster* cluster,std::vector<WCP::WCPointCloud<double>::WCPoint>& boundary_points, std::vector<WCP::WCPointCloud<double>::WCPoint>& independent_points, std::map<int,std::pair<double,double>>& dead_u_index, std::map<int,std::pair<double,double>>& dead_v_index, std::map<int,std::pair<double,double>>& dead_w_index, double length=0);
  

  std::vector<WCP::PR3DCluster*> Separate_2(WCP::PR3DCluster* cluster, double dis_cut = 5*units::cm);
  

  void Clustering_deghost(WCP::PR3DClusterSelection& live_clusters, std::map<WCP::PR3DCluster*,double>& cluster_length_map, std::map<int,std::pair<double,double>>& dead_u_index, std::map<int,std::pair<double,double>>& dead_v_index, std::map<int,std::pair<double,double>>& dead_w_index, double length_cut = 0);
  void Clustering_deghost(WCP::ToyCTPointCloud& ct_point_cloud, WCP::PR3DClusterSelection& live_clusters, std::map<WCP::PR3DCluster*,double>& cluster_length_map, std::map<int,std::pair<double,double>>& dead_u_index, std::map<int,std::pair<double,double>>& dead_v_index, std::map<int,std::pair<double,double>>& dead_w_index, double length_cut = 0);

  void Clustering_connect1(WCP::PR3DClusterSelection& live_clusters, std::map<WCP::PR3DCluster*,double>& cluster_length_map, WCP::DynamicToyPointCloud& global_point_cloud, std::map<int,std::pair<double,double>>& dead_u_index, std::map<int,std::pair<double,double>>& dead_v_index, std::map<int,std::pair<double,double>>& dead_w_index);
 

  void Clustering_examine_x_boundary(WCP::PR3DClusterSelection& live_clusters,std::map<WCP::PR3DCluster*,double>& cluster_length_map);

  void Clustering_protect_overclustering(WCP::PR3DClusterSelection& live_clusters,std::map<WCP::PR3DCluster*,double>& cluster_length_map, WCP::ToyCTPointCloud& ct_point_cloud);
  WCP::PR3DClusterSelection Examine_overclustering(WCP::PR3DCluster *cluster,  WCP::ToyCTPointCloud& ct_point_cloud);
  
  
  void Clustering_neutrino(WCP::PR3DClusterSelection& live_clusters,std::map<WCP::PR3DCluster*,double>& cluster_length_map, int num_try=0);
  
  void Clustering_dis(WCP::PR3DClusterSelection& live_clusters,std::map<WCP::PR3DCluster*,double>& cluster_length_map);
  

  std::map<WCP::PR3DCluster*,std::vector<std::pair<WCP::PR3DCluster*,double>>> Clustering_isolated(WCP::PR3DClusterSelection& live_clusters,std::map<WCP::PR3DCluster*,double>& cluster_length_map);
}

#endif
