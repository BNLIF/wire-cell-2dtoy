#ifndef WireCell2dToy_TOYCLUSTERING_H
#define WireCell2dToy_TOYCLUSTERING_H

#include "WireCellData/PR3DCluster.h"
//#include "WireCellSst/GeomDataSource.h"

namespace WireCell2dToy{
  double Find_Closeset_Points(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2,double length_1, double length_2, double length_cut, WireCell::SlimMergeGeomCell *mcell1, WireCell::SlimMergeGeomCell *mcell2, WireCell::Point& p1, WireCell::Point &p2);

  
  double cal_proj_angle_diff(TVector3& dir1, TVector3& dir2, double plane_angle);
  bool is_angle_consistent(TVector3& dir1, TVector3& dir2, bool same_direction, double angle_cut, double angle_u, double angle_v, double angle_w);
  
  
  void Clustering_jump_gap_cosmics(WireCell::PR3DClusterSelection& live_clusters, WireCell::PR3DClusterSelection& dead_clusters);

  
  // clustering live clusters associated with the same dead cluster
  void Clustering_live_dead(WireCell::PR3DClusterSelection& live_clusters, WireCell::PR3DClusterSelection& dead_clusters, std::map<WireCell::PR3DCluster*,double>& cluster_length_map, std::set<WireCell::PR3DCluster*>& cluster_connected_dead);
  // bool IsCrossing(WireCell::Point& p, double theta, double phi, WireCell::SlimMergeGeomCell *mcell, WireCellSst::GeomDataSource& gds);

  
  

  
  

  void Clustering_regular(WireCell::PR3DClusterSelection& live_clusters,  std::map<WireCell::PR3DCluster*,double>& cluster_length_map, std::set<WireCell::PR3DCluster*>& cluster_connected_dead, double length_cut = 45*units::cm, bool flag_enable_extend = true);
  bool Clustering_1st_round(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double length_1, double length_2, double length_cut = 45*units::cm, bool flag_enable_extend = true);
 
  void Clustering_parallel_prolong(WireCell::PR3DClusterSelection& live_clusters, std::map<WireCell::PR3DCluster*,double>& cluster_length_map, std::set<WireCell::PR3DCluster*>& cluster_connected_dead, double length_cut = 35*units::cm);
  bool Clustering_2nd_round(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double length_1, double length_2, double length_cut = 35*units::cm);

  void Clustering_close(WireCell::PR3DClusterSelection& live_clusters, std::map<WireCell::PR3DCluster*,double>& cluster_length_map, std::set<WireCell::PR3DCluster*>& cluster_connected_dead, double length_cut = 1*units::cm);
  bool Clustering_3rd_round(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double length_1, double length_2, double length_cut = 1*units::cm);

  void Clustering_extend(WireCell::PR3DClusterSelection& live_clusters, std::map<WireCell::PR3DCluster*,double>& cluster_length_map, std::set<WireCell::PR3DCluster*>& cluster_connected_dead, int flag, double length_cut = 150*units::cm, int num_try = 0);
  bool Clustering_4th_prol(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double length_2, WireCell::Point& earliest_p, TVector3& dir_earlp, double length_cut);
  bool Clustering_4th_para(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double length_2, WireCell::Point& earliest_p, TVector3& dir_earlp, double length_cut);
  bool Clustering_4th_reg(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double length_1, double length_2, WireCell::Point p1, double length_cut);
  bool Clustering_4th_dead(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double length_1, double length_2, double length_cut);
  
  
  //void Clustering_prolong(WireCell::PR3DClusterSelection& live_clusters, std::map<WireCell::PR3DCluster*,double>& cluster_length_map, double length_cut = 45*units::cm);


  
  
  

  std::map<WireCell::PR3DCluster*,std::vector<std::pair<WireCell::PR3DCluster*,double>>> Clustering_isolated(WireCell::PR3DClusterSelection& live_clusters);
}

#endif
