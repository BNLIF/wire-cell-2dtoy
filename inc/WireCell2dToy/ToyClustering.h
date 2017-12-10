#ifndef WireCell2dToy_TOYCLUSTERING_H
#define WireCell2dToy_TOYCLUSTERING_H

#include "WireCellData/PR3DCluster.h"
//#include "WireCellSst/GeomDataSource.h"

namespace WireCell2dToy{
  // clustering live clusters associated with the same dead cluster
  void Clustering_live_dead(WireCell::PR3DClusterSelection& live_clusters, WireCell::PR3DClusterSelection& dead_clusters);
  // bool IsCrossing(WireCell::Point& p, double theta, double phi, WireCell::SlimMergeGeomCell *mcell, WireCellSst::GeomDataSource& gds);

  
  void Clustering_jump_gap_cosmics(WireCell::PR3DClusterSelection& live_clusters);

  void Clustering_regular(WireCell::PR3DClusterSelection& live_clusters,  std::map<WireCell::PR3DCluster*,double>& cluster_length_map, double length_cut = 45*units::cm, bool flag_enable_extend = true);
  bool Clustering_1st_round(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double length_1, double length_2, double length_cut = 45*units::cm, bool flag_enable_extend = true);
 
  void Clustering_parallel_prolong(WireCell::PR3DClusterSelection& live_clusters, std::map<WireCell::PR3DCluster*,double>& cluster_length_map, double length_cut = 35*units::cm);
  bool Clustering_2nd_round(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double length_1, double length_2, double length_cut = 35*units::cm);

  void Clustering_close(WireCell::PR3DClusterSelection& live_clusters, std::map<WireCell::PR3DCluster*,double>& cluster_length_map, double length_cut = 1*units::cm);
  bool Clustering_3rd_round(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double length_1, double length_2, double length_cut = 1*units::cm);
  
  
  //void Clustering_prolong(WireCell::PR3DClusterSelection& live_clusters, std::map<WireCell::PR3DCluster*,double>& cluster_length_map, double length_cut = 45*units::cm);


  
  
  

  std::map<WireCell::PR3DCluster*,std::vector<std::pair<WireCell::PR3DCluster*,double>>> Clustering_isolated(WireCell::PR3DClusterSelection& live_clusters);
}

#endif
