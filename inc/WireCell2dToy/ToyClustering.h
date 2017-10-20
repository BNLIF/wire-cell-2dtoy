#ifndef WireCell2dToy_TOYCLUSTERING_H
#define WireCell2dToy_TOYCLUSTERING_H

#include "WireCellData/PR3DCluster.h"
//#include "WireCellSst/GeomDataSource.h"

namespace WireCell2dToy{
  void Clustering_live_dead(WireCell::PR3DClusterSelection& live_clusters, WireCell::PR3DClusterSelection& dead_clusters);
  // bool IsCrossing(WireCell::Point& p, double theta, double phi, WireCell::SlimMergeGeomCell *mcell, WireCellSst::GeomDataSource& gds);
  void Clustering_jump_gap_cosmics(WireCell::PR3DClusterSelection& live_clusters, double dis);
  std::pair<WireCell::SlimMergeGeomCell*, WireCell::SlimMergeGeomCell*> Get_Closest_MCells_Clusters(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2,double dis);
  
}

#endif
