#ifndef WireCell2dToy_TOYCLUSTERING_H
#define WireCell2dToy_TOYCLUSTERING_H

#include "WireCellData/PR3DCluster.h"
#include "WireCellSst/GeomDataSource.h"

namespace WireCell2dToy{
  void Clustering_live_dead(WireCell::PR3DClusterSelection& live_clusters, WireCell::PR3DClusterSelection& dead_clusters, WireCellSst::GeomDataSource& gds);
  // bool IsCrossing(WireCell::Point& p, double theta, double phi, WireCell::SlimMergeGeomCell *mcell, WireCellSst::GeomDataSource& gds);
}

#endif
