#ifndef WireCell2dToy_TOYCLUSTERING_H
#define WireCell2dToy_TOYCLUSTERING_H

#include "WireCellData/PR3DCluster.h"

namespace WireCell2dToy{
  void Clustering_live_dead(WireCell::PR3DClusterSelection& live_clusters, WireCell::PR3DClusterSelection& dead_clusters, std::map<WireCell::PR3DCluster*,std::vector<WireCell::PR3DCluster*>>& dead_live_cluster_mapping, std::map<WireCell::PR3DCluster*,std::vector<WireCell::SMGCSelection>>& dead_live_mcells_mapping);
}

#endif
