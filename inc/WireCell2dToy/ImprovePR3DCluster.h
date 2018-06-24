#ifndef WIRECELL_IMPROVEPR3DCluster_H
#define WIRECELL_IMPROVEPR3DCluster_H

#include "WireCellData/PR3DCluster.h"
#include "WireCellData/ToyCTPointCloud.h"
#include "WireCellSst/GeomDataSource.h"

namespace WireCell2dToy{
  WireCell::PR3DCluster* Improve_PR3DCluster(WireCell::PR3DCluster* cluster, WireCell::ToyCTPointCloud& ct_point_cloud,WireCellSst::GeomDataSource& gds);
  
}

#endif 
