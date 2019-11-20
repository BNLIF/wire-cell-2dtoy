#ifndef WIRECELL_IMPROVEPR3DCluster_H
#define WIRECELL_IMPROVEPR3DCluster_H

#include "WCPData/PR3DCluster.h"
#include "WCPData/ToyCTPointCloud.h"
#include "WCPSst/GeomDataSource.h"

namespace WCP2dToy{
  WCP::PR3DCluster* Improve_PR3DCluster(WCP::PR3DCluster* cluster, WCP::ToyCTPointCloud& ct_point_cloud,WCPSst::GeomDataSource& gds);
  
}

#endif 
