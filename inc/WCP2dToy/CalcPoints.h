#ifndef WIRECELL_CALC_POINTS_H
#define WIRECELL_CALC_POINTS_H

#include "WCPSst/GeomDataSource.h"
#include "WCPData/PR3DCluster.h"

namespace WCP2dToy{
  void calc_boundary_points_dead(WCP::GeomDataSource& gds, WCP::PR3DCluster* cluster);
  void calc_sampling_points(WCP::GeomDataSource& gds, WCP::PR3DCluster* cluster, int nrebin, int frame_length, double unit_dis);

  void calc_boundary_points_dead(WCP::GeomDataSource& gds, WCP::SlimMergeGeomCell* mcell);
  void calc_sampling_points(WCP::GeomDataSource& gds, WCP::SlimMergeGeomCell* mcell, int nrebin, int frame_length, double unit_dis);
  
}

#endif
