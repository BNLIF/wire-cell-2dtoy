#ifndef WIRECELL_CALC_POINTS_H
#define WIRECELL_CALC_POINTS_H

#include "WireCellSst/GeomDataSource.h"
#include "WireCellData/PR3DCluster.h"

namespace WireCell2dToy{
  void calc_boundary_points_dead(WireCell::GeomDataSource& gds, WireCell::PR3DCluster* cluster);
  void calc_sampling_points(WireCell::GeomDataSource& gds, WireCell::PR3DCluster* cluster, int nrebin, int frame_length, double unit_dis);

  void calc_boundary_points_dead(WireCell::GeomDataSource& gds, WireCell::SlimMergeGeomCell* mcell);
  void calc_sampling_points(WireCell::GeomDataSource& gds, WireCell::SlimMergeGeomCell* mcell, int nrebin, int frame_length, double unit_dis);
  
}

#endif
