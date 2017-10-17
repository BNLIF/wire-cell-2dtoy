#ifndef WIRECELL_CALC_POINTS_H
#define WIRECELL_CALC_POINTS_H

#include "WireCellSst/GeomDataSource.h"
#include "WireCellData/PR3DCluster.h"

namespace WireCell2dToy{
  void calc_boundary_points(WireCell::GeomDataSource& gds, WireCell::PR3DCluster* cluster);
  void calc_sampling_points(WireCell::GeomDataSource& gds, WireCell::PR3DCluster* cluster);

  void calc_boundary_points(WireCell::GeomDataSource& gds, WireCell::SlimMergeGeomCell* mcell);
  void calc_sampling_points(WireCell::GeomDataSource& gds, WireCell::SlimMergeGeomCell* mcell);
  
}

#endif
