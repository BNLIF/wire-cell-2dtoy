#ifndef WireCell2dToy_TOYMATCHING_H
#define WireCell2dToy_TOYMATCHING_H

#include "WireCellData/PR3DCluster.h"
//#include "WireCell2dToy/uBooNE_light_reco.h"
#include "WireCell2dToy/ToyLightReco.h"

namespace WireCell2dToy{
  // time_offset in us
  int convert_xyz_voxel_id(WireCell::Point& p);

  std::vector<std::tuple<WireCell::PR3DCluster*, WireCell::Opflash*, double, std::vector<double>>> tpc_light_match(int time_offset, int nrebin, std::map<WireCell::PR3DCluster*,std::vector<std::pair<WireCell::PR3DCluster*,double>>>& group_clusters, WireCell::OpflashSelection& flashes);
}

#endif
