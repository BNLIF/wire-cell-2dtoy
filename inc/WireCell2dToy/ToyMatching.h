#ifndef WireCell2dToy_TOYMATCHING_H
#define WireCell2dToy_TOYMATCHING_H

#include "WireCellData/PR3DCluster.h"
//#include "WireCell2dToy/uBooNE_light_reco.h"
#include "WireCell2dToy/ToyLightReco.h"
#include "WireCellData/FlashTPCBundle.h"

namespace WireCell2dToy{
  // time_offset in us
  int convert_xyz_voxel_id(WireCell::Point& p);

  WireCell::FlashTPCBundleSelection tpc_light_match(int time_offset, int nrebin, std::map<WireCell::PR3DCluster*,std::vector<std::pair<WireCell::PR3DCluster*,double>>>& group_clusters, WireCell::OpflashSelection& flashes, Int_t runno = 0, bool flag_data = true);

  // WireCell::FlashTPCBundleSelection tpc_light_match_ana(int time_offset, int nrebin, std::map<WireCell::PR3DCluster*,std::vector<std::pair<WireCell::PR3DCluster*,double>>>& group_clusters, WireCell::OpflashSelection& flashes);

  void organize_matched_bundles(WireCell::FlashTPCBundleSelection& results_bundles, Double_t *cos_pe_low, Double_t *cos_pe_mid, std::map<std::pair<WireCell::Opflash*,WireCell::PR3DCluster*>,WireCell::FlashTPCBundle*>& fc_bundles_map);
  
  
}

#endif
