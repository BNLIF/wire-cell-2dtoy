#ifndef WireCell2dToy_EXAMINECLUSTER_H
#define WireCell2dToy_EXAMINECLUSTER_H

#include "WireCellData/FlashTPCBundle.h"
#include "WireCellData/ToyCTPointCloud.h"

namespace WireCell2dToy{
  WireCell::FlashTPCBundleSelection ExamineBundles(WireCell::FlashTPCBundleSelection bundles, WireCell::ToyCTPointCloud& ct_point_cloud);
  WireCell::FlashTPCBundle* ExamineBundle(WireCell::FlashTPCBundle* bundle, std::set<int>& used_cluster_ids, WireCell::ToyCTPointCloud& ct_point_cloud);
  
}

#endif
