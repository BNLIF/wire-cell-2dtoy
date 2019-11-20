#ifndef WCP2dToy_EXAMINECLUSTER_H
#define WCP2dToy_EXAMINECLUSTER_H

#include "WCPData/FlashTPCBundle.h"
#include "WCPData/ToyCTPointCloud.h"

namespace WCP2dToy{
  WCP::FlashTPCBundleSelection ExamineBundles(WCP::FlashTPCBundleSelection bundles, WCP::ToyCTPointCloud& ct_point_cloud);
  WCP::FlashTPCBundle* ExamineBundle(WCP::FlashTPCBundle* bundle, std::set<int>& used_cluster_ids, WCP::ToyCTPointCloud& ct_point_cloud);
  
}

#endif
