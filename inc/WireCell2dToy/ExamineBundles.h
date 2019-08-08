#ifndef WireCell2dToy_EXAMINECLUSTER_H
#define WireCell2dToy_EXAMINECLUSTER_H

#include "WireCellData/FlashTPCBundle.h"

namespace WireCell2dToy{
  void ExamineBundles(WireCell::FlashTPCBundleSelection& bundles);
  void ExamineBundle(WireCell::FlashTPCBundle* bundle, std::set<int>& used_cluster_ids);
  
}

#endif
