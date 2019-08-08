#include "WireCell2dToy/ExamineBundles.h"

using namespace WireCell;

void WireCell2dToy::ExamineBundles(WireCell::FlashTPCBundleSelection& bundles){

  std::set<int> used_cluster_ids;
  
  for (auto it = bundles.begin(); it!= bundles.end(); it++){
    FlashTPCBundle *bundle = *it;
    ExamineBundle(bundle, used_cluster_ids);
  }
  
}

void WireCell2dToy::ExamineBundle(WireCell::FlashTPCBundle* bundle, std::set<int>& used_cluster_ids){
  int cluster_id;
  if (used_cluster_ids.size()==0){
    cluster_id = 1;
  }else{
    cluster_id = *used_cluster_ids.rbegin() + 1;
  }
  
  
  std::set<PR3DCluster*> orig_clusters;
  orig_clusters.insert(bundle->get_main_cluster());

  PR3DClusterSelection& other_clusters = bundle->get_other_clusters();
  PR3DClusterSelection& more_clusters = bundle->get_more_clusters();
  for (auto it = other_clusters.begin(); it!=other_clusters.end();it++){
    orig_clusters.insert(*it);
  }
  for (auto it = more_clusters.begin(); it!=more_clusters.end();it++){
    orig_clusters.insert(*it);
  }

  // create a new cluster
  PR3DCluster *new_cluster = new PR3DCluster(cluster_id);
  used_cluster_ids.insert(cluster_id);
  cluster_id++;
  std::set<SlimMergeGeomCell*> orig_mcells;
  for (auto it = orig_clusters.begin(); it!=orig_clusters.end();it++){
    SMGCSelection& mcells = (*it)->get_mcells();
    for (auto it1 = mcells.begin(); it1 != mcells.end();it1++){
      orig_mcells.insert(*it1);
    }
  }
  for (auto it = orig_mcells.begin(); it!=orig_mcells.end(); it++){
    new_cluster->AddCell(*it, (*it)->GetTimeSlice());
  }

  //set the new clusters ... 
  bundle->set_main_cluster(new_cluster);
  other_clusters.clear();
  more_clusters.clear();
  
  // delete original cluster
  for (auto it = orig_clusters.begin(); it!=orig_clusters.end(); it++){
    delete (*it);
  }
  
  
}
