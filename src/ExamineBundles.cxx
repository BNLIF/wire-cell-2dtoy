#include "WireCell2dToy/ExamineBundles.h"

using namespace WireCell;

void WireCell2dToy::ExamineBundles(WireCell::FlashTPCBundleSelection& bundles, WireCell::ToyCTPointCloud& ct_point_cloud){

  std::set<int> used_cluster_ids;
  
  for (auto it = bundles.begin(); it!= bundles.end(); it++){
    FlashTPCBundle *bundle = *it;
    ExamineBundle(bundle, used_cluster_ids, ct_point_cloud);
  }
  
}

void WireCell2dToy::ExamineBundle(WireCell::FlashTPCBundle* bundle, std::set<int>& used_cluster_ids, WireCell::ToyCTPointCloud& ct_point_cloud){
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
  //actual code to separate clusters ...
  std::vector<SMGCSelection> sep_mcells = new_cluster->Examine_graph(ct_point_cloud);
  delete new_cluster;

  std::vector<PR3DCluster*> new_clusters;
  for (auto it = sep_mcells.begin(); it!= sep_mcells.end(); it++){
    PR3DCluster *new_cluster = new PR3DCluster(cluster_id);
    used_cluster_ids.insert(cluster_id);
    cluster_id++;
    new_clusters.push_back(new_cluster);
    for (auto it1 = (*it).begin(); it1!= (*it).end(); it1++){
      new_cluster->AddCell(*it1, (*it1)->GetTimeSlice());
    }
  }

  
  //set the new clusters ... 
  bundle->set_main_cluster(new_clusters.at(0));
  other_clusters.clear();
  more_clusters.clear();
  for (size_t i=1; i<new_clusters.size();i++){
    other_clusters.push_back(new_clusters.at(i));
    more_clusters.push_back(new_clusters.at(i));
  }
  
  // delete original cluster
  for (auto it = orig_clusters.begin(); it!=orig_clusters.end(); it++){
    delete (*it);
  }
  
  
}
