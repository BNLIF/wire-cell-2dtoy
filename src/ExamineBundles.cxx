#include "WireCell2dToy/ExamineBundles.h"

using namespace WireCell;

FlashTPCBundleSelection WireCell2dToy::ExamineBundles(WireCell::FlashTPCBundleSelection bundles, WireCell::ToyCTPointCloud& ct_point_cloud){

  std::set<int> used_cluster_ids;
  FlashTPCBundleSelection new_bundles;
  
  for (auto it = bundles.begin(); it!= bundles.end(); it++){
    FlashTPCBundle *bundle = *it;
    FlashTPCBundle *new_bundle = ExamineBundle(bundle, used_cluster_ids, ct_point_cloud);
    new_bundles.push_back(new_bundle);  
  }

  std::set<PR3DCluster*> orig_clusters;
  for (size_t i=0;i!=bundles.size();i++){
    FlashTPCBundle *bundle = bundles.at(i);
    orig_clusters.insert(bundle->get_main_cluster());
    //delete bundle->get_main_cluster();
    for (size_t j=0;j!=bundle->get_other_clusters().size();j++){
      orig_clusters.insert(bundle->get_other_clusters().at(j));
      //delete bundle->get_other_clusters().at(j);
    }
    delete bundle;
  }
  for (auto it = orig_clusters.begin(); it!=orig_clusters.end(); it++){
    delete (*it);
  }
  
  return new_bundles;
  
}

WireCell::FlashTPCBundle* WireCell2dToy::ExamineBundle(WireCell::FlashTPCBundle* bundle, std::set<int>& used_cluster_ids, WireCell::ToyCTPointCloud& ct_point_cloud){
  
  int cluster_id;
  if (used_cluster_ids.size()==0){
    cluster_id = 1;
  }else{
    cluster_id = *used_cluster_ids.rbegin() + 1;
  }
  
  
  PR3DClusterSelection orig_clusters;
  orig_clusters.push_back(bundle->get_main_cluster());
  for (auto it1 = bundle->get_other_clusters().begin(); it1!=bundle->get_other_clusters().end();it1++){
    orig_clusters.push_back(*it1);
  }
  
  SMGCSelection orig_mcells;
  for (auto it = orig_clusters.begin(); it!=orig_clusters.end();it++){
    SMGCSelection& mcells = (*it)->get_mcells();
    //std::cout << (*it)->get_cluster_id() << " " << mcells.size() << std::endl;
    for (size_t i=0;i!=mcells.size();i++){
      orig_mcells.push_back(mcells.at(i));
    }
  }

  // create a new cluster
  PR3DCluster *new_cluster = new PR3DCluster(cluster_id);
   
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

  FlashTPCBundle *new_bundle =  new FlashTPCBundle(bundle->get_flash(), new_clusters.at(0), bundle->get_flash_index_id(), bundle->get_cluster_index_id());
  new_bundle->set_strength(bundle->get_strength());
  new_bundle->set_flag_close_to_PMT(bundle->get_flag_close_to_PMT());
  new_bundle->set_flag_at_x_boundary(bundle->get_flag_at_x_boundary());
  new_bundle->set_spec_end_flag(bundle->get_spec_end_flag());
  new_bundle->set_potential_bad_match_flag(bundle->get_potential_bad_match_flag());
  new_bundle->set_ks_dis(bundle->get_ks_dis());
  new_bundle->set_ndf(bundle->get_ndf());
  new_bundle->set_chi2(bundle->get_chi2());
  new_bundle->set_consistent_flag(new_bundle->get_consistent_flag());
  new_bundle->set_pred_pmt_light(bundle->get_pred_pmt_light());
  
  for (int i=1; i<new_clusters.size();i++){
    new_bundle->add_other_cluster(new_clusters.at(i));
  }
  
  return new_bundle;
}
