#include "WireCell2dToy/ToyClustering.h"

using namespace WireCell;

void WireCell2dToy::Clustering_live_dead(WireCell::PR3DClusterSelection& live_clusters, WireCell::PR3DClusterSelection& dead_clusters, std::map<WireCell::PR3DCluster*,std::vector<WireCell::PR3DCluster*>>& dead_live_cluster_mapping, std::map<WireCell::PR3DCluster*,std::vector<WireCell::SMGCSelection>>& dead_live_mcells_mapping){
  std::set<std::pair<PR3DCluster*, PR3DCluster*>> tested_pairs;
  std::set<std::pair<PR3DCluster*, PR3DCluster*>> to_be_merged_pairs;

  
  
  for (auto it = dead_live_cluster_mapping.begin(); it!= dead_live_cluster_mapping.end(); it++){
    PR3DCluster* the_dead_cluster = (*it).first;
    std::vector<PR3DCluster*> connected_live_clusters = (*it).second;
    std::vector<SMGCSelection> connected_live_mcells = dead_live_mcells_mapping[the_dead_cluster];
    if (connected_live_clusters.size()>1){
      //      std::cout << the_dead_cluster->get_cluster_id() << " " << connected_live_clusters.size() << std::endl;
      for (size_t i=0; i!= connected_live_clusters.size(); i++){
        PR3DCluster* cluster_1 = connected_live_clusters.at(i);
	SMGCSelection mcells_1 = connected_live_mcells.at(i);
	for (size_t j=i+1;j<connected_live_clusters.size(); j++){
	  PR3DCluster* cluster_2 = connected_live_clusters.at(j);
	  SMGCSelection mcells_2 = connected_live_mcells.at(j);

	  if (tested_pairs.find(std::make_pair(cluster_1,cluster_2))==tested_pairs.end()){
	    cluster_1->Create_point_cloud();
	    cluster_2->Create_point_cloud();
	    //std::cout << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << std::endl;
	    tested_pairs.insert(std::make_pair(cluster_1,cluster_2));
	    tested_pairs.insert(std::make_pair(cluster_2,cluster_1));
	    // starting the test ... 

	    bool flag_merge = false;
	    for (auto it1 = mcells_1.begin(); it1!= mcells_1.end(); it1++){
	      SlimMergeGeomCell* mcell_1 = (*it1);
	      Point mcell1_center = mcell_1->center();
	      for (auto it2 = mcells_2.begin(); it2!=mcells_2.end(); it2++){
		SlimMergeGeomCell* mcell_2 = (*it2);
		Point mcell2_center = mcell_2->center();
		//	std::cout << mcell1_center.x/units::cm << " " << mcell1_center.y/units::cm << " " << mcell1_center.z/units::cm << " " << mcell2_center.x/units::cm << " " << mcell2_center.y/units::cm << " " << mcell2_center.z/units::cm << std::endl;
		// test both sides
		

		
	      }
	      if (flag_merge) break;
	    }
	    if (flag_merge)
	      to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	  }
	}
      }
    }
  }

  
}
