#include "WireCellData/DynamicToyPointCloud.h"

using namespace WireCell;


void WireCell2dToy::Clustering_deghost(WireCell::PR3DClusterSelection& live_clusters, std::map<WireCell::PR3DCluster*,double>& cluster_length_map, std::set<int>& dead_u_index, std::set<int>& dead_v_index, std::set<int>& dead_w_index){
  
  // sort the clusters length ...
  {
    std::vector<std::pair<PR3DCluster*,double>> temp_pair_vec;
    for (auto it = cluster_length_map.begin(); it!=cluster_length_map.end();it++){
      temp_pair_vec.push_back(std::make_pair(it->first, it->second));
    }
    sort(temp_pair_vec.begin(), temp_pair_vec.end(), sortbysec);
    live_clusters.clear();
    for (auto it = temp_pair_vec.begin(); it!=temp_pair_vec.end();it++){
      live_clusters.push_back(it->first);
    }
  }

  // Create two point clouds ... 
  // One for the points ... point --> index --> cluster (vector) ...
  // The other for the skeleton of each track ...  point --> index --> cluster (vector)
  // Both cloud needs to be dynamic, keep adding things into it as we improve the knowledge

  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();

  
  DynamicToyPointCloud global_point_cloud(angle_u,angle_v,angle_w);
  DynamicToyPointCloud global_skeleton_cloud(angle_u,angle_v,angle_w);

  std::vector<PR3DCluster*> to_be_removed_clusters;
  std::set<std::pair<PR3DCluster*, PR3DCluster*>> to_be_merged_pairs;
  
  for (size_t i=0;i!=live_clusters.size();i++){
    if (i==0){
      // fill anyway ...
      live_clusters.at(i)->Create_point_cloud();
      global_point_cloud.AddPoints(live_clusters.at(i),0);
      if (cluster_length_map[live_clusters.at(i)]>30*units::cm){
	live_clusters.at(i)->Construct_skeleton();
	global_skeleton_cloud.AddPoints(live_clusters.at(i),1);
      }
    }else{ 
      // start the process to add things in and perform deghosting ... 
      
      // two cases, merge clusters or remove clusters
      
    }
  }



  //merge clusters
  std::vector<std::set<PR3DCluster*>> merge_clusters;
  for (auto it = to_be_merged_pairs.begin(); it!=to_be_merged_pairs.end(); it++){
    PR3DCluster *cluster1 = (*it).first;
    PR3DCluster *cluster2 = (*it).second;
    //  std::cout << cluster1 << " " << cluster2 << " " << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << std::endl;
    
    bool flag_new = true;
    std::vector<std::set<PR3DCluster*>> temp_set;
    for (auto it1 = merge_clusters.begin(); it1!=merge_clusters.end(); it1++){
      std::set<PR3DCluster*>& clusters = (*it1);
      if (clusters.find(cluster1)!=clusters.end() ||
	  clusters.find(cluster2)!=clusters.end()){
	clusters.insert(cluster1);
	clusters.insert(cluster2);
	flag_new = false;
	temp_set.push_back(clusters);
	//break;
      }
    }
    if (flag_new){
      std::set<PR3DCluster*> clusters;
      clusters.insert(cluster1);
      clusters.insert(cluster2);
      merge_clusters.push_back(clusters);
    }
    if (temp_set.size()>1){
      // merge them further ...
      std::set<PR3DCluster*> clusters;
      for (size_t i=0;i!=temp_set.size();i++){
	for (auto it1 = temp_set.at(i).begin(); it1!= temp_set.at(i).end(); it1++){
	  clusters.insert(*it1);
	}
	merge_clusters.erase(find(merge_clusters.begin(),merge_clusters.end(),temp_set.at(i)));
      }
      merge_clusters.push_back(clusters);
    }
  }

  // merge clusters into new clusters, delete old clusters 
  for (auto it = merge_clusters.begin(); it!=merge_clusters.end();it++){
    std::set<PR3DCluster*>& clusters = (*it);
    PR3DCluster *ncluster = new PR3DCluster((*clusters.begin())->get_cluster_id());
    live_clusters.push_back(ncluster);
    for (auto it1 = clusters.begin(); it1!=clusters.end(); it1++){
      PR3DCluster *ocluster = *(it1);
      // std::cout << ocluster->get_cluster_id() << " ";
      SMGCSelection& mcells = ocluster->get_mcells();
      for (auto it2 = mcells.begin(); it2!=mcells.end(); it2++){
  	SlimMergeGeomCell *mcell = (*it2);
  	//std::cout << ocluster->get_cluster_id() << " " << mcell << std::endl;
  	int time_slice = mcell->GetTimeSlice();
  	ncluster->AddCell(mcell,time_slice);
      }
      live_clusters.erase(find(live_clusters.begin(), live_clusters.end(), ocluster));
      cluster_length_map.erase(ocluster);
      delete ocluster;
    }
    std::vector<int> range_v1 = ncluster->get_uvwt_range();
    double length_1 = sqrt(2./3. * (pow(pitch_u*range_v1.at(0),2) + pow(pitch_v*range_v1.at(1),2) + pow(pitch_w*range_v1.at(2),2)) + pow(time_slice_width*range_v1.at(3),2));
    cluster_length_map[ncluster] = length_1;
    //std::cout << std::endl;
  }
  
  // delete clusters ... 
  for (auto it = to_be_removed_clusters.begin(); it!=to_be_removed_clusters.end(); it++){
    PR3DCluster *ocluster = *it;
    live_clusters.erase(find(live_clusters.begin(), live_clusters.end(), ocluster));
    cluster_length_map.erase(ocluster);
    delete ocluster;
  }

  

  
  /* Point p(0,0,0); */

  /* std::tuple<double, PR3DCluster*, size_t> test1 = global_point_cloud.get_closest_point_info(p); */
  /* std::tuple<double, PR3DCluster*, size_t> test2 = global_point_cloud.get_closest_2d_point_info(p,0); */
  /* std::tuple<double, PR3DCluster*, size_t> test3 = global_point_cloud.get_closest_2d_point_info(p,1); */
  /* std::tuple<double, PR3DCluster*, size_t> test4 = global_point_cloud.get_closest_2d_point_info(p,2); */

  /* std::cout << std::get<0>(test1)/units::cm << " " << std::get<1>(test1)->get_cluster_id() << " " */
  /* 	    << std::get<0>(test2)/units::cm << " " << std::get<1>(test2)->get_cluster_id() << " " */
  /* 	    << std::get<0>(test3)/units::cm << " " << std::get<1>(test3)->get_cluster_id() << " " */
  /* 	    << std::get<0>(test4)/units::cm << " " << std::get<1>(test4)->get_cluster_id() << " " */
  /* 	    << std::endl; */


  /* test1 = global_skeleton_cloud.get_closest_point_info(p); */
  /* test2 = global_skeleton_cloud.get_closest_2d_point_info(p,0); */
  /* test3 = global_skeleton_cloud.get_closest_2d_point_info(p,1); */
  /* test4 = global_skeleton_cloud.get_closest_2d_point_info(p,2); */
  
  /* std::cout << std::get<0>(test1)/units::cm << " " << std::get<1>(test1)->get_cluster_id() << " " */
  /* 	    << std::get<0>(test2)/units::cm << " " << std::get<1>(test2)->get_cluster_id() << " " */
  /* 	    << std::get<0>(test3)/units::cm << " " << std::get<1>(test3)->get_cluster_id() << " " */
  /* 	    << std::get<0>(test4)/units::cm << " " << std::get<1>(test4)->get_cluster_id() << " " */
  /* 	    << std::endl; */
}
