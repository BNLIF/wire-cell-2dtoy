std::map<PR3DCluster*,std::vector<std::pair<PR3DCluster*,double>>> WireCell2dToy::Clustering_isolated(WireCell::PR3DClusterSelection& live_clusters){
  PR3DClusterSelection big_clusters;
  PR3DClusterSelection small_clusters;
  for (size_t i=0;i!=live_clusters.size();i++){
    std::vector<int> ranges = live_clusters.at(i)->get_uvwt_range();
    int max = 0;
    for (int j=0;j!=4;j++){
      if (ranges.at(j)>max)
	max =ranges.at(j);
    }
    //    std::cout << ranges.at(0) << " " << ranges.at(1) << " " << ranges.at(2) << " " << ranges.at(3) << " " << max << std::endl;
    if (max < 80){
      small_clusters.push_back(live_clusters.at(i));
    }else{
      big_clusters.push_back(live_clusters.at(i));
    }
  }

  std::map<PR3DCluster*,std::vector<std::pair<PR3DCluster*,double>>> results;
  // Now for each small clusters, find the big clusters with the shortest distance and connect them
  for (auto it = small_clusters.begin(); it!= small_clusters.end(); it++){
    PR3DCluster* curr_cluster = (*it);
    curr_cluster->Create_point_cloud();
    ToyPointCloud *cloud1 = curr_cluster->get_point_cloud();
    double min_dis = 1e9; PR3DCluster* min_dis_cluster = 0;
    for (auto it1 = big_clusters.begin(); it1!=big_clusters.end(); it1++){
      PR3DCluster* big_cluster = (*it1);
      big_cluster->Create_point_cloud();
      ToyPointCloud *cloud2 = big_cluster->get_point_cloud();
      
      /* // pick any point and merged cell in cluster1, */
      /* SlimMergeGeomCell *prev_mcell1 = 0; */
      /* SlimMergeGeomCell *prev_mcell2 = 0; */
      /* SlimMergeGeomCell *mcell1 = curr_cluster->get_mcells().at(0); */
      /* Point p1 = mcell1->center(); */
      /* SlimMergeGeomCell *mcell2=0; */
      /* Point p2; */

      /* while(mcell1!=prev_mcell1 || mcell2!=prev_mcell2){ */
      /* 	prev_mcell1 = mcell1; */
      /* 	prev_mcell2 = mcell2; */
    
      /* 	// find the closest point and merged cell in cluster2 */
      /* 	std::pair<SlimMergeGeomCell*,Point> temp_results = big_cluster->get_closest_point_mcell(p1); */
      /* 	p2 = temp_results.second; */
      /* 	mcell2 = temp_results.first; */
      /* 	// find the closest point and merged cell in cluster1 */
      /* 	temp_results = curr_cluster->get_closest_point_mcell(p2); */
      /* 	p1 = temp_results.second; */
      /* 	mcell1 = temp_results.first; */
      /* } */
      /* double dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
 */

      std::tuple<int,int,double> results =  cloud2->get_closest_points(cloud1);
      double dis = std::get<2>(results);
      
      if (dis < min_dis){
	min_dis = dis;
	min_dis_cluster = big_cluster;
      }
    }
    if (min_dis < 0.8*units::m){
      if (results.find(min_dis_cluster)==results.end()){
	std::vector<std::pair<PR3DCluster*,double>> temp;
	temp.push_back(std::make_pair(curr_cluster,min_dis));
	results[min_dis_cluster]=temp;
      }else{
	results[min_dis_cluster].push_back(std::make_pair(curr_cluster,min_dis));
      }
    }else{
      big_clusters.push_back(curr_cluster);
      std::vector<std::pair<PR3DCluster*,double>> temp;
      results[curr_cluster] = temp;
    }
  }
  
  return results;
}
