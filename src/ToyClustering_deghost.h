

using namespace WCP;


void WCP2dToy::Clustering_deghost(WCP::PR3DClusterSelection& live_clusters, std::map<WCP::PR3DCluster*,double>& cluster_length_map, std::map<int,std::pair<double,double>>& dead_u_index, std::map<int,std::pair<double,double>>& dead_v_index, std::map<int,std::pair<double,double>>& dead_w_index, double length_cut){

  
  
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
      PR3DCluster* cluster = live_clusters.at(i);
      
      if (length_cut == 0 || cluster_length_map[live_clusters.at(i)] < length_cut)  {
	cluster->Create_point_cloud();
	WCP::WCPointCloud<double>& cloud = cluster->get_point_cloud()->get_cloud();
	int num_total_points = cloud.pts.size(); // total number of points
	int num_dead[3]={0,0,0}; // dead wires in each view 
	int num_unique[3]={0,0,0}; // points that are unique (not agree with any other clusters)
	std::map<PR3DCluster*, int> map_cluster_num[3];
	
	double dis_cut = 1.2*units::cm;
	
	
	
	for (size_t j=0;j!=num_total_points;j++){
	  Point test_point(cloud.pts.at(j).x,cloud.pts.at(j).y,cloud.pts.at(j).z);
	  
	  bool flag_dead = false;
	  
	  if (dead_u_index.find(cloud.pts.at(j).index_u)!=dead_u_index.end()){
	    if (cloud.pts.at(j).x >= dead_u_index[cloud.pts.at(j).index_u].first &&
		cloud.pts.at(j).x <= dead_u_index[cloud.pts.at(j).index_u].second){
	      flag_dead = true;
	    }
	  }
	  
	  if (!flag_dead){
	    std::tuple<double, PR3DCluster*, size_t> results = global_point_cloud.get_closest_2d_point_info(test_point, 0);
	    if (std::get<0>(results)<=dis_cut/3.){
	      if (map_cluster_num[0].find(std::get<1>(results))==map_cluster_num[0].end()){
		map_cluster_num[0][std::get<1>(results)] = 1;
	      }else{
		map_cluster_num[0][std::get<1>(results)] ++;
	      }  
	    }else{
	      results = global_skeleton_cloud.get_closest_2d_point_info(test_point, 0);
	      if (std::get<0>(results)<=dis_cut*2.0){
		if (map_cluster_num[0].find(std::get<1>(results))==map_cluster_num[0].end()){
		  map_cluster_num[0][std::get<1>(results)] = 1;
		}else{
		  map_cluster_num[0][std::get<1>(results)] ++;
		}
	      }else{
		num_unique[0]++;
	      }
	    }
	  }else{
	    num_dead[0]++;
	  }
	  
	  
	  flag_dead = false;
	  
	  if (dead_v_index.find(cloud.pts.at(j).index_v)!=dead_v_index.end()){
	    if (cloud.pts.at(j).x >= dead_v_index[cloud.pts.at(j).index_v].first &&
		cloud.pts.at(j).x <= dead_v_index[cloud.pts.at(j).index_v].second){
	      flag_dead = true;
	    }
	  }
	  
	  
	  if (!flag_dead){
	    std::tuple<double, PR3DCluster*, size_t> results = global_point_cloud.get_closest_2d_point_info(test_point, 1);
	    if (std::get<0>(results)<=dis_cut/3.){
	      if (map_cluster_num[1].find(std::get<1>(results))==map_cluster_num[1].end()){
		map_cluster_num[1][std::get<1>(results)] = 1;
	      }else{
		map_cluster_num[1][std::get<1>(results)] ++;
	      }
	    }else{
	      results = global_skeleton_cloud.get_closest_2d_point_info(test_point, 1);
	      if (std::get<0>(results)<=dis_cut*2.0){
		if (map_cluster_num[1].find(std::get<1>(results))==map_cluster_num[1].end()){
		  map_cluster_num[1][std::get<1>(results)] = 1;
		}else{
		  map_cluster_num[1][std::get<1>(results)] ++;
		}
	      }else{
		num_unique[1]++;
	      }
	    }
	  }else{
	    num_dead[1]++;
	  }
	  
	  flag_dead = false;
	  
	  if (dead_w_index.find(cloud.pts.at(j).index_w)!=dead_w_index.end()){
	    if (cloud.pts.at(j).x >= dead_w_index[cloud.pts.at(j).index_w].first &&
		cloud.pts.at(j).x <= dead_w_index[cloud.pts.at(j).index_w].second){
	      flag_dead = true;
	    }
	  }
	  
	  
	  if (!flag_dead){
	    std::tuple<double, PR3DCluster*, size_t> results = global_point_cloud.get_closest_2d_point_info(test_point, 2);
	    if (std::get<0>(results)<=dis_cut/3.){
	      if (map_cluster_num[2].find(std::get<1>(results))==map_cluster_num[2].end()){
		map_cluster_num[2][std::get<1>(results)] = 1;
	      }else{
		map_cluster_num[2][std::get<1>(results)] ++;
	      }
	    }else{
	      results = global_skeleton_cloud.get_closest_2d_point_info(test_point, 2);
	      if (std::get<0>(results)<=dis_cut*2.0){
		if (map_cluster_num[2].find(std::get<1>(results))==map_cluster_num[2].end()){
		  map_cluster_num[2][std::get<1>(results)] = 1;
		}else{
		  map_cluster_num[2][std::get<1>(results)] ++;
		}
	      }else{
		num_unique[2]++;
	      }
	    }
	  }else{
	    num_dead[2]++;
	  }
	  
	  
	  if ( (num_unique[1]+num_unique[0]+num_unique[2]) >= 0.24 * num_total_points &&
	       (num_unique[1]+num_unique[0]+num_unique[2]) > 25 )
	    break;
	}
	
	//std::cout << cluster->get_cluster_id() << " " << num_unique[0] << " " << num_unique[1] << " " << num_unique[2] << " " << num_dead[0] << " " << num_dead[1] << " " << num_dead[2] << " " << num_total_points << std::endl;
	
	bool flag_save = false;
	
	if (((num_unique[0] <= 0.1 * (num_total_points - num_dead[0]) || num_unique[0] <= 0.1 * num_total_points && num_unique[0] <= 8) &&
	     (num_unique[1] <= 0.1 * (num_total_points - num_dead[1]) || num_unique[1] <= 0.1 * num_total_points && num_unique[1] <= 8) &&
	     (num_unique[2] <= 0.1 * (num_total_points - num_dead[2]) || num_unique[2] <= 0.1 * num_total_points && num_unique[2] <= 8) &&
	     ((num_unique[0] + num_unique[1] + num_unique[2]) <= 0.05 * (num_total_points - num_dead[0] + num_total_points - num_dead[1] + num_total_points - num_dead[2]) || (num_unique[0] + num_unique[1] + num_unique[2]) < 0.15 * num_total_points && (num_unique[0] + num_unique[1] + num_unique[2]) <=8) ||
	     num_unique[0]==0 && num_unique[1]==0 && num_unique[2] < 0.24 * num_total_points ||
	     num_unique[0]==0 && num_unique[2]==0 && num_unique[1] < 0.24 * num_total_points ||
	     num_unique[2]==0 && num_unique[1]==0 && num_unique[0] < 0.24 * num_total_points ||
	     num_unique[0]==0 && (num_unique[1]+num_unique[2]) < 0.12 * num_total_points * 2 || 
	     num_unique[1]==0 && (num_unique[0]+num_unique[2]) < 0.12 * num_total_points * 2 || 
	     num_unique[2]==0 && (num_unique[1]+num_unique[0]) < 0.12 * num_total_points * 2 ||
	     (num_unique[1]+num_unique[0]+num_unique[2]) < 0.24 * num_total_points &&
	     (num_unique[1]<0.02 * num_total_points || num_unique[0]<0.02 * num_total_points
	      || num_unique[2]<0.02 * num_total_points)
	     ) &&
	    (num_unique[0] + num_unique[1] + num_unique[2])  <= 500 ){
	  flag_save = false;
	  
	  // now try to compare
	  // find the maximal for each map
	  PR3DCluster *max_cluster_u = 0, *max_cluster_v=0, *max_cluster_w=0;
	  int max_value_u = 0, max_value_v = 0, max_value_w = 0;
	  for (auto it = map_cluster_num[0].begin(); it!=map_cluster_num[0].end(); it++){
	    if (it->second > max_value_u){
	      max_value_u = it->second;
	      max_cluster_u = it->first;
	    }
	  }
	  for (auto it = map_cluster_num[1].begin(); it!=map_cluster_num[1].end(); it++){
	    if (it->second > max_value_v){
	      max_value_v = it->second;
	      max_cluster_v = it->first;
	    }
	  }
	  for (auto it = map_cluster_num[2].begin(); it!=map_cluster_num[2].end(); it++){
	    if (it->second > max_value_w){
	      max_value_w = it->second;
	      max_cluster_w = it->first;
	    }
	  }
	  bool flag_remove = true;
	  
	  /* std::cout << cluster->get_cluster_id() << " " << max_cluster_u << " " << max_value_u/(num_total_points-num_dead[0]+1e-9) << " " */
	  /* 	  << max_cluster_v << " " << max_value_v/(num_total_points-num_dead[1]+1e-9) << " " */
	  /* 	  << max_cluster_w << " " << max_value_w/(num_total_points-num_dead[2]+1e-9) << std::endl; */
	  
	  if (max_cluster_u==max_cluster_v && max_value_u > 0.8 * (num_total_points-num_dead[0]) && max_value_v > 0.8 * (num_total_points - num_dead[1]) ){
	    if (map_cluster_num[2].find(max_cluster_u)!=map_cluster_num[2].end()){
	      if (map_cluster_num[2][max_cluster_u]> 0.65 * (num_total_points - num_dead[2])){
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_u));
		  flag_remove = false;
		}
	      }
	    }else{
	      if (num_total_points == num_dead[2] && max_cluster_u!=0){
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_u));
		  flag_remove = false;
		}
	      }
	    }
	    
	  }else if (max_cluster_u==max_cluster_w && max_value_u > 0.8 * (num_total_points-num_dead[0]) && max_value_w > 0.8 * (num_total_points - num_dead[2])){
	    if (map_cluster_num[1].find(max_cluster_u)!=map_cluster_num[1].end()){
	      if (map_cluster_num[1][max_cluster_u]>0.65 * (num_total_points - num_dead[1])){
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_u));
		  flag_remove = false;
		}
	      }
	    }else{
	      if (num_total_points == num_dead[1] && max_cluster_u!=0){
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_u));
		  flag_remove = false;
		}
	      }
	    }
	  }else if (max_cluster_w==max_cluster_v && max_value_w > 0.8 * (num_total_points-num_dead[2]) && max_value_v > 0.8 * (num_total_points - num_dead[1]) ){
	    if (map_cluster_num[0].find(max_cluster_w)!=map_cluster_num[0].end()){
	      if (map_cluster_num[0][max_cluster_w]> 0.65 * (num_total_points - num_dead[0])){
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_w->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_w));
		  flag_remove = false;
		}
	      }
	    }else{
	      if (num_total_points == num_dead[0] && max_cluster_w!=0){
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_w->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_w));
		  flag_remove = false;
		}
	      }
	    }
	  }
	  if (flag_remove){
	    to_be_removed_clusters.push_back(cluster);
	  }
	  
	  //to_be_removed_clusters.push_back(cluster);
	  /* if (!flag_remove) */
	  //std::cout << cluster->get_cluster_id() << " " << num_dead[0] << " " << num_dead[1] << " " << num_dead[2] << " " << num_unique[0]/(num_total_points - num_dead[0]+1e-9) << " " << num_unique[1]/(num_total_points - num_dead[1]+1e-9) << " " << num_unique[2]/(num_total_points - num_dead[2]+1e-9) << " " << num_unique[0]+num_unique[1] + num_unique[2] << " " << (num_unique[0]+num_unique[1] + num_unique[2])/(num_total_points - num_dead[0] + num_total_points - num_dead[1] + num_total_points - num_dead[2]+1e-9) << " " << num_total_points << std::endl; 
	}else{
	  /* if (num_dead[0] + num_dead[1] + num_dead[2] >0) */
	  /* std::cout <<  cluster->get_cluster_id() << " " << num_dead[0] << " " << num_dead[1] << " " << num_dead[2] << " " << num_unique[0]/(num_total_points - num_dead[0]+1e-9) << " " << num_unique[1]/(num_total_points - num_dead[1]+1e-9) << " " << num_unique[2]/(num_total_points - num_dead[2]+1e-9) << " " << num_unique[0]+num_unique[1] + num_unique[2] << " " << (num_unique[0]+num_unique[1] + num_unique[2])/(num_total_points - num_dead[0] + num_total_points - num_dead[1] + num_total_points - num_dead[2]+1e-9) << " " << num_total_points << " " << num_unique[0] << " " << num_unique[1] << " " << num_unique[2] << " " << cluster_length_map[cluster]/units::cm << std::endl; */

	  flag_save = true;
	  if ((num_unique[0]+num_unique[1] + num_unique[2])/(num_total_points - num_dead[0] + num_total_points - num_dead[1] + num_total_points - num_dead[2]+1e-9)<0.15 && cluster_length_map[cluster] < 25*units::cm){
	    PR3DCluster *max_cluster_u = 0, *max_cluster_v=0, *max_cluster_w=0;
	    int max_value_u = 0, max_value_v = 0, max_value_w = 0;
	    for (auto it = map_cluster_num[0].begin(); it!=map_cluster_num[0].end(); it++){
	      if (it->second > max_value_u){
		max_value_u = it->second;
		max_cluster_u = it->first;
	      }
	    }
	    for (auto it = map_cluster_num[1].begin(); it!=map_cluster_num[1].end(); it++){
	      if (it->second > max_value_v){
		max_value_v = it->second;
		max_cluster_v = it->first;
	      }
	    }
	    for (auto it = map_cluster_num[2].begin(); it!=map_cluster_num[2].end(); it++){
	      if (it->second > max_value_w){
		max_value_w = it->second;
		max_cluster_w = it->first;
	      }
	    }
	    
	    /* std::cout << max_cluster_u << " " << max_value_u/(num_total_points-num_dead[0]+1e-9) << " " */
	    /* 	  << max_cluster_v << " " << max_value_v/(num_total_points-num_dead[1]+1e-9) << " " */
	    /* 	  << max_cluster_w << " " << max_value_w/(num_total_points-num_dead[2]+1e-9) << std::endl; */
	    
	    
	    if ( (max_cluster_u==max_cluster_v && max_cluster_v == max_cluster_w) ||
		 (max_cluster_u==max_cluster_v && max_cluster_w==0) ||
		 (max_cluster_w==max_cluster_v && max_cluster_u==0) ||
		 (max_cluster_u==max_cluster_w && max_cluster_v==0) ){
	      
	      //  std::cout << cluster->get_cluster_id() << " " << (num_unique[0]+num_unique[1] + num_unique[2])/(num_total_points - num_dead[0] + num_total_points - num_dead[1] + num_total_points - num_dead[2]+1e-9) << " " << (max_value_u+max_value_v+max_value_w)/(num_total_points  + num_total_points  + num_total_points +1e-9) << std::endl;
	      
	      if ((max_value_u+max_value_v+max_value_w)/(num_total_points  + num_total_points  + num_total_points +1e-9)>0.25 ){
		flag_save = false;
		if(max_cluster_u!=0){
		  ToyPointCloud *cloud1 = cluster->get_point_cloud();
		  ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
		  std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		  if (std::get<2>(temp_results)<20*units::cm){
		    to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_u));
		  }else{
		    to_be_removed_clusters.push_back(cluster);
		  }
		}else if (max_cluster_v!=0){
		  ToyPointCloud *cloud1 = cluster->get_point_cloud();
		  ToyPointCloud *cloud2 = max_cluster_v->get_point_cloud();
		  std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		  if (std::get<2>(temp_results)<20*units::cm){
		    to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_v));
		  }else{
		    to_be_removed_clusters.push_back(cluster);
		  }
		}else if (max_cluster_w!=0){
		  ToyPointCloud *cloud1 = cluster->get_point_cloud();
		  ToyPointCloud *cloud2 = max_cluster_w->get_point_cloud();
		  std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		  if (std::get<2>(temp_results)<20*units::cm){
		    to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_w));
		  }else{
		    to_be_removed_clusters.push_back(cluster);
		  }
		}
	      }
	    }else if (max_cluster_u==max_cluster_v && max_cluster_u!=0){
	      if ((max_value_u+max_value_v+map_cluster_num[2][max_cluster_u])/(num_total_points  + num_total_points  + num_total_points +1e-9)>0.25){
		flag_save = false;
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_u));
		}else{
		  to_be_removed_clusters.push_back(cluster);
		}
	      }
	    }else if (max_cluster_v==max_cluster_w && max_cluster_v!=0){
	      if ((map_cluster_num[0][max_cluster_v]+max_value_v+max_value_w)/(num_total_points  + num_total_points  + num_total_points +1e-9)>0.25){
		flag_save = false;
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_v->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_v));
		}else{
		  to_be_removed_clusters.push_back(cluster);
		}
	      }
	    }else if (max_cluster_u==max_cluster_w && max_cluster_w!=0){
	      if ((max_value_u+map_cluster_num[1][max_cluster_w]+max_value_w)/(num_total_points  + num_total_points  + num_total_points +1e-9)>0.25){
		flag_save = false;
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_w->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_w));
		}else{
		  to_be_removed_clusters.push_back(cluster);
		}
	      }
	    }
	  }
	  // two cases, merge clusters or remove clusters
	}
	
	
	if (flag_save){
	  live_clusters.at(i)->Create_point_cloud();
	  global_point_cloud.AddPoints(live_clusters.at(i),0);
	  if (cluster_length_map[live_clusters.at(i)]>30*units::cm){
	    live_clusters.at(i)->Construct_skeleton();
	    global_skeleton_cloud.AddPoints(live_clusters.at(i),1);
	  }
	}
      }else{
	live_clusters.at(i)->Create_point_cloud();
	global_point_cloud.AddPoints(live_clusters.at(i),0);
	if (cluster_length_map[live_clusters.at(i)]>30*units::cm){
	  live_clusters.at(i)->Construct_skeleton();
	  global_skeleton_cloud.AddPoints(live_clusters.at(i),1);
	}
      }
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


void WCP2dToy::Clustering_deghost(WCP::ToyCTPointCloud& ct_point_cloud, WCP::PR3DClusterSelection& live_clusters, std::map<WCP::PR3DCluster*,double>& cluster_length_map, std::map<int,std::pair<double,double>>& dead_u_index, std::map<int,std::pair<double,double>>& dead_v_index, std::map<int,std::pair<double,double>>& dead_w_index, double length_cut){

  
  
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
      if (cluster_length_map[live_clusters.at(i)]>30*units::cm){ // should be the default for most of them ...
	live_clusters.at(i)->Construct_skeleton(ct_point_cloud);
	global_skeleton_cloud.AddPoints(live_clusters.at(i),1);
      }else{
	global_skeleton_cloud.AddPoints(live_clusters.at(i),0);
      }
    }else{ 
      // start the process to add things in and perform deghosting ... 
      PR3DCluster* cluster = live_clusters.at(i);
      
      if (length_cut == 0 || cluster_length_map[live_clusters.at(i)] < length_cut)  {
	cluster->Create_point_cloud();
	WCP::WCPointCloud<double>& cloud = cluster->get_point_cloud()->get_cloud();
	int num_total_points = cloud.pts.size(); // total number of points
	int num_dead[3]={0,0,0}; // dead wires in each view 
	int num_unique[3]={0,0,0}; // points that are unique (not agree with any other clusters)
	std::map<PR3DCluster*, int> map_cluster_num[3];
	
	double dis_cut = 1.2*units::cm;
	

	
	for (size_t j=0;j!=num_total_points;j++){
	  Point test_point(cloud.pts.at(j).x,cloud.pts.at(j).y,cloud.pts.at(j).z);
	  
	  bool flag_dead = false;
	  
	  if (dead_u_index.find(cloud.pts.at(j).index_u)!=dead_u_index.end()){
	    if (cloud.pts.at(j).x >= dead_u_index[cloud.pts.at(j).index_u].first &&
		cloud.pts.at(j).x <= dead_u_index[cloud.pts.at(j).index_u].second){
	      flag_dead = true;
	    }
	  }
	 
	  if (!flag_dead){
	    std::tuple<double, PR3DCluster*, size_t> results = global_point_cloud.get_closest_2d_point_info(test_point, 0);
	    
	    if (std::get<0>(results)<=dis_cut/3.){
	      if (map_cluster_num[0].find(std::get<1>(results))==map_cluster_num[0].end()){
		map_cluster_num[0][std::get<1>(results)] = 1;
	      }else{
		map_cluster_num[0][std::get<1>(results)] ++;
	      }  
	    }else{
	      results = global_skeleton_cloud.get_closest_2d_point_info(test_point, 0);
	      if (std::get<0>(results)<=dis_cut*2.0){
		if (map_cluster_num[0].find(std::get<1>(results))==map_cluster_num[0].end()){
		  map_cluster_num[0][std::get<1>(results)] = 1;
		}else{
		  map_cluster_num[0][std::get<1>(results)] ++;
		}
	      }else{
		num_unique[0]++;
	      }
	    }
	  }else{
	    num_dead[0]++;
	  }

	  
	  
	  flag_dead = false;
	  
	  if (dead_v_index.find(cloud.pts.at(j).index_v)!=dead_v_index.end()){
	    if (cloud.pts.at(j).x >= dead_v_index[cloud.pts.at(j).index_v].first &&
		cloud.pts.at(j).x <= dead_v_index[cloud.pts.at(j).index_v].second){
	      flag_dead = true;
	    }
	  }
	  
	  
	  if (!flag_dead){
	    std::tuple<double, PR3DCluster*, size_t> results = global_point_cloud.get_closest_2d_point_info(test_point, 1);
	    if (std::get<0>(results)<=dis_cut/3.){
	      if (map_cluster_num[1].find(std::get<1>(results))==map_cluster_num[1].end()){
		map_cluster_num[1][std::get<1>(results)] = 1;
	      }else{
		map_cluster_num[1][std::get<1>(results)] ++;
	      }
	    }else{
	      results = global_skeleton_cloud.get_closest_2d_point_info(test_point, 1);
	      if (std::get<0>(results)<=dis_cut*2.0){
		if (map_cluster_num[1].find(std::get<1>(results))==map_cluster_num[1].end()){
		  map_cluster_num[1][std::get<1>(results)] = 1;
		}else{
		  map_cluster_num[1][std::get<1>(results)] ++;
		}
	      }else{
		num_unique[1]++;
	      }
	    }
	  }else{
	    num_dead[1]++;
	  }
	  
	  flag_dead = false;
	  
	  if (dead_w_index.find(cloud.pts.at(j).index_w)!=dead_w_index.end()){
	    if (cloud.pts.at(j).x >= dead_w_index[cloud.pts.at(j).index_w].first &&
		cloud.pts.at(j).x <= dead_w_index[cloud.pts.at(j).index_w].second){
	      flag_dead = true;
	    }
	  }
	  
	  
	  if (!flag_dead){
	    std::tuple<double, PR3DCluster*, size_t> results = global_point_cloud.get_closest_2d_point_info(test_point, 2);
	    if (std::get<0>(results)<=dis_cut/3.){
	      if (map_cluster_num[2].find(std::get<1>(results))==map_cluster_num[2].end()){
		map_cluster_num[2][std::get<1>(results)] = 1;
	      }else{
		map_cluster_num[2][std::get<1>(results)] ++;
	      }
	    }else{
	      results = global_skeleton_cloud.get_closest_2d_point_info(test_point, 2);
	      if (std::get<0>(results)<=dis_cut*2.0){
		if (map_cluster_num[2].find(std::get<1>(results))==map_cluster_num[2].end()){
		  map_cluster_num[2][std::get<1>(results)] = 1;
		}else{
		  map_cluster_num[2][std::get<1>(results)] ++;
		}
	      }else{
		num_unique[2]++;
	      }
	    }
	  }else{
	    num_dead[2]++;
	  }
	  
	  
	  if ( (num_unique[1]+num_unique[0]+num_unique[2]) >= 0.24 * num_total_points &&
	       (num_unique[1]+num_unique[0]+num_unique[2]) > 25 )
	    break;
	}
	
	//std::cout << cluster->get_cluster_id() << " " << num_unique[0] << " " << num_unique[1] << " " << num_unique[2] << " " << num_dead[0] << " " << num_dead[1] << " " << num_dead[2] << " " << num_total_points << std::endl;
	
	bool flag_save = false;
	
	if (((num_unique[0] <= 0.1 * (num_total_points - num_dead[0]) || num_unique[0] <= 0.1 * num_total_points && num_unique[0] <= 8) &&
	     (num_unique[1] <= 0.1 * (num_total_points - num_dead[1]) || num_unique[1] <= 0.1 * num_total_points && num_unique[1] <= 8) &&
	     (num_unique[2] <= 0.1 * (num_total_points - num_dead[2]) || num_unique[2] <= 0.1 * num_total_points && num_unique[2] <= 8) &&
	     ((num_unique[0] + num_unique[1] + num_unique[2]) <= 0.05 * (num_total_points - num_dead[0] + num_total_points - num_dead[1] + num_total_points - num_dead[2]) || (num_unique[0] + num_unique[1] + num_unique[2]) < 0.15 * num_total_points && (num_unique[0] + num_unique[1] + num_unique[2]) <=8) ||
	     num_unique[0]==0 && num_unique[1]==0 && num_unique[2] < 0.24 * num_total_points ||
	     num_unique[0]==0 && num_unique[2]==0 && num_unique[1] < 0.24 * num_total_points ||
	     num_unique[2]==0 && num_unique[1]==0 && num_unique[0] < 0.24 * num_total_points ||
	     num_unique[0]==0 && (num_unique[1]+num_unique[2]) < 0.12 * num_total_points * 2 || 
	     num_unique[1]==0 && (num_unique[0]+num_unique[2]) < 0.12 * num_total_points * 2 || 
	     num_unique[2]==0 && (num_unique[1]+num_unique[0]) < 0.12 * num_total_points * 2 ||
	     (num_unique[1]+num_unique[0]+num_unique[2]) < 0.24 * num_total_points &&
	     (num_unique[1]<0.02 * num_total_points || num_unique[0]<0.02 * num_total_points
	      || num_unique[2]<0.02 * num_total_points)
	     ) &&
	    (num_unique[0] + num_unique[1] + num_unique[2])  <= 500 ){
	  flag_save = false;
	  
	  // now try to compare
	  // find the maximal for each map
	  PR3DCluster *max_cluster_u = 0, *max_cluster_v=0, *max_cluster_w=0;
	  int max_value_u = 0, max_value_v = 0, max_value_w = 0;
	  for (auto it = map_cluster_num[0].begin(); it!=map_cluster_num[0].end(); it++){
	    if (it->second > max_value_u){
	      max_value_u = it->second;
	      max_cluster_u = it->first;
	    }
	  }
	  for (auto it = map_cluster_num[1].begin(); it!=map_cluster_num[1].end(); it++){
	    if (it->second > max_value_v){
	      max_value_v = it->second;
	      max_cluster_v = it->first;
	    }
	  }
	  for (auto it = map_cluster_num[2].begin(); it!=map_cluster_num[2].end(); it++){
	    if (it->second > max_value_w){
	      max_value_w = it->second;
	      max_cluster_w = it->first;
	    }
	  }
	  bool flag_remove = true;
	  
	  /* std::cout << cluster->get_cluster_id() << " " << max_cluster_u << " " << max_value_u/(num_total_points-num_dead[0]+1e-9) << " " */
	  /* 	  << max_cluster_v << " " << max_value_v/(num_total_points-num_dead[1]+1e-9) << " " */
	  /* 	  << max_cluster_w << " " << max_value_w/(num_total_points-num_dead[2]+1e-9) << std::endl; */
	  
	  if (max_cluster_u==max_cluster_v && max_value_u > 0.8 * (num_total_points-num_dead[0]) && max_value_v > 0.8 * (num_total_points - num_dead[1]) ){
	    if (map_cluster_num[2].find(max_cluster_u)!=map_cluster_num[2].end()){
	      if (map_cluster_num[2][max_cluster_u]> 0.65 * (num_total_points - num_dead[2])){
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_u));
		  flag_remove = false;
		}
	      }
	    }else{
	      if (num_total_points == num_dead[2] && max_cluster_u!=0){
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_u));
		  flag_remove = false;
		}
	      }
	    }
	    
	  }else if (max_cluster_u==max_cluster_w && max_value_u > 0.8 * (num_total_points-num_dead[0]) && max_value_w > 0.8 * (num_total_points - num_dead[2])){
	    if (map_cluster_num[1].find(max_cluster_u)!=map_cluster_num[1].end()){
	      if (map_cluster_num[1][max_cluster_u]>0.65 * (num_total_points - num_dead[1])){
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_u));
		  flag_remove = false;
		}
	      }
	    }else{
	      if (num_total_points == num_dead[1] && max_cluster_u!=0){
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_u));
		  flag_remove = false;
		}
	      }
	    }
	  }else if (max_cluster_w==max_cluster_v && max_value_w > 0.8 * (num_total_points-num_dead[2]) && max_value_v > 0.8 * (num_total_points - num_dead[1]) ){
	    if (map_cluster_num[0].find(max_cluster_w)!=map_cluster_num[0].end()){
	      if (map_cluster_num[0][max_cluster_w]> 0.65 * (num_total_points - num_dead[0])){
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_w->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_w));
		  flag_remove = false;
		}
	      }
	    }else{
	      if (num_total_points == num_dead[0] && max_cluster_w!=0){
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_w->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_w));
		  flag_remove = false;
		}
	      }
	    }
	  }
	  if (flag_remove){
	    to_be_removed_clusters.push_back(cluster);
	  }
	  
	  //to_be_removed_clusters.push_back(cluster);
	  /* if (!flag_remove) */
	  //std::cout << cluster->get_cluster_id() << " " << num_dead[0] << " " << num_dead[1] << " " << num_dead[2] << " " << num_unique[0]/(num_total_points - num_dead[0]+1e-9) << " " << num_unique[1]/(num_total_points - num_dead[1]+1e-9) << " " << num_unique[2]/(num_total_points - num_dead[2]+1e-9) << " " << num_unique[0]+num_unique[1] + num_unique[2] << " " << (num_unique[0]+num_unique[1] + num_unique[2])/(num_total_points - num_dead[0] + num_total_points - num_dead[1] + num_total_points - num_dead[2]+1e-9) << " " << num_total_points << std::endl; 
	}else{
	  /* if (num_dead[0] + num_dead[1] + num_dead[2] >0) */
	  /* std::cout <<  cluster->get_cluster_id() << " " << num_dead[0] << " " << num_dead[1] << " " << num_dead[2] << " " << num_unique[0]/(num_total_points - num_dead[0]+1e-9) << " " << num_unique[1]/(num_total_points - num_dead[1]+1e-9) << " " << num_unique[2]/(num_total_points - num_dead[2]+1e-9) << " " << num_unique[0]+num_unique[1] + num_unique[2] << " " << (num_unique[0]+num_unique[1] + num_unique[2])/(num_total_points - num_dead[0] + num_total_points - num_dead[1] + num_total_points - num_dead[2]+1e-9) << " " << num_total_points << " " << num_unique[0] << " " << num_unique[1] << " " << num_unique[2] << " " << cluster_length_map[cluster]/units::cm << std::endl; */

	  flag_save = true;
	  if ((num_unique[0]+num_unique[1] + num_unique[2])/(num_total_points - num_dead[0] + num_total_points - num_dead[1] + num_total_points - num_dead[2]+1e-9)<0.15 && cluster_length_map[cluster] < 25*units::cm){
	    PR3DCluster *max_cluster_u = 0, *max_cluster_v=0, *max_cluster_w=0;
	    int max_value_u = 0, max_value_v = 0, max_value_w = 0;
	    for (auto it = map_cluster_num[0].begin(); it!=map_cluster_num[0].end(); it++){
	      if (it->second > max_value_u){
		max_value_u = it->second;
		max_cluster_u = it->first;
	      }
	    }
	    for (auto it = map_cluster_num[1].begin(); it!=map_cluster_num[1].end(); it++){
	      if (it->second > max_value_v){
		max_value_v = it->second;
		max_cluster_v = it->first;
	      }
	    }
	    for (auto it = map_cluster_num[2].begin(); it!=map_cluster_num[2].end(); it++){
	      if (it->second > max_value_w){
		max_value_w = it->second;
		max_cluster_w = it->first;
	      }
	    }
	    
	    /* std::cout << max_cluster_u << " " << max_value_u/(num_total_points-num_dead[0]+1e-9) << " " */
	    /* 	  << max_cluster_v << " " << max_value_v/(num_total_points-num_dead[1]+1e-9) << " " */
	    /* 	  << max_cluster_w << " " << max_value_w/(num_total_points-num_dead[2]+1e-9) << std::endl; */
	    
	    
	    if ( (max_cluster_u==max_cluster_v && max_cluster_v == max_cluster_w) ||
		 (max_cluster_u==max_cluster_v && max_cluster_w==0) ||
		 (max_cluster_w==max_cluster_v && max_cluster_u==0) ||
		 (max_cluster_u==max_cluster_w && max_cluster_v==0) ){
	      
	      //  std::cout << cluster->get_cluster_id() << " " << (num_unique[0]+num_unique[1] + num_unique[2])/(num_total_points - num_dead[0] + num_total_points - num_dead[1] + num_total_points - num_dead[2]+1e-9) << " " << (max_value_u+max_value_v+max_value_w)/(num_total_points  + num_total_points  + num_total_points +1e-9) << std::endl;
	      
	      if ((max_value_u+max_value_v+max_value_w)/(num_total_points  + num_total_points  + num_total_points +1e-9)>0.25 ){
		flag_save = false;
		if(max_cluster_u!=0){
		  ToyPointCloud *cloud1 = cluster->get_point_cloud();
		  ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
		  std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		  if (std::get<2>(temp_results)<20*units::cm){
		    to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_u));
		  }else{
		    to_be_removed_clusters.push_back(cluster);
		  }
		}else if (max_cluster_v!=0){
		  ToyPointCloud *cloud1 = cluster->get_point_cloud();
		  ToyPointCloud *cloud2 = max_cluster_v->get_point_cloud();
		  std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		  if (std::get<2>(temp_results)<20*units::cm){
		    to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_v));
		  }else{
		    to_be_removed_clusters.push_back(cluster);
		  }
		}else if (max_cluster_w!=0){
		  ToyPointCloud *cloud1 = cluster->get_point_cloud();
		  ToyPointCloud *cloud2 = max_cluster_w->get_point_cloud();
		  std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		  if (std::get<2>(temp_results)<20*units::cm){
		    to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_w));
		  }else{
		    to_be_removed_clusters.push_back(cluster);
		  }
		}
	      }
	    }else if (max_cluster_u==max_cluster_v && max_cluster_u!=0){
	      if ((max_value_u+max_value_v+map_cluster_num[2][max_cluster_u])/(num_total_points  + num_total_points  + num_total_points +1e-9)>0.25){
		flag_save = false;
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_u));
		}else{
		  to_be_removed_clusters.push_back(cluster);
		}
	      }
	    }else if (max_cluster_v==max_cluster_w && max_cluster_v!=0){
	      if ((map_cluster_num[0][max_cluster_v]+max_value_v+max_value_w)/(num_total_points  + num_total_points  + num_total_points +1e-9)>0.25){
		flag_save = false;
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_v->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_v));
		}else{
		  to_be_removed_clusters.push_back(cluster);
		}
	      }
	    }else if (max_cluster_u==max_cluster_w && max_cluster_w!=0){
	      if ((max_value_u+map_cluster_num[1][max_cluster_w]+max_value_w)/(num_total_points  + num_total_points  + num_total_points +1e-9)>0.25){
		flag_save = false;
		ToyPointCloud *cloud1 = cluster->get_point_cloud();
		ToyPointCloud *cloud2 = max_cluster_w->get_point_cloud();
		std::tuple<int,int,double> temp_results = cloud1->get_closest_points(cloud2);
		if (std::get<2>(temp_results)<20*units::cm){
		  to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_w));
		}else{
		  to_be_removed_clusters.push_back(cluster);
		}
	      }
	    }
	  }
	  // two cases, merge clusters or remove clusters
	}
	
	
	if (flag_save){
	  live_clusters.at(i)->Create_point_cloud();
	  global_point_cloud.AddPoints(live_clusters.at(i),0);
	  if (cluster_length_map[live_clusters.at(i)]>30*units::cm){
	    live_clusters.at(i)->Construct_skeleton(ct_point_cloud);
	    global_skeleton_cloud.AddPoints(live_clusters.at(i),1);
	  }
	}
      }else{
	live_clusters.at(i)->Create_point_cloud();
	global_point_cloud.AddPoints(live_clusters.at(i),0);
	if (cluster_length_map[live_clusters.at(i)]>30*units::cm){
	  live_clusters.at(i)->Construct_skeleton(ct_point_cloud);
	  global_skeleton_cloud.AddPoints(live_clusters.at(i),1);
	}
      }
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
