
using namespace WireCell;

void WireCell2dToy::Clustering_connect1(WireCell::PR3DClusterSelection& live_clusters, std::map<WireCell::PR3DCluster*,double>& cluster_length_map, WireCell::DynamicToyPointCloud& global_point_cloud, std::map<int,std::pair<double,double>>& dead_u_index, std::map<int,std::pair<double,double>>& dead_v_index, std::map<int,std::pair<double,double>>& dead_w_index){
  
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

  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();

  DynamicToyPointCloud global_skeleton_cloud(angle_u,angle_v,angle_w);

  double extending_dis = 50*units::cm;
  double angle = 7.5;
  double loose_dis_cut = 7.5*units::cm;
  
  std::set<std::pair<PR3DCluster*, PR3DCluster*>> to_be_merged_pairs;

  WireCell::WCPointCloud<double>& global_cloud = global_skeleton_cloud.get_cloud();

  
  for (size_t i=0;i!=live_clusters.size();i++){
    PR3DCluster* cluster = live_clusters.at(i);
    cluster->Create_point_cloud();

  
    std::pair<Point,Point> extreme_points = cluster->get_two_extreme_points();
    TVector3 main_dir(extreme_points.second.x - extreme_points.first.x,
		      extreme_points.second.y - extreme_points.first.y,
		      extreme_points.second.z - extreme_points.first.z);
    TVector3 dir1, dir2;
    
    if (cluster_length_map[cluster] > 25*units::cm){
      dir1 = cluster->VHoughTrans(extreme_points.first,80*units::cm);
      dir2 = cluster->VHoughTrans(extreme_points.second,80*units::cm);
      if (dir1.Dot(main_dir)>0) dir1 *=-1;
      if (dir2.Dot(dir1)>0) dir2 *= -1;
    }else{
      dir1 = global_point_cloud.VHoughTrans(extreme_points.first,extending_dis);
      dir2 = global_point_cloud.VHoughTrans(extreme_points.second,extending_dis);
      if (dir1.Dot(main_dir)>0) dir1 *=-1;
      if (dir2.Dot(dir1)>0) dir2 *= -1;
    }
    
    if (fabs(extreme_points.first.x-97*units::cm)< 5*units::cm &&
	fabs(extreme_points.first.z-979*units::cm)< 5*units::cm &&
	fabs(extreme_points.first.y-17.5*units::cm)< 5*units::cm ){
      std::cout << cluster->get_cluster_id() << " " << extreme_points.first.x/units::cm << " " << extreme_points.first.y/units::cm << " " << extreme_points.first.z/units::cm << " " << extreme_points.second.x/units::cm << " " << extreme_points.second.y/units::cm << " " << extreme_points.second.z/units::cm << " "  << dir1.X() << " " << dir1.Y() << " " << dir1.Z() << " " << dir2.X() << " " << dir2.Y() << " " << dir2.Z() << std::endl;
    }
    
    
    if (i==0){
     
      
      global_skeleton_cloud.AddPoints(cluster,extreme_points.first,dir1,extending_dis,0.6*units::cm, angle);
      global_skeleton_cloud.AddPoints(cluster,extreme_points.second,dir2,extending_dis,0.6*units::cm, angle);
    }else{
      
      WireCell::WCPointCloud<double>& cloud = cluster->get_point_cloud()->get_cloud();
      int num_total_points = cloud.pts.size(); // total number of points
      int num_dead[3]={0,0,0}; // dead wires in each view
      int num_unique[3]={0,0,0}; // points that are unique (not agree with any other clusters)
      std::map<PR3DCluster*, int> map_cluster_num[3];
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
	  std::tuple<double, PR3DCluster*, size_t> results = global_skeleton_cloud.get_closest_2d_point_info(test_point, 0);
	  if (std::get<0>(results) <loose_dis_cut){
	    if (std::get<0>(results) < global_cloud.pts.at(std::get<2>(results)).index_u){
	      if (map_cluster_num[0].find(std::get<1>(results))==map_cluster_num[0].end()){
		map_cluster_num[0][std::get<1>(results)] = 1;
	      }else{
		map_cluster_num[0][std::get<1>(results)] ++;
	      }  
	    }else{
	      num_unique[0]++;
	    }
	  }else{
	    num_unique[0]++;
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
	  std::tuple<double, PR3DCluster*, size_t> results = global_skeleton_cloud.get_closest_2d_point_info(test_point, 1);
	  if (std::get<0>(results) <loose_dis_cut){
	    if (std::get<0>(results) < global_cloud.pts.at(std::get<2>(results)).index_v){
	      if (map_cluster_num[1].find(std::get<1>(results))==map_cluster_num[0].end()){
		map_cluster_num[1][std::get<1>(results)] = 1;
	      }else{
		map_cluster_num[1][std::get<1>(results)] ++;
	      }  
	    }else{
	      num_unique[1]++;
	    }
	  }else{
	    num_unique[1]++;
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
	  std::tuple<double, PR3DCluster*, size_t> results = global_skeleton_cloud.get_closest_2d_point_info(test_point, 2);
	  if (std::get<0>(results) <loose_dis_cut){
	    if (std::get<0>(results) < global_cloud.pts.at(std::get<2>(results)).index_w){
	      if (map_cluster_num[2].find(std::get<1>(results))==map_cluster_num[2].end()){
		map_cluster_num[2][std::get<1>(results)] = 1;
	      }else{
		map_cluster_num[2][std::get<1>(results)] ++;
	      }  
	    }else{
	      num_unique[2]++;
	    }
	  }else{
	    num_unique[2]++;
	  }
	}else{
	  num_dead[2]++;
	}
	
	
      }
      PR3DCluster *curr_cluster = cluster;


      if (fabs(extreme_points.first.x-97*units::cm)< 5*units::cm &&
	fabs(extreme_points.first.z-979*units::cm)< 5*units::cm &&
	fabs(extreme_points.first.y-17.5*units::cm)< 5*units::cm ){
	std::cout << cluster->get_cluster_id() << " " << num_dead[0] << " " << num_dead[1] << " " << num_dead[2] << " "
		  << num_unique[0] << " " << num_unique[1] << " " << num_unique[2] << " " << num_total_points << std::endl;
      }

      
      if ((num_unique[1]+num_unique[0]+num_unique[2]) < 0.24 * num_total_points ){
	

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

	if ( (max_cluster_u==max_cluster_v && max_cluster_v == max_cluster_w) ||
	     (max_cluster_u==max_cluster_v && max_cluster_w==0) ||
	     (max_cluster_w==max_cluster_v && max_cluster_u==0) ||
	     (max_cluster_u==max_cluster_w && max_cluster_v==0) ){
	  std::cout << cluster->get_cluster_id() << " " << (num_unique[0]+num_unique[1] + num_unique[2])/(num_total_points - num_dead[0] + num_total_points - num_dead[1] + num_total_points - num_dead[2]+1e-9) << " " << (max_value_u+max_value_v+max_value_w)/(num_total_points  + num_total_points  + num_total_points +1e-9) << " " << num_total_points <<std::endl;
	  if (max_cluster_u!=0) std::cout << 0 << " " << max_cluster_u->get_cluster_id() << std::endl;
	  if (max_cluster_v!=0) std::cout << 0 << " " << max_cluster_v->get_cluster_id() << std::endl;
	  if (max_cluster_w!=0) std::cout << 0 << " " << max_cluster_w->get_cluster_id() << std::endl;
	  
	  if ((max_value_u+max_value_v+max_value_w)/(num_total_points  + num_total_points  + num_total_points +1e-9)>0.25 ){
	    if (max_cluster_u!=0){
	      to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_u));
	      curr_cluster = max_cluster_u;
	    }else if (max_cluster_v!=0){
	      to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_v));
	      curr_cluster = max_cluster_v;
	    }else if (max_cluster_w!=0){
	      to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_w));
	      curr_cluster = max_cluster_w;
	    }
	  }else if (max_cluster_u==max_cluster_v && max_cluster_u!=0){
	    if ((max_value_u+max_value_v+map_cluster_num[2][max_cluster_u])/(num_total_points  + num_total_points  + num_total_points +1e-9)>0.25){
	       to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_u));
	       curr_cluster = max_cluster_u;
	    }
	  }else if (max_cluster_v==max_cluster_w && max_cluster_v!=0){
	    if ((map_cluster_num[0][max_cluster_v]+max_value_v+max_value_w)/(num_total_points  + num_total_points  + num_total_points +1e-9)>0.25){
	      to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_v));
	      curr_cluster = max_cluster_v;
	    }
	  }else if (max_cluster_u==max_cluster_w && max_cluster_w!=0){
	     if ((max_value_u+map_cluster_num[1][max_cluster_w]+max_value_w)/(num_total_points  + num_total_points  + num_total_points +1e-9)>0.25){
	       to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster_w));
	       curr_cluster = max_cluster_w;
	     }
	  }
	  
	}
	
      }
      
      // add extension points in ... 
      global_skeleton_cloud.AddPoints(curr_cluster,extreme_points.first,dir1,extending_dis,0.6*units::cm, angle);
      global_skeleton_cloud.AddPoints(curr_cluster,extreme_points.second,dir2,extending_dis,0.6*units::cm, angle);
      
    }

    
   


    
    //std::cout << extreme_points.first.x/units::cm << " " << extreme_points.first.y/units::cm << " " << extreme_points.first.z/units::cm << " " << extreme_points.second.x/units::cm << " " << extreme_points.second.y/units::cm << " " << extreme_points.second.z/units::cm << std::endl;
  }


  //to_be_merged_pairs.clear();


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
  
}
