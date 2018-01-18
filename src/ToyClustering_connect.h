
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
  
  std::set<std::pair<PR3DCluster*, PR3DCluster*>> to_be_merged_pairs;
  for (size_t i=0;i!=live_clusters.size();i++){
    live_clusters.at(i)->Create_point_cloud();
    
    std::pair<Point,Point> extreme_points = live_clusters.at(i)->get_two_extreme_points();
    TVector3 main_dir(extreme_points.second.x - extreme_points.first.x,
		      extreme_points.second.y - extreme_points.first.y,
		      extreme_points.second.z - extreme_points.first.z);
    TVector3 dir1 = global_point_cloud.VHoughTrans(extreme_points.first,80*units::cm);
    TVector3 dir2 = global_point_cloud.VHoughTrans(extreme_points.second,80*units::cm);
    if (dir1.Dot(main_dir)>0) dir1 *=-1;
    if (dir2.Dot(dir1)>0) dir2 *= -1;

    // add extension points in ... 
    
    
    //std::cout << extreme_points.first.x/units::cm << " " << extreme_points.first.y/units::cm << " " << extreme_points.first.z/units::cm << " " << extreme_points.second.x/units::cm << " " << extreme_points.second.y/units::cm << " " << extreme_points.second.z/units::cm << std::endl;
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
  
}
