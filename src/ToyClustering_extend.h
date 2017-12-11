void WireCell2dToy::Clustering_extend(WireCell::PR3DClusterSelection& live_clusters, std::map<PR3DCluster*,double>& cluster_length_map, double length_cut){
   // calculate the length ...
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  //std::vector<double> cluster_length_vec;
  //estimated length
  for (size_t i=0; i!= live_clusters.size(); i++){
    PR3DCluster* cluster_1 = live_clusters.at(i);
    std::vector<int> range_v1 = cluster_1->get_uvwt_range();
    double length_1 = sqrt(2./3. * (pow(pitch_u*range_v1.at(0),2) + pow(pitch_v*range_v1.at(1),2) + pow(pitch_w*range_v1.at(2),2)) + pow(time_slice_width*range_v1.at(3),2));
    //cluster_length_vec.push_back(length_1);
    cluster_length_map[cluster_1] = length_1;
  }
  

  for (int kk=0;kk!=1;kk++){
    std::set<std::pair<PR3DCluster*, PR3DCluster*>> to_be_merged_pairs;
    
    for (size_t i=0;i!=live_clusters.size();i++){
      PR3DCluster* cluster_1 = live_clusters.at(i);
      for (size_t j=i+1;j<live_clusters.size();j++){
	PR3DCluster* cluster_2 = live_clusters.at(j);
	//std::cout << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << std::endl;
	if (WireCell2dToy::Clustering_4th_round(cluster_1,cluster_2, cluster_length_map[cluster_1], cluster_length_map[cluster_2], length_cut))
	  to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	
      }
    }
    
    //    to_be_merged_pairs.clear();
    //std::cout << to_be_merged_pairs.size() << std::endl;
    
    
    std::vector<std::set<PR3DCluster*>> merge_clusters;
    for (auto it = to_be_merged_pairs.begin(); it!=to_be_merged_pairs.end(); it++){
      PR3DCluster *cluster1 = (*it).first;
      PR3DCluster *cluster2 = (*it).second;
      //std::cout << cluster1 << " " << cluster2 << " " << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << std::endl;
      
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
	    //std::cout << i << " " << (*it1) << " " << (*it1)->get_cluster_id()  << std::endl;
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
}

bool WireCell2dToy::Clustering_4th_round(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double length_1, double length_2, double length_cut){
  return false;
}
