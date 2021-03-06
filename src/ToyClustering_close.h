void WCP2dToy::Clustering_close(WCP::PR3DClusterSelection& live_clusters, WCP::map_pr3dcluster_double& cluster_length_map, WCP::PR3DClusterSet& cluster_connected_dead, double length_cut){
  // calculate the length ...
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  
  if (cluster_length_map.size()==0)
    for (size_t i=0; i!= live_clusters.size(); i++){
      PR3DCluster* cluster_1 = live_clusters.at(i);
      std::vector<int> range_v1 = cluster_1->get_uvwt_range();
      double length_1 = sqrt(2./3. * (pow(pitch_u*range_v1.at(0),2) + pow(pitch_v*range_v1.at(1),2) + pow(pitch_w*range_v1.at(2),2)) + pow(time_slice_width*range_v1.at(3),2));
      //cluster_length_vec.push_back(length_1);
      cluster_length_map[cluster_1] = length_1;
    }
  
  
  std::set<std::pair<PR3DCluster*, PR3DCluster*>> to_be_merged_pairs;
  std::set<PR3DCluster*> used_clusters;
  
  for (size_t i=0;i!=live_clusters.size();i++){
    PR3DCluster* cluster_1 = live_clusters.at(i);
    if (cluster_length_map[cluster_1] < 1.5*units::cm) continue;
    if (used_clusters.find(cluster_1)!=used_clusters.end()) continue;
    for (size_t j=i+1;j<live_clusters.size();j++){
      PR3DCluster* cluster_2 = live_clusters.at(j);
      if (used_clusters.find(cluster_2)!=used_clusters.end()) continue;
      if (cluster_length_map[cluster_2] < 1.5*units::cm) continue;
      if (WCP2dToy::Clustering_3rd_round(cluster_1,cluster_2, cluster_length_map[cluster_1], cluster_length_map[cluster_2], length_cut)){
	to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	if (cluster_length_map[cluster_1] < 5*units::cm){
	  used_clusters.insert(cluster_1);
	  break;
	}
	if (cluster_length_map[cluster_2] < 5*units::cm){
	  used_clusters.insert(cluster_2);
	}
      }
    }
  }
  
  

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
    bool flag_connected_to_dead = false;
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
      if (cluster_connected_dead.find(ocluster)!=cluster_connected_dead.end()){
	  
	flag_connected_to_dead = true;
	cluster_connected_dead.erase(ocluster);
      }
      delete ocluster;
    }
    if (flag_connected_to_dead)
      cluster_connected_dead.insert(ncluster);
    std::vector<int> range_v1 = ncluster->get_uvwt_range();
    double length_1 = sqrt(2./3. * (pow(pitch_u*range_v1.at(0),2) + pow(pitch_v*range_v1.at(1),2) + pow(pitch_w*range_v1.at(2),2)) + pow(time_slice_width*range_v1.at(3),2));
    cluster_length_map[ncluster] = length_1;
    
    //std::cout << std::endl;
  }
}


bool WCP2dToy::Clustering_3rd_round(WCP::PR3DCluster *cluster1, WCP::PR3DCluster *cluster2, double length_1, double length_2, double length_cut){
  cluster1->Create_point_cloud();
  cluster2->Create_point_cloud();

   // pick any point and merged cell in cluster1,
  SlimMergeGeomCell *prev_mcell1 = 0;
  SlimMergeGeomCell *prev_mcell2 = 0;
  SlimMergeGeomCell *mcell1 = 0; 
  Point p1;//
  SlimMergeGeomCell *mcell2=0;
  Point p2;
  double dis = Find_Closeset_Points(cluster1, cluster2, length_1, length_2, length_cut, mcell1, mcell2, p1,p2);

  
  TVector3 dir1, dir2;
  int num_p1, num_p2, num_tp1, num_tp2;

  // if very close merge anyway???
  if (dis < 0.5*units::cm){
    return true;
  }

  if (dis < 1.0*units::cm && length_2 < 12*units::cm && length_1 <12*units::cm)
    return true;

  
  if (dis < 2.0*units::cm && (length_2 >=12*units::cm || length_1 >=12*units::cm)){
    dir1 = cluster1->VHoughTrans(p1,50*units::cm); // cluster 1 direction based on hough
    dir2 = cluster2->VHoughTrans(p2,50*units::cm); // cluster 1 direction based on hough
    
    std::pair<int,int> num_ps_1 = cluster1->get_num_points(p1,dir1);
    std::pair<int,int> num_ps_2 = cluster2->get_num_points(p2,dir2);

    num_p1 = cluster1->get_num_points(p1, 10*units::cm);
    num_p2 = cluster2->get_num_points(p2, 10*units::cm);
    num_tp1 = cluster1->get_num_points();
    num_tp2 = cluster2->get_num_points();

    
    
    if (length_1 > 25*units::cm && length_2 > 25*units::cm){
      /* if (length_1 > 60*units::cm && length_2 > 60*units::cm ){ */
      /* 	//if (num_ps_1.second > num_ps_1.first * 0.05 || num_ps_2.second > num_ps_2.first * 0.05) */
      /* 	return false; */
      /* } */

      if ((num_ps_1.second < num_ps_1.first * 0.02 || num_ps_1.second <=3) &&
	  (num_ps_2.second < num_ps_2.first * 0.02 || num_ps_2.second <=3) ||
	  (num_ps_1.second < num_ps_1.first * 0.035 || num_ps_1.second <=6) &&
	  (num_ps_1.second <=1 || num_ps_1.second < num_ps_1.first * 0.005) ||
	  (num_ps_1.second <=1 || num_ps_2.second < num_ps_2.first * 0.005) &&
	  (num_ps_2.second < num_ps_2.first * 0.035 || num_ps_2.second <=6)
	  )
	return true;
    }
  }

  
  
  if (dis < length_cut && (length_2 >=12*units::cm || length_1 >=12*units::cm)){
    Point cluster1_ave_pos = cluster1->calc_ave_pos(p1,10*units::cm);
    Point cluster2_ave_pos = cluster2->calc_ave_pos(p2,10*units::cm);

    TVector3 tempV1(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
    TVector3 tempV2(cluster2_ave_pos.x - cluster1_ave_pos.x, cluster2_ave_pos.y - cluster1_ave_pos.y, cluster2_ave_pos.z - cluster1_ave_pos.z);

    /* if (length_1 > 150*units::cm || length_2 > 150*units::cm) */
    /*   std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << length_1/units::cm << " " << length_2/units::cm << " " << num_p1 << " " << num_p2 << " " << num_tp1 << " " << num_tp2 << std::endl; */
    /* return false; */
    
    // one small the other one is big 
    if (length_1 < 12 *units::cm && num_p1 > 0.5*num_tp1 && (num_p2> 50 || num_p2 > 0.25*num_tp2) ||
	length_2 < 12*units::cm && num_p2 > 0.5*num_tp2 && (num_p1>50 || num_p1 > 0.25*num_tp1) )
      return true;
    
    if (length_1 < 12*units::cm && num_p1 < 0.5*num_tp1) return false;
    if (length_2 < 12*units::cm && num_p2 < 0.5*num_tp2) return false;
        
    if ((num_p1 > 25 || num_p1 > 0.25*num_tp1 ) && (num_p2 > 25 || num_p2 > 0.25*num_tp2)){
      double angle5 = tempV1.Angle(tempV2);
              
      if (length_1 < 60*units::cm || length_2 < 60*units::cm){
	if (angle5 < 30/180.*3.1415926)
	  return true;
	if (angle5 < 90/180.*3.1415926 && (num_p1 > 50 && num_p2 > 50) && (num_p1>75 || num_p2>75))
	  return true;
      }

      
      
      
      if ((length_1 < 60*units::cm || num_p1 >40) && (length_2 < 60*units::cm || num_p2 > 40)){
	
	if ((3.1415926 - dir1.Angle(dir2))/3.1415926*180 < 30 &&
	    (3.1415926 - dir1.Angle(tempV1))/3.1415926*180. < 60 &&
	     dir2.Angle(tempV1)/3.1415926*180.<60 ||
	    (3.1415926 - dir1.Angle(dir2))/3.1415926*180 < 15)
	  return true;

	TVector3 dir3 = cluster1->VHoughTrans(cluster1_ave_pos,50*units::cm); // cluster 1 direction based on hough
	TVector3 dir4 = cluster2->VHoughTrans(cluster2_ave_pos,50*units::cm); // cluster 1 direction based on hough

	if ((3.1415926 - dir3.Angle(dir4))/3.1415926*180 < 25 &&
	    (3.1415926 - dir3.Angle(tempV2))/3.1415926*180. < 15 &&
	     dir4.Angle(tempV2)/3.1415926*180.<15 ||
	    (3.1415926 - dir3.Angle(dir4))/3.1415926*180 < 15)
	  return true;

	if (dis<0.6*units::cm && ((3.1415926 - dir3.Angle(tempV2))/3.1415926*180. < 45 && dir4.Angle(tempV2)/3.1415926*180. < 90 || (3.1415926 - dir3.Angle(tempV2))/3.1415926*180. < 90 && dir4.Angle(tempV2)/3.1415926*180. < 45))
	  return true;
	
      }
    }
  }
  
  return false;
}
