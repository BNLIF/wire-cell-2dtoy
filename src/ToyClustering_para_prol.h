void WireCell2dToy::Clustering_parallel_prolong(WireCell::PR3DClusterSelection& live_clusters, std::map<PR3DCluster*,double>& cluster_length_map, double length_cut){
  // calculate the length ...
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  
  
  std::set<std::pair<PR3DCluster*, PR3DCluster*>> to_be_merged_pairs;
  for (size_t i=0;i!=live_clusters.size();i++){
    PR3DCluster* cluster_1 = live_clusters.at(i);
    for (size_t j=i+1;j<live_clusters.size();j++){
      PR3DCluster* cluster_2 = live_clusters.at(j);
      if (WireCell2dToy::Clustering_2nd_round(cluster_1,cluster_2, cluster_length_map[cluster_1], cluster_length_map[cluster_2], length_cut))
	to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
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


bool WireCell2dToy::Clustering_2nd_round(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double length_1, double length_2, double length_cut){
  cluster1->Create_point_cloud();
  cluster2->Create_point_cloud();
  
  // pick any point and merged cell in cluster1,
  SlimMergeGeomCell *prev_mcell1 = 0;
  SlimMergeGeomCell *prev_mcell2 = 0;
  SlimMergeGeomCell *mcell1 = 0;//cluster1->get_mcells().at(0);
  Point p1;// = mcell1->center();
  SlimMergeGeomCell *mcell2=0;
  Point p2;

  double dis = Find_Closeset_Points(cluster1, cluster2, length_1, length_2, length_cut, mcell1, mcell2, p1,p2);

  /* while(mcell1!=prev_mcell1 || mcell2!=prev_mcell2){ */
  /*   prev_mcell1 = mcell1; */
  /*   prev_mcell2 = mcell2; */
    
  /*   // find the closest point and merged cell in cluster2 */
  /*   std::pair<SlimMergeGeomCell*,Point> temp_results = cluster2->get_closest_point_mcell(p1); */
  /*   p2 = temp_results.second; */
  /*   mcell2 = temp_results.first; */
  /*   // find the closest point and merged cell in cluster1 */
  /*   temp_results = cluster1->get_closest_point_mcell(p2); */
  /*   p1 = temp_results.second; */
  /*   mcell1 = temp_results.first; */
  /* } */

  /* double dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2)); */

  /* if (cluster1->get_cluster_id()==12 || cluster2->get_cluster_id()==12) */
  /*   //  if (length_1 > 50*units::cm && length_2 > 50*units::cm && dis < 5*units::cm) */
  /*   std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << */
  /*     dis/units::cm << " " << length_1/units::cm << " " << length_2/units::cm << " " << length_cut/units::cm << std::endl; */
  
  
  if (dis < length_cut || (dis < 80*units::cm && length_1 +length_2 > 50*units::cm && length_1>15*units::cm && length_2 > 15*units::cm)){
    Point cluster1_ave_pos = cluster1->calc_ave_pos(p1,10*units::cm);
    Point cluster2_ave_pos = cluster2->calc_ave_pos(p2,10*units::cm);

    bool flag_para = false;
    /* bool flag_prolonged_U = false; */
    /* bool flag_prolonged_V = false; */
    bool flag_para_U = false;
    bool flag_para_V = false;
    // parallel case 1 and perpendicular case 2 
    TVector3 drift_dir(1,0,0);
    // pronlonged case for U 3 and V 4 ...
    TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
    TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
    TVector3 W_dir(0,1,0);
    
    // deal the parallel case ...
    {
      TVector3 tempV1(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
      TVector3 tempV2(cluster2_ave_pos.x - cluster1_ave_pos.x, cluster2_ave_pos.y - cluster1_ave_pos.y, cluster2_ave_pos.z - cluster1_ave_pos.z);
      
      double angle1 = tempV1.Angle(drift_dir);
      double angle4 = tempV2.Angle(drift_dir);

      /* if ((length_1 > 50*units::cm && length_2 > 50*units::cm) && dis < 5*units::cm) */
      /* // if (cluster1->get_cluster_id()==7 || cluster2->get_cluster_id()==7 ||  */
      /* // 	  cluster1->get_cluster_id()==9 || cluster2->get_cluster_id()==9) */
      /* 	 std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << length_1/units::cm << " " << length_2/units::cm << " " << dis/units::cm << " " <<  angle1/3.1415926*180. << " " << angle4/3.1415926*180. <<  std::endl; */

      /* return false; */

      // looks like a parallel case
      if ( (fabs(angle1-3.1415926/2.)<10/180.*3.1415926 && dis > 10*units::cm ||
	    fabs(angle1-3.1415926/2.)<20/180.*3.1415926 && dis > 3*units::cm && dis <= 10*units::cm ||
	    fabs(angle1-3.1415926/2.)<45/180.*3.1415926 && dis <=3*units::cm)
	   && fabs(angle4-3.1415926/2.)<5/180.*3.1415926){

	TVector3 dir1 = cluster1->VHoughTrans(p1,60*units::cm); // cluster 1 direction based on hough
	TVector3 dir2 = cluster2->VHoughTrans(p2,60*units::cm); // cluster 2 direction based on hough

	double angle5 = dir1.Angle(drift_dir);
	double angle6 = dir2.Angle(drift_dir);

	if (fabs(angle5-3.1415926/2.)<5/180.*3.1415926 && fabs(angle6-3.1415926/2.)<20/180.*3.1415926 ||
	    fabs(angle5-3.1415926/2.)<20/180.*3.1415926 && fabs(angle6-3.1415926/2.)<5/180.*3.1415926){
	
	  flag_para = true;
	  
	  if (dis >= 3*length_1 && dis >= 3*length_2 && flag_para) return false;
	   
	  double angle2 = tempV1.Angle(U_dir);
	  double angle3 = tempV1.Angle(V_dir);

	  /* if ((length_1 > 100*units::cm || length_2 > 100*units::cm) && dis < 50*units::cm &&( cluster1->get_cluster_id()==50 || cluster2->get_cluster_id()==50)) */
	  /* // if (cluster1->get_cluster_id()==7 || cluster2->get_cluster_id()==7) */
	  /*   std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << (angle2-3.1415926/2.)/3.1415926*180. << " " << (angle3-3.1415926/2.)/3.1415926*180. << " " << (angle5-3.1415926/2.)/3.1415926*180. << " " << (angle6-3.1415926/2.)/3.1415926*180. << " " << dis/units::cm << " " << length_cut/units::cm << " " << length_1 / units::cm << " " << length_2 /units::cm << std::endl; */

	  // look at parallel U
	  if ((fabs(angle2-3.1415926/2.)<7.5/180.*3.1415926 || (fabs(angle2-3.1415926/2.)<15/180.*3.1415926)&&dis <6*units::cm)&& dis < length_cut){
	    flag_para_U = true;
	    //return true;
	    
	    if ((length_1 < 25*units::cm || length_2 < 25*units::cm) && fabs(angle2-3.1415926/2.)<5.0/180.*3.1415926  && dis < 15* units::cm || dis < 3*units::cm){
	      // for short or small distance one
	      return true;
	    }else if (fabs(angle6-3.1415926/2.)/3.1415926*180.<1 && fabs(angle5-3.1415926/2.)/3.1415926*180.<1 && fabs(angle2-3.1415926/2.)<5.0/180.*3.1415926 && dis < 20*units::cm && (length_1 < 50*units::cm || length_2 < 50*units::cm)){
	      return true;
	    }else if (dis < 15*units::cm && (length_1 < 60*units::cm || length_2 < 60*units::cm) &&
		      fabs(angle2-3.1415926/2.)<2.5/180.*3.1415926){
	      // parallel case for reasonably short distance one
	      return true;
	    }else if (fabs(angle2-3.1415926/2.)<2.5/180.*3.1415926 && fabs(angle5-3.1415926/2.)<5/180.*3.1415926 && fabs(angle6-3.1415926/2.)<5/180.*3.141592 ){
	      // parallel case, but exclude both very long tracks
	      if (length_1 < 60*units::cm || length_2 < 60*units::cm){
	      	return true;
	      }else if (dis <5*units::cm){
		return true;
	      }else{
	      	double angle7 = (3.1415926-dir1.Angle(dir2))/3.1415926*180.;
	      	double angle8 = (3.1415926-dir1.Angle(tempV1))/3.1415926*180.; // dir1 = -p1, tempV1 = p2 - p1
	      	double angle9 = dir2.Angle(tempV1)/3.1415926*180.; // dir2 = -p2
	      	if (angle7 < 30 && angle8 < 30 && angle9 < 30)
	      	  return true;
	      }
	      
	    }else{
	      // general case ... (not sure how useful though ...)
	      double angle7 = (3.1415926-dir1.Angle(dir2))/3.1415926*180.;
	      double angle8 = (3.1415926-dir1.Angle(tempV1))/3.1415926*180.; // dir1 = -p1, tempV1 = p2 - p1
	      double angle9 = dir2.Angle(tempV1)/3.1415926*180.; // dir2 = -p2
	      if ((angle7 < 30 && angle8 < 30 && angle9 < 30 ||
		  fabs(angle5-3.1415926/2.)<5/180.*3.1415926 && fabs(angle6-3.1415926/2.)<5/180.*3.141592 &&
		   angle7 < 45 && angle8 < 45 && angle9 < 45) && dis < 20*units::cm)
		return true;
	    }
	  }

	  // look at parallel V
	  if ((fabs(angle3-3.1415926/2.)<7.5/180.*3.1415926 || (fabs(angle3-3.1415926/2.)<15/180.*3.1415926)&&dis <6*units::cm )&& dis < length_cut){
	    flag_para_V = true;
	    //return true;
	    
	    if ((length_1 < 25*units::cm || length_2 < 25*units::cm) && fabs(angle3-3.1415926/2.)<5.0/180.*3.1415926 && dis < 15* units::cm || dis < 2*units::cm){
	      return true;
	    }else if (fabs(angle6-3.1415926/2.)/3.1415926*180.<1 && fabs(angle5-3.1415926/2.)/3.1415926*180.<1 && fabs(angle3-3.1415926/2.)<5.0/180.*3.1415926 && dis < 20*units::cm && (length_1 < 50*units::cm || length_2 < 50*units::cm)){
	      return true;
	    }else if (dis < 15*units::cm && fabs(angle3-3.1415926/2.)<2.5/180.*3.1415926 && (length_1 < 60*units::cm || length_2 < 60*units::cm) ){
	      return true;
	    }else if (fabs(angle3-3.1415926/2.)<2.5/180.*3.1415926 && fabs(angle5-3.1415926/2.)<5/180.*3.1415926 && fabs(angle6-3.1415926/2.)<5/180.*3.141592){
	      return true;
	    }else{
	      double angle7 = (3.1415926-dir1.Angle(dir2))/3.1415926*180.;
	      double angle8 = (3.1415926-dir1.Angle(tempV1))/3.1415926*180.; // dir1 = -p1, tempV1 = p2 - p1
	      double angle9 = dir2.Angle(tempV1)/3.1415926*180.; // dir2 = -p2
	      if (angle7 < 30 && angle8 < 30 && angle9 < 30||
		  fabs(angle5-3.1415926/2.)<5/180.*3.1415926 && fabs(angle6-3.1415926/2.)<5/180.*3.141592 &&
		  angle7 < 60 && angle8 < 60 && angle9 < 60)
		return true;
	    }
	  }
	}
      }
    }
    
    //return false;

    // look at prolonged case ... (add W case) 
    {
      TVector3 tempV1(0, p2.y - p1.y, p2.z - p1.z);
      TVector3 tempV5;
      double angle1 = tempV1.Angle(U_dir);
      tempV5.SetXYZ(fabs(p2.x-p1.x),sqrt(pow(p2.y - p1.y,2)+pow(p2.z - p1.z,2))*sin(angle1),0);
      angle1 = tempV5.Angle(drift_dir);
      
      double angle2 = tempV1.Angle(V_dir);
      tempV5.SetXYZ(fabs(p2.x-p1.x),sqrt(pow(p2.y - p1.y,2)+pow(p2.z - p1.z,2))*sin(angle2),0);
      angle2 = tempV5.Angle(drift_dir);

      double angle1p = tempV1.Angle(W_dir);
      tempV5.SetXYZ(fabs(p2.x-p1.x),sqrt(pow(p2.y - p1.y,2)+pow(p2.z - p1.z,2))*sin(angle1p),0);
      angle1p = tempV5.Angle(drift_dir);

      
      
      if (angle1<15/180.*3.1415926  ||
	  angle2<15/180.*3.1415926  ||
	  angle1p<15/180.*3.1415926 ){
	if (length_1 > 10*units::cm || length_2 > 10*units::cm){
	  TVector3 dir1 = cluster1->VHoughTrans(p1,60*units::cm); // cluster 1 direction based on hough
	  TVector3 dir2 = cluster2->VHoughTrans(p2,60*units::cm); // cluster 1 direction based on hough
	  TVector3 dir3(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z);
	  double angle3 = dir3.Angle(dir2);
	  double angle4 = 3.1415926-dir3.Angle(dir1);

	  /* if (cluster1->get_cluster_id()==12 || cluster2->get_cluster_id()==12) */
	  /*   std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << angle1*180./3.1415926 << " " << angle2*180./3.1415926 << " " << dis/units::cm << " " << length_1/units::cm << " " << length_2/units::cm << " " << angle3/3.1415826*180. << " " << angle4/3.1415926*180. << std::endl; */
	  // if (angle1p<15/180.*3.1415926 || (3.1415926-angle1p)<15/180.*3.1415926 )
	  //   std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << angle3/3.1415926*180. << " " << angle4/3.1415926*180. << " " << dis/units::cm << std::endl;
	  // if (cluster1->get_cluster_id()==52 || cluster2->get_cluster_id()==52)
	  //   std::cout << angle3/3.1415926*180 << " " << angle4/3.1415926*180 << std::endl;
	  if ((angle3<25/180.*3.1415926 || length_2<12*units::cm)&&(angle4<25/180.*3.1415926|| length_1<12*units::cm)&&dis<5*units::cm ||
	      (angle3<15/180.*3.1415926 || length_2<12*units::cm)&&(angle4<15/180.*3.1415926|| length_1<12*units::cm)&&dis<15*units::cm ||
	      (angle3<7.5/180.*3.1415926 || length_2<12*units::cm)&&(angle4<7.5/180.*3.1415926|| length_1<12*units::cm) ||
	      (angle3+angle4 < 15/180.*3.1415926 && angle3 < 10/180.*3.1415926 && angle4 < 10/180.*3.1415926)
	      )
	    return true;

	 
	      
	}
      }else{
	//regular cases (only for very short distance ... )
	if (dis < 5*units::cm){
	  if (length_1 > 10*units::cm || length_2 > 10*units::cm){
	    TVector3 dir1 = cluster1->VHoughTrans(p1,30*units::cm); // cluster 1 direction based on hough
	    TVector3 dir2 = cluster2->VHoughTrans(p2,30*units::cm); // cluster 1 direction based on hough
	    TVector3 dir3(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z);
	    double angle3 = dir3.Angle(dir2);
	    double angle4 = 3.1415926-dir3.Angle(dir1);

	    //std::cout << angle3/3.1415926*180. << " " << angle4/3.1415926*180. << std::endl;
	    if ((angle3<15/180.*3.1415926 || length_2<6*units::cm)&&(angle4<15/180.*3.1415926|| length_1<6*units::cm))
	      return true;
	  }
	}
	
      }
      
    }
  }
  
  return false;
}
