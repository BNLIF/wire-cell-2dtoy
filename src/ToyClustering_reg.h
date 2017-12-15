void WireCell2dToy::Clustering_regular(WireCell::PR3DClusterSelection& live_clusters, std::map<PR3DCluster*,double>& cluster_length_map, double length_cut, bool flag_enable_extend){
  
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
	if (WireCell2dToy::Clustering_1st_round(cluster_1,cluster_2, cluster_length_map[cluster_1], cluster_length_map[cluster_2], length_cut, flag_enable_extend))
	  to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	
      }
    }
    
    // to_be_merged_pairs.clear();
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
  //  return cluster_length_map;
}



bool WireCell2dToy::Clustering_1st_round(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double length_1, double length_2, double length_cut, bool flag_enable_extend){
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

  if (dis < length_cut){
    bool flag_para = false;
    bool flag_prolong_U = false;
    bool flag_prolong_V = false;
    bool flag_prolong_W= false;
    bool flag_para_U = false;
    bool flag_para_V = false;
    bool flag_para_W = false;
    bool flag_regular = false;
    bool flag_extend = false;
    bool flag_force_extend = false;

    // pronlonged case for U 3 and V 4 ...
    TVector3 drift_dir(1,0,0);
    TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
    TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
    TVector3 W_dir(0,1,0);

    // calculate average distance ... 
    Point cluster1_ave_pos = cluster1->calc_ave_pos(p1,5*units::cm);
    Point cluster2_ave_pos = cluster2->calc_ave_pos(p2,5*units::cm);

    TVector3 dir2_1(p2.x - p1.x+1e-9, p2.y - p1.y+1e-9, p2.z - p1.z+1e-9); // 2-1
    TVector3 dir2(cluster2_ave_pos.x - cluster1_ave_pos.x+1e-9, cluster2_ave_pos.y - cluster1_ave_pos.y+1e-9, cluster2_ave_pos.z - cluster1_ave_pos.z+1e-9); // 2-1
    dir2_1.SetMag(1);
    dir2.SetMag(1);
    
    
    // parallle case
    
    double angle1 = dir2_1.Angle(drift_dir);
    double angle2 = dir2.Angle(drift_dir);
    
    if (fabs(angle1-3.1415926/2.)<7.5/180.*3.1415926 ||
	fabs(angle2-3.1415926/2.)<7.5/180.*3.1415926){
      flag_para = true;
      
      double angle3 = dir2_1.Angle(U_dir);
      double angle4 = dir2_1.Angle(V_dir);
      double angle5 = dir2_1.Angle(W_dir);
      
      double angle3_1 = dir2.Angle(U_dir);
      double angle4_1 = dir2.Angle(V_dir);
      double angle5_1 = dir2.Angle(W_dir);
      
      if (fabs(angle3-3.1415926/2.)<7.5/180.*3.1415926 || fabs(angle3_1-3.1415926/2.)<7.5/180.*3.1415926 || ((fabs(angle3-3.1415926/2.)<15/180.*3.1415926 || fabs(angle3_1-3.1415926/2.)<15/180.*3.1415926)&&dis < 6*units::cm))
	flag_para_U = true;
      
      if (fabs(angle4-3.1415926/2.)<7.5/180.*3.1415926 || fabs(angle4_1-3.1415926/2.)<7.5/180.*3.1415926 || ((fabs(angle4-3.1415926/2.)<15/180.*3.1415926 || fabs(angle4_1-3.1415926/2.)<15/180.*3.1415926)&&dis < 6*units::cm))
	flag_para_V = true;
      
      if (fabs(angle5-3.1415926/2.)<7.5/180.*3.1415926 || fabs(angle5_1-3.1415926/2.)<7.5/180.*3.1415926 || ((fabs(angle5-3.1415926/2.)<15/180.*3.1415926 || fabs(angle5_1-3.1415926/2.)<15/180.*3.1415926)&&dis < 6*units::cm))
	flag_para_W = true;
    } 
    /*   if ( (cluster1->get_cluster_id()==46 || cluster2->get_cluster_id()==46) && (cluster1->get_cluster_id()==55 || cluster2->get_cluster_id()==55)) */
    /*   	std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << dis/units::cm << " " << length_1/units::cm << " " << length_2/units::cm << " " << flag_para << " " << flag_para_U << " " << flag_para_V << " " << flag_para_W << " " << */
    /*   	 fabs(angle3-3.1415926/2.)/3.1415926*180. << " " << fabs(angle3_1-3.1415926/2.)/3.1415926*180. << */
    /* 	   " " << fabs(angle4-3.1415926/2.)/3.1415926*180. << " " << fabs(angle4_1-3.1415926/2.)/3.1415926*180. <<std::endl; */
    /* return false; */
    
   
    
    
    // prolonged case
    
    TVector3 tempV3(0, p2.y - p1.y, p2.z - p1.z);
    TVector3 tempV4(0, cluster2_ave_pos.y - cluster1_ave_pos.y, cluster2_ave_pos.z - cluster1_ave_pos.z);
    TVector3 tempV5;
    
    double angle6 = tempV3.Angle(U_dir);
    tempV5.SetXYZ(fabs(p2.x-p1.x),sin(angle6),0);
    angle6 = tempV5.Angle(drift_dir);
    
    double angle7 = tempV3.Angle(V_dir);
    tempV5.SetXYZ(fabs(p2.x-p1.x),sin(angle7),0);
    angle7 = tempV5.Angle(drift_dir);
    
    double angle8 = tempV3.Angle(W_dir);
    tempV5.SetXYZ(fabs(p2.x-p1.x),sin(angle8),0);
    angle8 = tempV5.Angle(drift_dir);
    
    double angle6_1 = tempV4.Angle(U_dir);
    tempV5.SetXYZ(fabs(cluster2_ave_pos.x-cluster1_ave_pos.x),sin(angle6_1),0);
    angle6_1 = tempV5.Angle(drift_dir);
    double angle7_1 = tempV4.Angle(V_dir);
    tempV5.SetXYZ(fabs(cluster2_ave_pos.x-cluster1_ave_pos.x),sin(angle7_1),0);
    angle7_1 = tempV5.Angle(drift_dir);
    
    double angle8_1 = tempV4.Angle(W_dir);
    tempV5.SetXYZ(fabs(cluster2_ave_pos.x-cluster1_ave_pos.x),sin(angle8_1),0);
    angle8_1 = tempV5.Angle(drift_dir);
    
    if (angle6<7.5/180.*3.1415926  ||
	angle6_1<7.5/180.*3.1415926 )
      flag_prolong_U = true;
    
    if (angle7<7.5/180.*3.1415926  ||
	angle7_1<7.5/180.*3.1415926 )
      flag_prolong_V = true;
    
    if (angle8<7.5/180.*3.1415926  ||
	angle8_1<7.5/180.*3.1415926 )
      flag_prolong_W = true;
    
    
    
    // regular case
    
    if (dis <= 15*units::cm) flag_regular = true;
    
    if ((flag_para_U || flag_para_V || flag_para_W) && flag_para ||
	flag_prolong_U || flag_prolong_V || flag_prolong_W || flag_regular){
      double angle_cut=2.5;
      double para_angle_cut = 5.;
      double para_angle_cut_1 = 5;

      if (dis < 5*units::cm){
	angle_cut = 10;
      }else if (dis < 15*units::cm){
	angle_cut = 7.5;
      }else{
	angle_cut = 5;
      }

      if (dis > 45*units::cm){
	para_angle_cut = 5;
      }else if (dis > 15*units::cm){
	para_angle_cut = 15;
      }else if (dis > 5*units::cm){
	para_angle_cut = 30;
      }else{
	para_angle_cut = 60;
      }
      
      TVector3 dir1 = cluster1->VHoughTrans(cluster1_ave_pos,30*units::cm); // cluster 1 direction based on hough
      TVector3 dir1_1 = cluster1->VHoughTrans(p1,30*units::cm);
      
      TVector3 dir3 = cluster2->VHoughTrans(cluster2_ave_pos,30*units::cm);
      TVector3 dir3_1 = cluster2->VHoughTrans(p2,30*units::cm);

      dir1_1.SetMag(1);
      dir1.SetMag(1);
      dir3.SetMag(1);
      dir3_1.SetMag(1);
      
      double angle_diff1 = (3.1415926-dir1.Angle(dir2))/3.1415926*180.;
      double angle_diff1_1 = (3.1415926-dir1_1.Angle(dir2_1))/3.1415926*180.;

      double angle_diff2 = dir3.Angle(dir2)/3.1415926*180;
      double angle_diff2_1 = dir3_1.Angle(dir2_1)/3.1415926*180.;

      double angle_diff3 = (3.1415926 - dir1.Angle(dir3))/3.1415926*180;
      double angle_diff3_1 = (3.1415926 - dir1_1.Angle(dir3_1))/3.1415926*180.;

      if (dis<=3*units::cm){
	if ((angle_diff3 < angle_cut*1.5 || angle_diff3_1 < angle_cut*1.5) &&
	    ((angle_diff1 < angle_cut*4.5 || angle_diff1_1 < angle_cut*4.5 || length_1 < 12*units::cm) &&
	     (angle_diff2 < angle_cut*4.5 || angle_diff2_1 < angle_cut*4.5 || length_2 < 12*units::cm) ) )
	  return true;
	if ((angle_diff3 < angle_cut*1.5 || angle_diff3_1 < angle_cut*1.5) ||
	    (angle_diff1 < angle_cut*1.5 || angle_diff1_1 < angle_cut*1.5) ||
	    (angle_diff2 < angle_cut*1.5 || angle_diff2_1 < angle_cut*1.5))
	  flag_extend = true;
      }

      
      if (angle_diff1 < angle_cut || angle_diff1_1 < angle_cut){ // possible match 
	if (length_2 < 12*units::cm && dis < 2*units::cm) //small and very close
	  return true;
	
	if (angle_diff3 < angle_cut || angle_diff3_1 < angle_cut ){
	  return true;
	}
	 flag_extend = true;
      }

      if (angle_diff2 < angle_cut || angle_diff2_1 < angle_cut){
	if (length_1 < 12*units::cm && dis < 2*units::cm) // small and very close
	  return true;

	if (angle_diff3 < angle_cut || angle_diff3_1 < angle_cut)
	  return true;
	flag_extend = true;
      }

      /* if (flag_para && (!flag_para_U) && (!flag_para_V) && (!flag_para_W) && flag_regular && */
      /* 	  length_1 > 25*units::cm && length_2 > 25*units::cm) */
      /* 	flag_force_extend = true; */
      
      
      if (flag_para && (flag_para_U || flag_para_V || flag_para_W || flag_regular)){
	
	double dangle1 = (dir1.Angle(drift_dir)-3.1415926/2.)/3.1415926*180.;
	double dangle1_1 = (dir1_1.Angle(drift_dir)-3.1415926/2.)/3.1415926*180.;
	
	double dangle2 = (dir2.Angle(drift_dir)-3.1415926/2.)/3.1415926*180.;
	double dangle2_1 = (dir2_1.Angle(drift_dir)-3.1415926/2.)/3.1415926*180.;

	double dangle3 = (dir3.Angle(drift_dir)-3.1415926/2.)/3.1415926*180.;
	double dangle3_1 = (dir3_1.Angle(drift_dir)-3.1415926/2.)/3.1415926*180.;
	
	double dangle4 = dangle1 + dangle2;
	double dangle4_1 = dangle1_1 + dangle2_1;

	double dangle5 = dangle3 - dangle2;
	double dangle5_1 = dangle3_1 - dangle2_1;

	
	/* if (dis < 10*units::cm && (length_1 > 50*units::cm && length_2 > 50*units::cm))  */
	/*   std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << dis/units::cm << " " << length_1/units::cm << " " << length_2/units::cm << " " */
	/*     << dangle1 << " " << dangle1_1 << " " << dangle2 << " " << dangle2_1 << " " << dangle4 */
	/*     << " " << dangle4_1 << " " << angle_diff3 << " " << angle_diff3_1 << " " << */
	/*     angle_diff1 << " "<< angle_diff1_1 << " " << angle_diff2 << " " << angle_diff2_1 << " " */
	/* 	    << flag_para << " " << flag_para_U << " " << flag_para_V */
	/* 	    << std::endl; */

	
	if ((fabs(dangle1) < 2.5 || fabs(dangle1_1) < 2.5) &&
	    (fabs(dangle2) < 2.5 || fabs(dangle2_1) < 2.5) &&
	    (fabs(dangle3) < 2.5 || fabs(dangle3_1) < 2.5) &&
	    (length_1 > 25*units::cm && length_2 > 25*units::cm))
	  flag_force_extend = true;
	
	if (flag_para_U || flag_para_V || flag_para_W){
	  if ( (fabs(dangle1) < para_angle_cut_1 || fabs(dangle1_1) < para_angle_cut_1) &&
	       (fabs(dangle2) < para_angle_cut_1 || fabs(dangle2_1) < para_angle_cut_1) &&
	       (fabs(dangle4) < para_angle_cut_1/2. || fabs(dangle4_1) < para_angle_cut_1/2.) &&
	       (angle_diff1 < para_angle_cut || angle_diff1_1 < para_angle_cut)){
	    if (length_2 < 12*units::cm && dis < 2*units::cm) //small and very close
	      return true;
	    
	    if (angle_diff3 < para_angle_cut || angle_diff3_1 < para_angle_cut )
	      return true;
	    
	    flag_extend = true;
	    
	  }
	  
	  
	  if ( (fabs(dangle3) < para_angle_cut_1 || fabs(dangle3_1) < para_angle_cut_1) &&
	       (fabs(dangle2) < para_angle_cut_1 || fabs(dangle2_1) < para_angle_cut_1) &&
	       (fabs(dangle5) < para_angle_cut_1/2. || fabs(dangle5_1) < para_angle_cut_1/2.) &&
	       (angle_diff2 < para_angle_cut || angle_diff2_1 < para_angle_cut)){
	    if (length_1 < 12*units::cm && dis < 2*units::cm) //small and very close
	      return true;
	    
	    if (angle_diff3 < para_angle_cut || angle_diff3_1 < para_angle_cut )
	      return true;
	    
	    flag_extend = true;
	  }

	  if (dis<5*units::cm && length_1 > 25*units::cm && length_2 > 25*units::cm)
	    flag_extend = true;
	}
      }

      
      /* if (dis < 10*units::cm && (length_1 > 50*units::cm && length_2 > 50*units::cm)) */
      /*  	std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << dis/units::cm << " " << length_1/units::cm << " " << length_2/units::cm << " " << flag_extend << " " << flag_force_extend << std::endl; */
      
      
      if (flag_extend && flag_enable_extend || flag_force_extend){
      	// when to extend???
		
      	if ((flag_para && (flag_para_U || flag_para_V || flag_para_W || flag_regular)) ||
	    (!flag_para && (flag_prolong_U || flag_prolong_V || dis < 10*units::cm))){
      	  // look use 1 to predict 2
      	  // cluster1_ave_pos, dir1
      	  // calculate the average distance
      	  double ave_dis = sqrt(pow(cluster1_ave_pos.x-cluster2_ave_pos.x,2) + pow(cluster1_ave_pos.y-cluster2_ave_pos.y,2) + pow(cluster1_ave_pos.z-cluster2_ave_pos.z,2));
      	  Point test_point;
      	  double min_dis = 1e9, max_dis = -1e9;
      	  for (int i=-5;i!=6;i++){
      	    test_point.x = cluster1_ave_pos.x - dir1.X() * (ave_dis +i*2*units::cm);
      	    test_point.y = cluster1_ave_pos.y - dir1.Y() * (ave_dis +i*2*units::cm);
      	    test_point.z = cluster1_ave_pos.z - dir1.Z() * (ave_dis +i*2*units::cm);
	    
      	    std::pair<SlimMergeGeomCell*,Point> temp_results = cluster2->get_closest_point_mcell(test_point);
      	    //reuse this
      	    Point test_point1 = temp_results.second;
      	    if (sqrt(pow(test_point1.x-test_point.x,2)+pow(test_point1.y-test_point.y,2)+pow(test_point1.z-test_point.z,2))<1.5*units::cm){
      	      double temp_dis = (test_point1.x - cluster1_ave_pos.x)*dir1.X() + (test_point1.y - cluster1_ave_pos.y)*dir1.Y() + (test_point1.z - cluster1_ave_pos.z)*dir1.Z();
      	      temp_dis *=-1;
      	      if (temp_dis < min_dis) min_dis = temp_dis;
      	      if (temp_dis > max_dis) max_dis = temp_dis;
      	    }
      	  }
	  // std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << min_dis/units::cm << " " << max_dis/units::cm << " " << length_2/units::cm << std::endl;
	  
	  
	  
      	  if ((max_dis - min_dis)>2.5*units::cm) return true;
	  
      	  // look at the other side (repeat)
      	  // cluster2_ave_pos, dir2
      	  min_dis = 1e9;
      	  max_dis = -1e9;
      	  for (int i=-5;i!=6;i++){
      	    test_point.x = cluster2_ave_pos.x - dir3.X() * (ave_dis +i*2*units::cm);
      	    test_point.y = cluster2_ave_pos.y - dir3.Y() * (ave_dis +i*2*units::cm);
      	    test_point.z = cluster2_ave_pos.z - dir3.Z() * (ave_dis +i*2*units::cm);
	    
      	    std::pair<SlimMergeGeomCell*,Point> temp_results = cluster1->get_closest_point_mcell(test_point);
      	    //reuse this
      	    Point test_point1 = temp_results.second;
      	    if (sqrt(pow(test_point1.x-test_point.x,2)+pow(test_point1.y-test_point.y,2)+pow(test_point1.z-test_point.z,2))<1.5*units::cm){
      	      double temp_dis = (test_point1.x - cluster2_ave_pos.x)*dir3.X() + (test_point1.y - cluster2_ave_pos.y)*dir3.Y() + (test_point1.z - cluster2_ave_pos.z)*dir3.Z();
      	      temp_dis *=-1;
      	      if (temp_dis < min_dis) min_dis = temp_dis;
      	      if (temp_dis > max_dis) max_dis = temp_dis;
      	    }
      	  }
	  // std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << min_dis/units::cm << " " << max_dis/units::cm << " " << length_1/units::cm << std::endl;


      	  if ((max_dis - min_dis)>2.5*units::cm) return true;
	  
      	  // both side simutaneously? leave it for futhre
	  
      	}
      }
      
      
    }
  }

  
  return false;
}




  
