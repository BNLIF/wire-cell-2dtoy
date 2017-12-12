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
  //  return cluster_length_map;
}




  
bool WireCell2dToy::Clustering_1st_round(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double length_1, double length_2, double length_cut, bool flag_enable_extend){
  cluster1->Create_point_cloud();
  cluster2->Create_point_cloud();
  
  // pick any point and merged cell in cluster1,
  SlimMergeGeomCell *prev_mcell1 = 0;
  SlimMergeGeomCell *prev_mcell2 = 0;
  SlimMergeGeomCell *mcell1 = cluster1->get_mcells().at(0);
  Point p1 = mcell1->center();
  SlimMergeGeomCell *mcell2=0;
  Point p2;

  while(mcell1!=prev_mcell1 || mcell2!=prev_mcell2){
    prev_mcell1 = mcell1;
    prev_mcell2 = mcell2;
    
    // find the closest point and merged cell in cluster2
    std::pair<SlimMergeGeomCell*,Point> temp_results = cluster2->get_closest_point_mcell(p1);
    p2 = temp_results.second;
    mcell2 = temp_results.first;
    // find the closest point and merged cell in cluster1
    temp_results = cluster1->get_closest_point_mcell(p2);
    p1 = temp_results.second;
    mcell1 = temp_results.first;
  }
  
  //std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2))/units::cm<< std::endl;
  
  double dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));

  // if (cluster1->get_cluster_id()==7){
  //   std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << dis/units::cm << " " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm <<
  //     " " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm <<
  //     std::endl;
  // }
  
  
  if (dis < length_cut){
    Point mcell1_center = mcell1->center();
    Point mcell2_center = mcell2->center();
    Point cluster1_ave_pos = cluster1->calc_ave_pos(p1,5*units::cm);
    Point cluster2_ave_pos = cluster2->calc_ave_pos(p2,5*units::cm);

    bool flag_para = false;
    bool flag_para_U = false;
    bool flag_para_V = false;
    bool flag_perp = false;
    bool flag_prolonged_U = false;
    bool flag_prolonged_V = false;
    // parallel case 1 and perpendicular case 2 
    TVector3 drift_dir(1,0,0);
    // pronlonged case for U 3 and V 4 ...
    TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
    TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
    
    {
      TVector3 tempV1(mcell2_center.x - mcell1_center.x, mcell2_center.y - mcell1_center.y, mcell2_center.z - mcell1_center.z);
      TVector3 tempV2(cluster2_ave_pos.x - mcell1_center.x, cluster2_ave_pos.y - mcell1_center.y, cluster2_ave_pos.z - mcell1_center.z);
      TVector3 tempV3(mcell2_center.x - cluster1_ave_pos.x, mcell2_center.y - cluster1_ave_pos.y, mcell2_center.z - cluster1_ave_pos.z);
      TVector3 tempV4(cluster2_ave_pos.x - cluster1_ave_pos.x, cluster2_ave_pos.y - cluster1_ave_pos.y, cluster2_ave_pos.z - cluster1_ave_pos.z);
      TVector3 tempV5(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z);
      
      double angle1 = tempV1.Angle(drift_dir);
      double angle2 = tempV2.Angle(drift_dir);
      double angle3 = tempV3.Angle(drift_dir);
      double angle4 = tempV4.Angle(drift_dir);
      if (angle1<10/180.*3.1415926 || (3.1415926-angle1)<10/180.*3.1415926
	  ||angle2<10/180.*3.1415926 || (3.1415926-angle2)<10/180.*3.1415926
	  ||angle3<10/180.*3.1415926 || (3.1415926-angle3)<10/180.*3.1415926
	  ||angle4<10/180.*3.1415926 || (3.1415926-angle4)<10/180.*3.1415926)
	flag_perp = true;

      double angle5 = tempV5.Angle(U_dir);
      double angle6 = tempV5.Angle(V_dir);
      if (fabs(angle5-3.1415926/2.)<10/180.*3.1415926)
	flag_para_U = true;
      
      if (fabs(angle6-3.1415926/2.)<10/180.*3.1415926)
	flag_para_V = true;


      if (fabs(angle1-3.1415926/2.)<5/180.*3.1415926
	  ||fabs(angle2-3.1415926/2.)<5/180.*3.1415926
	  ||fabs(angle3-3.1415926/2.)<5/180.*3.1415926
	  ||fabs(angle4-3.1415926/2.)<5/180.*3.1415926)
	flag_para = true;
    }
    
    

    {
      TVector3 tempV1(0, mcell2_center.y - mcell1_center.y, mcell2_center.z - mcell1_center.z);
      TVector3 tempV2(0, cluster2_ave_pos.y - mcell1_center.y, cluster2_ave_pos.z - mcell1_center.z);
      TVector3 tempV3(0, mcell2_center.y - cluster1_ave_pos.y, mcell2_center.z - cluster1_ave_pos.z);
      TVector3 tempV4(0, cluster2_ave_pos.y - cluster1_ave_pos.y, cluster2_ave_pos.z - cluster1_ave_pos.z);

      TVector3 tempV5;
      double angle1 = tempV1.Angle(U_dir);
      tempV5.SetXYZ(fabs(mcell2_center.x - mcell1_center.x),sin(angle1),0);
      angle1 = tempV5.Angle(drift_dir);
      
      double angle2 = tempV2.Angle(U_dir);
      tempV5.SetXYZ(fabs(cluster2_ave_pos.x - mcell1_center.x),sin(angle2),0);
      angle2 = tempV5.Angle(drift_dir);
      
      double angle3 = tempV3.Angle(U_dir);
      tempV5.SetXYZ(fabs(mcell2_center.x - cluster1_ave_pos.x),sin(angle3),0);
      angle3 = tempV5.Angle(drift_dir);

      double angle4 = tempV4.Angle(U_dir);
      tempV5.SetXYZ(fabs(cluster2_ave_pos.x - cluster1_ave_pos.x),sin(angle4),0);
      angle4 = tempV5.Angle(drift_dir);

      if (angle1<5/180.*3.1415926 
	  ||angle2<5/180.*3.1415926 
	  ||angle3<5/180.*3.1415926
	  ||angle4<5/180.*3.1415926  )
	flag_prolonged_U = true;
      
      angle1 = tempV1.Angle(V_dir);
      tempV5.SetXYZ(fabs(mcell2_center.x - mcell1_center.x),sin(angle1),0);
      angle1 = tempV5.Angle(drift_dir);
      
      angle2 = tempV2.Angle(V_dir);
      tempV5.SetXYZ(fabs(cluster2_ave_pos.x - mcell1_center.x),sin(angle2),0);
      angle2 = tempV5.Angle(drift_dir);
      
      angle3 = tempV3.Angle(V_dir);
      tempV5.SetXYZ(fabs(mcell2_center.x - cluster1_ave_pos.x),sin(angle3),0);
      angle3 = tempV5.Angle(drift_dir);
      
      angle4 = tempV4.Angle(V_dir);
      tempV5.SetXYZ(fabs(cluster2_ave_pos.x - cluster1_ave_pos.x),sin(angle4),0);
      angle4 = tempV5.Angle(drift_dir);

      if (angle1<5/180.*3.1415926 
	  ||angle2<5/180.*3.1415926 
	  ||angle3<5/180.*3.1415926 
	  ||angle4<5/180.*3.1415926 )
	flag_prolonged_V = true;
    }
    
    
    bool flag_dir = false;
    double angle_cut=2.5;
    double para_angle_cut = 5.;
    double para_angle_cut_1 = 5;
    //    double point_angle_cut = 15;
    // if (dis < 2.5*units::cm){
    //   flag_dir = true;
    //   angle_cut = 20;
    // }else 
    if (dis < 5*units::cm ){
      // normal case ..., allow for 3 cm gap, and 10 degree cut ...
      flag_dir = true;
      angle_cut = 10;
    }else if (dis < 15*units::cm){
      flag_dir = true;
      angle_cut = 7.5;
    }else if (dis < 30*units::cm){
      flag_dir = true;
      angle_cut = 5.0;
    }
    
    if ((flag_prolonged_U||flag_prolonged_V||flag_perp)&&dis<45*units::cm){
      // normal case ..., allow for 45 cm gap
      flag_dir = true;
      // if (dis < 2.5*units::cm){
      // 	flag_dir = true;
      // 	angle_cut = 20;
      // // // }else if (dis < 2.5*units::cm){
      // // // 	flag_dir = true;
      // // // 	angle_cut = 20;
      // }else
      if (dis < 5*units::cm){
	angle_cut = 10;
      }else if (dis < 15*units::cm){
	angle_cut = 7.5;
      }else{
	angle_cut = 5;
      }
    }

    
    if (flag_para){
      flag_dir = true;
      // if (dis < 2.5*units::cm){
      // 	flag_dir = true;
      // 	angle_cut = 20;
      // // // }else if (dis < 2.5*units::cm){
      // // // 	flag_dir = true;
      // // // 	angle_cut = 20;
      // }else
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
      }else{
	para_angle_cut = 30;
      }
    }

    // {
    //   TVector3 tempV1(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
    //   TVector3 tempV2(cluster2_ave_pos.x - cluster1_ave_pos.x, cluster2_ave_pos.y - cluster1_ave_pos.y, cluster2_ave_pos.z - cluster1_ave_pos.z);
      
    // }

    
    
    if (flag_dir){
      bool flag_extend = false;
      
      // also with a complicated structures ...
      TVector3 dir1 = cluster1->VHoughTrans(cluster1_ave_pos,30*units::cm); // cluster 1 direction based on hough
      TVector3 dir1_1 = cluster1->VHoughTrans(p1,30*units::cm);
      
      TVector3 dir2 = cluster2->calc_dir(cluster1_ave_pos,cluster2_ave_pos,30*units::cm); // dir 2 --> 1
      TVector3 dir2_1 = cluster2->VHoughTrans(cluster1_ave_pos,30*units::cm); // 
      TVector3 dir2_2 = cluster2->calc_dir(cluster1_ave_pos,mcell2_center,10*units::cm); // dir 2 --> 1
      
      TVector3 dir3 = cluster2->VHoughTrans(cluster2_ave_pos,30*units::cm);
      TVector3 dir3_1 = cluster2->VHoughTrans(p2,30*units::cm);
      
      TVector3 dir4 = cluster1->calc_dir(cluster2_ave_pos,cluster1_ave_pos,30*units::cm); // dir 2 --> 1
      TVector3 dir4_1 = cluster1->VHoughTrans(cluster2_ave_pos,30*units::cm);//
      TVector3 dir4_2 = cluster1->calc_dir(cluster2_ave_pos,mcell1_center,10*units::cm); // dir 2 --> 1

	
      double angle_diff1 = (3.1415926-dir1.Angle(dir2))/3.1415926*180.; // 1 to 1->2
      double angle_diff1_1 = (3.1415926-dir1.Angle(dir2_1))/3.1415926*180.; // 1 to 2
      double angle_diff1_2 = (3.1415926-dir1.Angle(dir2_2))/3.1415926*180.; // 1 to 1->2

      double angle_diff2 = (3.1415926-dir1.Angle(dir3))/3.1415926*180.;  // 1 to 2
      double angle_diff2_1 = (3.1415926-dir1_1.Angle(dir3_1))/3.1415926*180.; // 1 to 2
      double angle_diff2_2 = (3.1415926-dir1.Angle(dir3_1))/3.1415926*180.; // 1 to 
      double angle_diff2_3 = (3.1415926-dir1_1.Angle(dir3))/3.1415926*180.;

      
      
      // possible match
      if ((angle_diff1 < angle_cut) || angle_diff1_1 < angle_cut || (angle_diff1_2 <angle_cut)){
	// direction match 
	if (angle_diff2 < angle_cut || angle_diff2_1 < angle_cut || (angle_diff2_2 < angle_cut && length_1 >= length_2) || (angle_diff2_3 < angle_cut && length_1 <= length_2)){
	  return true;
	}else{
	  if (length_2 < 12*units::cm && dis < 2*units::cm) //small and very close
	    return true;
	}
	flag_extend = true;
      }

      double angle_diff3 = (3.1415926-dir3.Angle(dir4))/3.1415926*180.;
      double angle_diff3_1 = (3.1415926-dir3.Angle(dir4_1))/3.1415926*180.;
      double angle_diff3_2 = (3.1415926-dir3.Angle(dir4_2))/3.1415926*180.;
      
      if ((angle_diff3 < angle_cut) || angle_diff3_1 < angle_cut || (angle_diff3_2 <angle_cut)){
	if (angle_diff2 < angle_cut || angle_diff2_1 < angle_cut || (angle_diff2_2 < angle_cut && length_1 >= length_2) || (angle_diff2_3 < angle_cut && length_1 <= length_2)) {
	  return true;
	}else{
	  if (length_1 < 12*units::cm && dis < 2*units::cm) // small and very close
	    return true;
	}
	flag_extend = true;
      }

      

      if (flag_para && (flag_para_U || flag_para_V)){
	TVector3 dir1_rot(dir1.Y(), dir1.Z(), dir1.X());
	TVector3 dir1_1_rot(dir1_1.Y(), dir1_1.Z(), dir1_1.X());
	TVector3 dir2_rot(dir2.Y(), dir2.Z(), dir2.X());
	
	double theta1_1 = (dir1_rot.Theta()-3.1415926/2.)/3.1415926*180.;
	double theta1_2 = (dir2_rot.Theta()-3.1415926/2.)/3.1415926*180.;
	double dphi1 = fabs(3.1415926 - fabs(dir1_rot.Phi()-dir2_rot.Phi()))/3.1415926*180.;
	
	TVector3 dir3_rot(dir3.Y(), dir3.Z(), dir3.X());
	TVector3 dir3_1_rot(dir3_1.Y(), dir3_1.Z(), dir3_1.X());
	TVector3 dir4_rot(dir4.Y(), dir4.Z(), dir4.X());
	
	double dphi3 = fabs(3.1415926 - fabs(dir1_rot.Phi()-dir3_rot.Phi()))/3.1415926*180.;
	double dphi3_1 = fabs(3.1415926 - fabs(dir1_1_rot.Phi()-dir3_1_rot.Phi()))/3.1415926*180.;
	double dphi3_2 = fabs(3.1415926 - fabs(dir1_1_rot.Phi()-dir3_rot.Phi()))/3.1415926*180.;
	double dphi3_3 = fabs(3.1415926 - fabs(dir1_rot.Phi()-dir3_1_rot.Phi()))/3.1415926*180.;
	
	if ( (fabs(theta1_1) < para_angle_cut_1 && fabs(theta1_2) < para_angle_cut_1 && fabs(theta1_1+theta1_2) < para_angle_cut_1/2. && dphi1 < para_angle_cut)){
	  if (dphi3 < para_angle_cut || dphi3_1 < para_angle_cut || (dphi3_2 < para_angle_cut && length_1 >= length_2) || (dphi3_3 < para_angle_cut && length_1 <= length_2)) {
	    return true;
	  }else{
	    if (length_2 < 12*units::cm && dis < 2 *units::cm) // very close and small
	      return true;
	  }
	  flag_extend = true;
	}

	double theta2_1 = (dir3_rot.Theta()-3.1415926/2.)/3.1415926*180.;
	double theta2_2 = (dir4_rot.Theta()-3.1415926/2.)/3.1415926*180.;
	double dphi2 = fabs(3.1415926 - fabs(dir3_rot.Phi()-dir4_rot.Phi()))/3.1415926*180.;
	
	if (fabs(theta2_1) < para_angle_cut_1 && fabs(theta2_2) < para_angle_cut_1 && fabs(theta2_1+theta2_2) < para_angle_cut_1/2. && dphi2 < para_angle_cut){
	  if (dphi3 < para_angle_cut || dphi3_1 < para_angle_cut || (dphi3_2 < para_angle_cut && length_1 >= length_2) || (dphi3_3 < para_angle_cut && length_1 <= length_2)) {
	    return true;
	  }else{
	    if (length_1 < 12*units::cm && dis < 2 *units::cm) // very close and small
	      return true;
	  }
	  flag_extend = true;
	}
      }

     
      
     // if (cluster1->get_cluster_id()==133 || cluster2->get_cluster_id()==133) 
     // 	    std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << dis/units::cm << " " << flag_para << " "<< flag_para_U << " " << flag_para_V << " " << flag_prolonged_U << " " << flag_prolonged_V << " " << flag_extend << " " << std::endl;

            
      if (flag_extend && flag_enable_extend ){
 	// when to extend???

	 
	/* if (cluster1->get_cluster_id()==12 || cluster2->get_cluster_id()==12 || */
	/*     cluster1->get_cluster_id()==13 || cluster2->get_cluster_id()==13 )  */
	/*   std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << dis/units::cm << " " << flag_para << " "<< flag_para_U << " " << flag_para_V << " " << flag_prolonged_U << " " << flag_prolonged_V << " " << flag_extend << " " << std::endl; */
	/* return false; */
	
	if ((flag_para && (flag_para_U || flag_para_V || dis < 10*units::cm)) ||
	    (!flag_para && (flag_prolonged_U || flag_prolonged_V || dis < 10*units::cm))){
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
	  //	std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << min_dis/units::cm << " " << max_dis/units::cm << " " << length_2/units::cm << std::endl;

	  /* if (cluster1->get_cluster_id()==12 || cluster2->get_cluster_id()==12 || */
	  /*   cluster1->get_cluster_id()==13 || cluster2->get_cluster_id()==13 )  */
	  /* std::cout << max_dis - min_dis << std::endl; */
	  
	  if ((max_dis - min_dis)>2.5*units::cm) return true;
	  
	  // look at the other side (repeat) 
	  // cluster2_ave_pos, dir2
	  min_dis = 1e9;
	  max_dis = -1e9;
	  for (int i=-5;i!=6;i++){
	    test_point.x = cluster2_ave_pos.x - dir2.X() * (ave_dis +i*2*units::cm);
	    test_point.y = cluster2_ave_pos.y - dir2.Y() * (ave_dis +i*2*units::cm);
	    test_point.z = cluster2_ave_pos.z - dir2.Z() * (ave_dis +i*2*units::cm);
	    
	    std::pair<SlimMergeGeomCell*,Point> temp_results = cluster1->get_closest_point_mcell(test_point);
	    //reuse this 
	    Point test_point1 = temp_results.second;
	    if (sqrt(pow(test_point1.x-test_point.x,2)+pow(test_point1.y-test_point.y,2)+pow(test_point1.z-test_point.z,2))<1.5*units::cm){
	      double temp_dis = (test_point1.x - cluster2_ave_pos.x)*dir1.X() + (test_point1.y - cluster2_ave_pos.y)*dir1.Y() + (test_point1.z - cluster2_ave_pos.z)*dir1.Z();
	      temp_dis *=-1;
	      if (temp_dis < min_dis) min_dis = temp_dis;
	      if (temp_dis > max_dis) max_dis = temp_dis;
	    }
	  }
	  //std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << min_dis/units::cm << " " << max_dis/units::cm << " " << length_1/units::cm << std::endl;

	  /* if (cluster1->get_cluster_id()==12 || cluster2->get_cluster_id()==12 || */
	  /*   cluster1->get_cluster_id()==13 || cluster2->get_cluster_id()==13 )  */
	  /* std::cout << max_dis - min_dis << std::endl; */
	  if ((max_dis - min_dis)>2.5*units::cm) return true;
	  
	  // both side simutaneously? leave it for futhre
	  
	}
      }
    }
  } // dis cut
  return false;
}
