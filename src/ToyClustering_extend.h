void WireCell2dToy::Clustering_extend(WireCell::PR3DClusterSelection& live_clusters, std::map<PR3DCluster*,double>& cluster_length_map, std::set<WireCell::PR3DCluster*>& cluster_connected_dead, int flag,  double length_cut, int num_try, double length_2_cut){
  
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


  TVector3 drift_dir(1,0,0);
  // pronlonged case for U 3 and V 4 ...
  TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
  TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
  TVector3 W_dir(0,1,0);

  std::set<PR3DCluster*> used_clusters;
  std::set<std::pair<PR3DCluster*, PR3DCluster*>> to_be_merged_pairs;
  
  for (size_t i=0;i!=live_clusters.size();i++){
    PR3DCluster* cluster_1 = live_clusters.at(i);
    cluster_1->Create_point_cloud();

    if (cluster_length_map[cluster_1] > 40*units::cm + num_try * 10*units::cm){
      Point highest_p, lowest_p, earliest_p, latest_p;
      bool flag_para = false;
      bool flag_prol = false;

      
      
      if (flag==1){// prolong case ... 

	std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> el_ps = cluster_1->get_earliest_latest_wcps();
	// find earliest point
	earliest_p.x = el_ps.first.x;
	earliest_p.y = el_ps.first.y;
	earliest_p.z = el_ps.first.z;
	
	TVector3 dir_earlp = cluster_1->VHoughTrans(earliest_p,60*units::cm);
	
	TVector3 tempV5,tempV1;
	tempV1.SetXYZ(0,dir_earlp.Y(),dir_earlp.Z());
	double angle1 = tempV1.Angle(U_dir);
	tempV5.SetXYZ(fabs(dir_earlp.X()),sqrt(pow(dir_earlp.Y(),2)+pow(dir_earlp.Z(),2))*sin(angle1),0);
	angle1 = tempV5.Angle(drift_dir);

	double angle2 = tempV1.Angle(V_dir);
	tempV5.SetXYZ(fabs(dir_earlp.X()),sqrt(pow(dir_earlp.Y(),2)+pow(dir_earlp.Z(),2))*sin(angle2),0);
	angle2 = tempV5.Angle(drift_dir);

	double angle3 = tempV1.Angle(W_dir);
	tempV5.SetXYZ(fabs(dir_earlp.X()),sqrt(pow(dir_earlp.Y(),2)+pow(dir_earlp.Z(),2))*sin(angle3),0);
	angle3 = tempV5.Angle(drift_dir);
	
	
	// find latest point 
	latest_p.x = el_ps.second.x;
	latest_p.y = el_ps.second.y;
	latest_p.z = el_ps.second.z;

	TVector3 dir_latep = cluster_1->VHoughTrans(latest_p, 60*units::cm);
	tempV1.SetXYZ(0,dir_latep.Y(),dir_latep.Z());
	double angle4 = tempV1.Angle(U_dir);
	tempV5.SetXYZ(fabs(dir_latep.X()),sqrt(pow(dir_latep.Y(),2)+pow(dir_latep.Z(),2))*sin(angle4),0);
	angle4 = tempV5.Angle(drift_dir);

	double angle5 = tempV1.Angle(V_dir);
	tempV5.SetXYZ(fabs(dir_latep.X()),sqrt(pow(dir_latep.Y(),2)+pow(dir_latep.Z(),2))*sin(angle5),0);
	angle5 = tempV5.Angle(drift_dir);

	double angle6 = tempV1.Angle(W_dir);
	tempV5.SetXYZ(fabs(dir_latep.X()),sqrt(pow(dir_latep.Y(),2)+pow(dir_latep.Z(),2))*sin(angle6),0);
	angle6 = tempV5.Angle(drift_dir);


	/* std::cout << cluster_1->get_cluster_id() << " " << angle1/3.1415926*180. << " " << angle2/3.1415926*180. << " " << angle3/3.1415926*180. << " " */
	/*   << angle4/3.1415926*180. << " " << angle5/3.1415926*180. << " " << angle6/3.1415926*180. << " " */
	/* 	  <<std::endl; */
	 

	
	if (angle1 <5./180.*3.1415926 || angle2 < 5./180.*3.1415926 || angle3 < 5./180.*3.1415926){
	  flag_prol = true;

	  for (size_t j=0;j!=live_clusters.size();j++){
	    PR3DCluster* cluster_2 = live_clusters.at(j);
	    if (cluster_2==cluster_1) continue;
	    if (Clustering_4th_prol(cluster_1,cluster_2,cluster_length_map[cluster_2],earliest_p,dir_earlp,length_cut))
	      to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	  }
	}
	
   

	if (angle4<5./180.*3.1415926 || angle5 < 5./180.*3.1415926 || angle6 < 5./180.*3.1415926){

	  flag_prol = true;
	  for (size_t j=0;j!=live_clusters.size();j++){
	    PR3DCluster* cluster_2 = live_clusters.at(j);
	    if (cluster_2==cluster_1) continue;
	    if (Clustering_4th_prol(cluster_1,cluster_2,cluster_length_map[cluster_2],latest_p,dir_latep,length_cut))
	      to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	  }
	}
      }else if (flag==2){ // parallel case
	
	std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> hl_ps = cluster_1->get_highest_lowest_wcps();
	// find highest point
	highest_p.x = hl_ps.first.x;
	highest_p.y = hl_ps.first.y;
	highest_p.z = hl_ps.first.z;

	highest_p = cluster_1->calc_ave_pos(highest_p,100);
	
	TVector3 dir_highp = cluster_1->VHoughTrans(highest_p,100*units::cm);
	//TVector3 dir_highp_1 = cluster_1->VHoughTrans(highest_p,25*units::cm);
	// find lowest point
	lowest_p.x = hl_ps.second.x;
	lowest_p.y = hl_ps.second.y;
	lowest_p.z = hl_ps.second.z;

	lowest_p = cluster_1->calc_ave_pos(lowest_p,100);
	
	TVector3 dir_lowp = cluster_1->VHoughTrans(lowest_p, 100*units::cm);
	//TVector3 dir_lowp_1 = cluster_1->VHoughTrans(lowest_p, 25*units::cm);

	//	std::cout << cluster_1->get_cluster_id() << " " << cluster_length_map[cluster_1]/units::cm << " " << fabs(dir_highp.Angle(drift_dir)-3.1415926/2.)/3.1415926*180. << " " << fabs(dir_lowp.Angle(drift_dir)-3.1415926/2.)/3.1415926*180. << " " << fabs(dir_highp_1.Angle(drift_dir)-3.1415926/2.)/3.1415926*180. << " " << fabs(dir_lowp_1.Angle(drift_dir)-3.1415926/2.)/3.1415926*180. <<std::endl; 
	
	 if (fabs(dir_highp.Angle(drift_dir)-3.1415926/2.)<5/180.*3.1415926){ 
	   flag_para = true; 

	   for (size_t j=0;j!=live_clusters.size();j++){
	     PR3DCluster* cluster_2 = live_clusters.at(j);
	     if (cluster_2==cluster_1) continue;
	     if (Clustering_4th_para(cluster_1,cluster_2,cluster_length_map[cluster_1],cluster_length_map[cluster_2],highest_p,dir_highp,length_cut))
	       to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	   }
	 }/* else if (fabs(dir_highp_1.Angle(drift_dir)-3.1415926/2.)<5/180.*3.1415926){ */
	 /*   for (size_t j=0;j!=live_clusters.size();j++){ */
	 /*     PR3DCluster* cluster_2 = live_clusters.at(j); */
	 /*     if (cluster_2 == cluster_1) continue; */
	 /*     if (Clustering_4th_para(cluster_1,cluster_2,cluster_length_map[cluster_2],highest_p,dir_highp_1,length_cut)) */
	 /*       to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2)); */
	 /*   } */
	 /* } */
	
	 if (fabs(dir_lowp.Angle(drift_dir)-3.1415926/2.)<5/180.*3.1415926 ){ 
	   flag_para = true; 

	   for (size_t j=0;j!=live_clusters.size();j++){
	     PR3DCluster* cluster_2 = live_clusters.at(j);
	     if (cluster_2==cluster_1) continue;
	     if (Clustering_4th_para(cluster_1,cluster_2,cluster_length_map[cluster_1],cluster_length_map[cluster_2],lowest_p,dir_lowp,length_cut))
	       to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	   }
	   
	 }/* else if (fabs(dir_lowp_1.Angle(drift_dir)-3.1415926/2.)<5/180.*3.1415926){ */
	 /*    for (size_t j=0;j!=live_clusters.size();j++){ */
	 /*     PR3DCluster* cluster_2 = live_clusters.at(j); */
	 /*     if (cluster_2 == cluster_1) continue; */
	 /*     if (Clustering_4th_para(cluster_1,cluster_2,cluster_length_map[cluster_2],lowest_p,dir_lowp_1,length_cut)) */
	 /*       to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2)); */
	 /*   } */
	 /* } */
	 
      }else if (flag==3){ // regular case ...
	std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> hl_ps = cluster_1->get_highest_lowest_wcps();
	std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> el_ps = cluster_1->get_earliest_latest_wcps();
	Point first_p, second_p;
	
	if (pow(hl_ps.first.x-hl_ps.second.x,2)+pow(hl_ps.first.y-hl_ps.second.y,2)+pow(hl_ps.first.z-hl_ps.second.z,2) > pow(el_ps.first.x-el_ps.second.x,2)+pow(el_ps.first.y-el_ps.second.y,2)+pow(el_ps.first.z-el_ps.second.z,2)){
	  first_p.x = hl_ps.first.x;
	  first_p.y = hl_ps.first.y;
	  first_p.z = hl_ps.first.z;

	  second_p.x = hl_ps.second.x;
	  second_p.y = hl_ps.second.y;
	  second_p.z = hl_ps.second.z;
	}else{
	  first_p.x = el_ps.first.x;
	  first_p.y = el_ps.first.y;
	  first_p.z = el_ps.first.z;

	  second_p.x = el_ps.second.x;
	  second_p.y = el_ps.second.y;
	  second_p.z = el_ps.second.z;
	}
	

	for (size_t j=0;j!=live_clusters.size();j++){
	  PR3DCluster* cluster_2 = live_clusters.at(j);
	  if (cluster_2==cluster_1) continue;

	  //  if (cluster_length_map[cluster_2] <40*units::cm) continue;
	  
	  if (Clustering_4th_reg(cluster_1,cluster_2,cluster_length_map[cluster_1],cluster_length_map[cluster_2],first_p,length_cut)){
	    to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	  }else if (Clustering_4th_reg(cluster_1,cluster_2,cluster_length_map[cluster_1],cluster_length_map[cluster_2],second_p,length_cut)){
	    to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	  }
		     	  
	}
      }else if (flag==4){
	if (cluster_connected_dead.find(cluster_1)!=cluster_connected_dead.end()){
	  used_clusters.insert(cluster_1);
	  for (size_t j=0;j!=live_clusters.size();j++){
	    PR3DCluster* cluster_2 = live_clusters.at(j);
	    if (cluster_length_map[cluster_2] < length_2_cut) continue;
	    if (used_clusters.find(cluster_2)!=used_clusters.end()) continue;
	    if (Clustering_4th_dead(cluster_1,cluster_2,cluster_length_map[cluster_1],cluster_length_map[cluster_2],length_cut)){
	      to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	    }
	  }
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


bool WireCell2dToy::Clustering_4th_dead(WireCell::PR3DCluster *cluster_1, WireCell::PR3DCluster *cluster_2, double length_1, double length_2, double length_cut){
  cluster_2->Create_point_cloud();

  TPCParams& mp = Singleton<TPCParams>::Instance();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  TVector3 drift_dir(1,0,0);

  
  SlimMergeGeomCell *mcell1 = 0;
  SlimMergeGeomCell *mcell2=0;
  Point p1;
  Point p2;
  double dis = Find_Closeset_Points(cluster_1, cluster_2, length_1, length_2, length_cut, mcell1, mcell2, p1,p2);

  //  if (length_2 < 10*units::cm && dis >20*units::cm) return false;
  
  if (dis < length_cut || (length_2 > 50*units::cm && dis < 80*units::cm)){

    Point cluster1_ave_pos_save;
    Point cluster2_ave_pos_save;
    TVector3 dir1_save;
    TVector3 dir3_save;
    
    for (int i=0;i!=3;i++){
      Point cluster1_ave_pos; 
      Point cluster2_ave_pos; 
      
      TVector3 dir1; 
      TVector3 dir3; 
      TVector3 dir2;
      
     

      
      
      if (i==0){
	cluster1_ave_pos = cluster_1->calc_ave_pos(p1,5*units::cm);
	cluster1_ave_pos_save = cluster1_ave_pos;
	cluster2_ave_pos = cluster_2->calc_ave_pos(p2,5*units::cm);
	cluster2_ave_pos_save = cluster2_ave_pos;
	
	dir1 = cluster_1->VHoughTrans(cluster1_ave_pos,80*units::cm);
	dir1_save = dir1;
	dir3 = cluster_2->VHoughTrans(cluster2_ave_pos,80*units::cm);
	dir3_save = dir3;
	dir2.SetXYZ(cluster2_ave_pos.x - cluster1_ave_pos.x+1e-9, cluster2_ave_pos.y - cluster1_ave_pos.y+1e-9, cluster2_ave_pos.z - cluster1_ave_pos.z+1e-9); // 2-1
      }else if (i==1 ){
	if (length_2 >= 15*units::cm){
	  cluster1_ave_pos = cluster1_ave_pos_save;//cluster_1->calc_ave_pos(p1,5*units::cm);
	  dir1 = dir1_save;//cluster_1->VHoughTrans(cluster1_ave_pos,80*units::cm);
	  
	  TVector3 dir_test(dir1);
	  dir_test.SetMag(1);
	  dir_test *= -1;
	  
	  std::pair<Point, double> temp_results = cluster_2->get_closest_point_along_vec(cluster1_ave_pos, dir_test, dis*2, 5*units::cm, 15, 10*units::cm);
	  
	  /* if ( length_1>200*units::cm)  */
	  /*   std::cout << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << " " << length_1/units::cm << " " << length_2/units::cm << " " <<  " a " << temp_results.second/units::cm << " " << cluster1_ave_pos.x/units::cm << " " << cluster1_ave_pos.y/units::cm << " " << cluster1_ave_pos.z/units::cm << " " << temp_results.first.x/units::cm << " " << temp_results.first.y/units::cm << " " << temp_results.first.z/units::cm << " " << dir_test.X() << " " << dir_test.Y() << " " << dir_test.Z() << " " << dis/units::cm << " " << std::endl;  */
	  
	  
	  if (temp_results.second < 100*units::cm){
	    cluster2_ave_pos = cluster_2->calc_ave_pos(temp_results.first,5*units::cm);
	    dir3 = cluster_2->VHoughTrans(cluster2_ave_pos,80*units::cm);
	    dir2.SetXYZ(cluster2_ave_pos.x - cluster1_ave_pos.x+1e-9, cluster2_ave_pos.y - cluster1_ave_pos.y+1e-9, cluster2_ave_pos.z - cluster1_ave_pos.z+1e-9); // 2-1
	    
	  }else{
	    continue;
	  }
	}else{
	  continue;
	}
      }else if (i==2){
	if (length_2 >=15*units::cm){
	  cluster2_ave_pos = cluster2_ave_pos_save;//cluster_2->calc_ave_pos(p2,5*units::cm);
	  dir3 = dir3_save;//cluster_2->VHoughTrans(cluster2_ave_pos,80*units::cm);
	  
	  TVector3 dir_test(dir3);
	  dir_test.SetMag(1);
	  dir_test *= -1;
	  
	  std::pair<Point, double> temp_results = cluster_1->get_closest_point_along_vec(cluster2_ave_pos, dir_test, dis*2, 5*units::cm, 15, 10*units::cm);
	  
	  if (temp_results.second < 100*units::cm){
	    cluster1_ave_pos = cluster_1->calc_ave_pos(temp_results.first,5*units::cm);
	    dir1 = cluster_1->VHoughTrans(cluster1_ave_pos,80*units::cm);
	    dir2.SetXYZ(cluster2_ave_pos.x - cluster1_ave_pos.x+1e-9, cluster2_ave_pos.y - cluster1_ave_pos.y+1e-9, cluster2_ave_pos.z - cluster1_ave_pos.z+1e-9); // 2-1
	    
	  }else{
	    continue;
	  }
	}else{
	  continue;
	}
      }
      
      TVector3 ave_dir(cluster1_ave_pos.x-cluster2_ave_pos.x,cluster1_ave_pos.y-cluster2_ave_pos.y,cluster1_ave_pos.z-cluster2_ave_pos.z);


     

      
      /* if (length_1 > 100*units::cm && i==0 && length_2 > 30*units::cm) */
      /* 	std::cout << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << " " << length_1/units::cm << " " << length_2/units::cm << " " <<  " c " << WireCell2dToy::is_angle_consistent(dir1,dir2,false,15,angle_u,angle_v,angle_w) << " " << WireCell2dToy::is_angle_consistent(dir3,dir2,true,15,angle_u,angle_v,angle_w) << " " << fabs(ave_dir.Angle(drift_dir)-3.1415926/2.)/3.1415926*180. << std::endl; */
      
      // use projection to deal with stuff ...
      if (fabs(ave_dir.Angle(drift_dir)-3.1415926/2.)/3.1415926*180.>7.5){
      	// non-parallel case ...
      	if (WireCell2dToy::is_angle_consistent(dir1,dir2,false,10,angle_u,angle_v,angle_w)){
      	  if (length_2 < 8*units::cm&& WireCell2dToy::is_angle_consistent(dir1,dir2,false,5,angle_u,angle_v,angle_w))
      	    return true;
      	  if (length_2 < 15*units::cm && WireCell2dToy::is_angle_consistent(dir1,dir2,false,7.5,angle_u,angle_v,angle_w))
      	    return true;
      	  if (WireCell2dToy::is_angle_consistent(dir3,dir2,true,10,angle_u,angle_v,angle_w)){
      	    return true;
      	  }
      	}
      }
      
      
      double angle1 = (3.1415926-dir1.Angle(dir2))/3.1415926*180.;
      double angle2 = dir3.Angle(dir2)/3.1415926*180.;
      double angle3 = (3.1415926-dir1.Angle(dir3))/3.1415926*180.;


      /* if ( length_2 > 100*units::cm && length_1>100*units::cm) */
      /* 	  std::cout << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << " " << length_1/units::cm << " " << length_2/units::cm << " " <<  angle1 << " " << angle2 << " " << angle3 << " "  << WireCell2dToy::is_angle_consistent(dir1,dir2,false,10,angle_u,angle_v,angle_w) << " " << WireCell2dToy::is_angle_consistent(dir3,dir2,true,10,angle_u,angle_v,angle_w) << std::endl; */


      
      if (length_2 <=10*units::cm){
	if (angle1 < 15 && (angle2 < 60 || length_2 < 5*units::cm)) return true;
      }else{
	if (angle1 < 15 && angle2 <15 && angle3 < 25)
	  return true;

	double ave_dis = sqrt(pow(cluster1_ave_pos.x-cluster2_ave_pos.x,2) + pow(cluster1_ave_pos.y-cluster2_ave_pos.y,2) + pow(cluster1_ave_pos.z-cluster2_ave_pos.z,2));
	Point test_point;
	double min_dis = 1e9, max_dis = -1e9;

	if (fabs(ave_dir.Angle(drift_dir)-3.1415926/2.)/3.1415926*180. > 7.5  && ave_dis < 30*units::cm){
	  
	  if (i==1){
	    for (int k=-5;k!=10;k++){
	      test_point.x = cluster1_ave_pos.x - dir1.X() * (ave_dis +k*2*units::cm);
	      test_point.y = cluster1_ave_pos.y - dir1.Y() * (ave_dis +k*2*units::cm);
	      test_point.z = cluster1_ave_pos.z - dir1.Z() * (ave_dis +k*2*units::cm);
	      
	      std::pair<SlimMergeGeomCell*,Point> temp_results = cluster_2->get_closest_point_mcell(test_point);
	      //reuse this
	      Point test_point1 = temp_results.second;
	      if (sqrt(pow(test_point1.x-test_point.x,2)+pow(test_point1.y-test_point.y,2)+pow(test_point1.z-test_point.z,2))<1.5*units::cm){
		double temp_dis = (test_point1.x - cluster1_ave_pos.x)*dir1.X() + (test_point1.y - cluster1_ave_pos.y)*dir1.Y() + (test_point1.z - cluster1_ave_pos.z)*dir1.Z();
		temp_dis *=-1;
		if (temp_dis < min_dis) min_dis = temp_dis;
		if (temp_dis > max_dis) max_dis = temp_dis;
	      }
	    }
	    if ((max_dis - min_dis)>2.5*units::cm) return true;
	  }else if (i==2){
	    for (int k=-5;k!=10;k++){
	      test_point.x = cluster2_ave_pos.x - dir3.X() * (ave_dis +k*2*units::cm);
	      test_point.y = cluster2_ave_pos.y - dir3.Y() * (ave_dis +k*2*units::cm);
	      test_point.z = cluster2_ave_pos.z - dir3.Z() * (ave_dis +k*2*units::cm);
	      
	      std::pair<SlimMergeGeomCell*,Point> temp_results = cluster_1->get_closest_point_mcell(test_point);
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
	  }
	}

      }
    }
  }
  
  return false;
}





bool WireCell2dToy::Clustering_4th_reg(WireCell::PR3DCluster *cluster_1, WireCell::PR3DCluster *cluster_2, double length_1, double length_2, WireCell::Point p1, double length_cut){
  cluster_2->Create_point_cloud();
  
  std::pair<SlimMergeGeomCell*,Point> temp_results = cluster_2->get_closest_point_mcell(p1);
  Point p2 = temp_results.second;

  temp_results = cluster_1->get_closest_point_mcell(p2);
  p1 = temp_results.second;
  /* temp_results = cluster_2->get_closest_point_mcell(p1); */
  /* p2 = temp_results.second; */
  
  double dis = sqrt(pow(p1.x-p2.x,2) + pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));

  TPCParams& mp = Singleton<TPCParams>::Instance();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  TVector3 drift_dir(1,0,0);
  
  if (dis < length_cut && (length_2 >= 40*units::cm || dis < 3*units::cm)){

    Point cluster1_ave_pos = cluster_1->calc_ave_pos(p1,5*units::cm);
    Point cluster2_ave_pos = cluster_2->calc_ave_pos(p2,5*units::cm);
    TVector3 dir1 = cluster_1->VHoughTrans(cluster1_ave_pos,30*units::cm);
    TVector3 dir3 = cluster_2->VHoughTrans(cluster2_ave_pos,30*units::cm);
    TVector3 dir2(cluster2_ave_pos.x - cluster1_ave_pos.x,
		  cluster2_ave_pos.y - cluster1_ave_pos.y,
		  cluster2_ave_pos.z - cluster1_ave_pos.z);

    /* if (fabs(dir2.Angle(drift_dir)-3.1415926/2.)/3.1415926*180.>7.5){ */
    /*   // non-parallel case ... */
    /*   if (WireCell2dToy::is_angle_consistent(dir1,dir2,false,7.5,angle_u,angle_v,angle_w)){ */
    /* 	if (WireCell2dToy::is_angle_consistent(dir3,dir2,true,7.5,angle_u,angle_v,angle_w)){ */
    /* 	  return true; */
    /* 	} */
    /*   } */
    /* } */

    
    double ave_dis = sqrt(pow(cluster1_ave_pos.x-cluster2_ave_pos.x,2) + pow(cluster1_ave_pos.y-cluster2_ave_pos.y,2) + pow(cluster1_ave_pos.z-cluster2_ave_pos.z,2));
      Point test_point;
      double min_dis = 1e9, max_dis = -1e9;
    
    if (dir2.Angle(dir1)>3.1415926/2. ){
          
      /* int num_p1 = cluster_1->get_num_points(p1, 10*units::cm); */
      /* int num_p2 = cluster_2->get_num_points(p2, 10*units::cm); */
      
      
      for (int i=-5;i!=10;i++){
	test_point.x = cluster1_ave_pos.x - dir1.X() * (ave_dis +i*2*units::cm);
	test_point.y = cluster1_ave_pos.y - dir1.Y() * (ave_dis +i*2*units::cm);
	test_point.z = cluster1_ave_pos.z - dir1.Z() * (ave_dis +i*2*units::cm);
	
	std::pair<SlimMergeGeomCell*,Point> temp_results = cluster_2->get_closest_point_mcell(test_point);
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
    }
    
    if (dir2.Angle(dir3)<3.1415926/2.){
      
      
      // look at the other side (repeat)
      // cluster2_ave_pos, dir2
      min_dis = 1e9;
      max_dis = -1e9;
      for (int i=-5;i!=10;i++){
	test_point.x = cluster2_ave_pos.x - dir3.X() * (ave_dis +i*2*units::cm);
	test_point.y = cluster2_ave_pos.y - dir3.Y() * (ave_dis +i*2*units::cm);
	test_point.z = cluster2_ave_pos.z - dir3.Z() * (ave_dis +i*2*units::cm);
	
	std::pair<SlimMergeGeomCell*,Point> temp_results = cluster_1->get_closest_point_mcell(test_point);
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
    }
  }else if (dis < 2 * length_cut && length_2 < 40*units::cm){
    TVector3 drift_dir(1,0,0);
    // pronlonged case for U 3 and V 4 ...
    TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
    TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
    
    TVector3 dir2(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z);
    bool flag_para = false, flag_prol =false, flag_reg = false;
    
    double angle1 = fabs(dir2.Angle(drift_dir)-3.1415926/2.)/3.1415926*180.;

     /* if (cluster_1->get_cluster_id()==261) */
     /*   std::cout << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << " " << angle1 << " " << length_2/units::cm << " " << dis/units::cm << " " << std::endl; */
     
    
    if (angle1 < 5 && dis < 2*length_cut || angle1 < 2 ){
      flag_para = true;
    }else if (dis < 2*length_cut){
      TVector3 tempV1(0, p2.y - p1.y, p2.z - p1.z);
      TVector3 tempV5;
      double angle2 = tempV1.Angle(U_dir);
      tempV5.SetXYZ(fabs(p2.x-p1.x),sqrt(pow(p2.y - p1.y,2)+pow(p2.z - p1.z,2))*sin(angle2),0);
      angle2 = tempV5.Angle(drift_dir)/3.1415926*180.;
      
      double angle3 = tempV1.Angle(V_dir);
      tempV5.SetXYZ(fabs(p2.x-p1.x),sqrt(pow(p2.y - p1.y,2)+pow(p2.z - p1.z,2))*sin(angle3),0);
      angle3 = tempV5.Angle(drift_dir)/3.1415926*180.;
      if (angle2<7.5 || angle3 < 7.5)
  	flag_prol = true;

      /* if (cluster_1->get_cluster_id()==448 || cluster_2->get_cluster_id()==85) */
      /* 	std::cout << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << " " << angle2 << " " << angle3 << " " << length_2/units::cm << " " << dis/units::cm << " " << std::endl; */
      
    }

    if (flag_para || flag_prol || flag_reg){
      TVector3 dir1 = cluster_1->VHoughTrans(p1,15*units::cm);
      TVector3 dir3 = cluster_2->VHoughTrans(p2,15*units::cm);
      
      double angle4 = (3.1415926-dir1.Angle(dir2))/3.1415926*180.;
      double angle5 = dir2.Angle(dir3)/3.1415926*180.;
      
      /* if (cluster_1->get_cluster_id()== 11 ) */
      /*   std::cout  << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << " " << angle4 << " " << angle5 << " " << dis/units::cm << " " << length_2/units::cm << std::endl; */

      // return false;
      
      if (flag_para){
  	if (angle4 < 30 && (length_2 < 12*units::cm || angle5 < 60))
  	  return true;
      }else if (flag_prol){
  	if (angle4 < 25 && (length_2 < 15*units::cm || angle5 < 25))
  	  return true;
      }

      if (fabs(dir2.Angle(drift_dir)-3.1415926/2.)/3.1415926*180.>7.5){
	// non-parallel case ... 
	if (WireCell2dToy::is_angle_consistent(dir1,dir2,false,10,angle_u,angle_v,angle_w)){
	  if (length_2 < 8*units::cm && WireCell2dToy::is_angle_consistent(dir1,dir2,false,5,angle_u,angle_v,angle_w)) 
	    return true; 
	  if (WireCell2dToy::is_angle_consistent(dir3,dir2,true,10,angle_u,angle_v,angle_w)){
	    return true;
	  }
	}
      }
      
      
    }
  }


  
  /* if (cluster_1->get_cluster_id()==11 || cluster_2->get_cluster_id()==55) */
  /*   std::cout << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << " " << length_2/units::cm << " " << dis/units::cm << " " << std::endl; */
  
  
  /* if (dis < 5 * length_cut){ */
  /*   // parallel case, prlonged case, additional dead channels in U plane ...  */
  /*   TVector3 drift_dir(1,0,0); */
  /*   // pronlonged case for U 3 and V 4 ... */
  /*   TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926)); */
  /*   TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926)); */
    
  /*   TVector3 dir2(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z); */
  /*   bool flag_para = false, flag_prol =false, flag_reg = false; */
    
  /*   double angle1 = fabs(dir2.Angle(drift_dir)-3.1415926/2.)/3.1415926*180.; */

  /*    /\* if (cluster_1->get_cluster_id()==261) *\/ */
  /*    /\*   std::cout << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << " " << angle1 << " " << length_2/units::cm << " " << dis/units::cm << " " << std::endl; *\/ */
     
    
  /*   if (angle1 < 5 && dis < 2*length_cut || angle1 < 2 ){ */
  /*     flag_para = true; */
  /*   }else if (dis < 2*length_cut){ */
  /*     TVector3 tempV1(0, p2.y - p1.y, p2.z - p1.z); */
  /*     TVector3 tempV5; */
  /*     double angle2 = tempV1.Angle(U_dir); */
  /*     tempV5.SetXYZ(fabs(p2.x-p1.x),sqrt(pow(p2.y - p1.y,2)+pow(p2.z - p1.z,2))*sin(angle2),0); */
  /*     angle2 = tempV5.Angle(drift_dir)/3.1415926*180.; */
      
  /*     double angle3 = tempV1.Angle(V_dir); */
  /*     tempV5.SetXYZ(fabs(p2.x-p1.x),sqrt(pow(p2.y - p1.y,2)+pow(p2.z - p1.z,2))*sin(angle3),0); */
  /*     angle3 = tempV5.Angle(drift_dir)/3.1415926*180.; */
  /*     if (angle2<7.5 || angle3 < 7.5) */
  /* 	flag_prol = true; */

  /*     /\* if (cluster_1->get_cluster_id()==448 || cluster_2->get_cluster_id()==85) *\/ */
  /*     /\* 	std::cout << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << " " << angle2 << " " << angle3 << " " << length_2/units::cm << " " << dis/units::cm << " " << std::endl; *\/ */
      
  /*   } */
  /*   if (dis < length_cut) */
  /*     flag_reg = true; */

  /*   if (flag_para || flag_prol || flag_reg){ */
  /*     TVector3 dir1 = cluster_1->VHoughTrans(p1,15*units::cm); */
  /*     TVector3 dir3 = cluster_2->VHoughTrans(p2,15*units::cm); */
      
  /*     double angle4 = (3.1415926-dir1.Angle(dir2))/3.1415926*180.; */
  /*     double angle5 = dir2.Angle(dir3)/3.1415926*180.; */
      
  /*     /\* if (cluster_1->get_cluster_id()== 11 ) *\/ */
  /*     /\*   std::cout  << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << " " << angle4 << " " << angle5 << " " << dis/units::cm << " " << length_2/units::cm << std::endl; *\/ */

  /*     // return false; */
      
  /*     if (flag_para){ */
  /* 	if (angle4 < 30 && (length_2 < 12*units::cm || angle5 < 60)) */
  /* 	  return true; */
  /*     }else if (flag_prol){ */
  /* 	if (angle4 < 25 && (length_2 < 15*units::cm || angle5 < 25)) */
  /* 	  return true; */
  /*     }else if (flag_reg){ */
  /* 	if ((angle4 < 20 || angle4 > 175 )&& (length_2 < 12*units::cm || angle5 < 20) || */
  /* 	    (dis< 10*units::cm && length_2 < 20*units::cm &&angle4 < 65 && angle5 < 90)) // cover the case missed by the live dead cases */
  /* 	  return true; */
  /*     } */
  /*   } */
  /* } */
  
  
  return false;
}



bool WireCell2dToy::Clustering_4th_para(WireCell::PR3DCluster *cluster_1, WireCell::PR3DCluster *cluster_2, double length_1, double length_2, WireCell::Point& earliest_p, TVector3& dir_earlp, double length_cut){
  cluster_2->Create_point_cloud();
 
  std::pair<SlimMergeGeomCell*,Point> temp_results = cluster_2->get_closest_point_mcell(earliest_p);
  Point p2 = temp_results.second;

  double dis = sqrt(pow(earliest_p.x-p2.x,2) + pow(earliest_p.y-p2.y,2)+pow(earliest_p.z-p2.z,2));
  
  
  if (dis < length_cut ){
    
    Point test_point; 
    double min_dis = 1e9, max_dis = -1e9;
    
    for (int i=-5;i!=10;i++){ 
 	test_point.x = earliest_p.x - dir_earlp.X() * (dis +i*2*units::cm); 
 	test_point.y = earliest_p.y - dir_earlp.Y() * (dis +i*2*units::cm); 
 	test_point.z = earliest_p.z - dir_earlp.Z() * (dis +i*2*units::cm); 
	
 	std::pair<SlimMergeGeomCell*,Point> temp_results = cluster_2->get_closest_point_mcell(test_point); 
 	//reuse this 
 	Point test_point1 = temp_results.second;

	/* std::cout << test_point.x/units::cm << " " << test_point.y/units::cm << " " << test_point.z/units::cm << " " << test_point1.x/units::cm << " " << test_point1.y/units::cm << " " << test_point1.z/units::cm << " " << sqrt(pow(test_point1.x-test_point.x,2)+pow(test_point1.y-test_point.y,2)+pow(test_point1.z-test_point.z,2))/units::cm << std::endl; */
	
 	if (sqrt(pow(test_point1.x-test_point.x,2)+pow(test_point1.y-test_point.y,2)+pow(test_point1.z-test_point.z,2))<1.5*units::cm){ 
 	  double temp_dis = (test_point1.x - earliest_p.x)*dir_earlp.X() + (test_point1.y - earliest_p.y)*dir_earlp.Y() + (test_point1.z - earliest_p.z)*dir_earlp.Z(); 
 	  temp_dis *=-1; 
 	  if (temp_dis < min_dis) min_dis = temp_dis; 
 	  if (temp_dis > max_dis) max_dis = temp_dis; 
 	} 
    }

    
    /* std::cout << cluster_1->get_cluster_id() << " " << length_2/units::cm << " " << max_dis/units::cm << " " << min_dis/units::cm << " " << dis/units::cm << " " << dir_earlp.X() << " " << dir_earlp.Y() << " " << dir_earlp.Z() << std::endl;  */

    
    if ((max_dis - min_dis)>2.5*units::cm) return true; 
    
  }
  return false;
}





bool WireCell2dToy::Clustering_4th_prol(WireCell::PR3DCluster *cluster_1, PR3DCluster *cluster_2, double length_2, Point& earliest_p, TVector3& dir_earlp, double length_cut){
  cluster_2->Create_point_cloud();
 
  std::pair<SlimMergeGeomCell*,Point> temp_results = cluster_2->get_closest_point_mcell(earliest_p);
  Point p2 = temp_results.second;
  double dis = sqrt(pow(earliest_p.x-p2.x,2) + pow(earliest_p.y-p2.y,2)+pow(earliest_p.z-p2.z,2));
  if (dis < length_cut){
    TVector3 dir_bp(p2.x-earliest_p.x,p2.y-earliest_p.y,p2.z-earliest_p.z);
    double angle_diff = (3.1415926-dir_bp.Angle( dir_earlp))/3.1415926*180.;
    if ( (angle_diff < 3 || angle_diff>177 || 
	  dis * sin(angle_diff/180.*3.1415926) < 6*units::cm)){
      if (length_2<10*units::cm){
	return true;
      }else{
	TVector3 dir = cluster_2->VHoughTrans(p2,60*units::cm);
	if ((3.14151926-dir.Angle(dir_earlp))/3.1415926*180.<5. ||
	    dir.Angle(dir_earlp)/3.1415926*180.<5.)
	  return true;
      }
    }
    
  }
    
  return false;
}






/* bool WireCell2dToy::Clustering_4th_round(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double length_1, double length_2, double length_cut){ */

/*   cluster1->Create_point_cloud(); */
/*   cluster2->Create_point_cloud(); */
  
/*   // pick any point and merged cell in cluster1, */
/*   SlimMergeGeomCell *prev_mcell1 = 0; */
/*   SlimMergeGeomCell *prev_mcell2 = 0; */
/*   SlimMergeGeomCell *mcell1 = cluster1->get_mcells().at(0); */
/*   Point p1 = mcell1->center(); */
/*   SlimMergeGeomCell *mcell2=0; */
/*   Point p2; */

/*   while(mcell1!=prev_mcell1 || mcell2!=prev_mcell2){ */
/*     prev_mcell1 = mcell1; */
/*     prev_mcell2 = mcell2; */
    
/*     // find the closest point and merged cell in cluster2 */
/*     std::pair<SlimMergeGeomCell*,Point> temp_results = cluster2->get_closest_point_mcell(p1); */
/*     p2 = temp_results.second; */
/*     mcell2 = temp_results.first; */
/*     // find the closest point and merged cell in cluster1 */
/*     temp_results = cluster1->get_closest_point_mcell(p2); */
/*     p1 = temp_results.second; */
/*     mcell1 = temp_results.first; */
/*   } */
  
/*   //std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2))/units::cm<< std::endl; */
  
/*   double dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2)); */

/*   // if (cluster1->get_cluster_id()==7){ */
/*   //   std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << dis/units::cm << " " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << */
/*   //     " " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << */
/*   //     std::endl; */
/*   // } */
  
  
/*   if (dis < length_cut){ */
/*     Point cluster1_ave_pos = cluster1->calc_ave_pos(p1,5*units::cm); */
/*     Point cluster2_ave_pos = cluster2->calc_ave_pos(p2,5*units::cm); */

/*     TVector3 dir1 = cluster1->VHoughTrans(cluster1_ave_pos,30*units::cm); */
/*     TVector3 dir2 = cluster2->VHoughTrans(cluster2_ave_pos,30*units::cm); */

/*     bool flag_para = false; */
/*     bool flag_prolong_U = false; */
/*     bool flag_prolong_V = false; */
/*     bool flag_prolong_W= false; */
/*     bool flag_para_U = false; */
/*     bool flag_para_V = false; */
/*     bool flag_para_W = false; */
    
/*     // pronlonged case for U 3 and V 4 ... */
/*     TVector3 drift_dir(1,0,0); */
/*     TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926)); */
/*     TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926)); */
/*     TVector3 W_dir(0,1,0); */
    
/*     TVector3 tempV1(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z); */
/*     TVector3 tempV2(cluster2_ave_pos.x - cluster1_ave_pos.x, cluster2_ave_pos.y - cluster1_ave_pos.y, cluster2_ave_pos.z - cluster1_ave_pos.z); */
      
/*     double angle1 = tempV1.Angle(drift_dir); */
/*     double angle2 = tempV2.Angle(drift_dir); */

/*     if (fabs(angle1-3.1415926/2.)<15/180.*3.1415926 && */
/*         fabs(angle2-3.1415926/2.)<10/180.*3.1415926){ */
/*       flag_para = true; */

/*       double angle3 = tempV1.Angle(U_dir); */
/*       double angle4 = tempV1.Angle(V_dir); */
/*       double angle5 = tempV1.Angle(W_dir); */

/*       double angle3_1 = tempV2.Angle(U_dir); */
/*       double angle4_1 = tempV2.Angle(V_dir); */
/*       double angle5_1 = tempV2.Angle(W_dir); */

/*       if (fabs(angle3-3.1415926/2.)<7.5/180.*3.1415926 || fabs(angle3_1-3.1415926/2.)<7.5/180.*3.1415926 ) */
/* 	flag_para_U = true; */

/*       if (fabs(angle4-3.1415926/2.)<7.5/180.*3.1415926 || fabs(angle4_1-3.1415926/2.)<7.5/180.*3.1415926 ) */
/* 	flag_para_V = true; */

/*       if (fabs(angle5-3.1415926/2.)<7.5/180.*3.1415926 || fabs(angle5_1-3.1415926/2.)<7.5/180.*3.1415926 ) */
/* 	flag_para_W = true; */
      
/*     } */

/*     { */
/*       TVector3 tempV3(0, p2.y - p1.y, p2.z - p1.z); */
/*       TVector3 tempV4(0, cluster2_ave_pos.y - cluster1_ave_pos.y, cluster2_ave_pos.z - cluster1_ave_pos.z); */

/*       TVector3 tempV5; */

/*       double angle6 = tempV3.Angle(U_dir); */
/*       tempV5.SetXYZ(fabs(p2.x-p1.x),sin(angle6),0); */
/*       angle6 = tempV5.Angle(drift_dir); */
				  
/*       double angle7 = tempV3.Angle(V_dir); */
/*       tempV5.SetXYZ(fabs(p2.x-p1.x),sin(angle7),0); */
/*       angle7 = tempV5.Angle(drift_dir); */

/*       double angle8 = tempV3.Angle(W_dir); */
/*       tempV5.SetXYZ(fabs(p2.x-p1.x),sin(angle8),0); */
/*       angle8 = tempV5.Angle(drift_dir); */
      
/*       double angle6_1 = tempV4.Angle(U_dir); */
/*       tempV5.SetXYZ(fabs(cluster2_ave_pos.x-cluster1_ave_pos.x),sin(angle6_1),0); */
/*       angle6_1 = tempV5.Angle(drift_dir); */
/*       double angle7_1 = tempV4.Angle(V_dir); */
/*       tempV5.SetXYZ(fabs(cluster2_ave_pos.x-cluster1_ave_pos.x),sin(angle7_1),0); */
/*       angle7_1 = tempV5.Angle(drift_dir); */
      
/*       double angle8_1 = tempV4.Angle(W_dir); */
/*       tempV5.SetXYZ(fabs(cluster2_ave_pos.x-cluster1_ave_pos.x),sin(angle8_1),0); */
/*       angle8_1 = tempV5.Angle(drift_dir); */
      
/*       if (angle6<7.5/180.*3.1415926  || */
/* 	  angle6_1<7.5/180.*3.1415926 ) */
/* 	flag_prolong_U = true; */

/*       if (angle7<7.5/180.*3.1415926  || */
/* 	  angle7_1<7.5/180.*3.1415926 ) */
/* 	flag_prolong_V = true; */

/*       if (angle8<7.5/180.*3.1415926  || */
/* 	  angle8_1<7.5/180.*3.1415926 ) */
/* 	flag_prolong_W = true; */
/*     } */


    
    
    
/*     if (flag_para_U || flag_para_V || flag_para_W || */
/* 	flag_prolong_U || flag_prolong_V || flag_prolong_W){ */
/*       TVector3 dir1_1 = cluster1->VHoughTrans(p1,30*units::cm); */
/*       TVector3 dir2_1 = cluster2->VHoughTrans(p2,30*units::cm); */
/*       double ave_dis = sqrt(pow(cluster1_ave_pos.x-cluster2_ave_pos.x,2) + pow(cluster1_ave_pos.y-cluster2_ave_pos.y,2) + pow(cluster1_ave_pos.z-cluster2_ave_pos.z,2)); */
      
/*       Point test_point; */
/*       double min_dis = 1e9, max_dis = -1e9; */
/*       for (int i=-5;i!=6;i++){ */
/* 	test_point.x = cluster1_ave_pos.x - dir1.X() * (ave_dis +i*2*units::cm); */
/* 	test_point.y = cluster1_ave_pos.y - dir1.Y() * (ave_dis +i*2*units::cm); */
/* 	test_point.z = cluster1_ave_pos.z - dir1.Z() * (ave_dis +i*2*units::cm); */
	
/* 	std::pair<SlimMergeGeomCell*,Point> temp_results = cluster2->get_closest_point_mcell(test_point); */
/* 	//reuse this  */
/* 	Point test_point1 = temp_results.second; */
/* 	if (sqrt(pow(test_point1.x-test_point.x,2)+pow(test_point1.y-test_point.y,2)+pow(test_point1.z-test_point.z,2))<1.5*units::cm){ */
/* 	  double temp_dis = (test_point1.x - cluster1_ave_pos.x)*dir1.X() + (test_point1.y - cluster1_ave_pos.y)*dir1.Y() + (test_point1.z - cluster1_ave_pos.z)*dir1.Z(); */
/* 	  temp_dis *=-1; */
/* 	  if (temp_dis < min_dis) min_dis = temp_dis; */
/* 	  if (temp_dis > max_dis) max_dis = temp_dis; */
/* 	} */
/*       } */
/*       if ((max_dis - min_dis)>2.5*units::cm) return true; */
      
      
/*       double ave_dis1 = sqrt(pow(p1.x-cluster2_ave_pos.x,2) + pow(p1.y-cluster2_ave_pos.y,2) + pow(p1.z-cluster2_ave_pos.z,2)); */
/*       min_dis = 1e9; */
/*       max_dis = -1e9; */
/*       for (int i=-5;i!=6;i++){ */
/* 	test_point.x = p1.x - dir1_1.X() * (ave_dis1 +i*2*units::cm); */
/* 	test_point.y = p1.y - dir1_1.Y() * (ave_dis1 +i*2*units::cm); */
/* 	test_point.z = p1.z - dir1_1.Z() * (ave_dis1 +i*2*units::cm); */
	
/* 	std::pair<SlimMergeGeomCell*,Point> temp_results = cluster2->get_closest_point_mcell(test_point); */
/* 	//reuse this  */
/* 	Point test_point1 = temp_results.second; */
/* 	if (sqrt(pow(test_point1.x-test_point.x,2)+pow(test_point1.y-test_point.y,2)+pow(test_point1.z-test_point.z,2))<1.5*units::cm){ */
/* 	  double temp_dis = (test_point1.x - p1.x)*dir1_1.X() + (test_point1.y - p1.y)*dir1_1.Y() + (test_point1.z - p1.z)*dir1_1.Z(); */
/* 	  temp_dis *=-1; */
/* 	  if (temp_dis < min_dis) min_dis = temp_dis; */
/* 	  if (temp_dis > max_dis) max_dis = temp_dis; */
/* 	} */
/*       } */
/*       if ((max_dis - min_dis)>2.5*units::cm) return true; */
      
      
/*       min_dis = 1e9; */
/*       max_dis = -1e9; */
/*       for (int i=-5;i!=6;i++){ */
/* 	test_point.x = cluster2_ave_pos.x - dir2.X() * (ave_dis +i*2*units::cm); */
/* 	test_point.y = cluster2_ave_pos.y - dir2.Y() * (ave_dis +i*2*units::cm); */
/* 	test_point.z = cluster2_ave_pos.z - dir2.Z() * (ave_dis +i*2*units::cm); */
	
/* 	std::pair<SlimMergeGeomCell*,Point> temp_results = cluster1->get_closest_point_mcell(test_point); */
/* 	//reuse this  */
/* 	Point test_point1 = temp_results.second; */
/* 	if (sqrt(pow(test_point1.x-test_point.x,2)+pow(test_point1.y-test_point.y,2)+pow(test_point1.z-test_point.z,2))<1.5*units::cm){ */
/* 	  double temp_dis = (test_point1.x - cluster2_ave_pos.x)*dir1.X() + (test_point1.y - cluster2_ave_pos.y)*dir1.Y() + (test_point1.z - cluster2_ave_pos.z)*dir1.Z(); */
/* 	  temp_dis *=-1; */
/* 	  if (temp_dis < min_dis) min_dis = temp_dis; */
/* 	  if (temp_dis > max_dis) max_dis = temp_dis; */
/* 	} */
/*       } */
/*       //std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << min_dis/units::cm << " " << max_dis/units::cm << " " << length_1/units::cm << std::endl; */
/*       if ((max_dis - min_dis)>2.5*units::cm) return true; */
      
      
/*       min_dis = 1e9; */
/*       max_dis = -1e9; */
/*       double ave_dis2 = sqrt(pow(cluster1_ave_pos.x-p2.x,2) + pow(cluster1_ave_pos.y-p2.y,2) + pow(cluster1_ave_pos.z-p2.z,2)); */
/*       for (int i=-5;i!=6;i++){ */
/* 	test_point.x = p2.x - dir2_1.X() * (ave_dis2 +i*2*units::cm); */
/* 	test_point.y = p2.y - dir2_1.Y() * (ave_dis2 +i*2*units::cm); */
/* 	test_point.z = p2.z - dir2_1.Z() * (ave_dis2 +i*2*units::cm); */
	
/* 	std::pair<SlimMergeGeomCell*,Point> temp_results = cluster1->get_closest_point_mcell(test_point); */
/* 	//reuse this  */
/* 	Point test_point1 = temp_results.second; */
/* 	if (sqrt(pow(test_point1.x-test_point.x,2)+pow(test_point1.y-test_point.y,2)+pow(test_point1.z-test_point.z,2))<1.5*units::cm){ */
/* 	  double temp_dis = (test_point1.x - p2.x)*dir1.X() + (test_point1.y - p2.y)*dir1.Y() + (test_point1.z - p2.z)*dir1.Z(); */
/* 	  temp_dis *=-1; */
/* 	  if (temp_dis < min_dis) min_dis = temp_dis; */
/* 	  if (temp_dis > max_dis) max_dis = temp_dis; */
/* 	} */
/*       } */
/*       //std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << min_dis/units::cm << " " << max_dis/units::cm << " " << length_1/units::cm << std::endl; */
/*       if ((max_dis - min_dis)>2.5*units::cm) return true; */
/*     } */
/*   } */


/*   return false; */
/* } */
