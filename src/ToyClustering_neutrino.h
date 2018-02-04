void WireCell2dToy::Clustering_examine_x_boundary(WireCell::PR3DClusterSelection& live_clusters,std::map<WireCell::PR3DCluster*,double>& cluster_length_map){
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();

  std::vector<PR3DCluster*> new_clusters;
  std::vector<PR3DCluster*> del_clusters;
  
  for (size_t i=0;i!=live_clusters.size();i++){
    PR3DCluster *cluster = live_clusters.at(i);
    // only examine big clusters ... 
    if (cluster_length_map[cluster] > 5*units::cm&&cluster_length_map[cluster]<150*units::cm){
      //cluster->Create_point_cloud();
      std::vector<PR3DCluster*> clusters = cluster->examine_x_boundary();
      if (clusters.size()!=0){
	del_clusters.push_back(cluster);
	std::copy(clusters.begin(),clusters.end(),std::back_inserter(new_clusters));
      }
    }
  }

  for (auto it=new_clusters.begin(); it!=new_clusters.end(); it++){
    PR3DCluster *ncluster = (*it);
    //ncluster->Create_point_cloud();
    std::vector<int> range_v1 = ncluster->get_uvwt_range();
    double length_1 = sqrt(2./3. * (pow(pitch_u*range_v1.at(0),2) + pow(pitch_v*range_v1.at(1),2) + pow(pitch_w*range_v1.at(2),2)) + pow(time_slice_width*range_v1.at(3),2));
    cluster_length_map[ncluster] = length_1;
    live_clusters.push_back(ncluster);
  }
  
  for (auto it=del_clusters.begin(); it!=del_clusters.end(); it++){
    PR3DCluster *ocluster = (*it);
    cluster_length_map.erase(ocluster);
    live_clusters.erase(find(live_clusters.begin(), live_clusters.end(), ocluster));
    delete ocluster;
  }
  
}


void WireCell2dToy::Clustering_neutrino(WireCell::PR3DClusterSelection& live_clusters,std::map<WireCell::PR3DCluster*,double>& cluster_length_map, int num_try){
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  TVector3 drift_dir(1,0,0);
  TVector3 vertical_dir(0,1,0);

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
  
  // find all the clusters that are inside the box ...
  WireCell::PR3DClusterSelection contained_clusters;
  WireCell::PR3DClusterSelection candidate_clusters;
  
  for (size_t i=0;i!=live_clusters.size();i++){
    PR3DCluster *cluster = live_clusters.at(i);
    cluster->Create_point_cloud();
    std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> hl_wcps = cluster->get_highest_lowest_wcps();
    std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> fb_wcps = cluster->get_front_back_wcps();
    std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> el_wcps = cluster->get_earliest_latest_wcps();

    //std::cout << cluster->get_cluster_id()  << " " << cluster_length_map[cluster] /units::cm << " " << el_wcps.first.x/units::cm << " " << el_wcps.second.x/units::cm << std::endl;
   
    
    if (el_wcps.first.x < -1*units::cm || el_wcps.second.x > 257*units::cm || cluster_length_map[cluster] < 6.0*units::cm) continue;

    
    bool flag_fy = false;
    bool flag_by = false;
    bool flag_fx = false;
    bool flag_bx = false;
    bool flag_fz = false;
    bool flag_bz = false;
   

    std::vector<WCPointCloud<double>::WCPoint> saved_wcps;
    if (hl_wcps.first.y > 101.5*units::cm){
      saved_wcps.push_back(hl_wcps.first);
      flag_fy = true;
    }

    if (hl_wcps.second.y < -99.5*units::cm){
      saved_wcps.push_back(hl_wcps.second);
      flag_by = true;
    }

    if (fb_wcps.first.z > 1022*units::cm){
      bool flag_save = true;
      for (size_t j=0;j!=saved_wcps.size();j++){
	double dis = sqrt(pow(saved_wcps.at(j).x - fb_wcps.first.x,2) + pow(saved_wcps.at(j).y - fb_wcps.first.y,2) + pow(saved_wcps.at(j).z - fb_wcps.first.z,2));
	if (dis < 15*units::cm){
	  flag_save = false;
	  break;
	}
      }
      if (flag_save){
	saved_wcps.push_back(fb_wcps.first);
	flag_bz = true;
      }
    }

    if (fb_wcps.second.z < 15*units::cm){
      bool flag_save = true;
      for (size_t j=0;j!=saved_wcps.size();j++){
	double dis = sqrt(pow(saved_wcps.at(j).x - fb_wcps.second.x,2) + pow(saved_wcps.at(j).y - fb_wcps.second.y,2) + pow(saved_wcps.at(j).z - fb_wcps.second.z,2));
	if (dis < 15*units::cm){
	  flag_save = false;
	  break;
	}
      }
      if ( flag_save){
	saved_wcps.push_back(fb_wcps.second);
	flag_fz = true;
      }
    }

    if (el_wcps.first.x < 1*units::cm){
      bool flag_save = true;
      for (size_t j=0;j!=saved_wcps.size();j++){
	double dis = sqrt(pow(saved_wcps.at(j).x - el_wcps.first.x,2) + pow(saved_wcps.at(j).y - el_wcps.first.y,2) + pow(saved_wcps.at(j).z - el_wcps.first.z,2));
	if (dis < 15*units::cm){
	  flag_save = false;
	  break;
	}
      }
      if (flag_save){
	saved_wcps.push_back(el_wcps.first);
	flag_fx = true;
      }
    }

    if (el_wcps.second.x > 255*units::cm){
      bool flag_save = true;
      for (size_t j=0;j!=saved_wcps.size();j++){
	double dis = sqrt(pow(saved_wcps.at(j).x - el_wcps.second.x,2) + pow(saved_wcps.at(j).y - el_wcps.second.y,2) + pow(saved_wcps.at(j).z - el_wcps.second.z,2));
	if (dis < 15*units::cm){
	  flag_save = false;
	  break;
	}
      }
      if (flag_save){
	saved_wcps.push_back(el_wcps.second);
	flag_bx = true;
      }
    }
    if (saved_wcps.size()<=1 ){
      candidate_clusters.push_back(cluster);
      contained_clusters.push_back(cluster);
    }

    

    if (saved_wcps.size()>=2 && (flag_fx && flag_bx || flag_fy && flag_by || flag_fz && flag_bz)){
    }else{
      contained_clusters.push_back(cluster); 
    } 


    

    /* if (cluster_length_map[cluster] > 60*units::cm) */
    //  std::cout << cluster->get_cluster_id() << " " << saved_wcps.size() << " " << cluster_length_map[cluster] /units::cm << std::endl;
   
    
    // std::cout << el_wcps.first.x/units::cm << " " << el_wcps.second.x/units::cm << " " << hl_wcps.first.y/units::cm << " " << hl_wcps.second.y/units::cm << " " << fb_wcps.first.z/units::cm << " " << fb_wcps.second.z/units::cm << " " << saved_wcps.size() << std::endl;
    
    //std::cout << el_wcps.first.x/units::cm << " " << el_wcps.second.x/units::cm << std::endl;
  }

  std::map<PR3DCluster*, std::pair<PR3DCluster*, double>> cluster_close_cluster_map;
  // calculate the closest distance??? ... 
  for (size_t i=0;i!=live_clusters.size();i++){
    PR3DCluster *cluster1 = live_clusters.at(i);
    ToyPointCloud *cloud1 = cluster1->get_point_cloud();
    for (size_t j=i+1;j!=live_clusters.size();j++){
      PR3DCluster *cluster2 = live_clusters.at(j);
      ToyPointCloud *cloud2 = cluster2->get_point_cloud();

      std::tuple<int,int,double> results =  cloud2->get_closest_points(cloud1);
      double dis = std::get<2>(results);

      if (cluster_close_cluster_map.find(cluster1)==cluster_close_cluster_map.end()){
	cluster_close_cluster_map[cluster1] = std::make_pair(cluster2, dis);
      }else{
	if (dis < cluster_close_cluster_map[cluster1].second)
	  cluster_close_cluster_map[cluster1] = std::make_pair(cluster2, dis);
      }

      if (cluster_close_cluster_map.find(cluster2)==cluster_close_cluster_map.end()){
	cluster_close_cluster_map[cluster2] = std::make_pair(cluster1, dis);
      }else{
	if (dis < cluster_close_cluster_map[cluster2].second)
	  cluster_close_cluster_map[cluster2] = std::make_pair(cluster1, dis);
      }
    }
  }
  
  
  //  std::cout << contained_clusters.size() << " " << candidate_clusters.size() << std::endl;

  std::set<std::pair<PR3DCluster*, PR3DCluster*>> to_be_merged_pairs;

  std::set< PR3DCluster*> used_clusters; 

  std::map<PR3DCluster*, ToyPointCloud*> cluster_cloud_map;

  for (auto it = candidate_clusters.begin(); it!=candidate_clusters.end(); it++){
    PR3DCluster *cluster1 = (*it);
    cluster1->Calc_PCA();
  }
  
  // ignore very small ones?
  // two short ones, NC pi0 case
  // one short one and one big one, CC pi0
  for (auto it = candidate_clusters.begin(); it!=candidate_clusters.end(); it++){
    PR3DCluster *cluster1 = (*it);
    cluster1->Create_point_cloud();
    ToyPointCloud *cloud1 = cluster1->get_point_cloud();
    for (auto it1 = contained_clusters.begin(); it1!=contained_clusters.end(); it1++){
      PR3DCluster *cluster2 = (*it1);
      //can not be the same
      if (cluster2==cluster1) continue;
      if (cluster_length_map[cluster2] > 150*units::cm){
	TVector3 dir1(cluster2->get_PCA_axis(0).x, cluster2->get_PCA_axis(0).y, cluster2->get_PCA_axis(0).z);
	if (fabs(dir1.Angle(vertical_dir)-3.1415926/2.)/3.1415926*180.>80)
	  continue;
      }
      cluster2->Create_point_cloud();
      ToyPointCloud *cloud2 = cluster2->get_point_cloud();

      std::tuple<int,int,double> results =  cloud2->get_closest_points(cloud1);
      double dis = std::get<2>(results);

      /* if (dis < 2.5*units::cm) */
      /* 	std::cout << dis/units::cm <<  " " << cluster_length_map[cluster1]/units::cm << " " << cluster_length_map[cluster2]/units::cm << " " << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << std::endl; */
      
      // cannot be too far ...
      /* if (num_try==0 && (dis > 80*units::cm || cluster_length_map[cluster1] > 60*units::cm   && dis > 10*units::cm)) continue; */
      if ( (dis > 80*units::cm || cluster_length_map[cluster1] > 80*units::cm   && dis > 10*units::cm)) continue;

      if (cluster_cloud_map.find(cluster1) == cluster_cloud_map.end()){
	cluster1->Calc_PCA();
	Point center = cluster1->get_center();
	TVector3 main_dir(cluster1->get_PCA_axis(0).x, cluster1->get_PCA_axis(0).y, cluster1->get_PCA_axis(0).z);
	main_dir.SetMag(1);
	
	ToyPointCloud* cloud1_ext = new ToyPointCloud(angle_u, angle_v, angle_w);
	cluster_cloud_map[cluster1] = cloud1_ext;
	WireCell::PointVector pts;
	std::pair<Point,Point> extreme_pts = cluster1->get_two_extreme_points();

	TVector3 dir1 = cluster1->VHoughTrans(extreme_pts.first,30*units::cm);
	TVector3 dir2 = cluster1->VHoughTrans(extreme_pts.second,30*units::cm);

	bool flag_enable_temp = false;
	std::pair<Point,Point> temp_extreme_pts;
	TVector3 temp_dir1;
	TVector3 temp_dir2;
	int num_clusters = 0;

	//std::cout << cluster1->get_PCA_value(1) << " " << cluster1->get_PCA_value(0) << " " << cluster_length_map[cluster1]/units::cm << std::endl;
	
	if (cluster1->get_num_points(extreme_pts.first,15*units::cm) <= 75 && cluster1->get_num_points() >75 ||
	    cluster1->get_num_points(extreme_pts.second,15*units::cm) <= 75 && cluster1->get_num_points() > 75 ||
	    cluster1->get_PCA_value(1) > 0.022 * cluster1->get_PCA_value(0) && cluster_length_map[cluster1] > 45*units::cm){
	  std::vector<PR3DCluster*> sep_clusters = Separate_2(cluster1,2.5*units::cm);
	  PR3DCluster *largest_cluster = 0;
	  int max_num_points = 0;
	  for (size_t j=0;j!=sep_clusters.size();j++){
	    sep_clusters.at(j)->Create_point_cloud();
	    if (sep_clusters.at(j)->get_num_points() > max_num_points){
	      max_num_points = sep_clusters.at(j)->get_num_points();
	      largest_cluster = sep_clusters.at(j);
	    }
	  }
	  
	  temp_extreme_pts = largest_cluster->get_two_extreme_points();
	  center = largest_cluster->get_center();
	  main_dir.SetXYZ(largest_cluster->get_PCA_axis(0).x,largest_cluster->get_PCA_axis(0).y,largest_cluster->get_PCA_axis(0).z );
	  num_clusters = sep_clusters.size(); 

	  // std::cout << cluster1->get_cluster_id() << " " << num_clusters << std::endl;
	  
	  for (size_t j=0;j!=sep_clusters.size();j++){
	    delete sep_clusters.at(j);
	  }

	  // if (fabs(dir1.Angle(main_dir)-3.1415926/2.)<75/180.*3.1415926 ){
	  flag_enable_temp = true;
	    /* if (pow(extreme_pts.first.x-temp_extreme_pts.first.x,2)+pow(extreme_pts.first.y-temp_extreme_pts.first.y,2)+pow(extreme_pts.first.z-temp_extreme_pts.first.z,2) < pow(extreme_pts.first.x-temp_extreme_pts.second.x,2)+pow(extreme_pts.first.y-temp_extreme_pts.second.y,2)+pow(extreme_pts.first.z-temp_extreme_pts.second.z,2)){ */
	  // extreme_pts.first = temp_extreme_pts.first;
	      /* }else{ */
	    /*   extreme_pts.first = temp_extreme_pts.second; */
	    /* } */
	  temp_dir1 = cluster1->VHoughTrans(temp_extreme_pts.first,30*units::cm);
	    // }

	   // if (fabs(dir2.Angle(main_dir)-3.1415926/2.)<75/180.*3.1415926){
	   //  extreme_pts.second = temp_extreme_pts.second;
	    
	  temp_dir2 = cluster1->VHoughTrans(temp_extreme_pts.second,30*units::cm);
	    //}

	  
	    /* if (pow(extreme_pts.second.x-temp_extreme_pts.second.x,2)+pow(extreme_pts.second.y-temp_extreme_pts.second.y,2)+pow(extreme_pts.second.z-temp_extreme_pts.second.z,2) < pow(extreme_pts.second.x-temp_extreme_pts.first.x,2)+pow(extreme_pts.second.y-temp_extreme_pts.first.y,2)+pow(extreme_pts.second.z-temp_extreme_pts.first.z,2)){ */
	    //extreme_pts.second = temp_extreme_pts.second;
	    /* }else{ */
	    /*   extreme_pts.second = temp_extreme_pts.first; */
	    /* } */
	    
	    
	    
	 
	}


	/* if (cluster_length_map[cluster2]>275*units::cm){  */
	//  std::cout << cluster1->get_cluster_id()<< " " << extreme_pts.first.x/units::cm << " " << extreme_pts.first.y/units::cm << " " << extreme_pts.first.z/units::cm << " " << extreme_pts.second.x/units::cm << " " << extreme_pts.second.y/units::cm << " " << extreme_pts.second.z/units::cm << std::endl;
	/* }  */
	
	/* std::cout << fabs(dir1.Angle(main_dir)-3.1415926/2.)/3.1415926*180. << " " << fabs(dir2.Angle(main_dir)-3.1415926/2.)/3.1415926*180. << std::endl; */
	
	dir1 *=-1;
	dir2 *=-1;
	dir1.SetMag(1);
	dir2.SetMag(1);
	if (flag_enable_temp){
	  temp_dir1 *=-1;
	  temp_dir2 *=-1;
	  temp_dir1.SetMag(1);
	  temp_dir2.SetMag(1);
	}

	bool flag_add1 = true;
	if (cluster1->get_num_points(extreme_pts.first,15*units::cm) <= 75 && cluster_length_map[cluster1] > 60*units::cm
	    || flag_enable_temp && num_clusters >=4  )
	  flag_add1 = false;
	
	bool flag_add2 = true;
	if (cluster1->get_num_points(extreme_pts.second,15*units::cm) <= 75 && cluster_length_map[cluster1] > 60*units::cm
	    || flag_enable_temp && num_clusters >=4  )
	  flag_add2 = false;
	
	for (size_t j=0;j!=150;j++){
	  if (flag_add1){
	    Point pt1(extreme_pts.first.x + dir1.X() * (j+1) * 0.5*units::cm,
		      extreme_pts.first.y + dir1.Y() * (j+1) * 0.5*units::cm,
		      extreme_pts.first.z + dir1.Z() * (j+1) * 0.5*units::cm);
	    pts.push_back(pt1);
	  }
	  if (flag_add2){
	    Point pt2(extreme_pts.second.x + dir2.X() * (j+1) * 0.5*units::cm,
		      extreme_pts.second.y + dir2.Y() * (j+1) * 0.5*units::cm,
		      extreme_pts.second.z + dir2.Z() * (j+1) * 0.5*units::cm);
	    pts.push_back(pt2);
	  }
	  if (flag_enable_temp){
	    Point pt1(temp_extreme_pts.first.x + temp_dir1.X() * (j+1) * 0.5*units::cm,
	  	      temp_extreme_pts.first.y + temp_dir1.Y() * (j+1) * 0.5*units::cm,
	  	      temp_extreme_pts.first.z + temp_dir1.Z() * (j+1) * 0.5*units::cm);
	    pts.push_back(pt1);
	    Point pt2(temp_extreme_pts.second.x + temp_dir2.X() * (j+1) * 0.5*units::cm,
	  	      temp_extreme_pts.second.y + temp_dir2.Y() * (j+1) * 0.5*units::cm,
	  	      temp_extreme_pts.second.z + temp_dir2.Z() * (j+1) * 0.5*units::cm);
	    pts.push_back(pt2);
	  }
	  
	  if ((!flag_add1) && (!flag_add2)){
	    pts.push_back(extreme_pts.first);
	    pts.push_back(extreme_pts.second);
	  }
	   
	    
	  if (cluster_length_map[cluster1] < 60*units::cm ){
	    TVector3 temp1(extreme_pts.first.x - center.x, extreme_pts.first.y - center.y, extreme_pts.first.z - center.z);
	    double length1 = temp1.Dot(main_dir);
	    if (length1 > 0){
	      Point pt3(center.x + main_dir.X() * (length1 + (j+1) * 0.5*units::cm),
	  		center.y + main_dir.Y() * (length1 + (j+1) * 0.5*units::cm),
	  		center.z + main_dir.Z() * (length1 + (j+1) * 0.5*units::cm));
	      pts.push_back(pt3);
	    }else{
	      Point pt3(center.x + main_dir.X() * (length1 - (j+1) * 0.5*units::cm) ,
	  		center.y + main_dir.Y() * (length1 - (j+1) * 0.5*units::cm),
	  		center.z + main_dir.Z() * (length1 - (j+1) * 0.5*units::cm));
	      pts.push_back(pt3);
	    }
	    
	    
	    TVector3 temp2(extreme_pts.second.x - center.x, extreme_pts.second.y - center.y, extreme_pts.second.z -center.z);
	    double length2 = temp2.Dot(main_dir);

	    if (length2 > 0 ){
	      Point pt4(center.x + main_dir.X() * (length2 + (j+1) * 0.5*units::cm),
	  		center.y + main_dir.Y() * (length2 + (j+1) * 0.5*units::cm),
	  		center.z + main_dir.Z() * (length2 + (j+1) * 0.5*units::cm));
	      pts.push_back(pt4);
	    }else{
	      Point pt4(center.x + main_dir.X() * (length2 - (j+1) * 0.5*units::cm),
	  		center.y + main_dir.Y() * (length2 - (j+1) * 0.5*units::cm),
	  		center.z + main_dir.Z() * (length2 - (j+1) * 0.5*units::cm));
	      pts.push_back(pt4);
	    }
	    
	  }
	}
	cloud1_ext->AddPoints(pts);
	cloud1_ext->build_kdtree_index();
	
      }
      
      if (cluster_cloud_map.find(cluster2) == cluster_cloud_map.end()){
	cluster2->Calc_PCA();
	Point center = cluster2->get_center();
	TVector3 main_dir(cluster2->get_PCA_axis(0).x, cluster2->get_PCA_axis(0).y, cluster2->get_PCA_axis(0).z);
	main_dir.SetMag(1);
	
	ToyPointCloud* cloud2_ext = new ToyPointCloud(angle_u, angle_v, angle_w);
	cluster_cloud_map[cluster2] = cloud2_ext;
	WireCell::PointVector pts;
	std::pair<Point,Point> extreme_pts = cluster2->get_two_extreme_points();
	
	TVector3 dir1 = cluster2->VHoughTrans(extreme_pts.first,30*units::cm);
	TVector3 dir2 = cluster2->VHoughTrans(extreme_pts.second,30*units::cm);

	bool flag_enable_temp = false;
	std::pair<Point,Point> temp_extreme_pts;
	TVector3 temp_dir1;
	TVector3 temp_dir2;
	int num_clusters = 0;
	
	if (cluster2->get_num_points(extreme_pts.first,15*units::cm) <= 75 && cluster2->get_num_points() >75 ||
	    cluster2->get_num_points(extreme_pts.second,15*units::cm) <= 75 && cluster2->get_num_points() > 75 ||
	    cluster2->get_PCA_value(1) > 0.022 * cluster2->get_PCA_value(0) && cluster_length_map[cluster2] > 45*units::cm){
	  std::vector<PR3DCluster*> sep_clusters = Separate_2(cluster2,2.5*units::cm);
	  PR3DCluster *largest_cluster = 0;
	  int max_num_points = 0;
	  for (size_t j=0;j!=sep_clusters.size();j++){
	    sep_clusters.at(j)->Create_point_cloud();
	    if (sep_clusters.at(j)->get_num_points() > max_num_points){
	      max_num_points = sep_clusters.at(j)->get_num_points();
	      largest_cluster = sep_clusters.at(j);
	    }
	  }
	  temp_extreme_pts = largest_cluster->get_two_extreme_points();
	  center = largest_cluster->get_center();
	  main_dir.SetXYZ(largest_cluster->get_PCA_axis(0).x,largest_cluster->get_PCA_axis(0).y,largest_cluster->get_PCA_axis(0).z );
	  num_clusters = sep_clusters.size();

	  //	  std::cout << cluster2->get_cluster_id() << " " << num_clusters << std::endl;

	  for (size_t j=0;j!=sep_clusters.size();j++){
	    delete sep_clusters.at(j);
	  }
	  temp_dir1 = cluster2->VHoughTrans(temp_extreme_pts.first,30*units::cm);
	  temp_dir2 = cluster2->VHoughTrans(temp_extreme_pts.second,30*units::cm);
	  
	  //	  if (fabs(dir1.Angle(main_dir)-3.1415926/2.)<75/180.*3.1415926 ){
	    /* if (pow(extreme_pts.first.x-temp_extreme_pts.first.x,2)+pow(extreme_pts.first.y-temp_extreme_pts.first.y,2)+pow(extreme_pts.first.z-temp_extreme_pts.first.z,2) < pow(extreme_pts.first.x-temp_extreme_pts.second.x,2)+pow(extreme_pts.first.y-temp_extreme_pts.second.y,2)+pow(extreme_pts.first.z-temp_extreme_pts.second.z,2)){ */
	  flag_enable_temp = true;
	    //extreme_pts.first = temp_extreme_pts.first; 
	    /* }else{ */
	    /*   extreme_pts.first = temp_extreme_pts.second; */
	    /* } */
	  
	    //	  }
	  
	    //	  if (fabs(dir2.Angle(main_dir)-3.1415926/2.)<75/180.*3.1415926){
	    //extreme_pts.second = temp_extreme_pts.second;
	  
	    //}
	  
	  /* if (){ */
	  /*   /\* if (pow(extreme_pts.second.x-temp_extreme_pts.second.x,2)+pow(extreme_pts.second.y-temp_extreme_pts.second.y,2)+pow(extreme_pts.second.z-temp_extreme_pts.second.z,2) < pow(extreme_pts.second.x-temp_extreme_pts.first.x,2)+pow(extreme_pts.second.y-temp_extreme_pts.first.y,2)+pow(extreme_pts.second.z-temp_extreme_pts.first.z,2)){ *\/ */
	  /*   extreme_pts.second = temp_extreme_pts.second; */
	  /*   /\* }else{ *\/ */
	  /*   /\*   extreme_pts.second = temp_extreme_pts.first; *\/ */
	  /*   /\* } *\/ */
	  /*   dir2 = cluster2->VHoughTrans(extreme_pts.second,30*units::cm); */
	  /* } */
	}

	/* if (cluster2->get_cluster_id()==19 || cluster2->get_cluster_id()==20){ */
	//std::cout << cluster2->get_cluster_id()<< " " << extreme_pts.first.x/units::cm << " " << extreme_pts.first.y/units::cm << " " << extreme_pts.first.z/units::cm << " " << extreme_pts.second.x/units::cm << " " << extreme_pts.second.y/units::cm << " " << extreme_pts.second.z/units::cm << " " << cluster2->get_num_points(extreme_pts.first,15*units::cm)<< " " << cluster2->get_num_points(extreme_pts.second,15*units::cm) << std::endl; 
	/* } */


	
	
	dir1 *=-1;
	dir2 *=-1;
	dir1.SetMag(1);
	dir2.SetMag(1);
	if (flag_enable_temp){
	  temp_dir1 *= -1;
	  temp_dir2 *= -1;
	  temp_dir1.SetMag(1);
	  temp_dir2.SetMag(1);
	}

	bool flag_add1 = true;
	if (cluster2->get_num_points(extreme_pts.first,15*units::cm) <= 75 && cluster_length_map[cluster2] > 60*units::cm
	    || flag_enable_temp && num_clusters >=4 )
	  flag_add1 = false;
	bool flag_add2 = true;
	if (cluster2->get_num_points(extreme_pts.second,15*units::cm) <= 75 && cluster_length_map[cluster2] > 60*units::cm
	    || flag_enable_temp && num_clusters >=4  )
	  flag_add2 = false;

	
	for (size_t j=0;j!=150;j++){
	  if (flag_add1){
	    Point pt1(extreme_pts.first.x + dir1.X() * (j+1) * 0.5*units::cm,
		      extreme_pts.first.y + dir1.Y() * (j+1) * 0.5*units::cm,
		      extreme_pts.first.z + dir1.Z() * (j+1) * 0.5*units::cm);
	    pts.push_back(pt1);
	  }
	  if (flag_add2){
	    Point pt2(extreme_pts.second.x + dir2.X() * (j+1) * 0.5*units::cm,
		      extreme_pts.second.y + dir2.Y() * (j+1) * 0.5*units::cm,
		      extreme_pts.second.z + dir2.Z() * (j+1) * 0.5*units::cm);
	    pts.push_back(pt2);
	  }

	  if (flag_enable_temp){
	    Point pt1(temp_extreme_pts.first.x + temp_dir1.X() * (j+1) * 0.5*units::cm,
	  	      temp_extreme_pts.first.y + temp_dir1.Y() * (j+1) * 0.5*units::cm,
	  	      temp_extreme_pts.first.z + temp_dir1.Z() * (j+1) * 0.5*units::cm);
	    pts.push_back(pt1);
	    Point pt2(temp_extreme_pts.second.x + temp_dir2.X() * (j+1) * 0.5*units::cm,
	  	      temp_extreme_pts.second.y + temp_dir2.Y() * (j+1) * 0.5*units::cm,
	  	      temp_extreme_pts.second.z + temp_dir2.Z() * (j+1) * 0.5*units::cm);
	    pts.push_back(pt2);
	  }
	  
	  if ((!flag_add1) && (!flag_add2)){
	    pts.push_back(extreme_pts.first);
	    pts.push_back(extreme_pts.second);
	  }
	  
	  if (cluster_length_map[cluster2] < 60*units::cm){
	    TVector3 temp1(extreme_pts.first.x - center.x, extreme_pts.first.y - center.y, extreme_pts.first.z - center.z);
	    double length1 = temp1.Dot(main_dir);
	    if (length1 > 0){
	      Point pt3(center.x + main_dir.X() * (length1 + (j+1) * 0.5*units::cm),
	  		center.y + main_dir.Y() * (length1 + (j+1) * 0.5*units::cm),
	  		center.z + main_dir.Z() * (length1 + (j+1) * 0.5*units::cm));
	      pts.push_back(pt3);
	    }else{
	      Point pt3(center.x + main_dir.X() * (length1 - (j+1) * 0.5*units::cm) ,
	  		center.y + main_dir.Y() * (length1 - (j+1) * 0.5*units::cm),
	  		center.z + main_dir.Z() * (length1 - (j+1) * 0.5*units::cm));
	      pts.push_back(pt3);
	    }
	    
	    
	    TVector3 temp2(extreme_pts.second.x - center.x, extreme_pts.second.y - center.y, extreme_pts.second.z -center.z);
	    double length2 = temp2.Dot(main_dir);

	    if (length2 > 0 ){
	      Point pt4(center.x + main_dir.X() * (length2 + (j+1) * 0.5*units::cm),
	  		center.y + main_dir.Y() * (length2 + (j+1) * 0.5*units::cm),
	  		center.z + main_dir.Z() * (length2 + (j+1) * 0.5*units::cm));
	      pts.push_back(pt4);
	    }else{
	      Point pt4(center.x + main_dir.X() * (length2 - (j+1) * 0.5*units::cm),
	  		center.y + main_dir.Y() * (length2 - (j+1) * 0.5*units::cm),
	  		center.z + main_dir.Z() * (length2 - (j+1) * 0.5*units::cm));
	      pts.push_back(pt4);
	    }
	  }
	  
	}
	cloud2_ext->AddPoints(pts);
	cloud2_ext->build_kdtree_index();
      }

      ToyPointCloud *cloud1_ext = cluster_cloud_map[cluster1];
      ToyPointCloud *cloud2_ext = cluster_cloud_map[cluster2];

      
      int merge_type = 0;
      bool flag_merge = false;
      {
	std::tuple<int,int,double> results_1 =  cloud1_ext->get_closest_points(cloud2);
	Point test_pt(cloud1_ext->get_cloud().pts.at(std::get<0>(results_1)).x,
		      cloud1_ext->get_cloud().pts.at(std::get<0>(results_1)).y,
		      cloud1_ext->get_cloud().pts.at(std::get<0>(results_1)).z);
	
	Point test_pt1(cloud2->get_cloud().pts.at(std::get<1>(results_1)).x,
		       cloud2->get_cloud().pts.at(std::get<1>(results_1)).y,
		       cloud2->get_cloud().pts.at(std::get<1>(results_1)).z);
	
	
	double dis1 = std::get<2>(results_1);
	double dis2 = cloud1->get_closest_dis(test_pt);

	//	std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << flag_merge << std::endl;
	
	/* if (dis1 < std::max(3.5*units::cm,dis2*sin(15/180.*3.1415926)) && num_try == 0)  */
	/*   flag_merge = true; */
	if (dis1 < std::min(std::max(4.5*units::cm,dis2*sin(15/180.*3.1415926)),12*units::cm) &&
	    (cluster_length_map[cluster2]>25*units::cm ||
	     cluster_length_map[cluster1] <= cluster_length_map[cluster2]) ||
	    dis1 < std::min(std::max(2.5*units::cm,dis2*sin(10/180.*3.1415926)),10*units::cm)) {

	  if (cluster_length_map[cluster1]>25*units::cm && cluster1->get_PCA_value(1) < 0.0015 * cluster1->get_PCA_value(0)){
	    flag_merge = false;
	  }else{
	    if (cluster_length_map[cluster2]<30*units::cm){ 
	      flag_merge = true;
	      if (cluster_length_map[cluster1]>15*units::cm && cluster1->get_PCA_value(1) < 0.012 * cluster1->get_PCA_value(0)){
		if (dis1 > std::max(2.5*units::cm, dis2*sin(7.5/180.*3.1415926)))
		  flag_merge = false;
		//	std::cout << JudgeSeparateDec_1(cluster1,drift_dir) << " A1 " << cluster1->get_PCA_value(0) << " " << cluster1->get_PCA_value(1) << std::endl; 
	      }
	      //std::cout << JudgeSeparateDec_1(cluster1,drift_dir) << " A2 " << cluster1->get_PCA_value(0) << " " << cluster1->get_PCA_value(1) << std::endl; 
	    }else if (JudgeSeparateDec_1(cluster2,drift_dir)){
	      flag_merge = cluster2->judge_vertex(test_pt1);
	      if (cluster_length_map[cluster1]>15*units::cm && cluster1->get_PCA_value(1) < 0.012 * cluster1->get_PCA_value(0)){
		if (dis1 > std::max(2.5*units::cm, dis2*sin(7.5/180.*3.1415926)))
		  flag_merge = false;
		//std::cout << JudgeSeparateDec_1(cluster1,drift_dir) << " B1 " << cluster1->get_PCA_value(0) << " " << cluster1->get_PCA_value(1) << std::endl; 
	      }
	      //  std::cout << JudgeSeparateDec_1(cluster1,drift_dir) << " B2 " << cluster1->get_PCA_value(0) << " " << cluster1->get_PCA_value(1) << " " << test_pt1.x/units::cm << " " << test_pt1.y/units::cm << " " << test_pt1.z/units::cm << std::endl; 
	    }else{
	      //std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << std::endl;
	      if (dis2 < 5*units::cm){
		flag_merge = cluster2->judge_vertex(test_pt1,2./3.);
	      }else if (dis < 0.5*units::cm){
		flag_merge = cluster2->judge_vertex(test_pt1,0.5,0.6);
	      }else{
		flag_merge = cluster2->judge_vertex(test_pt1);
	      }

	      if (cluster_length_map[cluster1]>15*units::cm && cluster1->get_PCA_value(1) < 0.012 * cluster1->get_PCA_value(0)){
		if (dis1 > std::max(2.5*units::cm, dis2*sin(7.5/180.*3.1415926)))
		  flag_merge = false;
	      }
	      
	    }
	    merge_type = 1;
	  }
	  // std::cout << flag_merge << " " << merge_type << std::endl;

	  if (cluster_close_cluster_map[cluster1].second < 1.2*units::cm && cluster_close_cluster_map[cluster1].first != cluster2 && cluster_length_map[cluster_close_cluster_map[cluster1].first] > 60*units::cm){
	    flag_merge = false;
	  }
	  
	  if (flag_merge && cluster_length_map[cluster1]>150*units::cm && cluster_length_map[cluster2]>150*units::cm){
	    // protect against two long tracks ... 
	    cluster1->Calc_PCA();
	    cluster2->Calc_PCA();
	    TVector3 temp_dir1(cluster1->get_PCA_axis(0).x,cluster1->get_PCA_axis(0).y,cluster1->get_PCA_axis(0).z);
	    TVector3 temp_dir2(cluster2->get_PCA_axis(0).x,cluster2->get_PCA_axis(0).y,cluster2->get_PCA_axis(0).z);
	    if (fabs(temp_dir1.Angle(temp_dir2)-3.1415926/2.)<60/180.*3.1415926)
	      flag_merge = false;
	  }
	}

	/* if (flag_merge)  */
	/* /\*   // std::cout << cluster_length_map[cluster1]/units::cm << " " << cluster1->get_PCA_value(0) << " " << cluster1->get_PCA_value(1) << " " << cluster_length_map[cluster2]/units::cm << " " << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << cluster_close_cluster_map[cluster1].first->get_cluster_id() << " " << cluster_close_cluster_map[cluster1].second/units::cm << std::endl; *\/ */
	/* /\* /\\*   /\\\* if (dis < 3*units::cm && dis1 < 3*units::cm) *\\\/ *\\/ *\/ */
	/*   std::cout << dis1/units::cm << " " << dis2/units::cm << " " << dis/units::cm << " " << */
	/*     cluster_length_map[cluster1]/units::cm << " " << cluster_length_map[cluster2]/units::cm << " " << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << flag_merge << " " << cluster1->get_PCA_value(0) << " " << cluster1->get_PCA_value(1) << " " << cluster2->get_PCA_value(0) << " " << cluster2->get_PCA_value(1) << std::endl; */

	
	if (dis < 1.8*units::cm && cluster_length_map[cluster1] < 75*units::cm && cluster_length_map[cluster2] < 75*units::cm && (cluster_length_map[cluster1]+cluster_length_map[cluster2]) < 120*units::cm){
	  flag_merge = true;
	  merge_type = 2;
	}
      }

      if (!flag_merge){
      	std::tuple<int,int,double> results_2 =  cloud1_ext->get_closest_points(cloud2_ext);
      	Point test_pt(cloud1_ext->get_cloud().pts.at(std::get<0>(results_2)).x,
      		      cloud1_ext->get_cloud().pts.at(std::get<0>(results_2)).y,
      		      cloud1_ext->get_cloud().pts.at(std::get<0>(results_2)).z);
      	Point test_pt1(cloud2_ext->get_cloud().pts.at(std::get<1>(results_2)).x,
      		       cloud2_ext->get_cloud().pts.at(std::get<1>(results_2)).y,
      		       cloud2_ext->get_cloud().pts.at(std::get<1>(results_2)).z);
      	double dis1 = std::get<2>(results_2);
      	double dis2 = cloud1->get_closest_dis(test_pt);
      	double dis3 = cloud2->get_closest_dis(test_pt1);
      	if (dis1 < std::min(std::max(4.5*units::cm, (dis2+dis3)/2.*sin(15/180.*3.1415926)),12*units::cm) &&
      	    dis1 < std::min(std::max(4.5*units::cm, (dis3+dis2)/2.*sin(15/180.*3.1415926)),12*units::cm) &&
	    dis2 + dis3 < 72*units::cm  && cluster_length_map[cluster2] < 60*units::cm
	    && cluster_length_map[cluster1] < 60*units::cm){
	  flag_merge = true;
	  merge_type = 3;
	}else if (dis2 + dis3 < 90*units::cm && dis1 < 2.7*units::cm && dis < 20*units::cm &&
		  cluster_length_map[cluster2] > 30*units::cm &&
		  cluster_length_map[cluster1] > 30*units::cm){
	  cluster1->Calc_PCA();
	  cluster2->Calc_PCA();
	  if (cluster1->get_PCA_value(1) > 0.0015 *cluster1->get_PCA_value(0) && cluster2->get_PCA_value(1)>0.0015 * cluster2->get_PCA_value(0)) {
	    flag_merge = true;
	    merge_type = 3;

	    /* std::cout << dis/units::cm << " " << dis1/units::cm << " A " << dis2/units::cm << " " << dis3/units::cm << " " << */
	    /*   cluster_length_map[cluster1]/units::cm << " " << cluster_length_map[cluster2]/units::cm << " " << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << flag_merge << std::endl; */
	  }
	}


	
      }

      

      if (flag_merge){
	bool flag_proceed = true;
	if (merge_type==1){
	  if (used_clusters.find(cluster1) != used_clusters.end())
	    flag_proceed = false;
	}else if (merge_type == 3){
	  if (used_clusters.find(cluster1) != used_clusters.end() &&
	      used_clusters.find(cluster2) != used_clusters.end())
	    flag_proceed = false;
	}

	if (flag_proceed){
	  // std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << merge_type << " " <<std::endl;

	  to_be_merged_pairs.insert(std::make_pair(cluster1,cluster2));

	  //  if (cluster_length_map[cluster1] < 20*units::cm )
	  if (merge_type==1){
	    used_clusters.insert(cluster1);
	  }else if (merge_type==3){
	    used_clusters.insert(cluster1);
	    used_clusters.insert(cluster2);
	  }
	  /* if (cluster_length_map[cluster2] < 20*units::cm ) */
	  /*   used_clusters.insert(cluster2); */
	  
	}
      }
    }
  }
  
  
  for (auto it = cluster_cloud_map.begin(); it!= cluster_cloud_map.end(); it++){
    delete it->second;
  }

  // to_be_merged_pairs.clear();


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


void WireCell2dToy::Clustering_dis(WireCell::PR3DClusterSelection& live_clusters,std::map<WireCell::PR3DCluster*,double>& cluster_length_map){
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  TVector3 drift_dir(1,0,0);

  int range_cut = 80;
  
  
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
    if (max < range_cut&&cluster_length_map[live_clusters.at(i)] < 15*units::cm){
      //std::cout << cluster_length_map[live_clusters.at(i)] << " " << live_clusters.at(i)->get_cluster_id() << std::endl;
      small_clusters.push_back(live_clusters.at(i));
    }else{
      if (cluster_length_map[live_clusters.at(i)] < 60*units::cm){
	if (JudgeSeparateDec_1(live_clusters.at(i),drift_dir)){
      	  std::vector<PR3DCluster*> sep_clusters = Separate_2(live_clusters.at(i),2.5*units::cm);
      	  int max = 0;
	  double max_length = 0;
      	  for (auto it = sep_clusters.begin(); it!= sep_clusters.end(); it++){
      	    std::vector<int> ranges = (*it)->get_uvwt_range();
	    double length_1 = sqrt(2./3. * (pow(pitch_u*ranges.at(0),2) + pow(pitch_v*ranges.at(1),2) + pow(pitch_w*ranges.at(2),2)) + pow(time_slice_width*ranges.at(3),2));
      	    for (int j=0;j!=4;j++){
      	      if (ranges.at(j)>max){
      		max =ranges.at(j);
		max_length = length_1;
	      }
	    }
      	    if (max >= range_cut || max_length >= 15*units::cm) break;
      	  }

      	  for (size_t j=0;j!=sep_clusters.size();j++){
      	    delete sep_clusters.at(j);
      	  }
	  
      	  if (max < range_cut && max_length < 15*units::cm){
      	    small_clusters.push_back(live_clusters.at(i));
      	  }else{
      	    big_clusters.push_back(live_clusters.at(i));
      	  }
      	}else{ 
	  big_clusters.push_back(live_clusters.at(i)); 
	} 
      }else{
	big_clusters.push_back(live_clusters.at(i));
      }
    }
  }

  std::set<std::pair<PR3DCluster*, PR3DCluster*>> to_be_merged_pairs;

    
  // clustering small with big ones ... 
  double small_big_dis_cut = 80*units::cm;
  std::set<PR3DCluster*> used_small_clusters;
  for (auto it = small_clusters.begin(); it!= small_clusters.end(); it++){
    PR3DCluster* curr_cluster = (*it);
    curr_cluster->Create_point_cloud();
    ToyPointCloud *cloud1 = curr_cluster->get_point_cloud();
    double min_dis = 1e9; PR3DCluster* min_dis_cluster = 0;
    
    for (auto it1 = big_clusters.begin(); it1!=big_clusters.end(); it1++){
      PR3DCluster* big_cluster = (*it1);
      big_cluster->Create_point_cloud();
      ToyPointCloud *cloud2 = big_cluster->get_point_cloud();
      std::tuple<int,int,double> results =  cloud2->get_closest_points(cloud1);
      double dis = std::get<2>(results);
      if (dis < min_dis){
	min_dis = dis;
	min_dis_cluster = big_cluster;
      }
    }
    if (min_dis < small_big_dis_cut){
      to_be_merged_pairs.insert(std::make_pair(min_dis_cluster,curr_cluster));
      used_small_clusters.insert(curr_cluster);
    }
  }

  // small distance ... 
  double small_small_dis_cut = 5*units::cm;
  for (size_t i=0;i!=small_clusters.size(); i++){
    PR3DCluster *cluster1 = small_clusters.at(i);
    ToyPointCloud *cloud1 = cluster1->get_point_cloud();
    for (size_t j=i+1; j!=small_clusters.size(); j++){
      PR3DCluster *cluster2 = small_clusters.at(j);
      ToyPointCloud *cloud2 = cluster2->get_point_cloud();
      std::tuple<int,int,double> results =  cloud2->get_closest_points(cloud1);
      double dis = std::get<2>(results);
      if (dis < small_small_dis_cut){
	if (used_small_clusters.find(cluster1)!=used_small_clusters.end() && used_small_clusters.find(cluster2)==used_small_clusters.end()||
	    used_small_clusters.find(cluster2)!=used_small_clusters.end() && used_small_clusters.find(cluster1)==used_small_clusters.end()){
	  to_be_merged_pairs.insert(std::make_pair(cluster1,cluster2));
	  used_small_clusters.insert(cluster1);
	  used_small_clusters.insert(cluster2);
	}
      }
    }
  }
  
  
  std::vector<PR3DCluster*> remaining_small_clusters;
  for (auto it = small_clusters.begin(); it!= small_clusters.end(); it++){
    PR3DCluster* curr_cluster = (*it);
    if (used_small_clusters.find(curr_cluster)==used_small_clusters.end())
      remaining_small_clusters.push_back(curr_cluster);
  }

  
  // clustering small with small ones ... 
  small_small_dis_cut = 50*units::cm;
  for (size_t i=0;i!=remaining_small_clusters.size(); i++){
    PR3DCluster *cluster1 = remaining_small_clusters.at(i);
    ToyPointCloud *cloud1 = cluster1->get_point_cloud();
    for (size_t j=i+1; j!=remaining_small_clusters.size(); j++){
      PR3DCluster *cluster2 = remaining_small_clusters.at(j);
      ToyPointCloud *cloud2 = cluster2->get_point_cloud();
      std::tuple<int,int,double> results =  cloud2->get_closest_points(cloud1);
      double dis = std::get<2>(results);
      if (dis < small_small_dis_cut){
	to_be_merged_pairs.insert(std::make_pair(cluster1,cluster2));
      }
    }
  }

  // clustering big ones ...
  std::set<PR3DCluster*> used_big_clusters;
  double big_dis_cut = 3*units::cm;
  double big_dis_range_cut = 16*units::cm;
  for (size_t i=0;i!=big_clusters.size();i++){
    PR3DCluster *cluster1 = big_clusters.at(i);
    for (size_t j=i+1; j!=big_clusters.size();j++){
      PR3DCluster *cluster2 = big_clusters.at(j);
      ToyPointCloud *cloud1, *cloud2;
      if (cluster_length_map[cluster1] > cluster_length_map[cluster2]){
	cloud1 = cluster1->get_point_cloud();
	cloud2 = cluster2->get_point_cloud();
	if (used_big_clusters.find(cluster2)!=used_big_clusters.end()) continue;
      }else{
	cloud1 = cluster2->get_point_cloud();
	cloud2 = cluster1->get_point_cloud();
	if (used_big_clusters.find(cluster1)!=used_big_clusters.end()) continue;
      }
      std::tuple<int,int,double> results =  cloud2->get_closest_points(cloud1);
      double min_dis = std::get<2>(results);
      if (min_dis < big_dis_cut){
	bool flag_merge = true;

	int num_outside_points = 0;
	/* int num_close_points = 0; */
	const int N = cloud2->get_num_points();
	WireCell::WCPointCloud<double>& cloud = cloud2->get_cloud();
	for (int k=0;k!=N;k++){
	  Point test_p1(cloud.pts[k].x,cloud.pts[k].y,cloud.pts[k].z);
	  double close_dis = cloud1->get_closest_dis(test_p1);
	  if (close_dis > big_dis_range_cut){
	    num_outside_points ++ ;
	    if (num_outside_points > 400)
	      break;
	  }
	  /* if (close_dis <= big_dis_cut) */
	  /*   num_close_points ++;  */
	}

	//std::cout << num_outside_points << " " << N << " " << num_close_points << std::endl;
	
	if (num_outside_points > 0.125 * N  || num_outside_points > 400   )
	  flag_merge = false;
	
	if (flag_merge) {
	  to_be_merged_pairs.insert(std::make_pair(cluster1,cluster2));
	  if (cluster_length_map[cluster1] < cluster_length_map[cluster2]){
	    used_big_clusters.insert(cluster1);
	  }else{
	    used_big_clusters.insert(cluster2);
	  }	
	}
      }
    }
  }
  
  
  //merge clusters
  std::vector<std::set<PR3DCluster*>> merge_clusters;
  for (auto it = to_be_merged_pairs.begin(); it!=to_be_merged_pairs.end(); it++){
    PR3DCluster *cluster1 = (*it).first;
    PR3DCluster *cluster2 = (*it).second;
    /* std::cout << cluster1 << " " << cluster2 << " " << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << cluster_length_map[cluster1]/units::cm << " " << cluster_length_map[cluster2]/units::cm << std::endl; */
    
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
