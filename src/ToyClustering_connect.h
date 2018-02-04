#include "WireCellData/Line.h"

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

  std::map<PR3DCluster*, TVector3 > map_cluster_dir1;
  std::map<PR3DCluster*, TVector3 > map_cluster_dir2;


  TVector3 drift_dir(1,0,0);
  TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
  TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
  TVector3 W_dir(0,1,0);
  
  for (size_t i=0;i!=live_clusters.size();i++){
    PR3DCluster* cluster = live_clusters.at(i);
    cluster->Create_point_cloud();
  
    std::pair<Point,Point> extreme_points = cluster->get_two_extreme_points();
    TVector3 main_dir(extreme_points.second.x - extreme_points.first.x,
		      extreme_points.second.y - extreme_points.first.y,
		      extreme_points.second.z - extreme_points.first.z);
    TVector3 dir1, dir2;
    

    if (main_dir.Mag() > 10*units::cm && fabs(main_dir.Angle(drift_dir)-3.1415926/2.) < 5 * 3.1415926 / 180.){
      dir1 = main_dir;
      dir1 *= -1;
      dir2 = main_dir;
    }else if (cluster_length_map[cluster] > 25*units::cm){
      dir1 = cluster->VHoughTrans(extreme_points.first,80*units::cm);
      if (dir1.Mag()!=0) dir1.SetMag(1);
      if (fabs(dir1.Angle(drift_dir) - 3.1415926/2.) < 5*3.1415926/180.){
	dir1.SetXYZ(dir1.X(),(extreme_points.second.y - extreme_points.first.y)/main_dir.Mag(),
		    (extreme_points.second.z - extreme_points.first.z)/main_dir.Mag());
	dir1 *= -1;
      }
      dir2 = cluster->VHoughTrans(extreme_points.second,80*units::cm);
      if (dir2.Mag()!=0) dir2.SetMag(1);
      if (fabs(dir2.Angle(drift_dir) - 3.1415926/2.) < 5*3.1415926/180.){
	dir2.SetXYZ(dir2.X(),(extreme_points.second.y - extreme_points.first.y)/main_dir.Mag(),
		    (extreme_points.second.z - extreme_points.first.z)/main_dir.Mag());
      }
      if (dir1.Dot(main_dir)>0) dir1 *=-1;
      if (dir2.Dot(dir1)>0) dir2 *= -1;
    }else{
      dir1 = global_point_cloud.VHoughTrans(extreme_points.first,extending_dis);
      dir2 = global_point_cloud.VHoughTrans(extreme_points.second,extending_dis);
      if (dir1.Dot(main_dir)>0) dir1 *=-1;
      if (dir2.Dot(dir1)>0) dir2 *= -1;
    }

    bool flag_add_dir1 = true;
    bool flag_add_dir2 = true;
    map_cluster_dir1[cluster] = dir1;
    map_cluster_dir2[cluster] = dir2;

    // judge if something is good ...
    if (cluster_length_map[cluster] < 15*units::cm && i!=0){
      flag_add_dir1 = false;
      flag_add_dir2 = false;

      if (fabs(dir1.Angle(drift_dir)-3.1415926/2.)<7.5*3.1415926/180. ){
    	flag_add_dir1 = true;
      }else{
    	TVector3 tempV1(0,dir1.Y(), dir1.Z());
    	TVector3 tempV5;
    	double angle1 = tempV1.Angle(U_dir);
    	tempV5.SetXYZ(fabs(dir1.X()),sqrt(pow(dir1.Y(),2) + pow(dir1.Z(),2)) * sin(angle1),0);
    	angle1 = tempV5.Angle(drift_dir);
	
    	if (angle1 < 7.5/180.*3.1415926 ){
    	  flag_add_dir1 = true;
	}else{
    	  angle1 = tempV1.Angle(V_dir);
    	  tempV5.SetXYZ(fabs(dir1.X()),sqrt(pow(dir1.Y(),2) + pow(dir1.Z(),2)) * sin(angle1),0);
    	  angle1 = tempV5.Angle(drift_dir);

	  /* if (extreme_points.first.z<20*units::cm) */
	  /*   std::cout << angle1/3.1415926*180. << " " << cluster_length_map[cluster]/units::cm << std::endl; */
	  
    	  if (angle1 < 7.5/180.*3.1415926 ){
    	    flag_add_dir1 = true;
	  }else{
    	    angle1 = tempV1.Angle(W_dir);
    	    tempV5.SetXYZ(fabs(dir1.X()),sqrt(pow(dir1.Y(),2) + pow(dir1.Z(),2)) * sin(angle1),0);
    	    angle1 = tempV5.Angle(drift_dir);

    	  
	    
    	    if (angle1 < 7.5/180.*3.1415926  ){
    	      flag_add_dir1 = true;
    	      
    	    }
    	  }
    	}
      }

      if (fabs(dir2.Angle(drift_dir)-3.1415926/2.)<7.5*3.1415926/180. ){
	flag_add_dir2 = true;
      }else{
	TVector3 tempV2(0,dir2.Y(), dir2.Z());
    	TVector3 tempV6;
    	double angle2 = tempV2.Angle(U_dir);
    	tempV6.SetXYZ(fabs(dir2.X()),sqrt(pow(dir2.Y(),2) + pow(dir2.Z(),2)) * sin(angle2),0);
    	angle2 = tempV6.Angle(drift_dir);
	if (angle2 < 7.5/180.*3.1415926){
	  flag_add_dir2 = true;
	}else{
	  angle2 = tempV2.Angle(V_dir);
    	  tempV6.SetXYZ(fabs(dir2.X()),sqrt(pow(dir2.Y(),2) + pow(dir2.Z(),2)) * sin(angle2),0);
    	  angle2 = tempV6.Angle(drift_dir);
	  if (angle2 < 7.5/180.*3.1415926){
	    flag_add_dir2 = true;
	  }else{
	    angle2 = tempV2.Angle(W_dir);
    	    tempV6.SetXYZ(fabs(dir2.X()),sqrt(pow(dir2.Y(),2) + pow(dir2.Z(),2)) * sin(angle2),0);
    	    angle2 = tempV6.Angle(drift_dir);
	    if (angle2 < 7.5 /180.*3.1415926){
	      flag_add_dir2 = true;
	    }
	  }
	}
      }

    }
    
    
    if (i==0){
      if (fabs(dir1.Angle(drift_dir) - 3.1415926/2.) < 5*3.1415926/180.){
	global_skeleton_cloud.AddPoints(cluster,extreme_points.first,dir1,extending_dis*2,1.2*units::cm, angle/2.);
	dir1 *= -1;
	global_skeleton_cloud.AddPoints(cluster,extreme_points.first,dir1,extending_dis*2,1.2*units::cm, angle/2.);
      }else{
	global_skeleton_cloud.AddPoints(cluster,extreme_points.first,dir1,extending_dis,1.2*units::cm, angle);
	dir1 *= -1;
	global_skeleton_cloud.AddPoints(cluster,extreme_points.first,dir1,extending_dis,1.2*units::cm, angle);
      }

      if (fabs(dir2.Angle(drift_dir) - 3.1415926/2.) < 5*3.1415926/180.){
	global_skeleton_cloud.AddPoints(cluster,extreme_points.second,dir2,extending_dis*2.0,1.2*units::cm, angle/2.);
	dir2 *= -1;
	global_skeleton_cloud.AddPoints(cluster,extreme_points.second,dir2,extending_dis*2.0,1.2*units::cm, angle/2.);
      }else{
	global_skeleton_cloud.AddPoints(cluster,extreme_points.second,dir2,extending_dis,1.2*units::cm, angle);
	dir2 *= -1;
	global_skeleton_cloud.AddPoints(cluster,extreme_points.second,dir2,extending_dis,1.2*units::cm, angle);
      }

      /* if (cluster_length_map[cluster] > extending_dis * 0.8){ */
      /* 	live_clusters.at(i)->Construct_skeleton(); */
      /* 	global_skeleton_cloud.AddPoints(live_clusters.at(i),1); */
      /* } */
    }else{

      if (cluster_length_map[cluster] < 100*units::cm  ||
	  fabs(dir2.Angle(drift_dir) - 3.1415926/2.) < 5*3.1415926/180. && fabs(dir1.Angle(drift_dir) - 3.1415926/2.) < 5*3.1415926/180.  && cluster_length_map[cluster] < 200*units::cm){
	WireCell::WCPointCloud<double>& cloud = cluster->get_point_cloud()->get_cloud();
	int num_total_points = cloud.pts.size(); // total number of points
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
	    std::vector<std::tuple<double, PR3DCluster*, size_t>> results = global_skeleton_cloud.get_2d_points_info(test_point, loose_dis_cut, 0);
	    bool flag_unique = true;
	    if (results.size()>0){
	      std::set<PR3DCluster*> temp_clusters;
	      for (size_t k = 0; k!= results.size(); k++){
		if (std::get<0>(results.at(k)) < global_cloud.pts.at(std::get<2>(results.at(k))).index_u){
		  flag_unique = false;
		  temp_clusters.insert(std::get<1>(results.at(k)));
		}
	      }
	      for (auto it = temp_clusters.begin(); it!= temp_clusters.end(); it++){
		if (map_cluster_num[0].find(*it)==map_cluster_num[0].end()){
		  map_cluster_num[0][*it] = 1;
		}else{
		  map_cluster_num[0][*it] ++;
		}
	      }
	    }
	    if (flag_unique)
	      num_unique[0]++;
	  }else{
	    std::vector<std::tuple<double, PR3DCluster*, size_t>> results = global_skeleton_cloud.get_2d_points_info(test_point, loose_dis_cut, 0);
	    bool flag_unique = true;
	    if (results.size()>0){
	      std::set<PR3DCluster*> temp_clusters;
	      for (size_t k = 0; k!= results.size(); k++){
		if (std::get<0>(results.at(k)) < loose_dis_cut/3.*2.){
		  flag_unique = false;
		  temp_clusters.insert(std::get<1>(results.at(k)));
		}
	      }
	      for (auto it = temp_clusters.begin(); it!= temp_clusters.end(); it++){
		if (map_cluster_num[0].find(*it)==map_cluster_num[0].end()){
		  map_cluster_num[0][*it] = 1;
		}else{
		  map_cluster_num[0][*it] ++;
		}
	      }
	    }
	    if (flag_unique)
	      num_unique[0]++;
	  }
	  
	  
	  flag_dead = false;
	  if (dead_v_index.find(cloud.pts.at(j).index_v)!=dead_v_index.end()){
	    if (cloud.pts.at(j).x >= dead_v_index[cloud.pts.at(j).index_v].first &&
		cloud.pts.at(j).x <= dead_v_index[cloud.pts.at(j).index_v].second){
	      flag_dead = true;
	    }
	  }
	  
	  if (!flag_dead){
	    std::vector<std::tuple<double, PR3DCluster*, size_t>> results = global_skeleton_cloud.get_2d_points_info(test_point, loose_dis_cut, 1);
	    bool flag_unique = true;
	    if (results.size()>0){
	      std::set<PR3DCluster*> temp_clusters;
	      for (size_t k = 0; k!= results.size(); k++){
		if (std::get<0>(results.at(k)) < global_cloud.pts.at(std::get<2>(results.at(k))).index_v){
		  flag_unique = false;
		  temp_clusters.insert(std::get<1>(results.at(k)));
		}
	      }
	      for (auto it = temp_clusters.begin(); it!= temp_clusters.end(); it++){
		if (map_cluster_num[1].find(*it)==map_cluster_num[1].end()){
		  map_cluster_num[1][*it] = 1;
		}else{
		  map_cluster_num[1][*it] ++;
		}
	      }
	    }
	    if (flag_unique)
	      num_unique[1]++;
	  }else{
	    std::vector<std::tuple<double, PR3DCluster*, size_t>> results = global_skeleton_cloud.get_2d_points_info(test_point, loose_dis_cut, 1);
	    bool flag_unique = true;
	    if (results.size()>0){
	      std::set<PR3DCluster*> temp_clusters;
	      for (size_t k = 0; k!= results.size(); k++){
		if (std::get<0>(results.at(k)) < loose_dis_cut/3.*2.){
		  flag_unique = false;
		  temp_clusters.insert(std::get<1>(results.at(k)));
		}
	      }
	      for (auto it = temp_clusters.begin(); it!= temp_clusters.end(); it++){
		if (map_cluster_num[1].find(*it)==map_cluster_num[1].end()){
		  map_cluster_num[1][*it] = 1;
		}else{
		  map_cluster_num[1][*it] ++;
		}
	      }
	    }
	    if (flag_unique)
	      num_unique[1]++;
	  }
	  
	  
	  
	  
	  flag_dead = false;
	  if (dead_w_index.find(cloud.pts.at(j).index_w)!=dead_w_index.end()){
	    if (cloud.pts.at(j).x >= dead_w_index[cloud.pts.at(j).index_w].first &&
		cloud.pts.at(j).x <= dead_w_index[cloud.pts.at(j).index_w].second){
	      flag_dead = true;
	    }
	  }
	  
	  if (!flag_dead){
	    std::vector<std::tuple<double, PR3DCluster*, size_t>> results = global_skeleton_cloud.get_2d_points_info(test_point, loose_dis_cut, 2);
	    bool flag_unique = true;
	    if (results.size()>0){
	      std::set<PR3DCluster*> temp_clusters;
	      for (size_t k = 0; k!= results.size(); k++){
		if (std::get<0>(results.at(k)) < global_cloud.pts.at(std::get<2>(results.at(k))).index_w){
		  flag_unique = false;
		  temp_clusters.insert(std::get<1>(results.at(k)));
		}
	      }
	      for (auto it = temp_clusters.begin(); it!= temp_clusters.end(); it++){
		if (map_cluster_num[2].find(*it)==map_cluster_num[2].end()){
		  map_cluster_num[2][*it] = 1;
		}else{
		  map_cluster_num[2][*it] ++;
		}
	      }
	    }
	    if (flag_unique)
	      num_unique[2]++;
	  }else{
	    std::vector<std::tuple<double, PR3DCluster*, size_t>> results = global_skeleton_cloud.get_2d_points_info(test_point, loose_dis_cut, 2);
	    bool flag_unique = true;
	    if (results.size()>0){
	      std::set<PR3DCluster*> temp_clusters;
	      for (size_t k = 0; k!= results.size(); k++){
		if (std::get<0>(results.at(k)) < loose_dis_cut/3.*2.){
		  flag_unique = false;
		  temp_clusters.insert(std::get<1>(results.at(k)));
		}
	      }
	      for (auto it = temp_clusters.begin(); it!= temp_clusters.end(); it++){
		if (map_cluster_num[2].find(*it)==map_cluster_num[2].end()){
		  map_cluster_num[2][*it] = 1;
		}else{
		  map_cluster_num[2][*it] ++;
		}
	      }
	    }
	    if (flag_unique)
	      num_unique[2]++;
	  }
	  
	  
	}
	//      PR3DCluster *curr_cluster = cluster;
	
	
	bool flag_merge = false;
	
	
	/* if (fabs(extreme_points.first.x-140.7*units::cm)< 15*units::cm && */
	/* 	  fabs(extreme_points.first.y+72.3*units::cm)< 15*units::cm && */
	/* 	  fabs(extreme_points.first.z-30.1*units::cm)< 15*units::cm || cluster->get_cluster_id()==155) */
	/* 	std::cout << cluster->get_cluster_id()  << " A " */
	/* 		  << num_unique[0] << " " << num_unique[1] << " " << num_unique[2] << " " << num_total_points << " " << extreme_points.first.x/units::cm << " " << extreme_points.first.y/units::cm << " " << extreme_points.first.z/units::cm << " " << extreme_points.second.x/units::cm << " " << extreme_points.second.y/units::cm << " " << extreme_points.second.z/units::cm << " " */
	/* 		  << dir1.X() << " " << dir1.Y() << " " << dir1.Z() << " " << dir2.X() << " " << dir2.Y() << " " << dir2.Z() << " " << cluster_length_map[cluster] /units::cm << std::endl; */
	
	
	
	{
	  PR3DCluster *max_cluster_u = 0, *max_cluster_v=0, *max_cluster_w=0;
	  int max_value_u[3] = {0,0,0};
	  int max_value_v[3] = {0,0,0};
	  int max_value_w[3] = {0,0,0};
	  
	  int max_value[3]={0,0,0};
	  PR3DCluster *max_cluster=0;
	  
	  for (auto it = map_cluster_num[0].begin(); it!=map_cluster_num[0].end(); it++){
	    if (it->second > max_value_u[0]){
	      max_value_u[0] = it->second;
	      max_cluster_u = it->first;
	      
	      if (map_cluster_num[1].find(max_cluster_u)!=map_cluster_num[1].end()){
		max_value_u[1] = map_cluster_num[1][max_cluster_u];
	      }else{
		max_value_u[1] = 0;
	      }
	      
	      if (map_cluster_num[2].find(max_cluster_u)!=map_cluster_num[2].end()){
		max_value_u[2] = map_cluster_num[2][max_cluster_u];
	      }else{
		max_value_u[2] = 0;
	      }
	      
	    }
	  }
	  for (auto it = map_cluster_num[1].begin(); it!=map_cluster_num[1].end(); it++){
	    if (it->second > max_value_v[1]){
	      max_value_v[1] = it->second;
	      max_cluster_v = it->first;
	      
	      if (map_cluster_num[0].find(max_cluster_v)!=map_cluster_num[0].end()){
		max_value_v[0] = map_cluster_num[0][max_cluster_v];
	      }else{
		max_value_v[0] = 0;
	      }
	      if (map_cluster_num[2].find(max_cluster_v)!=map_cluster_num[2].end()){
		max_value_v[2] = map_cluster_num[2][max_cluster_v];
	      }else{
		max_value_v[2] = 0;
	      }
	    }
	  }
	  for (auto it = map_cluster_num[2].begin(); it!=map_cluster_num[2].end(); it++){
	    if (it->second > max_value_w[2]){
	      max_value_w[2] = it->second;
	      max_cluster_w = it->first;
	      
	      if (map_cluster_num[1].find(max_cluster_w)!=map_cluster_num[1].end()){
		max_value_w[1] = map_cluster_num[1][max_cluster_w];
	      }else{
		max_value_w[1] = 0;
	      }
	      if (map_cluster_num[0].find(max_cluster_w)!=map_cluster_num[0].end()){
		max_value_w[0] = map_cluster_num[0][max_cluster_w];
	      }else{
		max_value_w[0] = 0;
	      }
	      
	    }
	  }

	  /* std::cout << max_value_u[0] << " " << max_value_u[1] << " " << max_value_u[2] << " " */
	  /* 	    << max_value_v[0] << " " << max_value_v[1] << " " << max_value_v[2] << " " */
	  /* 	    << max_value_w[0] << " " << max_value_w[1] << " " << max_value_w[2] << " " << std::endl; */
	  
	  if ((max_value_u[0] > 0.33 * num_total_points || max_value_u[0] > 100) &&
	      (max_value_u[1] > 0.33 * num_total_points || max_value_u[1] > 100) &&
	      (max_value_u[2] > 0.33 * num_total_points || max_value_u[2] > 100)){
	    if (max_value_u[0] + max_value_u[1] + max_value_u[2] > max_value[0] + max_value[1] + max_value[2]){
	      max_value[0] = max_value_u[0];
	      max_value[1] = max_value_u[1];
	      max_value[2] = max_value_u[2];
	      max_cluster = max_cluster_u;
	    }
	  }
	  if ((max_value_v[0] > 0.33 * num_total_points || max_value_v[0] > 100) &&
	      (max_value_v[1] > 0.33 * num_total_points || max_value_v[1] > 100) &&
	      (max_value_v[2] > 0.33 * num_total_points || max_value_v[2] > 100)){
	    if (max_value_v[0] + max_value_v[1] + max_value_v[2] > max_value[0] + max_value[1] + max_value[2]){
	      max_value[0] = max_value_v[0];
	      max_value[1] = max_value_v[1];
	      max_value[2] = max_value_v[2];
	      max_cluster = max_cluster_v;
	    }
	  }
	  if ((max_value_w[0] > 0.33 * num_total_points || max_value_w[0] > 100) &&
	      (max_value_w[1] > 0.33 * num_total_points || max_value_w[1] > 100) &&
	      (max_value_w[2] > 0.33 * num_total_points || max_value_w[2] > 100)){
	    if (max_value_w[0] + max_value_w[1] + max_value_w[2] > max_value[0] + max_value[1] + max_value[2]){
	      max_value[0] = max_value_w[0];
	      max_value[1] = max_value_w[1];
	      max_value[2] = max_value_w[2];
	      max_cluster = max_cluster_w;
	    }
	  }


	  /* if (extreme_points.first.x < 117*units::cm && extreme_points.first.x > 105*units::cm && */
	  /*     extreme_points.first.z < 280*units::cm){ */
	  /*   std::cout << extreme_points.first.x /units::cm << " " << extreme_points.first.y/units::cm << " " << extreme_points.first.z/units::cm << " " << cluster->get_cluster_id() << " B " << max_value[0] << " " << max_value[1] << " " << max_value[2] << " " << num_total_points << " " << cluster_length_map[cluster]/units::cm << std::endl; */
	  /* } */
	  
	 
	  // if overlap a lot merge
	  if ((max_value[0]+max_value[1]+max_value[2]) > 0.75 *(num_total_points  + num_total_points  + num_total_points) &&
	      ((num_unique[1]+num_unique[0]+num_unique[2]) < 0.24 * num_total_points ||
	       ((num_unique[1]+num_unique[0]+num_unique[2]) < 0.45 * num_total_points &&
		(num_unique[1]+num_unique[0]+num_unique[2]) < 25))){

	    if (fabs(dir1.Angle(map_cluster_dir1[max_cluster])-3.1415926/2.) >= 70*3.1415926/180. ||
		fabs(dir1.Angle(map_cluster_dir2[max_cluster])-3.1415926/2.) >= 70*3.1415926/180. ||
		fabs(dir2.Angle(map_cluster_dir1[max_cluster])-3.1415926/2.) >= 70*3.1415926/180. ||
		fabs(dir2.Angle(map_cluster_dir2[max_cluster])-3.1415926/2.) >= 70*3.1415926/180.){
	      flag_merge = true;
	      to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster)); 
	    //curr_cluster = max_cluster;
	    }
	    
	    if (fabs(dir1.Angle(map_cluster_dir1[max_cluster])-3.1415926/2.) < 75*3.1415926/180. &&
		fabs(dir1.Angle(map_cluster_dir2[max_cluster])-3.1415926/2.) < 75*3.1415926/180. ){
	      flag_add_dir1 = false;
	    }else{
	      flag_add_dir1 = true;
	    }
	    if (fabs(dir2.Angle(map_cluster_dir1[max_cluster])-3.1415926/2.) < 75*3.1415926/180. &&
		fabs(dir2.Angle(map_cluster_dir2[max_cluster])-3.1415926/2.) < 75*3.1415926/180. ){
	      flag_add_dir2 = false;
	    }else{
	      flag_add_dir2 = true;
	    }
	    
	    /* max_cluster->Calc_PCA(); */
	    /* TVector3 p2_dir(max_cluster->get_PCA_axis(0).x, max_cluster->get_PCA_axis(0).y, max_cluster->get_PCA_axis(0).z); */
	    /* if (fabs(p2_dir.Angle(dir1)-3.1415926/2.) < 75 || fabs(p2_dir.Angle(dir2)-3.1415926/2.) < 75) */
	    /*   flag_add = false; */
	    /* if ( flag_merge ) */
	    /*   std::cout << extreme_points.first.x /units::cm << " " << extreme_points.first.y/units::cm << " " << extreme_points.first.z/units::cm << " " << cluster->get_cluster_id() << " B " << max_cluster->get_cluster_id() << " " << cluster_length_map[cluster]/units::cm << " " << cluster_length_map[max_cluster]/units::cm << " " << max_value[0] << " " << max_value[1] << " " << max_value[2] << " " << num_total_points << " " << num_unique[0] << " " << num_unique[1] << " " << num_unique[2] << std::endl; */
	    
	  }

	  if ((max_value[0]+max_value[1]+max_value[2]) > 300 && !flag_merge){
	    if (cluster_length_map[cluster]> 25*units::cm || cluster_length_map[max_cluster]> 25*units::cm){
	      // if overlap significant, compare the PCA
	      cluster->Calc_PCA();
	      Point p1_c = cluster->get_center();
	      TVector3 p1_dir(cluster->get_PCA_axis(0).x, cluster->get_PCA_axis(0).y, cluster->get_PCA_axis(0).z);
	      max_cluster->Calc_PCA();
	      Point p2_c = max_cluster->get_center();
	      TVector3 p2_dir(max_cluster->get_PCA_axis(0).x, max_cluster->get_PCA_axis(0).y, max_cluster->get_PCA_axis(0).z);
	      
	      double angle_diff = p1_dir.Angle(p2_dir)/3.1415926*180.;
	      double angle1_drift = p1_dir.Angle(drift_dir)/3.1415926*180.;
	      double angle2_drift = p2_dir.Angle(drift_dir)/3.1415926*180.;
	      Line l1(p1_c,p1_dir);
	      Line l2(p2_c,p2_dir);
	      double dis = l1.closest_dis(l2);
	      double dis1 = sqrt(pow(p1_c.x - p2_c.x,2) + pow(p1_c.y - p2_c.y,2) + pow(p1_c.z - p2_c.z,2));
	      
	    
	      if ( (angle_diff < 5 || angle_diff > 175 || fabs(angle1_drift-90) < 5 && fabs(angle2_drift-90) < 5 && fabs(angle1_drift-90)+fabs(angle2_drift-90) < 6 && (angle_diff < 30 || angle_diff > 150)) && dis < 1.5*units::cm ||
		   (angle_diff < 10 || angle_diff > 170) && dis < 0.9*units::cm &&
		   dis1 > (cluster_length_map[cluster] + cluster_length_map[max_cluster])/3.){
		to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster));
		//curr_cluster = max_cluster;
		flag_merge = true;
	      }else if (((angle_diff < 5 || angle_diff > 175) && dis < 2.5*units::cm ||
			 (angle_diff < 10 || angle_diff > 170) && dis < 1.2*units::cm) &&
			dis1 > (cluster_length_map[cluster] + cluster_length_map[max_cluster])/3.){
		to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster));
		flag_merge = true;
	      }

	      if (fabs(dir2.Angle(drift_dir) - 3.1415926/2.) < 5*3.1415926/180. && fabs(dir1.Angle(drift_dir) - 3.1415926/2.) < 5*3.1415926/180. && (max_value[0]+max_value[1]+max_value[2]) > 0.7 *(num_total_points  + num_total_points  + num_total_points)){
		to_be_merged_pairs.insert(std::make_pair(cluster,max_cluster));
		flag_merge = true;
	      }

	      std::cout << fabs(dir2.Angle(drift_dir) - 3.1415926/2.) /3.1415928*180. << " " << fabs(dir1.Angle(drift_dir) - 3.1415926/2.)/3.1415926*180. << " " << (max_value[0]+max_value[1]+max_value[2]) << " " << num_total_points << std::endl;
	      
	      //std::cout << extreme_points.first.x /units::cm << " " << extreme_points.first.y/units::cm << " " << extreme_points.first.z/units::cm << " " << cluster->get_cluster_id() << " " << max_cluster->get_cluster_id() << " " << cluster_length_map[cluster]/units::cm << " " << cluster_length_map[max_cluster]/units::cm << " " << angle_diff << " " << angle1_drift << " " << angle2_drift << " " << dis/units::cm << " " << dis1/units::cm << " " << flag_merge << std::endl;
	    }

	   

	  }

	  
	 
	  
	  
	  


	  
	  /* if (fabs(extreme_points.first.x-228.6*units::cm)< 5*units::cm && */
	  /*     fabs(extreme_points.first.y+5.2*units::cm)< 5*units::cm && */
	  /*     fabs(extreme_points.first.z-448.4*units::cm)< 5*units::cm ){ */
	  /*   // if (max_cluster!=0) */
	  /*   std::cout << cluster->get_cluster_id() << " " << max_cluster << " " << max_value[0] << " " << max_value[1] << " " << max_value[2] << " " << num_total_points << std::endl; */
	  /*   std::cout << max_value_u[0] << " " << max_value_u[1] << " " << max_value_u[2] << " " << max_cluster_u->get_cluster_id() << " "<< max_cluster_u << " " << max_value_v[0] << " " << max_value_v[1] << " " << max_value_v[2] << " " << max_cluster_v->get_cluster_id() << " " << max_cluster_v << " " << max_value_w[0] << " " << max_value_w[1] << " " << max_value_w[2] << " " << max_cluster_w->get_cluster_id() << " " << max_cluster_w << std::endl; */
	  /*   std::cout << map_cluster_num[0][max_cluster_w] << " " << map_cluster_num[1][max_cluster_w] << " " << map_cluster_num[2][max_cluster_w] << " " << std::endl; */
	  /* } */
	}


	// when added points in
	// if overlap a lot merge 
	// if overlap significant, compare the PCA
	

	
	  
	  //if (max_value[0] + max_value[1] + max_value[2]>0)
	/* if (extreme_points.first.z <=1020*units::cm && extreme_points.first.z >= 960*units::cm && */
	/*     extreme_points.first.y < 0) */
	/*   std::cout << extreme_points.first.x/units::cm << " " << extreme_points.first.y/units::cm << " " << extreme_points.first.z/units::cm << " " << max_value[0] << " " << max_value[1] << " " << max_value[2] << " " << num_total_points << " " << num_unique[0] << " " << num_unique[1] << " " << num_unique[2] << " " << flag_merge << " " << flag_add_dir1 << " " << flag_add_dir2 << std::endl; */
	
	
	
      } // length cut ... 
     
      if (flag_add_dir1){
	// add extension points in ... 
	if (fabs(dir1.Angle(drift_dir) - 3.1415926/2.) < 5*3.1415926/180.){
	  global_skeleton_cloud.AddPoints(cluster,extreme_points.first,dir1,extending_dis*2,1.2*units::cm, angle/2.);
	  dir1 *= -1;
	  global_skeleton_cloud.AddPoints(cluster,extreme_points.first,dir1,extending_dis*2,1.2*units::cm, angle/2.);
	}else{
	  global_skeleton_cloud.AddPoints(cluster,extreme_points.first,dir1,extending_dis,1.2*units::cm, angle);
	  dir1 *= -1;
	  global_skeleton_cloud.AddPoints(cluster,extreme_points.first,dir1,extending_dis,1.2*units::cm, angle);
	}
      }
      
      if (flag_add_dir2){
	if (fabs(dir2.Angle(drift_dir) - 3.1415926/2.) < 5*3.1415926/180.){
	  global_skeleton_cloud.AddPoints(cluster,extreme_points.second,dir2,extending_dis*2.0,1.2*units::cm, angle/2.);
	  dir2 *= -1;
	  global_skeleton_cloud.AddPoints(cluster,extreme_points.second,dir2,extending_dis*2.0,1.2*units::cm, angle/2.);
	}else{
	  global_skeleton_cloud.AddPoints(cluster,extreme_points.second,dir2,extending_dis,1.2*units::cm, angle);
	  dir2 *= -1;
	  global_skeleton_cloud.AddPoints(cluster,extreme_points.second,dir2,extending_dis,1.2*units::cm, angle);
	}
      }
    } // not the first cluster ... 


    
    //std::cout << extreme_points.first.x/units::cm << " " << extreme_points.first.y/units::cm << " " << extreme_points.first.z/units::cm << " " << extreme_points.second.x/units::cm << " " << extreme_points.second.y/units::cm << " " << extreme_points.second.z/units::cm << std::endl;
  } // loop over clusters ... 


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



  WireCell::PR3DClusterSelection new_clusters;
  
  // merge clusters into new clusters, delete old clusters 
  for (auto it = merge_clusters.begin(); it!=merge_clusters.end();it++){
    std::set<PR3DCluster*>& clusters = (*it);
    PR3DCluster *ncluster = new PR3DCluster((*clusters.begin())->get_cluster_id());
    live_clusters.push_back(ncluster);
    
    new_clusters.push_back(ncluster);

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


  to_be_merged_pairs.clear();
  for (auto it = new_clusters.begin(); it!= new_clusters.end(); it++){
    PR3DCluster *cluster_1 = (*it);
    cluster_1->Calc_PCA();
    Point p1_c = cluster_1->get_center();
    TVector3 p1_dir(cluster_1->get_PCA_axis(0).x, cluster_1->get_PCA_axis(0).y, cluster_1->get_PCA_axis(0).z);
    Line l1(p1_c,p1_dir);
    for (auto it1 = live_clusters.begin(); it1 != live_clusters.end(); it1++){
      PR3DCluster *cluster_2 = (*it1);
      if (cluster_length_map[cluster_2] < 3*units::cm) continue;
      if (cluster_2 == cluster_1) continue;

      if (cluster_length_map[cluster_1]> 25*units::cm || cluster_length_map[cluster_2]> 25*units::cm ||
	  (cluster_length_map[cluster_1]+ cluster_length_map[cluster_2]) > 30*units::cm){
	cluster_2->Calc_PCA();
	Point p2_c = cluster_2->get_center();
	TVector3 p2_dir(cluster_2->get_PCA_axis(0).x, cluster_2->get_PCA_axis(0).y, cluster_2->get_PCA_axis(0).z);

	TVector3 cc_dir(p2_c.x - p1_c.x, p2_c.y-p1_c.y, p2_c.z - p1_c.z);
	
	double angle_diff = fabs(p1_dir.Angle(p2_dir)-3.1415926/2.)/3.1415926*180.;
	double angle_diff1 = fabs(cc_dir.Angle(p1_dir)-3.1415926/2.)/3.1415926*180;
	double angle_diff2 = fabs(cc_dir.Angle(p2_dir)-3.1415926/2.)/3.1415926*180;
	  
	Line l2(p2_c,p2_dir);
	double dis = l1.closest_dis(l2);
	
	double dis1 = sqrt(pow(p1_c.x - p2_c.x,2) + pow(p1_c.y - p2_c.y,2) + pow(p1_c.z - p2_c.z,2));
	
	if (p1_dir.Mag()!=0) p1_dir.SetMag(1);
	if (p2_dir.Mag()!=0) p2_dir.SetMag(1);

	//if (cluster_2->get_cluster_id()==431 || cluster_1->get_cluster_id()==431)
	/* if ((fabs(p2_c.z/units::cm-420) <20 && fabs(p2_c.x/units::cm-250)<20 || */
	/*      fabs(p1_c.z/units::cm-420) <20 && fabs(p1_c.x/units::cm-250)<20 */
	/*      ) && cluster_length_map[cluster_1]/units::cm> 5 && cluster_length_map[cluster_2]/units::cm > 5)  */
	
	bool flag_merge = false;
	
	if (((angle_diff >85) && (angle_diff1 > 90 - 1.5 * (90-angle_diff)) &&
	     (angle_diff2 > 90 - 1.5 * (90-angle_diff)) && dis < 2.5*units::cm ||
	     (angle_diff >80) && angle_diff1 > 80 && angle_diff2 > 80 && dis < 1.2*units::cm) &&
	    dis1 > (cluster_length_map[cluster_2] + cluster_length_map[cluster_1])/3. ){
	  to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	  flag_merge = true;
	}else if ((angle_diff >87) && (angle_diff1 > 90 - 1.5 * (90-angle_diff) ) &&
		  (angle_diff2 > 90 - 1.5 * (90-angle_diff) ) && dis < 4.0*units::cm &&
		  dis1 > (cluster_length_map[cluster_2] + cluster_length_map[cluster_1])/2. &&
		  cluster_length_map[cluster_2] > 15*units::cm && cluster_length_map[cluster_1] > 15*units::cm &&
		  cluster_length_map[cluster_2]+cluster_length_map[cluster_1] > 45*units::cm){
	  to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	  flag_merge = true;
	}

	/* if (flag_merge) */
	/*   std::cout << cluster_2->get_cluster_id() << " " << cluster_1->get_cluster_id() << " " << angle_diff << " " << angle_diff1 << " " << dis/units::cm << " " << dis1/units::cm <<  " " << cluster_length_map[cluster_1]/units::cm << " " << cluster_length_map[cluster_2]/units::cm << " " << p1_dir.X() << " " << p1_dir.Y() << " " <<  p1_dir.Z() << " " << p2_dir.X() << " " << p2_dir.Y() << " " << p2_dir.Z() << " " << p1_c.x/units::cm << " " << p1_c.y/units::cm << " " << p1_c.z/units::cm << " " << p2_c.x/units::cm << " " << p2_c.y/units::cm << " " << p2_c.z/units::cm << std::endl;  */
	
      }
    }
  }


  merge_clusters.clear();
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
