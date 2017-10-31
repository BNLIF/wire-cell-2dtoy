#include "WireCell2dToy/ToyClustering.h"

using namespace WireCell;

std::map<PR3DCluster*,std::vector<std::pair<PR3DCluster*,double>>> WireCell2dToy::Clustering_isolated(WireCell::PR3DClusterSelection& live_clusters){
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
    if (max < 40){
      small_clusters.push_back(live_clusters.at(i));
    }else{
      big_clusters.push_back(live_clusters.at(i));
    }
  }

  std::map<PR3DCluster*,std::vector<std::pair<PR3DCluster*,double>>> results;
  // Now for each small clusters, find the big clusters with the shortest distance and connect them
  for (auto it = small_clusters.begin(); it!= small_clusters.end(); it++){
    PR3DCluster* curr_cluster = (*it);
    curr_cluster->Create_point_cloud();
    double min_dis = 1e9; PR3DCluster* min_dis_cluster = 0;
    for (auto it1 = big_clusters.begin(); it1!=big_clusters.end(); it1++){
      PR3DCluster* big_cluster = (*it1);
      big_cluster->Create_point_cloud();
  
      // pick any point and merged cell in cluster1,
      SlimMergeGeomCell *prev_mcell1 = 0;
      SlimMergeGeomCell *prev_mcell2 = 0;
      SlimMergeGeomCell *mcell1 = curr_cluster->get_mcells().at(0);
      Point p1 = mcell1->center();
      SlimMergeGeomCell *mcell2=0;
      Point p2;

      while(mcell1!=prev_mcell1 || mcell2!=prev_mcell2){
	prev_mcell1 = mcell1;
	prev_mcell2 = mcell2;
    
	// find the closest point and merged cell in cluster2
	std::pair<SlimMergeGeomCell*,Point> temp_results = big_cluster->get_closest_point_mcell(p1);
	p2 = temp_results.second;
	mcell2 = temp_results.first;
	// find the closest point and merged cell in cluster1
	temp_results = curr_cluster->get_closest_point_mcell(p2);
	p1 = temp_results.second;
	mcell1 = temp_results.first;
      }
      double dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
      if (dis < min_dis){
	min_dis = dis;
	min_dis_cluster = big_cluster;
      }
    }
    if (min_dis < 0.8*units::m){
      if (results.find(min_dis_cluster)==results.end()){
	std::vector<std::pair<PR3DCluster*,double>> temp;
	temp.push_back(std::make_pair(curr_cluster,min_dis));
	results[min_dis_cluster]=temp;
      }else{
	results[min_dis_cluster].push_back(std::make_pair(curr_cluster,min_dis));
      }
    }else{
      big_clusters.push_back(curr_cluster);
      std::vector<std::pair<PR3DCluster*,double>> temp;
      results[curr_cluster] = temp;
    }
  }
  
  return results;
}

void WireCell2dToy::Clustering_jump_gap_cosmics(WireCell::PR3DClusterSelection& live_clusters){

  for (int kk=0;kk!=1;kk++){
    std::set<std::pair<PR3DCluster*, PR3DCluster*>> to_be_merged_pairs;
    
    for (size_t i=0;i!=live_clusters.size();i++){
      PR3DCluster* cluster_1 = live_clusters.at(i);
      for (size_t j=i+1;j<live_clusters.size();j++){
	PR3DCluster* cluster_2 = live_clusters.at(j);
	//std::cout << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << std::endl;
	if (WireCell2dToy::Clustering_jump_gap_cosmics(cluster_1,cluster_2))
	  to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	
      }
    }
    
    //to_be_merged_pairs.clear();
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
	delete ocluster;
      }
      //std::cout << std::endl;
    }
    
  }  
}
  
bool WireCell2dToy::Clustering_jump_gap_cosmics(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2){
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
    
    //  std::cout << mcell1 << " " << mcell2 << " " << sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2))/units::cm<< std::endl;
  }
  
  //std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2))/units::cm<< std::endl;


  
  double dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));

  // if (cluster1->get_cluster_id()==106 || cluster1->get_cluster_id()==143 ||
  //     cluster1->get_cluster_id()==214 || cluster1->get_cluster_id()==251){
  //   std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << dis/units::cm << " " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm <<
  //     " " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm <<
  //     std::endl;
  // }
  
  
  if (dis < 45*units::cm){
    Point mcell1_center = mcell1->center();
    Point mcell2_center = mcell2->center();
    Point cluster1_ave_pos = cluster1->calc_ave_pos(p1,5*units::cm);
    Point cluster2_ave_pos = cluster2->calc_ave_pos(p2,5*units::cm);

    bool flag_para = false;
    bool flag_perp = false;
    bool flag_prolonged_U = false;
    bool flag_prolonged_V = false;
    // parallel case 1 and perpendicular case 2 
    TVector3 drift_dir(1,0,0);
    {
      TVector3 tempV1(mcell2_center.x - mcell1_center.x, mcell2_center.y - mcell1_center.y, mcell2_center.z - mcell1_center.z);
      TVector3 tempV2(cluster2_ave_pos.x - mcell1_center.x, cluster2_ave_pos.y - mcell1_center.y, cluster2_ave_pos.z - mcell1_center.z);
      TVector3 tempV3(mcell2_center.x - cluster1_ave_pos.x, mcell2_center.y - cluster1_ave_pos.y, mcell2_center.z - cluster1_ave_pos.z);
      TVector3 tempV4(cluster2_ave_pos.x - cluster1_ave_pos.x, cluster2_ave_pos.y - cluster1_ave_pos.y, cluster2_ave_pos.z - cluster1_ave_pos.z);
      
      double angle1 = tempV1.Angle(drift_dir);
      double angle2 = tempV2.Angle(drift_dir);
      double angle3 = tempV3.Angle(drift_dir);
      double angle4 = tempV4.Angle(drift_dir);
      if (angle1<10/180.*3.1415926 || (3.1415926-angle1)<10/180.*3.1415926
	  ||angle2<10/180.*3.1415926 || (3.1415926-angle2)<10/180.*3.1415926
	  ||angle3<10/180.*3.1415926 || (3.1415926-angle3)<10/180.*3.1415926
	  ||angle4<10/180.*3.1415926 || (3.1415926-angle4)<10/180.*3.1415926)
	flag_perp = true;

      if (fabs(angle1-3.1415926/2.)<5/180.*3.1415926
	  ||fabs(angle2-3.1415926/2.)<5/180.*3.1415926
	  ||fabs(angle3-3.1415926/2.)<5/180.*3.1415926
	  ||fabs(angle4-3.1415926/2.)<5/180.*3.1415926)
	flag_para = true;
    }
    
    // pronlonged case for U 3 and V 4 ...
    TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
    TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));

    {
      TVector3 tempV1(0, mcell2_center.y - mcell1_center.y, mcell2_center.z - mcell1_center.z);
      TVector3 tempV2(0, cluster2_ave_pos.y - mcell1_center.y, cluster2_ave_pos.z - mcell1_center.z);
      TVector3 tempV3(0, mcell2_center.y - cluster1_ave_pos.y, mcell2_center.z - cluster1_ave_pos.z);
      TVector3 tempV4(0, cluster2_ave_pos.y - cluster1_ave_pos.y, cluster2_ave_pos.z - cluster1_ave_pos.z);

      double angle1 = tempV1.Angle(U_dir);
      double angle2 = tempV2.Angle(U_dir);
      double angle3 = tempV3.Angle(U_dir);
      double angle4 = tempV4.Angle(U_dir);
      if (angle1<5/180.*3.1415926 || (3.1415926-angle1)<5/180.*3.1415926
	  ||angle2<5/180.*3.1415926 || (3.1415926-angle2)<5/180.*3.1415926
	  ||angle3<5/180.*3.1415926 || (3.1415926-angle3)<5/180.*3.1415926
	  ||angle4<5/180.*3.1415926 || (3.1415926-angle4)<5/180.*3.1415926)
	flag_prolonged_U = true;
      
      angle1 = tempV1.Angle(V_dir);
      angle2 = tempV2.Angle(V_dir);
      angle3 = tempV3.Angle(V_dir);
      angle4 = tempV4.Angle(V_dir);
      if (angle1<5/180.*3.1415926 || (3.1415926-angle1)<5/180.*3.1415926
	  ||angle2<5/180.*3.1415926 || (3.1415926-angle2)<5/180.*3.1415926
	  ||angle3<5/180.*3.1415926 || (3.1415926-angle3)<5/180.*3.1415926
	  ||angle4<5/180.*3.1415926 || (3.1415926-angle4)<5/180.*3.1415926)
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
    
    if (flag_dir){
      if (cluster1->get_point_cloud()->get_num_points()>60 || cluster2->get_mcells().size()>3){
	// also with a complicated structures ...
	TVector3 dir1 = cluster1->VHoughTrans(cluster1_ave_pos,30*units::cm); // cluster 1 direction
	// need to be careful if this is reversed ...
	// TVector3 temp_dir1(cluster1_ave_pos.x-p1.x,cluster1_ave_pos.y-p1.y,cluster1_ave_pos.z-p1.z);
	// if (temp_dir1.Angle(dir1)>135/180.*3.1415926)
	//   dir1.SetXYZ(-dir1.X(),-dir1.Y(),-dir1.Z());
	
	TVector3 dir2 = cluster2->calc_dir(cluster1_ave_pos,cluster2_ave_pos,30*units::cm); // dir 2 --> 1
	TVector3 dir2_1 = cluster2->VHoughTrans(cluster1_ave_pos,30*units::cm);// 
	TVector3 dir2_2 = cluster2->calc_dir(cluster1_ave_pos,mcell2_center,10*units::cm); // dir 2 --> 1
	
	TVector3 dir1_rot(dir1.Y(), dir1.Z(), dir1.X());
	TVector3 dir2_rot(dir2.Y(), dir2.Z(), dir2.X());
	
	double angle_diff1 = (3.1415926-dir1.Angle(dir2))/3.1415926*180.;
	double angle_diff1_1 = (3.1415926-dir1.Angle(dir2_1))/3.1415926*180.;
	double angle_diff1_2 = (3.1415926-dir1.Angle(dir2_2))/3.1415926*180.;
	
	double theta1 = (dir1_rot.Theta()-3.1415926/2.)/3.1415926*180.;
	double theta2 = (dir2_rot.Theta()-3.1415926/2.)/3.1415926*180.;
	double dphi = fabs(3.1415926 - fabs(dir1_rot.Phi()-dir2_rot.Phi()))/3.1415926*180.;

	// if (cluster1->get_cluster_id()==142|| cluster2->get_cluster_id()==142) 
	//   std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << angle_diff1 << " " << angle_diff1_1 << " " << angle_diff1_2 << " " << theta1 << " " << theta2 << " " << dphi << " " << dis/units::cm << " " << cluster1->get_mcells().size() << " " << cluster2->get_mcells().size() << std::endl;
	
	if ((angle_diff1 < angle_cut) || angle_diff1_1 < angle_cut || (angle_diff1_2 <angle_cut )
	    || (flag_para && fabs(theta1) < para_angle_cut_1 && fabs(theta2) < para_angle_cut_1 && fabs(theta1+theta2) < para_angle_cut_1/2. && dphi < para_angle_cut)
	    //   || ((angle_diff1_1<point_angle_cut || angle_diff1_2 < point_angle_cut ) && dis < 3*units::cm)
	    ){
	  return true;
	}
      }
      // reverse the test ...
      if (cluster2->get_point_cloud()->get_num_points()>60 || cluster2->get_mcells().size()>3){
	TVector3 dir1 = cluster2->VHoughTrans(cluster2_ave_pos,30*units::cm); // cluster 1 direction
	// need to be careful if this is reversed ...
	// TVector3 temp_dir1(cluster2_ave_pos.x-p2.x,cluster2_ave_pos.y-p2.y,cluster2_ave_pos.z-p2.z);
	// if (temp_dir1.Angle(dir1)>135/180.*3.1415926)
	//   dir1.SetXYZ(-dir1.X(),-dir1.Y(),-dir1.Z());
	
	TVector3 dir2 = cluster1->calc_dir(cluster2_ave_pos,cluster1_ave_pos,30*units::cm); // dir 2 --> 1
	TVector3 dir2_1 = cluster1->VHoughTrans(cluster2_ave_pos,30*units::cm);//
	TVector3 dir2_2 = cluster1->calc_dir(cluster2_ave_pos,mcell1_center,10*units::cm); // dir 2 --> 1
	
	TVector3 dir1_rot(dir1.Y(), dir1.Z(), dir1.X());
	TVector3 dir2_rot(dir2.Y(), dir2.Z(), dir2.X());
	
	double angle_diff1 = (3.1415926-dir1.Angle(dir2))/3.1415926*180.;
	double angle_diff1_1 = (3.1415926-dir1.Angle(dir2_1))/3.1415926*180.;
	double angle_diff1_2 = (3.1415926-dir1.Angle(dir2_2))/3.1415926*180.;
	
	double theta1 = (dir1_rot.Theta()-3.1415926/2.)/3.1415926*180.;
	double theta2 = (dir2_rot.Theta()-3.1415926/2.)/3.1415926*180.;
	double dphi = fabs(3.1415926 - fabs(dir1_rot.Phi()-dir2_rot.Phi()))/3.1415926*180.;

	// if (cluster1->get_cluster_id()==142 || cluster2->get_cluster_id()==142 ) 
	//   std::cout << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << " " << angle_diff1 << " " << angle_diff1_1 << " " << angle_diff1_2 << " " << theta1 << " " << theta2 << " " << dphi << " " << dis/units::cm << " " << cluster1->get_mcells().size() << " " << cluster2->get_mcells().size() << std::endl;
	
	
	if ((angle_diff1 < angle_cut ) || angle_diff1_1 < angle_cut || (angle_diff1_2 < angle_cut )
	    || (flag_para && fabs(theta1) < para_angle_cut_1 && fabs(theta2) < para_angle_cut_1 && fabs(theta1+theta2) < para_angle_cut_1/2. && dphi < para_angle_cut)
	    //  || ((angle_diff1_1<point_angle_cut || angle_diff1_2 < point_angle_cut) && dis < 3*units::cm)
	    ){
	  return true;
	}
      }

    }
  } // dis cut
  return false;
}


void WireCell2dToy::Clustering_live_dead(WireCell::PR3DClusterSelection& live_clusters, WireCell::PR3DClusterSelection& dead_clusters){

   // dead to live clusters mapping ... 
   std::map<PR3DCluster*,std::vector<PR3DCluster*>> dead_live_cluster_mapping;
   std::map<PR3DCluster*,std::vector<SMGCSelection>> dead_live_mcells_mapping;
   
   // form map between live and dead clusters ... 
   for (size_t i = 0; i!=live_clusters.size(); i++){
     for (size_t j = 0; j!= dead_clusters.size(); j++){
       SMGCSelection mcells = (live_clusters.at(i))->Is_Connected(dead_clusters.at(j),2);
       int live_cluster_id = live_clusters.at(i)->get_cluster_id();
       int dead_cluster_id = dead_clusters.at(j)->get_cluster_id();
       // if (live_cluster_id==5 || live_cluster_id == 21 || live_cluster_id == 76){
       // 	 if (dead_cluster_id==8 || dead_cluster_id==100)
       // 	   std::cout << live_cluster_id << " " << dead_cluster_id << " " << mcells.size() << std::endl;
       // }
       if ( mcells.size()>0 ){
	 if (dead_live_cluster_mapping.find(dead_clusters.at(j))==dead_live_cluster_mapping.end() ){
	   std::vector<PR3DCluster*> temp_clusters;
	   temp_clusters.push_back(live_clusters.at(i));
	   dead_live_cluster_mapping[dead_clusters.at(j)] = temp_clusters;
	   std::vector<SMGCSelection> temp_mcells;
	   temp_mcells.push_back(mcells);
	   dead_live_mcells_mapping[dead_clusters.at(j)] = temp_mcells;
	 }else{
	   dead_live_cluster_mapping[dead_clusters.at(j)].push_back(live_clusters.at(i));
	   dead_live_mcells_mapping[dead_clusters.at(j)].push_back(mcells);
	 }
       }
       // if (mcells.size()>0)
       // 	 std::cout << mcells.size() << " " <<
       // 	   live_clusters.at(i)->get_cluster_id() << " "
       // 		   << live_clusters.at(i)->get_num_mcells() << " "
       // 		   << live_clusters.at(i)->get_num_time_slices() << " " 
       // 		   << dead_clusters.at(j)->get_cluster_id() << " "
       // 		   << dead_clusters.at(j)->get_num_mcells() << " "
       // 		   << dead_clusters.at(j)->get_num_time_slices() << std::endl;
     }
   }
  
  std::set<std::pair<PR3DCluster*, PR3DCluster*>> tested_pairs;
  std::set<std::pair<PR3DCluster*, PR3DCluster*>> to_be_merged_pairs;

  
  for (auto it = dead_live_cluster_mapping.begin(); it!= dead_live_cluster_mapping.end(); it++){
    PR3DCluster* the_dead_cluster = (*it).first;
    std::vector<PR3DCluster*> connected_live_clusters = (*it).second;
    std::vector<SMGCSelection> connected_live_mcells = dead_live_mcells_mapping[the_dead_cluster];
    if (connected_live_clusters.size()>1){
      //      std::cout << the_dead_cluster->get_cluster_id() << " " << connected_live_clusters.size() << std::endl;

      // std::cout << "Candidates: ";
      // for (size_t i=0; i!= connected_live_clusters.size(); i++){
      //   PR3DCluster* cluster_1 = connected_live_clusters.at(i);
      // 	std::cout << cluster_1->get_cluster_id() << " ";
      // }
      // std::cout << std::endl;
      
      for (size_t i=0; i!= connected_live_clusters.size(); i++){
        PR3DCluster* cluster_1 = connected_live_clusters.at(i);
	SMGCSelection mcells_1 = connected_live_mcells.at(i);
	for (size_t j=i+1;j<connected_live_clusters.size(); j++){
	  PR3DCluster* cluster_2 = connected_live_clusters.at(j);
	  SMGCSelection mcells_2 = connected_live_mcells.at(j);	  
	  if (tested_pairs.find(std::make_pair(cluster_1,cluster_2))==tested_pairs.end()){
	    cluster_1->Create_point_cloud();
	    cluster_2->Create_point_cloud();

	    tested_pairs.insert(std::make_pair(cluster_1,cluster_2));
	    tested_pairs.insert(std::make_pair(cluster_2,cluster_1));
	    // starting the test ... 

	    bool flag_merge = false;
	    
	    // pick any point and merged cell in cluster1,
	    SlimMergeGeomCell *prev_mcell1 = 0;
	    SlimMergeGeomCell *prev_mcell2 = 0;
	    SlimMergeGeomCell *mcell1 = mcells_1.at(0);
	    Point p1 = mcell1->center();
	    SlimMergeGeomCell *mcell2=0;
	    Point p2;
	    
	    while(mcell1!=prev_mcell1 || mcell2!=prev_mcell2){
	      prev_mcell1 = mcell1;
	      prev_mcell2 = mcell2;
	      
	      // find the closest point and merged cell in cluster2
	      std::pair<SlimMergeGeomCell*,Point> temp_results = cluster_2->get_closest_point_mcell(p1);
	      p2 = temp_results.second;
	      mcell2 = temp_results.first;
	      // find the closest point and merged cell in cluster1
	      temp_results = cluster_1->get_closest_point_mcell(p2);
	      p1 = temp_results.second;
	      mcell1 = temp_results.first;
	      //  std::cout << mcell1 << " " << mcell2 << " " << sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2))/units::cm<< std::endl;
	    }

	    Point mcell1_center = cluster_1->calc_ave_pos(p1,5*units::cm);
	    TVector3 dir1 = cluster_1->VHoughTrans(mcell1_center,30*units::cm);
	    // test reversed direction
	    // TVector3 temp_dir1(mcell1_center.x-p1.x,mcell1_center.y-p1.y,mcell1_center.z-p1.z);
	    // if (temp_dir1.Angle(dir1)>135/180.*3.1415926)
	    //   dir1.SetXYZ(-dir1.X(),-dir1.Y(),-dir1.Z());
	    
	    Point mcell2_center = cluster_2->calc_ave_pos(p2,5*units::cm);
	    TVector3 dir2 = cluster_2->calc_dir(mcell1_center,mcell2_center,30*units::cm);
	    double angle_diff = (3.1415926-dir1.Angle(dir2))/3.1415926*180.;
	    double dis = sqrt(pow(p1.x-p2.x,2)
			      +pow(p1.y-p2.y,2)
			      +pow(p1.z-p2.z,2));
	    // if (angle_diff < 30)
	    //std::cout << "Xin1: " << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << " " << angle_diff << " " << dis/units::cm << " " << mcells_1.size() << " " << mcells_2.size() << " " << dis/units::cm << " " << fabs(p1.x-p2.x)/units::cm<<
	    // " " << p1.y/units::cm << " " << p1.z/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << std::endl;
	    if ( (dis <= 3*units::cm|| fabs(p1.x-p2.x)<0.5*units::cm)  && angle_diff <= 45 
		|| dis <= 10*units::cm && angle_diff <=15 
		|| angle_diff<5 
		//  || mcells_1.size()==1 && mcells_2.size()==1 && dis < 15*units::cm
		){
	      flag_merge = true;
	    }

	    if (!flag_merge){
	      dir1 = cluster_2->VHoughTrans(mcell2_center,30*units::cm);
	      // test reversed direction
	      // TVector3 temp_dir1(mcell2_center.x-p2.x,mcell2_center.y-p2.y,mcell2_center.z-p2.z);
	      // if (temp_dir1.Angle(dir1)>135/180.*3.1415926)
	      // 	dir1.SetXYZ(-dir1.X(),-dir1.Y(),-dir1.Z());
	      
	     
	      dir2 = cluster_1->calc_dir(mcell2_center,mcell1_center,30*units::cm);
	      angle_diff = (3.1415926-dir1.Angle(dir2))/3.1415926*180.;
	      dis = sqrt(pow(p1.x-p2.x,2)
			 +pow(p1.y-p2.y,2)
			 +pow(p1.z-p2.z,2));
	      // if (angle_diff < 30)
	      //std::cout << "Xin2: " << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << " " << angle_diff << " " << dis/units::cm << " " << mcells_1.size() << " " << mcells_2.size() << " " << dis/units::cm << " " << fabs(p1.x-p2.x)/units::cm << std::endl;
	      if ( (dis <= 3*units::cm || fabs(p1.x-p2.x)<0.5*units::cm) && angle_diff <= 45 
		  || dis <= 10*units::cm && angle_diff <=15 
		  || angle_diff<5  
		  //  || mcells_1.size()==1 && mcells_2.size()==1 && dis < 15*units::cm
		  ){
		flag_merge = true;
	      }
	    }
	    
	    
	    if (flag_merge)
	      to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	  }
	}
      }
      //  std::cout << std::endl;
    }
  }
  
  //to_be_merged_pairs.clear();
  
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
      delete ocluster;
    }
    //std::cout << std::endl;
  }
  
  //std::cout << merge_clusters.size() << std::endl;
  
  //std::cout << to_be_merged_pairs.size() << std::endl;
  // 
  //   std::cout << "Pair: " << (*it).first->get_cluster_id() << " " << (*it).second->get_cluster_id() << std::endl;
  // }


  
  
}


// bool WireCell2dToy::IsCrossing(WireCell::Point& p, double theta, double phi, WireCell::SlimMergeGeomCell *mcell, WireCellSst::GeomDataSource& gds){
//   // dx = sintheta*cosphi*dl
//   // dy = sintheta*sinphi*dl
//   // dz = costheta*dl
//   Point pc = mcell->center();
//   Point p_test;
//   if (sin(theta)*cos(phi)!=0){
//     double dl = (pc.x - p.x)/(sin(theta)*cos(phi));
//     if (dl<0) return false;
    
//     p_test.x = pc.x;
//     p_test.y = p.y + sin(theta)*sin(phi)*dl;
//     p_test.z = p.z + cos(theta)*dl;

//     // std::cout << pc.x/units::cm << " " << p.x/units::cm << " " << dl/units::cm << " " << sqrt(pow(p_test.x - pc.x,2)+pow(p_test.y - pc.y,2)+pow(p_test.z - pc.z,2))/units::cm << std::endl;
    
//   }else{
//     p_test.x = pc.x;
//     double dl = (pc.y-p.y)*sin(theta)*sin(phi) + (pc.z-p.z)*cos(theta);
//     if (dl<0) return false;
//     p_test.y = p.y + sin(theta)*sin(phi)*dl;
//     p_test.z = p.z + cos(theta)*dl;
//   }
//   const GeomWire *uwire = gds.closest(p_test,WirePlaneType_t(0));
//   const GeomWire *vwire = gds.closest(p_test,WirePlaneType_t(1));
//   const GeomWire *wwire = gds.closest(p_test,WirePlaneType_t(2));

//   GeomWireSelection uwires = mcell->get_uwires();
//   GeomWireSelection vwires = mcell->get_vwires();
//   GeomWireSelection wwires = mcell->get_wwires();
  
//   // std::cout << uwires.front()->index() << " " << uwire->index() << " " << uwires.back()->index() << " "
//   // 	    << vwires.front()->index() << " " << vwire->index() << " " << vwires.back()->index() << " "
//   // 	    << wwires.front()->index() << " " << wwire->index() << " " << wwires.back()->index() << " " << std::endl;
    
  
  

//   if (uwire->index() >= uwires.front()->index() && uwire->index() <= uwires.back()->index() &&
//       vwire->index() >= vwires.front()->index() && vwire->index() <= vwires.back()->index() &&
//       wwire->index() >= wwires.front()->index() && wwire->index() <= wwires.back()->index()){
//     return true;
//   }else{
//     return false;
//   }
  
// }
