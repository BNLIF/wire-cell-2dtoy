#include "WireCellData/DynamicToyPointCloud.h"

using namespace WireCell;


void WireCell2dToy::Clustering_deghost(WireCell::PR3DClusterSelection& live_clusters, std::map<WireCell::PR3DCluster*,double>& cluster_length_map){
  
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

  // Create two point clouds ... 
  // One for the points ... point --> index --> cluster (vector) ...
  // The other for the skeleton of each track ...  point --> index --> cluster (vector)
  // Both cloud needs to be dynamic, keep adding things into it as we improve the knowledge

  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();

  
  DynamicToyPointCloud global_point_cloud(angle_u,angle_v,angle_w);
  DynamicToyPointCloud global_skeleton_cloud(angle_u,angle_v,angle_w);
  
  for (size_t i=0;i!=live_clusters.size();i++){
    if (i==0){
      // fill anyway ...
      live_clusters.at(i)->Create_point_cloud();
      global_point_cloud.AddPoints(live_clusters.at(i),0);
      if (cluster_length_map[live_clusters.at(i)]>30*units::cm){
	live_clusters.at(i)->Construct_skeleton();
	global_skeleton_cloud.AddPoints(live_clusters.at(i),1);
      }
    }else{ 
      // start the process to add things in ... 
    } 
  }

  Point p(0,0,0);

  std::tuple<double, PR3DCluster*, size_t> test1 = global_point_cloud.get_closest_point_info(p);
  std::tuple<double, PR3DCluster*, size_t> test2 = global_point_cloud.get_closest_2d_point_info(p,0);
  std::tuple<double, PR3DCluster*, size_t> test3 = global_point_cloud.get_closest_2d_point_info(p,1);
  std::tuple<double, PR3DCluster*, size_t> test4 = global_point_cloud.get_closest_2d_point_info(p,2);

  std::cout << std::get<0>(test1)/units::cm << " " << std::get<1>(test1)->get_cluster_id() << " "
	    << std::get<0>(test2)/units::cm << " " << std::get<1>(test2)->get_cluster_id() << " "
	    << std::get<0>(test3)/units::cm << " " << std::get<1>(test3)->get_cluster_id() << " "
	    << std::get<0>(test4)/units::cm << " " << std::get<1>(test4)->get_cluster_id() << " "
	    << std::endl;


  test1 = global_skeleton_cloud.get_closest_point_info(p);
  test2 = global_skeleton_cloud.get_closest_2d_point_info(p,0);
  test3 = global_skeleton_cloud.get_closest_2d_point_info(p,1);
  test4 = global_skeleton_cloud.get_closest_2d_point_info(p,2);
  
  std::cout << std::get<0>(test1)/units::cm << " " << std::get<1>(test1)->get_cluster_id() << " "
	    << std::get<0>(test2)/units::cm << " " << std::get<1>(test2)->get_cluster_id() << " "
	    << std::get<0>(test3)/units::cm << " " << std::get<1>(test3)->get_cluster_id() << " "
	    << std::get<0>(test4)/units::cm << " " << std::get<1>(test4)->get_cluster_id() << " "
	    << std::endl;
}
