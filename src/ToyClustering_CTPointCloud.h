using namespace WCP;

void WCP2dToy::Clustering_CTPointCloud(WCP::ToyCTPointCloud& ct_point_cloud, WCP::PR3DClusterSelection& live_clusters, WCP::map_pr3dcluster_double& cluster_length_map, std::map<int,std::pair<double,double>>& dead_u_index, std::map<int,std::pair<double,double>>& dead_v_index, std::map<int,std::pair<double,double>>& dead_w_index){

    // test a few different functions and then print out ...

    // is_good_point

    Point p(0*units::cm, 0*units::cm, 0*units::cm);
    auto p1 = live_clusters[0]->get_point_cloud()->get_closest_point(p);
    p = p1.second;
    std::cout << p << std::endl;

    bool flag = ct_point_cloud.is_good_point(p, 0.6*units::cm, 1,1);
    std::cout << "is_good_point: " << flag << std::endl;

    // get_closest_points
    auto closest_points_u= ct_point_cloud.get_closest_points(p, 0.6*units::cm, 0);
    auto closest_points_v= ct_point_cloud.get_closest_points(p, 0.6*units::cm, 1);
    auto closest_points_w= ct_point_cloud.get_closest_points(p, 0.6*units::cm, 2);


    std::cout << "get_closest_points: " << closest_points_u.pts.size() << " " << closest_points_v.pts.size() << " " << closest_points_w.pts.size() << std::endl;

    // get_closest_dead_chs
    bool flag_dead_chs_u = ct_point_cloud.get_closest_dead_chs(p, 0, 1);
    bool flag_dead_chs_v = ct_point_cloud.get_closest_dead_chs(p, 1, 1);
    bool flag_dead_chs_w = ct_point_cloud.get_closest_dead_chs(p, 2, 1);
    std::cout << "get_closest_dead_chs: " << flag_dead_chs_u << " " << flag_dead_chs_v << " " << flag_dead_chs_w << std::endl;

    // Convert_3Dpoint_time_ch
    auto time_ch = ct_point_cloud.convert_3Dpoint_time_ch(p);
    std::cout << "Convert_3Dpoint_time_ch: " << time_ch.at(0) << " " << time_ch.at(1) << " " << time_ch.at(2) << " " << time_ch.at(3) << std::endl;


}
