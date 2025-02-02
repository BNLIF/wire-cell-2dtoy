using namespace WCP;

void WCP2dToy::Clustering_CTPointCloud(WCP::ToyCTPointCloud& ct_point_cloud, WCP::PR3DClusterSelection& live_clusters, WCP::map_pr3dcluster_double& cluster_length_map, std::map<int,std::pair<double,double>>& dead_u_index, std::map<int,std::pair<double,double>>& dead_v_index, std::map<int,std::pair<double,double>>& dead_w_index){

    // test a few different functions and then print out ...

    // is_good_point
    Point p(0*units::cm, 0*units::cm, 0*units::cm);
    // for (size_t i=0; i!=live_clusters.size(); i++){
    //     auto p1 = live_clusters[i]->get_point_cloud()->get_closest_point(p);
    //     if (cluster_length_map[live_clusters.at(i)]/units::cm  > 200)
    //     std::cout << cluster_length_map[live_clusters.at(i)]/units::cm << " " << p1.second << std::endl;
    // }
    auto p1 = live_clusters[0]->get_point_cloud()->get_closest_point(p);
    p = p1.second;
    std::cout << p << std::endl;

    bool flag = ct_point_cloud.is_good_point(p, 0.6*units::cm, 1,1);
    bool flag_wc = ct_point_cloud.is_good_point_wc(p, 0.6*units::cm, 1,1);
    std::cout << "is_good_point: " << flag << " " << flag_wc << std::endl;

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

    std::cout << "Number of Points: " << ct_point_cloud.get_num_points(0)<< " " << ct_point_cloud.get_num_points(1) << " " << ct_point_cloud.get_num_points(2) << std::endl;

    auto num_planes = ct_point_cloud.test_good_point(p, 0.6*units::cm, 1);
    std::cout << "test_good_point: " << num_planes.at(0) << " " << num_planes.at(1) << " " << num_planes.at(2) << " " << num_planes.at(3) << " " << num_planes.at(4) << " " << num_planes.at(5) << std::endl;

    std::cout << "test ave charge " << ct_point_cloud.get_ave_3d_charge(p, 1.0*units::cm) << " " << ct_point_cloud.get_ave_charge(p, 1.0*units::cm, 0) << " " << ct_point_cloud.get_ave_charge(p, 1.0*units::cm, 1) << " " << ct_point_cloud.get_ave_charge(p, 1.0*units::cm, 2) << std::endl;

    auto point1 = ct_point_cloud.convert_time_ch_2Dpoint(10, 10, 0);
    auto point2 = ct_point_cloud.convert_time_ch_2Dpoint(10, 10+2400, 1);
    auto point3 = ct_point_cloud.convert_time_ch_2Dpoint(10, 10+4800, 2);

    std::cout << "test 2D conversion " 
              << point1.first << ", " << point1.second << " "
              << point2.first << ", " << point2.second << " "
              << point3.first << ", " << point3.second << std::endl;

    auto dead_chs_u = ct_point_cloud.get_overlap_dead_chs(10, 1000, 0, 2400, 0);
    auto dead_chs_v = ct_point_cloud.get_overlap_dead_chs(10, 1000, 2400, 4800, 1);
    auto dead_chs_w = ct_point_cloud.get_overlap_dead_chs(10, 1000, 4800, 8256, 2);
    std::cout << "test Overlap dead chs: " << dead_chs_u.size() << " " << dead_chs_v.size() << " " << dead_chs_w.size() << std::endl;

    std::cout << "test all dead chs " << ct_point_cloud.get_all_dead_chs().size() << std::endl; 

    std::cout << "test good chs " << ct_point_cloud.get_overlap_good_ch_charge(10,1000,0,2400,0).size() << " " << ct_point_cloud.get_overlap_good_ch_charge(10,1000,2400,4800,1).size() << " " << ct_point_cloud.get_overlap_good_ch_charge(10,1000,4800,8256,2).size() << " " << std::endl;

    std::cout << "Test new functions in Cluster" << std::endl;
    for (size_t i=0; i!=live_clusters.size(); i++){
        auto p1 = live_clusters[i]->get_point_cloud()->get_closest_point(p);
        if (cluster_length_map[live_clusters.at(i)]/units::cm  > 239){
            std::cout << cluster_length_map[live_clusters.at(i)]/units::cm << " " << p1.second << std::endl;
            auto points = live_clusters[i]->get_main_axis_wcps();
            std::cout << "(" << points.first.x << ", "  << points.first.y << ", " << points.first.z << ") " 
            << "(" << points.second.x << ", " << points.second.y << ", " << points.second.z << ")" << std::endl;

            Point p1(points.first.x, points.first.y, points.first.z);
            Point p2(points.second.x, points.second.y, points.second.z);
            auto dir1 = live_clusters[i]->calc_dir(p1, p2, 10*units::cm);
            std::cout << dir1.X() << " " << dir1.Y() << " "<< dir1.Z() << std::endl;

            PointVector points1;
            points1.push_back(p1);
            points1.push_back(p2);
            Point p3(-1204.49*units::mm, -57.85*units::mm, 5635*units::mm);
            points1.push_back(p3);
            live_clusters[i]->Calc_PCA(points1);
            std::cout << live_clusters[i]->get_center() << " " << live_clusters[i]->get_PCA_axis(0) << " " << live_clusters[i]->get_PCA_axis(1) << " " << live_clusters[i]->get_PCA_axis(2) << std::endl;

            Point p4(0,0,0);
            auto dir2 = live_clusters[i]->calc_PCA_dir(p4, points1);
            std::cout << dir2.X() << " " << dir2.Y() << " "<< dir2.Z() << std::endl;

            auto p5 = live_clusters[i]->calc_ave_pos(p1, 10);
            std::cout << p5 << std::endl;

            // test shortest path ...
            live_clusters[i]->dijkstra_shortest_paths(points.first, ct_point_cloud);
            live_clusters[i]->cal_shortest_path(points.second);



            auto path_wcps = live_clusters[i]->get_path_wcps();
            std::cout << path_wcps.size() << " (" 
                      << path_wcps.front().x << ", " << path_wcps.front().y << ", " << path_wcps.front().z << ") ("
                      << path_wcps.back().x << ", " << path_wcps.back().y << ", " << path_wcps.back().z << ")" << std::endl;
            
            // for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
            //     std::cout << "{geo_point_t temp_p(" << it->x << ", " << it->y << ", " << it->z << "); points6.push_back(temp_p);}" << std::endl;
            // }

            std::vector<WCPointCloud<double>::WCPoint> path_wcps_vec;//(path_wcps.begin(), path_wcps.end());
            live_clusters[i]->organize_wcps_path(path_wcps_vec, 0.6*units::cm);
            std::cout << path_wcps_vec.size() << " (" 
                      << path_wcps_vec.front().x << ", " << path_wcps_vec.front().y << ", " << path_wcps_vec.front().z << ") ("
                      << path_wcps_vec.back().x << ", " << path_wcps_vec.back().y << ", " << path_wcps_vec.back().z << ")" << std::endl;

            std::vector<WCPointCloud<double>::WCPoint> path_wcps_vec1(path_wcps.begin(), path_wcps.end());
            live_clusters[i]->organize_wcps_vec_path(path_wcps_vec1, 0.6*units::cm);
            std::cout << path_wcps_vec1.size() << " (" 
                      << path_wcps_vec1.front().x << ", " << path_wcps_vec1.front().y << ", " << path_wcps_vec1.front().z << ") ("
                      << path_wcps_vec1.back().x << ", " << path_wcps_vec1.back().y << ", " << path_wcps_vec1.back().z << ")" << std::endl;
        }
    }
}
