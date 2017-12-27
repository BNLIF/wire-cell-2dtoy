void WireCell2dToy::Clustering_separate(WireCell::PR3DClusterSelection& live_clusters, std::map<WireCell::PR3DCluster*,double>& cluster_length_map){
  for (size_t i=0;i!=live_clusters.size();i++){
    PR3DCluster *cluster = live_clusters.at(i);
    if (cluster_length_map[cluster]> 100*units::cm){
      ToyPointCloud* cloud = cluster->get_point_cloud();
      //std::vector<int> boundary_indices = cloud->get_hull();
      std::vector<WCPointCloud<double>::WCPoint> boundary_points = cloud->get_hull();
      
      cluster->Calc_PCA();
      
      //  std::cout << i << " " << cluster_length_map[cluster]/units::cm << " " << results.size() << " " << cluster->get_PCA_value(0) << " " << cluster->get_PCA_axis(0) << " " << cluster->get_PCA_value(1) << " " << cluster->get_PCA_axis(1) << " " << cluster->get_PCA_value(2) << " " << cluster->get_PCA_axis(2) << " " << std::endl;
      std::cout << i << " " << cluster->get_cluster_id() << " " << cluster_length_map[cluster]/units::cm << " " << boundary_points.size() << " " << cluster->get_PCA_value(1)/cluster->get_PCA_value(0) << std::endl;
      
    }
  }
}
