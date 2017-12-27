void WireCell2dToy::Clustering_separate(WireCell::PR3DClusterSelection& live_clusters, std::map<WireCell::PR3DCluster*,double>& cluster_length_map){
  for (size_t i=0;i!=live_clusters.size();i++){
    PR3DCluster *cluster = live_clusters.at(i);
    ToyPointCloud* cloud = cluster->get_point_cloud();
    std::vector<std::pair<WireCell::SlimMergeGeomCell*, WireCell::Point>> results = cloud->get_hull();
    std::cout << i << " " << cluster_length_map[cluster]/units::cm << " " << results.size() << std::endl;
  }
}
