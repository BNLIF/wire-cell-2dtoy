
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
  
  
}
