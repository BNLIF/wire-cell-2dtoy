 
void WireCell2dToy::Clustering_protect_overclustering(WireCell::PR3DClusterSelection& live_clusters,std::map<WireCell::PR3DCluster*,double>& cluster_length_map, WireCell::ToyCTPointCloud& ct_point_cloud){
  
  // can follow ToyClustering_separate to add clusters ...

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
    PR3DClusterSelection clusters = Examine_overclustering(cluster, ct_point_cloud);
    if (clusters.size()>0){
      del_clusters.push_back(cluster);
      std::copy(clusters.begin(), clusters.end(), std::back_inserter(new_clusters));
    }
  }
  
  // at the end of it ... 
  for (auto it=new_clusters.begin(); it!=new_clusters.end(); it++){
    PR3DCluster *ncluster = (*it);
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


PR3DClusterSelection WireCell2dToy::Examine_overclustering(PR3DCluster *cluster,  WireCell::ToyCTPointCloud& ct_point_cloud){
  //copy the graph operations here ... 
}
