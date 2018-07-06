
#include <boost/graph/connected_components.hpp>

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
  //copy the create_graph from the PR3D Cluster ... 

  PR3DClusterSelection new_clusters;
  
  cluster->Create_point_cloud();
  const int N = cluster->get_num_points();
  MCUGraph *graph = new MCUGraph(N);
  
  ToyPointCloud* point_cloud = cluster->get_point_cloud();
  SMGCSelection& mcells = cluster->get_mcells();
  std::map<int,SMGCSet>& time_cells_set_map = cluster->get_time_cells_set_map();

  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  WireCell::WC2DPointCloud<double>& cloud_u = point_cloud->get_cloud_u();
  WireCell::WC2DPointCloud<double>& cloud_v = point_cloud->get_cloud_v();
  WireCell::WC2DPointCloud<double>& cloud_w = point_cloud->get_cloud_w();
  
  std::map<SlimMergeGeomCell*, std::map<int, std::set<int>>> map_mcell_uindex_wcps;
  std::map<SlimMergeGeomCell*, std::map<int, std::set<int>>> map_mcell_vindex_wcps;
  std::map<SlimMergeGeomCell*, std::map<int, std::set<int>>> map_mcell_windex_wcps;

   for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    std::map<int, std::set<int>> map_uindex_wcps;
    std::map<int, std::set<int>> map_vindex_wcps;
    std::map<int, std::set<int>> map_windex_wcps;
    std::vector<int>& wcps = point_cloud->get_mcell_indices(mcell);
    for (auto it1 = wcps.begin(); it1!=wcps.end(); it1++){
      WCPointCloud<double>::WCPoint& wcp = cloud.pts[*it1]; 
      // int index = wcp.index;
      // std::cout << index << " " << wcp.x << " " << wcp.y << " " << wcp.z << " " << wcp.index_u << " " << wcp.index_v << " " << wcp.index_w << std::endl;
      
      auto v = vertex(wcp.index, *graph); // retrieve vertex descriptor
      (*graph)[v].index = wcp.index;
      if (map_uindex_wcps.find(wcp.index_u)==map_uindex_wcps.end()){
   	std::set<int> wcps;
   	wcps.insert(wcp.index);
   	map_uindex_wcps[wcp.index_u] = wcps;
      }else{
   	map_uindex_wcps[wcp.index_u].insert(wcp.index);
      }
      
      if (map_vindex_wcps.find(wcp.index_v)==map_vindex_wcps.end()){
  	std::set<int> wcps;
  	wcps.insert(wcp.index);
  	map_vindex_wcps[wcp.index_v] = wcps;
      }else{
  	map_vindex_wcps[wcp.index_v].insert(wcp.index);
      }

      if (map_windex_wcps.find(wcp.index_w)==map_windex_wcps.end()){
  	std::set<int> wcps;
  	wcps.insert(wcp.index);
  	map_windex_wcps[wcp.index_w] = wcps;
      }else{
  	map_windex_wcps[wcp.index_w].insert(wcp.index);
      }
      
      
    }
    map_mcell_uindex_wcps[mcell] = map_uindex_wcps;
    map_mcell_vindex_wcps[mcell] = map_vindex_wcps;
    map_mcell_windex_wcps[mcell] = map_windex_wcps;
  }

    int num_edges = 0;
  
  // create graph for points inside the same mcell
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    std::vector<int>& wcps = point_cloud->get_mcell_indices(mcell);
    int max_wire_interval = mcell->get_max_wire_interval();
    int min_wire_interval = mcell->get_min_wire_interval();
    std::map<int, std::set<int>>* map_max_index_wcps;
    std::map<int, std::set<int>>* map_min_index_wcps;
    if (mcell->get_max_wire_type()==WirePlaneType_t(0)){
      map_max_index_wcps = &map_mcell_uindex_wcps[mcell];
    }else if (mcell->get_max_wire_type()==WirePlaneType_t(1)){
      map_max_index_wcps = &map_mcell_vindex_wcps[mcell];
    }else{
      map_max_index_wcps = &map_mcell_windex_wcps[mcell];
    }
    if (mcell->get_min_wire_type()==WirePlaneType_t(0)){
      map_min_index_wcps = &map_mcell_uindex_wcps[mcell];
    }else if (mcell->get_min_wire_type()==WirePlaneType_t(1)){
      map_min_index_wcps = &map_mcell_vindex_wcps[mcell];
    }else{
      map_min_index_wcps = &map_mcell_windex_wcps[mcell];
    }
    
    for (auto it1 = wcps.begin(); it1!=wcps.end(); it1++){
      WCPointCloud<double>::WCPoint& wcp1 = cloud.pts[*it1];
      int index1 = wcp1.index;
      int index_max_wire;
      int index_min_wire;
      if (mcell->get_max_wire_type()==WirePlaneType_t(0)){
  	index_max_wire = wcp1.index_u;
      }else if (mcell->get_max_wire_type()==WirePlaneType_t(1)){
  	index_max_wire = wcp1.index_v;
      }else{
  	index_max_wire = wcp1.index_w;
      }
      if (mcell->get_min_wire_type()==WirePlaneType_t(0)){
  	index_min_wire = wcp1.index_u;
      }else if (mcell->get_min_wire_type()==WirePlaneType_t(1)){
  	index_min_wire = wcp1.index_v;
      }else{
  	index_min_wire = wcp1.index_w;
      }

      std::vector<std::set<int>*> max_wcps_set;
      std::vector<std::set<int>*> min_wcps_set;
      
      // go through the first map and find the ones satisfying the condition
      for (auto it2 = map_max_index_wcps->begin(); it2!=map_max_index_wcps->end(); it2++){
   	if (fabs(it2->first - index_max_wire)<=max_wire_interval){
   	  max_wcps_set.push_back(&(it2->second));
   	}
      }
      // go through the second map and find the ones satisfying the condition
      for (auto it2 = map_min_index_wcps->begin(); it2!=map_min_index_wcps->end(); it2++){
   	if (fabs(it2->first - index_min_wire)<=min_wire_interval){
   	  min_wcps_set.push_back(&(it2->second));
   	}
      }

      std::set<int> wcps_set1;
      std::set<int> wcps_set2;

      for (auto it2 = max_wcps_set.begin(); it2!=max_wcps_set.end(); it2++){
	wcps_set1.insert((*it2)->begin(), (*it2)->end());
      }
      for (auto it3 = min_wcps_set.begin(); it3!=min_wcps_set.end(); it3++){
	wcps_set2.insert((*it3)->begin(), (*it3)->end());
      }
      
      
      // std::cout << max_wcps_set.size() << " " << min_wcps_set.size() << std::endl;
      // for (auto it2 = max_wcps_set.begin(); it2!=max_wcps_set.end(); it2++){
      //	for (auto it3 = min_wcps_set.begin(); it3!=min_wcps_set.end(); it3++){
      {
	std::set<int> common_set;
	set_intersection(wcps_set1.begin(), wcps_set1.end(), wcps_set2.begin(), wcps_set2.end(),std::inserter(common_set,common_set.begin()));

	//	std::cout << "S0: " << common_set.size() << std::endl;

	for (auto it4 = common_set.begin(); it4!=common_set.end(); it4++){
	  WCPointCloud<double>::WCPoint& wcp2 = cloud.pts[*it4];
	  if (wcp2.index != wcp1.index){
	    int index2 = wcp2.index;
	    //  std::cout << index1 << " " << index2 << std::endl;
	    // add edge ...
	    auto edge = add_edge(index1,index2,*graph);
	    if (edge.second){
	      (*graph)[edge.first].dist = sqrt(pow(wcp1.x-wcp2.x,2)+pow(wcp1.y-wcp2.y,2)+pow(wcp1.z-wcp2.z,2));
	      num_edges ++;
	      // std::cout << wcp1.x << " " << wcp1.y << " " << wcp1.z << " " << wcp1.index_u << " " << wcp1.index_v << " " << wcp1.index_w << " " << wcp2.index_u << " " << wcp2.index_v << " " << wcp2.index_w << std::endl;
	    }
	  }
	}
	//}
      }
    }
  }


  //  std::cout << "Xin: " << num_edges << " " << N << std::endl;
  
  
  
  std::vector<int> time_slices;
  for (auto it1 = time_cells_set_map.begin(); it1!=time_cells_set_map.end(); it1++){
    time_slices.push_back((*it1).first);
  }

  std::vector<std::pair<SlimMergeGeomCell*,SlimMergeGeomCell*>> connected_mcells;
  
  for (size_t i=0; i!= time_slices.size(); i++){
    SMGCSet& mcells_set = time_cells_set_map[time_slices.at(i)];
    
    // create graph for points in mcell inside the same time slice
    if (mcells_set.size()>=2){
      for (auto it2 = mcells_set.begin(); it2!=mcells_set.end();it2++){
  	SlimMergeGeomCell *mcell1 = *it2;
  	auto it2p = it2;
  	if (it2p!=mcells_set.end()){
  	  it2p++;
  	  for (auto it3 = it2p; it3!=mcells_set.end(); it3++){
  	    SlimMergeGeomCell *mcell2 = *(it3);
  	    //std::cout << mcell1 << " " << mcell2 << " " << mcell1->Overlap_fast(mcell2,2) << std::endl;
	    if (mcell1->Overlap_fast(mcell2,2))
	      connected_mcells.push_back(std::make_pair(mcell1,mcell2));
  	  }
  	}
      }
    }
    // create graph for points between connected mcells in adjacent time slices + 1, if not, + 2
    std::vector<SMGCSet> vec_mcells_set;
    if (i+1 < time_slices.size()){
      if (time_slices.at(i+1)-time_slices.at(i)==1){
	vec_mcells_set.push_back(time_cells_set_map[time_slices.at(i+1)]);
	if (i+2 < time_slices.size())
	  if (time_slices.at(i+2)-time_slices.at(i)==2)
	    vec_mcells_set.push_back(time_cells_set_map[time_slices.at(i+2)]);
      }else if (time_slices.at(i+1) - time_slices.at(i)==2){
	vec_mcells_set.push_back(time_cells_set_map[time_slices.at(i+1)]);
      }
    }
    //    bool flag = false;
    for (size_t j=0; j!=vec_mcells_set.size(); j++){
      //      if (flag) break;
      SMGCSet& next_mcells_set = vec_mcells_set.at(j);
      for (auto it1 = mcells_set.begin(); it1!= mcells_set.end(); it1++){
	SlimMergeGeomCell *mcell1 = (*it1);
	for (auto it2 = next_mcells_set.begin(); it2!=next_mcells_set.end(); it2++){
	  SlimMergeGeomCell *mcell2 = (*it2);
	  if (mcell1->Overlap_fast(mcell2,2)){
	    //	    flag = true; // correct???
	    connected_mcells.push_back(std::make_pair(mcell1,mcell2));
	  }
	}
      }
    }
  }
  
  // establish edge ... 
  std::map<std::pair<int,int>,std::pair<int,double> > closest_index;

  // std::cout << connected_mcells.size() << std::endl;
  for (auto it = connected_mcells.begin(); it!= connected_mcells.end(); it++){
    SlimMergeGeomCell *mcell1 = (*it).first;
    SlimMergeGeomCell *mcell2 = (*it).second;

    std::vector<int>& wcps1 = point_cloud->get_mcell_indices(mcell1);
    std::vector<int>& wcps2 = point_cloud->get_mcell_indices(mcell2);

    // test 2 against 1 ... 
    int max_wire_interval = mcell1->get_max_wire_interval();
    int min_wire_interval = mcell1->get_min_wire_interval();
    std::map<int, std::set<int>>* map_max_index_wcps;
    std::map<int, std::set<int>>* map_min_index_wcps;
    
    if (mcell1->get_max_wire_type()==WirePlaneType_t(0)){
      map_max_index_wcps = &map_mcell_uindex_wcps[mcell2];
    }else if (mcell1->get_max_wire_type()==WirePlaneType_t(1)){
      map_max_index_wcps = &map_mcell_vindex_wcps[mcell2];
    }else{
      map_max_index_wcps = &map_mcell_windex_wcps[mcell2];
    }
    if (mcell1->get_min_wire_type()==WirePlaneType_t(0)){
      map_min_index_wcps = &map_mcell_uindex_wcps[mcell2];
    }else if (mcell1->get_min_wire_type()==WirePlaneType_t(1)){
      map_min_index_wcps = &map_mcell_vindex_wcps[mcell2];
    }else{
      map_min_index_wcps = &map_mcell_windex_wcps[mcell2];
    }

    for (auto it1 = wcps1.begin(); it1!=wcps1.end(); it1++){
      WCPointCloud<double>::WCPoint& wcp1 = cloud.pts[*it1];
      int index1 = wcp1.index;
      int index_max_wire;
      int index_min_wire;
      if (mcell1->get_max_wire_type()==WirePlaneType_t(0)){
  	index_max_wire = wcp1.index_u;
      }else if (mcell1->get_max_wire_type()==WirePlaneType_t(1)){
  	index_max_wire = wcp1.index_v;
      }else{
  	index_max_wire = wcp1.index_w;
      }
      if (mcell1->get_min_wire_type()==WirePlaneType_t(0)){
  	index_min_wire = wcp1.index_u;
      }else if (mcell1->get_min_wire_type()==WirePlaneType_t(1)){
  	index_min_wire = wcp1.index_v;
      }else{
  	index_min_wire = wcp1.index_w;
      }
      std::vector<std::set<int>*> max_wcps_set;
      std::vector<std::set<int>*> min_wcps_set;
      // go through the first map and find the ones satisfying the condition
      for (auto it2 = map_max_index_wcps->begin(); it2!=map_max_index_wcps->end(); it2++){
   	if (fabs(it2->first - index_max_wire)<=max_wire_interval){
   	  max_wcps_set.push_back(&(it2->second));
   	}
      }
      // go through the second map and find the ones satisfying the condition
      for (auto it2 = map_min_index_wcps->begin(); it2!=map_min_index_wcps->end(); it2++){
   	if (fabs(it2->first - index_min_wire)<=min_wire_interval){
   	  min_wcps_set.push_back(&(it2->second));
   	}
      }

      std::set<int> wcps_set1;
      std::set<int> wcps_set2;

      for (auto it2 = max_wcps_set.begin(); it2!=max_wcps_set.end(); it2++){
	wcps_set1.insert((*it2)->begin(), (*it2)->end());
      }
      for (auto it3 = min_wcps_set.begin(); it3!=min_wcps_set.end(); it3++){
	wcps_set2.insert((*it3)->begin(), (*it3)->end());
      }

      
      //   for (auto it2 = max_wcps_set.begin(); it2!=max_wcps_set.end(); it2++){
      //	for (auto it3 = min_wcps_set.begin(); it3!=min_wcps_set.end(); it3++){
      {
	std::set<int> common_set;
	set_intersection(wcps_set1.begin(), wcps_set1.end(), wcps_set2.begin(), wcps_set2.end(),std::inserter(common_set,common_set.begin()));

	//	std::cout << "S1: " << common_set.size() << std::endl;
	//	  std::cout << common_set.size() << std::endl;

	//	std::map<int,std::pair<int,double> > closest_index;
	
	for (auto it4 = common_set.begin(); it4!=common_set.end(); it4++){
	  WCPointCloud<double>::WCPoint& wcp2 = cloud.pts[*it4];
	  if (wcp2.index != wcp1.index){
	    int index2 = wcp2.index;
	    double dis = sqrt(pow(wcp1.x-wcp2.x,2)+pow(wcp1.y-wcp2.y,2)+pow(wcp1.z-wcp2.z,2));
	    
	    if (closest_index.find(std::make_pair(index1,wcp2.mcell->GetTimeSlice()))==closest_index.end()){
	      closest_index[std::make_pair(index1,wcp2.mcell->GetTimeSlice())] = std::make_pair(index2,dis);
	    }else{
	      if (dis < closest_index[std::make_pair(index1,wcp2.mcell->GetTimeSlice())].second)
		closest_index[std::make_pair(index1,wcp2.mcell->GetTimeSlice())] = std::make_pair(index2,dis);
	    }
	  }
	}

	//	std::cout << closest_index.size() << std::endl;
	// for (auto it4 = closest_index.begin(); it4!=closest_index.end(); it4++){
	//   int index2 = it4->second.first;
	//   double dis = it4->second.second;
	//   auto edge = add_edge(index1,index2,*graph);
	//   if (edge.second){
	//     (*graph)[edge.first].dist = dis;
	//     num_edges ++;
	//   }
	// }
      }
      //}
      
    }


    // test 1 against 2 ...
    max_wire_interval = mcell2->get_max_wire_interval();
    min_wire_interval = mcell2->get_min_wire_interval();
    if (mcell2->get_max_wire_type()==WirePlaneType_t(0)){
      map_max_index_wcps = &map_mcell_uindex_wcps[mcell1];
    }else if (mcell2->get_max_wire_type()==WirePlaneType_t(1)){
      map_max_index_wcps = &map_mcell_vindex_wcps[mcell1];
    }else{
      map_max_index_wcps = &map_mcell_windex_wcps[mcell1];
    }
    if (mcell2->get_min_wire_type()==WirePlaneType_t(0)){
      map_min_index_wcps = &map_mcell_uindex_wcps[mcell1];
    }else if (mcell2->get_min_wire_type()==WirePlaneType_t(1)){
      map_min_index_wcps = &map_mcell_vindex_wcps[mcell1];
    }else{
      map_min_index_wcps = &map_mcell_windex_wcps[mcell1];
    }
    for (auto it1 = wcps2.begin(); it1!=wcps2.end(); it1++){
      WCPointCloud<double>::WCPoint& wcp1 = cloud.pts[*it1];
      int index1 = wcp1.index;
      int index_max_wire;
      int index_min_wire;
      if (mcell2->get_max_wire_type()==WirePlaneType_t(0)){
  	index_max_wire = wcp1.index_u;
      }else if (mcell2->get_max_wire_type()==WirePlaneType_t(1)){
  	index_max_wire = wcp1.index_v;
      }else{
  	index_max_wire = wcp1.index_w;
      }
      if (mcell2->get_min_wire_type()==WirePlaneType_t(0)){
  	index_min_wire = wcp1.index_u;
      }else if (mcell2->get_min_wire_type()==WirePlaneType_t(1)){
  	index_min_wire = wcp1.index_v;
      }else{
  	index_min_wire = wcp1.index_w;
      }
      std::vector<std::set<int>*> max_wcps_set;
      std::vector<std::set<int>*> min_wcps_set;
      // go through the first map and find the ones satisfying the condition
      for (auto it2 = map_max_index_wcps->begin(); it2!=map_max_index_wcps->end(); it2++){
   	if (fabs(it2->first - index_max_wire)<=max_wire_interval){
   	  max_wcps_set.push_back(&(it2->second));
   	}
      }
      // go through the second map and find the ones satisfying the condition
      for (auto it2 = map_min_index_wcps->begin(); it2!=map_min_index_wcps->end(); it2++){
   	if (fabs(it2->first - index_min_wire)<=min_wire_interval){
   	  min_wcps_set.push_back(&(it2->second));
   	}
      }

      std::set<int> wcps_set1;
      std::set<int> wcps_set2;

      for (auto it2 = max_wcps_set.begin(); it2!=max_wcps_set.end(); it2++){
	wcps_set1.insert((*it2)->begin(), (*it2)->end());
      }
      for (auto it3 = min_wcps_set.begin(); it3!=min_wcps_set.end(); it3++){
	wcps_set2.insert((*it3)->begin(), (*it3)->end());
      }

      
      // for (auto it2 = max_wcps_set.begin(); it2!=max_wcps_set.end(); it2++){
      // 	for (auto it3 = min_wcps_set.begin(); it3!=min_wcps_set.end(); it3++){
      {
	std::set<int> common_set;
	set_intersection(wcps_set1.begin(), wcps_set1.end(), wcps_set2.begin(), wcps_set2.end(),std::inserter(common_set,common_set.begin()));

	//	std::cout << "S2: " << common_set.size() << std::endl;

	//	std::map<int,std::pair<int,double> > closest_index;
	
	for (auto it4 = common_set.begin(); it4!=common_set.end(); it4++){
	  WCPointCloud<double>::WCPoint& wcp2 = cloud.pts[*it4];
	  if (wcp2.index != wcp1.index){
	    int index2 = wcp2.index;
	    double dis = sqrt(pow(wcp1.x-wcp2.x,2)+pow(wcp1.y-wcp2.y,2)+pow(wcp1.z-wcp2.z,2));
	    
	    if (closest_index.find(std::make_pair(index1,wcp2.mcell->GetTimeSlice()))==closest_index.end()){
	      closest_index[std::make_pair(index1,wcp2.mcell->GetTimeSlice())] = std::make_pair(index2,dis);
	    }else{
	      if (dis < closest_index[std::make_pair(index1,wcp2.mcell->GetTimeSlice())].second)
		closest_index[std::make_pair(index1,wcp2.mcell->GetTimeSlice())] = std::make_pair(index2,dis);
	    }
	  }
	}

	//std::cout << closest_index.size() << std::endl;
	// for (auto it4 = closest_index.begin(); it4!=closest_index.end(); it4++){
	//   int index2 = it4->second.first;
	//   double dis = it4->second.second;
	//   auto edge = add_edge(index1,index2,*graph);
	//   if (edge.second){
	//     (*graph)[edge.first].dist = dis;
	//     num_edges ++;
	//   }
	// }

	
	// for (auto it4 = common_set.begin(); it4!=common_set.end(); it4++){
	//   WCPointCloud<double>::WCPoint& wcp2 = cloud.pts[*it4];
	//   if (wcp2.index != wcp1.index){
	//     int index2 = wcp2.index;
	//     auto edge = add_edge(index1,index2,*graph);
	//     if (edge.second){
	//       (*graph)[edge.first].dist = sqrt(pow(wcp1.x-wcp2.x,2)+pow(wcp1.y-wcp2.y,2)+pow(wcp1.z-wcp2.z,2));
	//       num_edges ++;
	//     }
	//   }
	// }
      }
      //      }
    }
  }

  for (auto it4 = closest_index.begin(); it4!=closest_index.end(); it4++){
    int index1 = it4->first.first;
    int index2 = it4->second.first;
    double dis = it4->second.second;
    auto edge = add_edge(index1,index2,*graph);
    if (edge.second){
      (*graph)[edge.first].dist = dis;
      num_edges ++;
    }
  }
  // end of copying ... 


  // now form the connected components
  std::vector<int> component(num_vertices(*graph));
  const int num = connected_components(*graph,&component[0]);
  if (num > 1){
    std::vector<ToyPointCloud*> pt_clouds;
    for (int j=0;j!=num;j++){
      ToyPointCloud *pt_cloud = new ToyPointCloud();
      pt_clouds.push_back(pt_cloud);
    }
    
    std::vector<int>::size_type i;
    for (i=0;i!=component.size(); ++i){
      pt_clouds.at(component[i])->AddPoint(cloud.pts[i],cloud_u.pts[i],cloud_v.pts[i],cloud_w.pts[i]);
      //   std::cout << "Vertex " << i << " " << cloud.pts[i].x << " " << cloud.pts[i].y << " " << cloud.pts[i].z << " " << cloud.pts[i].index_u << " " << cloud.pts[i].index_v << " " << cloud.pts[i].index_w << " " << cloud.pts[i].mcell << " " << cloud.pts[i].mcell->GetTimeSlice()  << " is in component " << component[i] << std::endl;
    }
    for (int j=0;j!=num;j++){
      pt_clouds.at(j)->build_kdtree_index();
    }
    
    std::vector< std::vector< std::tuple<int,int,double> > > index_index_dis(num, std::vector< std::tuple<int,int,double> >(num));
    //std::vector< std::vector< std::tuple<int,int,double> > > index_index_dis_mst(num, std::vector< std::tuple<int,int,double> >(num));
    std::vector< std::vector< std::tuple<int,int,double> > > index_index_dis_dir1(num, std::vector< std::tuple<int,int,double> >(num));
    std::vector< std::vector< std::tuple<int,int,double> > > index_index_dis_dir2(num, std::vector< std::tuple<int,int,double> >(num));
    
    for (int j=0;j!=num;j++){
      for (int k=0;k!=num;k++){
	index_index_dis[j][k] = std::make_tuple(-1,-1,1e9);
	//index_index_dis_mst[j][k] = std::make_tuple(-1,-1,1e9);
	index_index_dis_dir1[j][k] = std::make_tuple(-1,-1,1e9);
	index_index_dis_dir2[j][k] = std::make_tuple(-1,-1,1e9);
      }
    }
    
    // check against the closest distance ...
    // no need to have MST ... 
    for (int j=0;j!=num;j++){
      for (int k=j+1;k!=num;k++){
	index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(pt_clouds.at(k));
	
	if (num < 100 && pt_clouds.at(j)->get_num_points()>100 && pt_clouds.at(k)->get_num_points()>100 &&
	    (pt_clouds.at(j)->get_num_points()+pt_clouds.at(k)->get_num_points()) > 400 ||
	    pt_clouds.at(j)->get_num_points()>500 && pt_clouds.at(k)->get_num_points()>500){
	  WCPointCloud<double>::WCPoint wp1 = cloud.pts.at(std::get<0>(index_index_dis[j][k]));
	  WCPointCloud<double>::WCPoint wp2 = cloud.pts.at(std::get<1>(index_index_dis[j][k]));
	  Point p1(wp1.x,wp1.y,wp1.z);
	  Point p2(wp2.x,wp2.y,wp2.z);
	  
	  TVector3 dir1 = cluster->VHoughTrans(p1, 30*units::cm);
	  TVector3 dir2 = cluster->VHoughTrans(p2, 30*units::cm);
	  dir1 *= -1;
	  dir2 *= -1;
	  
	  std::pair<int,double> result1 = pt_clouds.at(k)->get_closest_point_along_vec(p1, dir1, 80*units::cm, 5*units::cm, 7.5, 3*units::cm);
	  
	  if (result1.first >=0){
	    Point test_p1(cloud.pts.at(std::get<0>(index_index_dis[j][k])).x,cloud.pts.at(std::get<0>(index_index_dis[j][k])).y,cloud.pts.at(std::get<0>(index_index_dis[j][k])).z);
	    Point test_p2(cloud.pts.at(result1.first).x,cloud.pts.at(result1.first).y,cloud.pts.at(result1.first).z);
	    double dis = sqrt(pow(test_p2.x-test_p1.x,2)+pow(test_p2.y-test_p1.y,2)+pow(test_p2.z-test_p1.z,2));
	    int num_points = dis/(1.5*units::cm)+1;
	    int num_cut_points = 0;
	    for (size_t k1=0; k1!=num_points-1; k1++){
	      Point test_p3(test_p1.x + (test_p2.x-test_p1.x) * (k1+1)/num_points ,
			    test_p1.y + (test_p2.y-test_p1.y) * (k1+1)/num_points ,
			    test_p1.z + (test_p2.z-test_p1.z) * (k1+1)/num_points );
	      double dis1 = point_cloud->get_closest_dis(test_p3);
	      if (dis1 < 1*units::cm)
		num_cut_points ++;
	    }
	    if (num_cut_points <=8 && num_cut_points< 0.25 * num_points + 2 && dis > 5*units::cm)
	      index_index_dis_dir1[j][k] = std::make_tuple(std::get<0>(index_index_dis[j][k]), result1.first, result1.second);
	  }
	  
	  std::pair<int,double> result2 = pt_clouds.at(j)->get_closest_point_along_vec(p2, dir2, 80*units::cm, 5*units::cm, 7.5, 3*units::cm);
	  
	  if (result2.first >=0){
	    
	    Point test_p1(cloud.pts.at(std::get<1>(index_index_dis[j][k])).x,cloud.pts.at(std::get<1>(index_index_dis[j][k])).y,cloud.pts.at(std::get<1>(index_index_dis[j][k])).z);
	    Point test_p2(cloud.pts.at(result2.first).x,cloud.pts.at(result2.first).y,cloud.pts.at(result2.first).z);
	    double dis = sqrt(pow(test_p2.x-test_p1.x,2)+pow(test_p2.y-test_p1.y,2)+pow(test_p2.z-test_p1.z,2));
	    int num_points = dis/(1.5*units::cm)+1;
	    int num_cut_points = 0;
	    for (size_t k1=0; k1!=num_points-1; k1++){
	      Point test_p3(test_p1.x + (test_p2.x-test_p1.x) * (k1+1)/num_points ,
			    test_p1.y + (test_p2.y-test_p1.y) * (k1+1)/num_points ,
			    test_p1.z + (test_p2.z-test_p1.z) * (k1+1)/num_points );
	      double dis1 = point_cloud->get_closest_dis(test_p3);
	      if ( dis1 < 1*units::cm)
		num_cut_points ++;
	    }
	    
	    if (num_cut_points <=8 && num_cut_points < 0.25 * num_points + 2 && dis > 5*units::cm)
	      index_index_dis_dir2[j][k] = std::make_tuple(result2.first, std::get<1>(index_index_dis[j][k]), result2.second);
	  }
	}
	
	// Now check the path ... 
	{
	  WCPointCloud<double>::WCPoint wp1 = cloud.pts.at(std::get<0>(index_index_dis[j][k]));
	  WCPointCloud<double>::WCPoint wp2 = cloud.pts.at(std::get<1>(index_index_dis[j][k]));
	  Point p1(wp1.x,wp1.y,wp1.z);
	  Point p2(wp2.x,wp2.y,wp2.z);
	  
	  double dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
	  double step_dis = 1.0*units::cm;
	  int num_steps = dis/step_dis + 1;
	  int num_bad = 0;
	  for (int ii=0;ii!=num_steps;ii++){
	    Point test_p;
	    test_p.x = p1.x + (p2.x-p1.x)/num_steps*(ii+1);
	    test_p.y = p1.y + (p2.y-p1.y)/num_steps*(ii+1);
	    test_p.z = p1.z + (p2.z-p1.z)/num_steps*(ii+1);
	    if (!ct_point_cloud.is_good_point(test_p)){
	      num_bad ++;
	      /* if (cluster->get_cluster_id()==11) */
	      /* 	std::cout << test_p.x/units::cm << " " << test_p.y/units::cm << " " << test_p.z/units::cm << std::endl; */
	    }
	  }
	  
	  //  std::cout << cluster->get_cluster_id() << " " << j << " " << k << " " << num_bad << " " << num_steps << std::endl;
	   
	  if (num_bad > 7 ||
	      num_bad > 2 && num_bad >=0.75*num_steps){
	    index_index_dis[j][k] = std::make_tuple(-1,-1,1e9);
	  }
	  //	    std::cout << j << " " << k << " " << num_bad << " " << num_steps << std::endl;
	}
	
	// Now check the path ... 
	if (std::get<0>(index_index_dis_dir1[j][k])>=0){
	  WCPointCloud<double>::WCPoint wp1 = cloud.pts.at(std::get<0>(index_index_dis_dir1[j][k]));
	  WCPointCloud<double>::WCPoint wp2 = cloud.pts.at(std::get<1>(index_index_dis_dir1[j][k]));
	  Point p1(wp1.x,wp1.y,wp1.z);
	  Point p2(wp2.x,wp2.y,wp2.z);
	  
	  double dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
	  double step_dis = 1.0*units::cm;
	  int num_steps = dis/step_dis + 1;
	  int num_bad = 0;
	  for (int ii=0;ii!=num_steps;ii++){
	    Point test_p;
	    test_p.x = p1.x + (p2.x-p1.x)/num_steps*(ii+1);
	    test_p.y = p1.y + (p2.y-p1.y)/num_steps*(ii+1);
	    test_p.z = p1.z + (p2.z-p1.z)/num_steps*(ii+1);
	    if (!ct_point_cloud.is_good_point(test_p))
	      num_bad ++;
	  }
	  
	  
	  if (num_bad > 7 ||
	      num_bad > 2 && num_bad >=0.75*num_steps){
	    index_index_dis_dir1[j][k] = std::make_tuple(-1,-1,1e9);
	  }
	  //	    std::cout << j << " " << k << " " << num_bad << " " << num_steps << std::endl;
	}
	
	
	// Now check the path ... 
	if (std::get<0>(index_index_dis_dir2[j][k])>=0){
	  WCPointCloud<double>::WCPoint wp1 = cloud.pts.at(std::get<0>(index_index_dis_dir2[j][k]));
	  WCPointCloud<double>::WCPoint wp2 = cloud.pts.at(std::get<1>(index_index_dis_dir2[j][k]));
	  Point p1(wp1.x,wp1.y,wp1.z);
	  Point p2(wp2.x,wp2.y,wp2.z);
	  
	  double dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
	  double step_dis = 1.0*units::cm;
	  int num_steps = dis/step_dis + 1;
	  int num_bad = 0;
	  for (int ii=0;ii!=num_steps;ii++){
	    Point test_p;
	    test_p.x = p1.x + (p2.x-p1.x)/num_steps*(ii+1);
	    test_p.y = p1.y + (p2.y-p1.y)/num_steps*(ii+1);
	    test_p.z = p1.z + (p2.z-p1.z)/num_steps*(ii+1);
	    if (!ct_point_cloud.is_good_point(test_p))
	      num_bad ++;
	  }
	 
	  
	  if (num_bad > 7 ||
	      num_bad > 2 && num_bad >=0.75*num_steps){
	    index_index_dis_dir2[j][k] = std::make_tuple(-1,-1,1e9);
	  }
	  //	    std::cout << j << " " << k << " " << num_bad << " " << num_steps << std::endl;
	}
	
	
	// establish the path ... 
	if (std::get<0>(index_index_dis[j][k])>=0){
	  auto edge = add_edge(std::get<0>(index_index_dis[j][k]),std::get<1>(index_index_dis[j][k]),*graph);
	  if (edge.second){
	    if (std::get<2>(index_index_dis[j][k])>5*units::cm){
	      (*graph)[edge.first].dist = std::get<2>(index_index_dis[j][k]);
	    }else{
	      (*graph)[edge.first].dist = std::get<2>(index_index_dis[j][k]);
	    }
	  }
	}
	
	if (std::get<0>(index_index_dis_dir1[j][k])>=0){
	  auto edge = add_edge(std::get<0>(index_index_dis_dir1[j][k]),std::get<1>(index_index_dis_dir1[j][k]),*graph);
	  if (edge.second){
	    if (std::get<2>(index_index_dis_dir1[j][k])>5*units::cm){
	      (*graph)[edge.first].dist = std::get<2>(index_index_dis_dir1[j][k]);
	    }else{
	      (*graph)[edge.first].dist = std::get<2>(index_index_dis_dir1[j][k]);
	    }
	  }
	}
	if (std::get<0>(index_index_dis_dir2[j][k])>=0){
	  auto edge = add_edge(std::get<0>(index_index_dis_dir2[j][k]),std::get<1>(index_index_dis_dir2[j][k]),*graph);
	  if (edge.second){
	    if (std::get<2>(index_index_dis_dir2[j][k])>5*units::cm){
	      (*graph)[edge.first].dist = std::get<2>(index_index_dis_dir2[j][k]);
	    }else{
	      (*graph)[edge.first].dist = std::get<2>(index_index_dis_dir2[j][k]);
	    }
	  }
	}
	// end check ... 
      }
    }
    
    
    // study the independent component again ... 
    {
      std::vector<int> component1(num_vertices(*graph));
      const int num1 = connected_components(*graph,&component1[0]);
      
      std::cout << cluster->get_cluster_id() << " " << num << " " << num1 << std::endl;
      if (num1 >1){
	// form new clusters ...
	std::vector<SMGCSet> vec_mcells_set;
	for (int ii=0;ii!=num1;ii++){
	  SMGCSet mcells_set;
	  vec_mcells_set.push_back(mcells_set);
	}
	
	std::vector<int>::size_type i;
	for (i=0;i!=component1.size(); ++i){
	  vec_mcells_set.at(component1[i]).insert(cloud.pts[i].mcell);
	  //      pt_clouds.at(component[i])->AddPoint(cloud.pts[i],cloud_u.pts[i],cloud_v.pts[i],cloud_w.pts[i]);
	  //   std::cout << "Vertex " << i << " " << cloud.pts[i].x << " " << cloud.pts[i].y << " " << cloud.pts[i].z << " " << cloud.pts[i].index_u << " " << cloud.pts[i].index_v << " " << cloud.pts[i].index_w << " " << cloud.pts[i].mcell << " " << cloud.pts[i].mcell->GetTimeSlice()  << " is in component " << component[i] << std::endl;
	}
	for (int ii=0;ii!=num1;ii++){
	  PR3DCluster *cluster1 = new PR3DCluster(1);
	  for (auto it = vec_mcells_set.at(ii).begin(); it!=vec_mcells_set.at(ii).end(); it++){
	    SlimMergeGeomCell *mcell = (*it);
	    //	  std::cout << vec_mcells_set.at(ii).size() << std::endl;
	    cluster1->AddCell(mcell,mcell->GetTimeSlice());
	  }
	  new_clusters.push_back(cluster1);
	}
	
      }
    }
    
    
    
    for (int i=0;i!=num;i++){
      delete pt_clouds.at(i);
    }
    
  } // if (num > 1) 
  
  
  delete graph;
  return new_clusters;
}
