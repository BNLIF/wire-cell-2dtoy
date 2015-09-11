#include "WireCell2dToy/ToyTracking.h"
#include "TVector3.h"

using namespace WireCell;

WireCell2dToy::ToyTracking::ToyTracking(WireCell2dToy::ToyCrawler& toycrawler){
  CreateVertices(toycrawler);
  RemoveSame(); // get rid of duplicated ones
  
  MergeVertices();   // merge things that are common ...
  OrganizeTracks();  // improve the end points nearby and break things
  Associate();    //associate the rest
   
  MergeVertices(0); //merge things, but require the center must be at the end of common track
  //CheckVertices(toycrawler);
  OrganizeTracks(); //improve the end points nearby and break things
  Associate(); //associate the rest
  
  MergeVertices(2);  // merge some vertices with distance cut only 
  CheckVertices(toycrawler);
  Crawl();
  OrganizeTracks(); //associate the rest
  Associate(); //associate the rest
    
  MergeVertices(2);  // merge some vertices together, allow single track
  CheckVertices(toycrawler);
  BreakTracks();    //improve the end points and break things
  OrganizeTracks(); //associate the rest

  // MergeVertices(2);  // merge some vertices together, allow single track
  // CheckVertices(toycrawler);
  // OrganizeTracks(); //associate the rest

  Associate();  //associate the rest .. 
  CleanUpVertex(); //do some final clean up ... 
  

  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    
    bool success = vertex->FindVertex();
    //std::cout << i << " " << vertices.size() << std::endl;
    if (success){
      if (!ExamineVertex(vertex,toycrawler)){
    	bool success1 = vertex->FindVertex();
    	if (success1){
    	  ExamineVertex(vertex,toycrawler);
    	}
      }
    } else{
      // to be improved later ...
      ExamineVertex(vertex,toycrawler);
      success = vertex->FindVertex(1);
      if (success){
       	ExamineVertex(vertex,toycrawler);
      }else{
    	vertex->reset_center();
      }
    }

    //std::cout << i << " Vertex " << vertex->get_ntracks() << " " << success << " " << vertex->Center().x/units::cm << " " << vertex->Center().y/units::cm << " " << vertex->Center().z/units::cm << " " << vertex->get_msc()->Get_Center().x/units::cm << std::endl;
    
    // for (int j=0;j!=vertex->get_ntracks();j++){
    //   std::cout << vertex->get_tracks().at(j)->get_end_scells().at(0)->Get_Center().x/units::cm << " " 
    // 		<< vertex->get_tracks().at(j)->get_end_scells().at(1)->Get_Center().x/units::cm << " " 
    // 		<< vertex->get_tracks().at(j)->get_all_cells().front()->Get_Center().x/units::cm << " " 
    // 		<< vertex->get_tracks().at(j)->get_all_cells().back()->Get_Center().x/units::cm << " " 
    // 		<< std::endl;
    // }
  }
  CheckVertices(toycrawler); // redefine all the tracks ... 
  // for (int i=0;i!=vertices.size();i++){
  //   WCVertex *vertex = vertices.at(i);
  //   std::cout << vertex->Center().x/units::cm << " " << vertex->get_msc()->Get_Center().x/units::cm <<  std::endl;
  //   for (int j=0;j!=vertex->get_ntracks();j++){
  //     WCTrack *track = vertex->get_tracks().at(j);
  //     std::cout << "abc " << track->get_end_scells().at(0)->Get_Center().x/units::cm << " " <<  track->get_end_scells().at(1)->Get_Center().x/units::cm << " " <<  std::endl;
  //   }
  // }
  // no need any more
  //OrganizeTracks(); 

  

  // for (int i=0;i!=tracks.size();i++){
  //   for (int j=0;j!=tracks.at(i)->get_all_cells().size();j++){
  //     std::cout << i << " " << j << " " << tracks.at(i)->get_all_cells().at(j)->Get_Center().x/units::cm << " " << tracks.at(i)->get_all_cells().at(j)->Get_Center().y/units::cm << " " << tracks.at(i)->get_all_cells().at(j)->Get_Center().z/units::cm << std::endl;
  //   }
  // }


  //deal with wiggle tracks
  //std::cout << "Deal with Wiggle Tracks " << std::endl;
  deal_wiggle_tracks();
  CheckVertices(toycrawler); // redefine all the tracks ... 
  update_maps();
  
  
  //CheckVertices(toycrawler); // redefine all the tracks ... 
  
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    
    bool success = vertex->get_fit_success();
    std::cout << i << " Vertex " << vertex->get_ntracks() << " " << success << " " << vertex->Center().x/units::cm << " " << vertex->Center().y/units::cm << " " << vertex->Center().z/units::cm << " " << vertex->get_msc()->Get_Center().x/units::cm << std::endl;
    // for (int j=0;j!=vertex->get_ntracks();j++){
    //   std::cout << j << " " << vertex->get_tracks().at(j)->get_end_scells().at(0)->Get_Center().x/units::cm 
    // 		<< " " << vertex->get_tracks().at(j)->get_end_scells().at(1)->Get_Center().x/units::cm  << std::endl;
    // }
  }

  if (vertices.size() > 0 ) {
    
    //Now do fine tracking for existing tracks
    std::cout << "FineTracking " << std::endl; 
    update_maps();
    fine_tracking();
    cleanup_bad_tracks();
    
    update_maps1();
    
    
    
    
    // Now need to figure out how to judge whether this is a shower or track ... 
    // just from the number of tracks and the connectivities?
    bool shower_flag = IsThisShower(toycrawler);
    
    std::cout << "Shower? " << shower_flag << std::endl;
    
    // separate various issues ... 
    
    if (!shower_flag){
      //not a shower
      std::cout << "Grown single track " << std::endl;
      if (grow_track_fill_gap(toycrawler)){
	std::cout << "FineTracking Again" << std::endl; 
	update_maps();
	fine_tracking();
      }
      
      form_parallel_tiny_tracks(toycrawler);
      
    }else{
      //is a shower  
      
    }
  }
}

void WireCell2dToy::ToyTracking::form_parallel_tiny_tracks(WireCell2dToy::ToyCrawler& toycrawler){
  
  //form all the cells
  MergeSpaceCellMap& mcells_map = toycrawler.Get_mcells_map();
  MergeSpaceCellSelection mcells;
  for (auto it = mcells_map.begin(); it != mcells_map.end(); it++){
    mcells.push_back(it->first);
  }
  
  //examine all the cells to find out outlier 
  for (int i=0;i!=mcells.size();i++){
    MergeSpaceCell *mcell = mcells.at(i);
    std::vector<int> track_no;
    
    int flag = -1;
    for (int j=0;j!=good_tracks.size();j++){
      if (good_tracks.at(j)->IsContained(mcell)){
	flag = 1;
	track_no.push_back(j);
      }
    }
    
    if (flag == -1){
      //entire thing ... 
      if (mcell->Get_all_spacecell().size()>0){
	//orig_mcells.push_back(mcell);
	new_mcells.push_back(mcell);
	new_mcells_map[mcell] = mcell;
      }
    }else{
      MergeSpaceCell *nmcell = new MergeSpaceCell();
      
      int flag2 = 0;

      for (int j=0;j!=mcell->Get_all_spacecell().size();j++){
	SpaceCell *cell = mcell->Get_all_spacecell().at(j);
	int flag1 = 0;

	for (int k=0;k!=track_no.size();k++){
	  double dist = good_tracks.at(track_no.at(k))->dist_proj(mcell,cell)/units::mm;
	  if (dist < 6.0){
	    flag1 = 1;
	    break;
	  }
	}
	
	if (flag1 == 0){
	  flag2 ++;
	  nmcell->AddSpaceCell(cell);
	  // partial thing
	}
      }
      // if (flag2 > 0)
      // 	std::cout << flag2 << std::endl;
      if (flag2 >0){
	//orig_mcells.push_back(mcell);
	new_mcells.push_back(nmcell);
	new_mcells_map[nmcell] = mcell;
      }else{
	delete nmcell;
      }
    }
  }


  std::vector<MergeSpaceCellSelection> cluster_msc;
  // Need to cluster these by whether they are connected ...
  for (int i=0;i!=new_mcells.size();i++){
    MergeSpaceCell *mcell = new_mcells.at(i);

    if (cluster_msc.size() == 0 ){
      MergeSpaceCellSelection mscs;
      mscs.push_back(mcell);
      cluster_msc.push_back(mscs);
    }else{
      int flag = 0;
      for (int j=0;j!=cluster_msc.size();j++){
   	MergeSpaceCellSelection& mscs = cluster_msc.at(j);
   	for (int k=0;k!=mscs.size();k++){
   	  MergeSpaceCell *mcell1 = mscs.at(k);
	  
   	  if (fabs(mcell1->Get_Center().x - mcell->Get_Center().x) < mcell1->thickness() + 0.2*units::mm &&
	      mcell1->Overlap(*mcell)){
   	    mscs.push_back(mcell);
   	    flag = 1;
	    break;
   	  }
   	}
	if (flag == 1)
	  break;
      }
      if (flag == 0){
  	MergeSpaceCellSelection mscs;
  	mscs.push_back(mcell);
  	cluster_msc.push_back(mscs);
      }
    }
  }
  
  int flag = 1;
  while(flag){
    flag = 0;
    for (int i=0;i!=cluster_msc.size();i++){
      MergeSpaceCellSelection& mscs_1 = cluster_msc.at(i);
      for (int j=i+1;j< cluster_msc.size();j++){
  	MergeSpaceCellSelection& mscs_2 = cluster_msc.at(j);
	
  	for (int k1 = 0; k1 != mscs_1.size(); k1++){
  	  MergeSpaceCell *mcell1 = mscs_1.at(k1);
  	  for (int k2 = 0; k2!= mscs_2.size(); k2++){
  	    MergeSpaceCell * mcell2 = mscs_2.at(k2);
	    
  	    if (fabs(mcell1->Get_Center().x - mcell2->Get_Center().x) < mcell1->thickness() + 0.2*units::mm &&
  		mcell1->Overlap(*mcell2)){

  	      cluster_msc.at(i).insert(cluster_msc.at(i).end(),cluster_msc.at(j).begin(),cluster_msc.at(j).end());
  	      cluster_msc.erase(cluster_msc.begin() + j);
  	      
	      //std::cout << flag << std::endl;
  	      flag = 1;
  	      break;
  	    }
  	  }
  	  if (flag == 1) break;
  	}
  	if (flag == 1) break;
      }
      if (flag == 1) break;
    }
  }
  
  flag = 1;
  while (flag){
    flag = 0;
    for (int i=0;i!=cluster_msc.size();i++){
      MergeSpaceCellSelection& mscs_1 = cluster_msc.at(i);
      int sum = 0;
      for (int j=0;j!=mscs_1.size();j++){
  	sum += mscs_1.at(j)->Get_all_spacecell().size();
      }
      if (sum < 5){
  	cluster_msc.erase(cluster_msc.begin() + i);
  	flag = 1;
  	break;
      }
    }
  }


  // std::cout << vertices.size() << " " << new_mcells.size() << " " << cluster_msc.size() << std::endl;
  // for (int i=0;i!=cluster_msc.size();i++){
  //   std::cout << cluster_msc.at(i).at(0)->Get_Center().x/units::cm << " " << cluster_msc.at(i).size() << std::endl;
  // }
  
  // Need to judge if the cluster is around a vertex or not (what criteria?)
  for (int i=0;i!=cluster_msc.size();i++){
    WCVertex *vertex;
    int flag = 0;
    for (int k=0;k!=cluster_msc.at(i).size();k++){
      MergeSpaceCell *mcell1 = cluster_msc.at(i).at(k);
      MergeSpaceCell *mcell2 = new_mcells_map[mcell1];
      for (int j=0;j!=vertices.size();j++){
	vertex = vertices.at(j);
	MergeSpaceCell *vertex_cell = vertex->get_msc();
	
	if (fabs(vertex_cell->Get_Center().x - mcell2->Get_Center().x) < 2*vertex_cell->thickness() + 0.5*units::mm && vertex_cell->Overlap(*mcell2)){
	  flag = 1;
	  break;
	}
      }
      if (flag == 1)
	break;
    }
    if (flag == 1){
      //std::cout << flag << " " << vertex->Center().x/units::cm << " " << std::endl;
    }else{
      WCTrack *track = new WCTrack(cluster_msc.at(i));
      short_tracks.push_back(track);
    }
  }

  

  // If not, do short track, just save the nmcell etc
  // if yes, do parallel track finder ... 
  


  //  
  


}


bool WireCell2dToy::ToyTracking::grow_track_fill_gap(WireCell2dToy::ToyCrawler& toycrawler){
  bool result = false;
  MergeSpaceCellMap& mcells_map = toycrawler.Get_mcells_map();
  
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    if (vertex->get_ntracks() == 1){
      WCTrack *track = vertex->get_tracks().at(0);
      MergeSpaceCell *vertex_cell = vertex->get_msc();
      MergeSpaceCellSelection cells = mcells_map[vertex_cell];

      MergeSpaceCellSelection saved_cells;

      if(cells.size() >1){
	TVector3 vec(1,vertex->get_ky(track),vertex->get_kz(track));
	float theta = vec.Theta();
	float phi  = vec.Phi();
	Point p = vertex->Center();

	MergeSpaceCell *prev_cell; //previous cell
	MergeSpaceCellSelection test_cells; // hold cells to be grow
	MergeSpaceCell *final_cell = vertex_cell;
	
	//std::cout << track->get_all_cells().size() << std::endl;

	for (int j=0;j!=cells.size();j++){
	  MergeSpaceCell *cell = cells.at(j);
	  auto it = find(track->get_all_cells().begin(), track->get_all_cells().end(), cell);
	  if (it == track->get_all_cells().end() && fabs(cell->Get_Center().x - final_cell->Get_Center().x)>0.1*units::mm ){
	    test_cells.push_back(cell);
	  }else{
	    prev_cell = cell;
	  }
	}

	//start to crawl ... 
	int flag = 1;
	while(flag){
	  flag = 0;
	  if (test_cells.size()==1){
	    prev_cell = final_cell;
	    final_cell = test_cells.at(0);
	    flag = 1;
	  }else{
	    for (int j=0;j!=test_cells.size();j++){
	      if (test_cells.at(j)->CrossCell(p,theta,phi)){
		prev_cell = final_cell;
		final_cell = test_cells.at(j);
		flag = 1;
		break;
	      }
	    }
	  }
	  
	  if (flag == 1){
	    //track->replace_end_scells(final_cell);
	    //	    track->get_all_cells().insert(track->get_all_cells().begin(),final_cell);
	    saved_cells.push_back(final_cell);

	    //test if final cell belong to another track, if so end
	    int flag1 = 0;
	    for (auto it1 = wct_wcv_map.begin(); it1 != wct_wcv_map.end(); it1++){
	      WCTrack *ntrack = it1->first;
	      auto it2 = find(ntrack->get_all_cells().begin(),ntrack->get_all_cells().end(), final_cell);
	      if (it2 != ntrack->get_all_cells().end()){
		flag1 = 1;
		break;
	      }
	    }
	    if (flag1 == 1){
	      flag = 0;
	    }else{
	      test_cells.clear();
	      for (int j=0;j!=mcells_map[final_cell].size();j++){
		if ((mcells_map[final_cell].at(j)->Get_Center().x - final_cell->Get_Center().x) * 
		    (final_cell->Get_Center().x - prev_cell->Get_Center().x) >0 
		    && fabs(mcells_map[final_cell].at(j)->Get_Center().x - final_cell->Get_Center().x) > 0.1*units::mm){
		  test_cells.push_back(mcells_map[final_cell].at(j));
		}
	      }
	    }
	  }
	}

	float dis1 = fabs(vertex_cell->Get_Center().x - track->get_end_scells().at(0)->Get_Center().x);
	float dis2 = fabs(vertex_cell->Get_Center().x - track->get_end_scells().at(1)->Get_Center().x);
	
	if (dis1 < dis2){
	  for (int qx = 0; qx!=saved_cells.size();qx++){
	    track->get_all_cells().insert(track->get_all_cells().begin(),saved_cells.at(qx));
	  }
	  //track->ReplaceEndCell(saved_cells.back(),track->get_end_scells().at(1));
	    // track->get_end_scells().clear();
	  // track->get_end_scells().push_back(center);
	  // track->get_end_scells().push_back(cells.back());
	}else{
	  for (int qx = 0; qx!=saved_cells.size();qx++){
	    track->get_all_cells().push_back(saved_cells.at(qx));
	  }
	  //	  track->ReplaceEndCell(saved_cells.back(),track->get_end_scells().at(0));
	  // cells.push_back(center);
	  // track->get_end_scells().clear();
	  // track->get_end_scells().push_back(cells.front());
	  // track->get_end_scells().push_back(center);
	
	}
	
	if (final_cell != vertex_cell){
	  vertex->set_msc(final_cell);
	  
	  result = true;
	  // std::cout << track->get_all_cells().size() << std::endl;
	  // std::cout << final_cell->Get_Center().x/units::cm << std::endl;
	  
	  
	  
	  Point p1;
	  p1.x = final_cell->Get_Center().x;
	  p1.y = vertex_cell->Get_Center().y + vertex->get_ky(track) * (p1.x - vertex_cell->Get_Center().x);
	  p1.z = vertex_cell->Get_Center().z + vertex->get_kz(track) * (p1.x - vertex_cell->Get_Center().x);
	  
	  // std::cout << vertex_cell << " " << final_cell << " " << vertex->Center().x/units::cm << " " << vertex->Center().y/units::cm << " " << vertex->Center().z/units::cm << " " 
	  // 	    << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << std::endl;
	  
	  vertex->set_center(p1);
	//	vertex->reset_center();
	// vertex->FindVertex();
	// std::cout << final_cell->Get_Center().x/units::cm << " " << 
	//       final_cell->Get_Center().y/units::cm << " " << 
	//       final_cell->Get_Center().z/units::cm << " " << std::endl;
	}
      }
    }
  }
  return result;
}


void WireCell2dToy::ToyTracking::deal_wiggle_tracks(){
  
  int flag1 = 1;
  WCTrackSelection used_tracks;
  
  while(flag1 ==1){
    flag1 = 0;
    for (auto it = type3_tracks.begin();it!=type3_tracks.end();it++){
      if (it->second == 2 ){
	WCTrack *track = it->first;
	
	if (find(used_tracks.begin(),used_tracks.end(),track) != used_tracks.end())
	  continue;

	Point end_p1 = track->get_end_scells().at(0)->Get_Center();
	Point end_p2 = track->get_end_scells().at(1)->Get_Center();
	
	MergeSpaceCell *vertex_cell1;
	MergeSpaceCell *vertex_cell2;
	
	float dis;
	
	int flag = 0;
	
	
	for (int i=0;i!=vertices.size();i++){
	  Point p = vertices.at(i)->get_msc()->Get_Center();
	  auto it = find(track->get_all_cells().begin(), track->get_all_cells().end(),vertices.at(i)->get_msc());
	  if (it != track->get_all_cells().end()){
	    float dis1 = sqrt(pow(end_p1.x - p.x,2) + pow(end_p1.y - p.y,2) + pow(end_p1.z - p.z,2));
	    float dis2 = sqrt(pow(end_p2.x - p.x,2) + pow(end_p2.y - p.y,2) + pow(end_p2.z - p.z,2));
	    if (dis1<dis2){
	      vertex_cell1 = vertices.at(i)->get_msc();
	      vertex_cell2 = track->get_end_scells().at(1);
	    }else{
	      vertex_cell1 = vertices.at(i)->get_msc();
	      vertex_cell2 = track->get_end_scells().at(0);
	    }
	    flag = 1;
	    break;
	  }
	  
	  
	  dis = sqrt(pow(end_p1.x - p.x,2) + pow(end_p1.y - p.y,2) + pow(end_p1.z - p.z,2));
	  //	  std::cout << dis/units::cm << std::endl;
	  if (dis < 2* units::cm && fabs(end_p1.x-p.x) < 0.65*units::cm ){
	    vertex_cell1 = track->get_end_scells().at(0);
	    vertex_cell2 = track->get_end_scells().at(1);
	    flag = 1;
	    break;
	  }
	  
	  dis = sqrt(pow(end_p2.x - p.x,2) + pow(end_p2.y - p.y,2) + pow(end_p2.z - p.z,2));
	  //std::cout << dis/units::cm << std::endl;
	  if (dis < 2* units::cm && fabs(end_p2.x-p.x) < 0.65*units::cm){
	    vertex_cell1 = track->get_end_scells().at(0);
	    vertex_cell2 = track->get_end_scells().at(1);
	    //vertex_cell = track->get_end_scells().at(1);
	    flag = 1;
	    break;
	  }
	}
	
	//std::cout << flag << " " << vertices.size() << std::endl;
	if (flag == 1){
	  flag1 = 1;
	  used_tracks.push_back(track);
	  
	  WCVertex *vertex = new WCVertex(*vertex_cell1);
	  vertex->get_tracks().push_back(track);
	  vertices.push_back(vertex);
	  vertex->FindVertex();
	  
	  WCVertex *vertex1 = new WCVertex(*vertex_cell2);
	  vertex1->get_tracks().push_back(track);
	  vertices.push_back(vertex1);
	  vertex1->FindVertex();
	}
	//std::cout << vertices.size() << std::endl;
	
	//      std::cout << track->get_end_scells().at(0)->Get_Center().x/units::cm << " " << 
	//	track->get_end_scells().at(1)->Get_Center().x/units::cm << std::endl;
      }
    }
  }


}

bool WireCell2dToy::ToyTracking::IsThisShower(WireCell2dToy::ToyCrawler& toycrawler){
  // find out how many separated tracks are there? 
  // 1. # of tracks
  int ntracks = 0;
  MergeSpaceCellSelection all_track_mscs;
  for (int i=0;i!=tracks.size();i++){
    auto it = find(good_tracks.begin(),good_tracks.end(),tracks.at(i));
    if (it != good_tracks.end()){
    //    if (tracks.at(i)->get_centerVP().size() > 0 ){
      all_track_mscs.insert(all_track_mscs.end(),tracks.at(i)->get_centerVP_cells().begin(),tracks.at(i)->get_centerVP_cells().end()) ; 
      ntracks ++;
    }
  }

  // find out how many blobs are not contained in the tracks
  // 1. How many blobs are not accounted in tracks (need a cut on size)
  int nmcells = 0;
  std::set<int> time_mcells_set;
  std::set<int> time_mcells_set1;
  std::set<int> time_mcells_set2;
  MergeSpaceCellMap mcells_map = toycrawler.Get_mcells_map();
  for (auto it = mcells_map.begin(); it!= mcells_map.end();it ++ ){
    MergeSpaceCell *mscell = it->first;
    auto it1 = find(all_track_mscs.begin(),all_track_mscs.end(),mscell);
    int time = mscell->Get_Center().x/mscell->thickness();
    if (it1 == all_track_mscs.end()){
      if (mscell->Get_all_spacecell().size() > 50){
	time_mcells_set.insert(time);
	nmcells ++;
      }
      time_mcells_set1.insert(time);
    }else{
      time_mcells_set2.insert(time);
    }
  }
  
  // 2. the connectivity of tracks  (how many big clusters)
  std::vector<WCTrackSelection> track_cluster;
  WCTrackSelection used_tracks;
  WCVertexSelection used_vertices;

  for (auto it = wct_wcv_map.begin(); it!= wct_wcv_map.end(); it++){
    WCTrack *track = it->first;
    auto it1 = find(used_tracks.begin(),used_tracks.end(),track);
    if (it1 == used_tracks.end()){
      used_tracks.push_back(track);
      WCTrackSelection cur_tracks;
      cur_tracks.push_back(track);
    
      
      // fill in vertices
      WCVertexSelection vertices = it->second;
      WCVertexSelection cur_vertices;
      cur_vertices.insert(cur_vertices.begin(),vertices.begin(),vertices.end());
      
      //std::cout << vertices.size() << " " << used_tracks.size() << " " << used_vertices.size() << std::endl;

      int flag = 1;
      while(flag){
  	flag = 0;
  	if (cur_vertices.size() > 0){
  	  WCVertex *vertex1 = cur_vertices.back();
	  cur_vertices.pop_back();
  	  flag = 1;
	  
	  //std::cout << wcv_wct_map[vertex1].size() << " " << used_vertices.size() << std::endl;

  	  auto it2 = find(used_vertices.begin(),used_vertices.end(),vertex1);
  	  if (it2 == used_vertices.end()){
  	    used_vertices.push_back(vertex1);
  	    // find all the tracks associated with vertex
	    // add them if not used, and not in current one
	    for (int j=0; j!=wcv_wct_map[vertex1].size();j++){
	      WCTrack *ntrack = wcv_wct_map[vertex1].at(j);
	      auto it3 = find(used_tracks.begin(),used_tracks.end(),ntrack);
	      if (it3 == used_tracks.end()){
		used_tracks.push_back(ntrack);
		cur_tracks.push_back(ntrack);
		// find all the vertices associated with these tracks
		for (int k=0;k!=wct_wcv_map[ntrack].size();k++){
		  WCVertex *nvertex = wct_wcv_map[ntrack].at(k);
		  auto it4 = find(used_vertices.begin(),used_vertices.end(),nvertex);
		  if (it4 == used_vertices.end()){
		    //used_vertices.push_back(nvertex);
		    cur_vertices.push_back(nvertex);
		  }
		}
	      }else{
		//	std::cout << ntrack << std::endl;
	      }
	    }
  	  }else{
	    //std::cout << vertex1 << std::endl;
	  }
  	}
	//std::cout << cur_vertices.size() << " a " << cur_tracks.size() << " " << used_tracks.size() << std::endl;
      } // end while

      track_cluster.push_back(cur_tracks);
    } // end if
    //    WCVertexSelection vertices = it->second;
  }
  
  // for (auto it =  wcv_wct_map.begin(); it!= wcv_wct_map.end();it++){
  //   std::cout << "Vertex " << " " << it->first->Center().x/units::cm << " " << it->second.size() << std::endl; 
  // }
  // for (auto it =  wct_wcv_map.begin(); it!= wct_wcv_map.end();it++){
  //   std::cout << "Track " << it->second.size() << std::endl; 
  // }


  std::cout << "Number of Tracks " << ntracks << " "  << time_mcells_set.size() << " " << track_cluster.size() << " " << time_mcells_set1.size() << " " << time_mcells_set2.size() << std::endl; 
  
  
  //need better tuning .... 
  // no tracks 
  if (ntracks == 0 && time_mcells_set.size() >= 4){
    return true;
  }else{
    //unacounted ones spam more than 10
    if (time_mcells_set.size() >=10){
      return true; 
    }else{
      if (track_cluster.size() >=2 && time_mcells_set.size() >=4){
	return true;
      }else{
	if (time_mcells_set1.size() > time_mcells_set2.size() * 0.75 && time_mcells_set1.size() > 5){
	  return true;
	}
      }
    }
  }
  


  return false;
  


}


void WireCell2dToy::ToyTracking::cleanup_bad_tracks(){
  WCVertexSelection to_be_removed;

  for (int i=0;i!=good_tracks.size();i++){
    if (good_tracks.at(i)->IsBadTrack()){      
      bad_tracks.push_back(good_tracks.at(i));
      //bad_tracks.at(i)->reset_fine_tracking();
      
      for (int j=0;j!=vertices.size();j++){
        WCVertex *vertex = vertices.at(j);
	WCTrackSelection& vtracks = vertex->get_tracks();
	auto it = find(vtracks.begin(),vtracks.end(),good_tracks.at(i));
	if (it != vtracks.end()){
	  vtracks.erase(it);
	}
	if (vtracks.size()==0){
	  //eliminate this vertex
	  auto it1 = find(to_be_removed.begin(),to_be_removed.end(),vertices.at(j));
	  if (it1 == to_be_removed.end())
	    to_be_removed.push_back(vertices.at(j));
	}
      }
    }  
  }

  for (int i=0;i!=bad_tracks.size();i++){
    auto it = find(good_tracks.begin(),good_tracks.end(),bad_tracks.at(i));
    good_tracks.erase(it);
  }

  // for (int i=0;i!=vertices.size();i++){
  //   auto it1 = find(to_be_removed.begin(),to_be_removed.end(),vertices.at(i));
  //   if (vertices.at(i)->get_ntracks()==0 && it1==to_be_removed.end())
  //     to_be_removed.push_back(vertices.at(i));
  // }

  for (int i=0;i!=to_be_removed.size();i++){
    WCVertex *vertex = to_be_removed.at(i);
    auto it = find(vertices.begin(),vertices.end(),vertex);
    if (it != vertices.end()){
      vertices.erase(it);
      delete vertex;
    }
  }

}

void WireCell2dToy::ToyTracking::update_maps1(){
  // update maps with good tracks ... 
  wct_wcv_map.clear();
  wcv_wct_map.clear();

  for (int i=0; i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    WCTrackSelection tracks_1 = vertex->get_tracks();

    // push good tracks in
    WCTrackSelection tracks1;
    for (int j=0;j!=tracks_1.size();j++){
      auto it =  find(good_tracks.begin(),good_tracks.end(),tracks_1.at(j));
      if (it != good_tracks.end()){
      //      if (tracks_1.at(j)->get_centerVP().size() > 0){
	tracks1.push_back(tracks_1.at(j));
      }
    }

    // judge the rest of tracks
    for (int j=0;j!=tracks.size();j++){
      auto it3 = find(good_tracks.begin(),good_tracks.end(),tracks.at(j));
      //  if (tracks.at(j)->get_centerVP().size()>0){
      if (it3 != good_tracks.end()){
	auto it = find(tracks1.begin(),tracks1.end(),tracks.at(j));
	if (it == tracks1.end()){
	  auto it1 = find(tracks.at(j)->get_centerVP_cells().begin(),
			  tracks.at(j)->get_centerVP_cells().end(),
			  vertex->get_msc());
	  if (it1 !=  tracks.at(j)->get_centerVP_cells().end()){
	    tracks1.push_back(tracks.at(j));
	  }
	}
      }

    }

    
    if (tracks1.size()>0){
      wcv_wct_map[vertex] = tracks1;
    
      for (int j=0;j!=tracks1.size();j++){
	WCTrack *track = tracks1.at(j);
	if (wct_wcv_map.find(track)==wct_wcv_map.end()){
	  WCVertexSelection vertexes;
	  vertexes.push_back(vertex);
	  wct_wcv_map[track] = vertexes;
	}else{
	  wct_wcv_map[track].push_back(vertex);
	}
      }
    }
  }
  


}


void WireCell2dToy::ToyTracking::update_maps(){
  wct_wcv_map.clear();
  wcv_wct_map.clear();

  for (int i=0; i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    WCTrackSelection tracks = vertex->get_tracks();
    wcv_wct_map[vertex] = tracks;
    
    for (int j=0;j!=tracks.size();j++){
      WCTrack *track = tracks.at(j);
      if (wct_wcv_map.find(track)==wct_wcv_map.end()){
	WCVertexSelection vertexes;
	vertexes.push_back(vertex);
	wct_wcv_map[track] = vertexes;
      }else{
	wct_wcv_map[track].push_back(vertex);
      }
    }
  }
  
  //protect against the fine_tracking
  for (auto it = wct_wcv_map.begin(); it!= wct_wcv_map.end();it++){
    WCTrack *track = it->first;
    WCVertexSelection& vertices1 = it->second;
  //   if (vertices1.size() == 1){
  //     MergeSpaceCell *esell1 = track->get_end_scells().at(0);
  //     MergeSpaceCell *esell2 = track->get_end_scells().at(1);
  //     MergeSpaceCell *esell;
  //     float dis1 = pow(esell1->Get_Center().x-vertices1.at(0)->get_msc()->Get_Center().x,2)
  // 	+pow(esell1->Get_Center().y-vertices1.at(0)->get_msc()->Get_Center().y,2)
  // 	+pow(esell1->Get_Center().z-vertices1.at(0)->get_msc()->Get_Center().z,2);
  //     float dis2 = pow(esell2->Get_Center().x-vertices1.at(0)->get_msc()->Get_Center().x,2)
  // 	+pow(esell2->Get_Center().y-vertices1.at(0)->get_msc()->Get_Center().y,2)
  // 	+pow(esell2->Get_Center().z-vertices1.at(0)->get_msc()->Get_Center().z,2);

  //     if (dis1 < dis2){
  // 	esell = esell2;
  //     }else{
  // 	esell = esell1;
  //     }
  //     //loop through all vertices
  //     WCVertex *min_vertex = vertices.at(0);
  //     float min_dis = pow(esell->Get_Center().x-min_vertex->get_msc()->Get_Center().x,2)
  // 	+pow(esell->Get_Center().y-min_vertex->get_msc()->Get_Center().y,2)
  // 	+pow(esell->Get_Center().z-min_vertex->get_msc()->Get_Center().z,2);
  //     for (int i = 0;i!=vertices.size();i++){
  // 	float dis3 = pow(esell->Get_Center().x-vertices.at(i)->get_msc()->Get_Center().x,2)
  // 	  +pow(esell->Get_Center().y-vertices.at(i)->get_msc()->Get_Center().y,2)
  // 	  +pow(esell->Get_Center().z-vertices.at(i)->get_msc()->Get_Center().z,2);
  // 	if (dis3 < min_dis){
  // 	  min_dis = dis3;
  // 	  min_vertex = vertices.at(i);
  // 	}
  //     }
  //     if (min_vertex != vertices1.at(0)){
  // 	vertices1.push_back(min_vertex);
  //     }
  //     if (wcv_wct_map.find(min_vertex) == wcv_wct_map.end()){
  // 	WCTrackSelection temp_tracks;
  // 	temp_tracks.push_back(track);
  // 	wcv_wct_map[min_vertex] = temp_tracks;
  //     }else{
  // 	wcv_wct_map[min_vertex].push_back(track);
  //     }

  //   }else 
    if (vertices1.size() > 2){
      MergeSpaceCell *esell1 = track->get_end_scells().at(0);
      MergeSpaceCell *esell2 = track->get_end_scells().at(1);

      WCVertex *min_vertex1 = vertices1.at(0);
      float min_dis1 = pow(esell1->Get_Center().x-min_vertex1->get_msc()->Get_Center().x,2)
  	+pow(esell1->Get_Center().y-min_vertex1->get_msc()->Get_Center().y,2)
  	+pow(esell1->Get_Center().z-min_vertex1->get_msc()->Get_Center().z,2);
      
      for (int i=0;i!=vertices1.size();i++){
  	float dis3 = pow(esell1->Get_Center().x-vertices1.at(i)->get_msc()->Get_Center().x,2)
  	  +pow(esell1->Get_Center().y-vertices1.at(i)->get_msc()->Get_Center().y,2)
  	  +pow(esell1->Get_Center().z-vertices1.at(i)->get_msc()->Get_Center().z,2);
  	if (dis3 < min_dis1){
  	  min_dis1 = dis3;
  	  min_vertex1 = vertices1.at(i);
  	}
      } 


      WCVertex *min_vertex2 = vertices1.at(0);
      float min_dis2 = pow(esell2->Get_Center().x-min_vertex2->get_msc()->Get_Center().x,2)
  	+pow(esell2->Get_Center().y-min_vertex2->get_msc()->Get_Center().y,2)
  	+pow(esell2->Get_Center().z-min_vertex2->get_msc()->Get_Center().z,2);
      
      for (int i=0;i!=vertices1.size();i++){
  	float dis3 = pow(esell2->Get_Center().x-vertices1.at(i)->get_msc()->Get_Center().x,2)
  	  +pow(esell2->Get_Center().y-vertices1.at(i)->get_msc()->Get_Center().y,2)
  	  +pow(esell2->Get_Center().z-vertices1.at(i)->get_msc()->Get_Center().z,2);
  	if (dis3 < min_dis2){
  	  min_dis2 = dis3;
  	  min_vertex2 = vertices1.at(i);
  	}
      } 

      
      //hack for now .... Not a good map ... 
      //for (auto it1 = vertices1.begin();it1!=vertices1.end();it1++){
      //	if (*it1!=min_vertex1 && *it1 !=min_vertex2)
      //	  it1 = vertices1.erase(it1);
      vertices1.clear();
      vertices1.push_back(min_vertex1);
      vertices1.push_back(min_vertex2);
      // }
    }
  }


}

void WireCell2dToy::ToyTracking::fine_tracking(){
  for (int i=0;i!=tracks.size();i++){
    // std::cout << i << " " << tracks.size() << std::endl;
    WCTrack *track = tracks.at(i);
    if (wct_wcv_map.find(track)!=wct_wcv_map.end()){
      //std::cout << wct_wcv_map[track].size() << std::endl;
      if (wct_wcv_map[track].size()==2){
	//std::cout << i << "abc1 " << std::endl;
	WCVertex *vertex1 = wct_wcv_map[track].at(0);
	//std::cout << i << "abc2 " << std::endl;
	WCVertex *vertex2 = wct_wcv_map[track].at(1);
	//std::cout << i << "abc3 " << std::endl;
	if (vertex1 == vertex2) continue;
	Point p1 = vertex1->Center();
	Point p2 = vertex2->Center();

	// std::cout << p1.x/units::cm << " " << p1.y/units::cm << " " <<p1.z/units::cm 
	// 	  << " " << p2.x/units::cm << " " << p2.y/units::cm << " " <<p2.z/units::cm 
	// 	  << std::endl;

	// find the directions 
	double ky1, kz1, ky2, kz2;
	ky1 = vertex1->get_ky(track);
	kz1 = vertex1->get_kz(track);
	
	ky2 = vertex2->get_ky(track);
	kz2 = vertex2->get_kz(track);
	int np1 = vertex1->get_ntracks();
	int np2 = vertex2->get_ntracks();
	//std::cout << i << "abc4 " << std::endl;
	track->reset_fine_tracking();
	track->fine_tracking(np1,p1,ky1,kz1,np2,p2,ky2,kz2);
	good_tracks.push_back(track);
	//std::cout << i << "abc5 " << std::endl;
      }
    }
  }
}

bool WireCell2dToy::ToyTracking::ExamineVertex(WCVertex* vertex, WireCell2dToy::ToyCrawler& toycrawler){
  Point vertex_location = vertex->Center();
  
  if (fabs(vertex_location.x - vertex->get_msc()->Get_Center().x) < 0.32*5*units::cm){
  }else{
    vertex->reset_center();
    return false;
  }

  MergeSpaceCellSelection cells1;

  for (int i=0;i!=vertex->get_tracks().size();i++){
    WCTrack *track = vertex->get_tracks().at(i);
    for (int j=0;j!=track->get_all_cells().size();j++){
      MergeSpaceCell *cell1 = track->get_all_cells().at(j);
      cells1.push_back(cell1);
    }
  }
  
  MergeSpaceCell *vertex_cell = toycrawler.GetClosestMSC(vertex_location,cells1);
  //std::cout << vertex_location.x/units::cm << " " << vertex_cell->Get_Center().x/units::cm << std::endl;
  
  // std::cout << vertex_cell << " " << cells1.size() << std::endl;

  if (vertex_cell!=0){
    vertex->set_msc(vertex_cell);
    
    // Now need to examine the vertex ... 
    MergeSpaceCell *center = vertex->get_msc();
    WCTrackSelection tracks = vertex->get_tracks();
    MergeSpaceCellMap& mcells_map = toycrawler.Get_mcells_map();
    
    WCTrackSelection tracks_to_be_removed;
    
    for (int j=0;j!=tracks.size();j++){
      WCTrack *track = tracks.at(j);
      MergeSpaceCellSelection& cells = track->get_all_cells();
      
      auto it = find(cells.begin(),cells.end(),center);
      if (it == cells.end()){  //does not contain the cell
	// need to fill in ... 
	float dis1 = fabs(center->Get_Center().x - cells.front()->Get_Center().x);
	float dis2 = fabs(center->Get_Center().x - cells.back()->Get_Center().x);
	MergeSpaceCell *cell1; // the one ... 
	if (dis1 < dis2){
	  cell1 = cells.front();
	}else{
	  cell1 = cells.back();
	}
	
	double x1 = cell1->Get_Center().x;
	double y1 = cell1->Get_Center().y;
	double z1 = cell1->Get_Center().z;
	double ky = (center->Get_Center().y - y1)/(center->Get_Center().x - x1);
	double kz = (center->Get_Center().z - z1)/(center->Get_Center().x - x1);
	double dis = center->Get_Center().x - x1;
	
	MergeSpaceCell *prev_cell2 = 0;
	
	MergeSpaceCell *cell2 = cell1;
	while( (cell2->Get_Center().x - cell1->Get_Center().x) * (cell2->Get_Center().x - center->Get_Center().x) < 0 || cell2 == cell1){
	  
	  
	  //	std::cout << cell2->Get_Center().x/units::cm << " " << cell2 << std::endl;
	  
	  // start to crawl ... 
	  MergeSpaceCellSelection cells2 = mcells_map[cell2];
	  
	  MergeSpaceCell *min_cell3;
	  double min = 1e9;
	  
	  for (int k=0;k!=cells2.size();k++){
	    MergeSpaceCell *cell3 = cells2.at(k);
	    if (dis * (cell3->Get_Center().x - cell2->Get_Center().x) < 0)
	      continue;
	    
	    double dy = cell3->get_dy();
	    double dz = cell3->get_dz();
	    
	    if (dy == 0) dy = 0.3/2 * units::cm;
	    if (dz == 0) dz = 0.3/2 * units::cm;
	    
	    double yp = y1 + ky * (cell3->Get_Center().x-x1);
	    double zp = z1 + kz * (cell3->Get_Center().x-x1);
	    double dis_sigma = pow(cell3->Get_Center().y - yp,2)/pow(dy,2)
	      + pow(cell3->Get_Center().z - zp,2)/pow(dz,2);
	    
	    if (dis_sigma < min){
	      min = dis_sigma;
	      min_cell3 = cell3;
	    }
	  }
	  prev_cell2 = cell2;
	  cell2 = min_cell3;
	  
	  if (cell2 == 0){
	    tracks_to_be_removed.push_back(track);
	    break;
	  }else{
	    if (fabs(cell2->Get_Center().x - center->Get_Center().x) < 0.1 *units::mm && cell2 != center){
	      tracks_to_be_removed.push_back(track);
	      break;
	    }
	  }
	  if (cell2 == prev_cell2) break;
	} //while loop
	
      }
    }
    
    
    
    if (tracks_to_be_removed.size()>0){
      std::cout << "Remove Tracks!!! " << std::endl;
      for (int i = 0 ;i!=tracks_to_be_removed.size();i++){
	auto it = find(vertex->get_tracks().begin(),vertex->get_tracks().end(),tracks_to_be_removed.at(i));
	vertex->get_tracks().erase(it);
      }
      return false;
    }else{
      return true;
    }
  }
  return false;
}



void WireCell2dToy::ToyTracking::CheckVertices(WireCell2dToy::ToyCrawler& toycrawler){
  
  MergeSpaceCellMap& mcells_map = toycrawler.Get_mcells_map();
  
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    MergeSpaceCell *center = vertex->get_msc();
    Point pcenter = vertex->Center();
    WCTrackSelection tracks = vertex->get_tracks();

    // std::cout << center->Get_Center().x/units::cm << " " << tracks.size() << std::endl;

    //deal the case where center is vertex ...
    // if (center->Get_Center().x == pcenter.x &&
    // 	center->Get_Center().y == pcenter.y &&
    // 	center->Get_Center().z == pcenter.z){
    for (int j=0;j!=tracks.size();j++){
      WCTrack *track = tracks.at(j);
      MergeSpaceCellSelection& cells = track->get_all_cells();
      auto it = find(cells.begin(),cells.end(),center);

      // std::cout << center->Get_Center().x/units::cm << " " << it - cells.end() << std::endl;

      if (it == cells.end() ){
	//(cells.at(0)->Get_Center().x - center->Get_Center().x) * (cells.at(cells.size()-1)->Get_Center().x - center->Get_Center().x)>0){
	// need to fill in ... 
	float dis1 = fabs(center->Get_Center().x - cells.front()->Get_Center().x);
	float dis2 = fabs(center->Get_Center().x - cells.back()->Get_Center().x);
	MergeSpaceCell *cell1; // the one ... 
	if (dis1 < dis2){
	  cell1 = cells.front();
	}else{
	  cell1 = cells.back();
	}
	
	double x1 = cell1->Get_Center().x;
	double y1 = cell1->Get_Center().y;
	double z1 = cell1->Get_Center().z;
	
	double ky = (center->Get_Center().y - y1)/(center->Get_Center().x - x1);
	double kz = (center->Get_Center().z - z1)/(center->Get_Center().x - x1);
	
	double dis = center->Get_Center().x - x1;
	
	MergeSpaceCell *prev_cell2 = 0;
	MergeSpaceCell *cell2 = cell1;
	
	
	while( ((cell2->Get_Center().x - cell1->Get_Center().x) * (cell2->Get_Center().x - center->Get_Center().x) < 0 || cell2 == cell1) && (cell2 !=0)){
	  // start to do the surgery ... 
	  MergeSpaceCellSelection cells2 = mcells_map[cell2];
	  MergeSpaceCell *min_cell3 = 0;
	  
	  double min = 1e9;
	  
	  for (int k=0;k!=cells2.size();k++){
	    MergeSpaceCell *cell3 = cells2.at(k);
	    if (dis * (cell3->Get_Center().x - cell2->Get_Center().x) < 0)
	      continue;
	    
	    double dy = cell3->get_dy();
	    double dz = cell3->get_dz();
	    
	    if (dy == 0) dy = 0.3/2 * units::cm;
	    if (dz == 0) dz = 0.3/2 * units::cm;
	    
	    double yp = y1 + ky * (cell3->Get_Center().x-x1);
	    double zp = z1 + kz * (cell3->Get_Center().x-x1);
	    double dis_sigma = pow(cell3->Get_Center().y - yp,2)/pow(dy,2)
	      + pow(cell3->Get_Center().z - zp,2)/pow(dz,2);
	    
	    if (dis_sigma < min){
	      min = dis_sigma;
	      min_cell3 = cell3;
	    }
	    
	  }
	  prev_cell2 = cell2;
	  cell2 = min_cell3;
	 
	  if (cell2 !=0){
	    // put it in  ...
	    if ((cell2->Get_Center().x - cell1->Get_Center().x) * (cell2->Get_Center().x - center->Get_Center().x) < 0){
	      if (dis1 < dis2){
		cells.insert(cells.begin(),cell2);
	      }else{
		cells.push_back(cell2);
	      }
	    }
	  }
	  if (cell2 == prev_cell2 || cell2 == 0) break;
	}
	
	// put the final center in ...
	if (dis1 < dis2){
	  cells.insert(cells.begin(),center);
	  track->get_end_scells().clear();
	  track->get_end_scells().push_back(center);
	  track->get_end_scells().push_back(cells.back());
	}else{
	  cells.push_back(center);
	  track->get_end_scells().clear();
	  track->get_end_scells().push_back(cells.front());
	  track->get_end_scells().push_back(center);
	}
	  
      }else{
	
	

      }
    }
    // }else{
    //   //after fit .... 
    // }
    

  }
}


void WireCell2dToy::ToyTracking::CleanUpVertex(){
  WCVertexSelection temp;
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    if (vertex->get_ntracks()==0){
      temp.push_back(vertex);
    }else if (vertex->get_ntracks()==1){
      MergeSpaceCell *cell = vertex->get_msc();
      auto it = find(vertex->get_tracks().at(0)->get_all_cells().begin(),
		     vertex->get_tracks().at(0)->get_all_cells().end(),
		     cell);
      if (it == vertex->get_tracks().at(0)->get_all_cells().end()){
	temp.push_back(vertex);
      }
    }
  }

  for (int i=0;i!=temp.size();i++){
    auto it = find(vertices.begin(),vertices.end(),temp.at(i));
    vertices.erase(it);
  }

}


void WireCell2dToy::ToyTracking::OrganizeTracks(int flag){
  MergeSpaceCellSelection saved_cells;
  for (int i=0;i!=vertices.size();i++){
    saved_cells.push_back(vertices.at(i)->get_msc());
  } 

  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    vertex->OrganizeEnds(saved_cells, flag);
    //std::cout << vertex->get_ntracks() << " " << vertex->Center().x/units::cm << " " << vertex->Center().y/units::cm << " " << vertex->Center().z/units::cm << " " << std::endl;
  }
}


void WireCell2dToy::ToyTracking::Crawl(){
  // first shift the vertex locations
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);

    //if (i==1){
    vertex->OrganizeTracks();
      //std::cout << vertex->get_ntracks() << " " << vertex->Center().x/units::cm << " " << vertex->Center().y/units::cm << " " << vertex->Center().z/units::cm << std::endl;
      // }
  }
}

void WireCell2dToy::ToyTracking::BreakTracks(){
  // Now need to break the track?
   

  // first break the track if a vertice is in the middle
  WCTrackSelection break_tracks;
  int flag = 1;
  while(flag){
    flag = 0;
    for(int i=0;i!=vertices.size();i++){
      WCVertex *vertex1 = vertices.at(i);
      break_tracks = vertex1->BreakTracks();
      if (break_tracks.size()==2){
       	flag = 1;
       	tracks.push_back(break_tracks.at(1));
       	for (int j=0;j!=vertices.size();j++){
       	  WCVertex *vertex2 = vertices.at(j);
       	  if (vertex2 != vertex1){
       	    vertex2->ProcessTracks(break_tracks);
       	  }
       	}
       	break;
      }
    }
  }


  // second break the track if there is a direction change ... 
  break_tracks.clear();
  WCTrackSelection finished_tracks;
  WCVertexSelection NewVertices;
  
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    break_tracks = vertex->BreakTracksAngle(finished_tracks);
    if (break_tracks.size()>0){

      for (int j=1;j!=break_tracks.size();j++){
      	tracks.push_back(break_tracks.at(j));
      }
      // Form new vertices
      for (int j=0;j!=break_tracks.size()-1;j++){
      	//WCVertex *vertex2 = FormNewVertex(break_tracks.at(j),break_tracks.at(j+1));
      	WCVertex *vertex2 = new WCVertex(*(break_tracks.at(j)->get_all_cells().back()));
      	vertex2->get_tracks().push_back(break_tracks.at(j));
      	vertex2->get_tracks().push_back(break_tracks.at(j+1));
      	NewVertices.push_back(vertex2);
      }


      // need to deal with other vertices ...
      WCTrackSelection break_tracks1;
      break_tracks1.push_back(break_tracks.at(0));
      break_tracks1.push_back(break_tracks.at(break_tracks.size()-1));
      for (int j=0;j!=vertices.size();j++){
      	WCVertex *vertex2 = vertices.at(j);
      	if (vertex2 != vertex){
      	  vertex2->ProcessTracks(break_tracks1);
      	}
      }
      }
  }


  for (int i=0;i!=NewVertices.size();i++){
    vertices.push_back(NewVertices.at(i));
  }


}


void WireCell2dToy::ToyTracking::MergeVertices(int flag){
  WCVertexSelection to_be_removed;
  //prepare to merge vertices 
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex1 = vertices.at(i);

    if (flag!=2)
      if (vertex1->get_ntracks()==1) continue;

    auto it1 = find(to_be_removed.begin(),to_be_removed.end(),vertex1);
    if (it1 == to_be_removed.end()){
      // std::cout << i << " " << vertex1->get_ntracks() << std::endl;
      for (int j=0;j!=vertices.size();j++){
  	WCVertex *vertex2 = vertices.at(j);

	if (flag !=2)
	  if (vertex2->get_ntracks()==1 ) continue;

	if (vertex2 == vertex1) continue;
  	auto it2 = find(to_be_removed.begin(),to_be_removed.end(),vertex2);
  	if (it2 == to_be_removed.end()){
  	  //std::cout << i << " " << j << " " << vertex2->get_ntracks() << std::endl;
  	  if (vertex1->AddVertex(vertex2, flag)){
	    
	    // std::cout << "remove2 " << vertex1->Center().x/units::cm << " " <<
	    //   vertex1->Center().y/units::cm << " " << vertex1->Center().z/units::cm << " " <<
	    //   vertex2->Center().x/units::cm << " " <<
	    //   vertex2->Center().y/units::cm << " " << vertex2->Center().z/units::cm << " " <<std::endl;
	    
  	    to_be_removed.push_back(vertex2);
  	  }
	  //	  if (vertex1->AddVertex(vertex2, 3)){
	  //  to_be_removed.push_back(vertex2);
  	  //}

  	}
      }
    }
  }
  
  
  for (int i=0;i!=to_be_removed.size();i++){
    WCVertex *vertex = to_be_removed.at(i);
    auto it = find(vertices.begin(),vertices.end(),vertex);
    vertices.erase(it);
  }
}


void WireCell2dToy::ToyTracking::RemoveSame(){
  // remove the contained one ... 
  WCVertexSelection to_be_removed;
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex1 = vertices.at(i);
    if (vertex1->get_ntracks()==1) continue;
    for (int j=0;j!=vertices.size();j++){
      WCVertex * vertex2 = vertices.at(j);
      if (vertex2->get_ntracks()== 1 || vertex2 == vertex1) continue;
      //if (vertex2 == vertex1) continue;
      auto it1 = find(to_be_removed.begin(),to_be_removed.end(),vertex2);
      if (it1 == to_be_removed.end()){
  	if (vertex1->IsInside(vertex2) >=0){
  	  // std::cout << "remove " << vertex1->Center().x/units::cm << " " <<
  	  //   vertex1->Center().y/units::cm << " " << vertex1->Center().z/units::cm << " " <<
  	  //   vertex2->Center().x/units::cm << " " <<
  	  //   vertex2->Center().y/units::cm << " " << vertex2->Center().z/units::cm << " " <<std::endl;
  	  auto it = find(to_be_removed.begin(),to_be_removed.end(),vertex1);   
  	  if (it == to_be_removed.end()){
  	    to_be_removed.push_back(vertex1);
  	  }
  	}
      }
    }
  }
  for (int i=0;i!=to_be_removed.size();i++){
    WCVertex *vertex = to_be_removed.at(i);
    auto it = find(vertices.begin(),vertices.end(),vertex);
    vertices.erase(it);
  }
}



void WireCell2dToy::ToyTracking::CreateVertices(ToyCrawler& toycrawler){
  std::map<MergeSpaceCell*,WCTrackSelection> msc_wct_map;

   // fill in tracks ... 
  for (int i=0;i!=toycrawler.Get_allMCT().size();i++){
    toycrawler.Get_allMCT().at(i)->Organize();

    WCTrack *track = new WCTrack(*toycrawler.Get_allMCT().at(i));
    tracks.push_back(track);
    for (int j=0;j!=track->get_all_cells().size();j++){
      MergeSpaceCell *cell = track->get_all_cells().at(j);
      cell->CalMinMax();
      //std::cout << cell->get_dy() << " " <<cell->get_dz() << std::endl;
    }


    for (int j=0;j!=track->get_end_scells().size();j++){
      MergeSpaceCell *cell = track->get_end_scells().at(j);
      

      auto it = msc_wct_map.find(cell);
      if (it == msc_wct_map.end()){
       	WCTrackSelection temp_tracks;
       	temp_tracks.push_back(track);
       	msc_wct_map[cell] = temp_tracks;
      }else{
       	msc_wct_map[cell].push_back(track);
      }
    }
    // std::cout << "Xin " << track->get_end_scells().size()<<std::endl;
  }

  // MergeSpaceCellSelection vertex_candidate;
  WCTrackSelection saved_tracks;

  //start to judge vertex ... 
  for (auto it = msc_wct_map.begin();it!=msc_wct_map.end();it++){
    MergeSpaceCell *cell = it->first;
    WCTrackSelection temp_tracks = it->second;
    WCVertex *vertex = new WCVertex(*cell);
    for (int i=0;i!=temp_tracks.size();i++){
      WCTrack *temp_track = temp_tracks.at(i);
      int type = temp_track->TrackType(*cell);
      //std::cout << "Xin " << type << " " << cell->Get_Center().x/units::cm << " " << cell->Get_Center().y/units::cm << " " << cell->Get_Center().z/units::cm << " " << temp_track->get_mct().Get_TimeLength() << std::endl;

      // save all the wiggle tracks to be dealt later ... 
      if (type==3){
	if (type3_tracks.find(temp_track)==type3_tracks.end()){
	  type3_tracks[temp_track] = 1;
	}else{
	  type3_tracks[temp_track] ++;
	}
      }
      //

      if (type==2){
      	vertex->Add(temp_track);
	auto it = find(saved_tracks.begin(),saved_tracks.end(),temp_track);
	if (it == saved_tracks.end()){
	  saved_tracks.push_back(temp_track);
	}
      }
    }
    if (vertex->get_ntracks() > 0){
      vertices.push_back(vertex);
      //     std::cout << vertex->get_ntracks() << std::endl;
    }else{
      delete vertex;
    }
    //  	break;
  }

  //complete the other side of the track? ... 
  MergeSpaceCellSelection existing_cells;
  for (int i = 0;i!=vertices.size();i++){
    existing_cells.push_back(vertices.at(i)->get_msc());
  }
  for (int i=0;i!=saved_tracks.size();i++){
     WCTrack *temp_track = saved_tracks.at(i);
     for (int j=0;j!=temp_track->get_end_scells().size();j++){
       MergeSpaceCell *cell = temp_track->get_end_scells().at(j);
       auto it = find(existing_cells.begin(),existing_cells.end(),cell);
       if (it == existing_cells.end()){
  	 WCVertex *vertex = new WCVertex(*cell);
  	 vertex->Add(temp_track);
  	 vertices.push_back(vertex);
       }
     }
  }



  // Now need to match tracks with vertices ... 
  for (int i= 0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    MergeSpaceCell *cell = vertex->get_msc();
    WCTrackSelection& tracks = vertex->get_tracks();
    for (int j=0;j!=saved_tracks.size();j++){
      WCTrack *temp_track = saved_tracks.at(j);
      auto it = find(tracks.begin(),tracks.end(),temp_track);
      if (it == tracks.end()){
	//	MergeClusterTrack& mct = temp_track->get_mct();
	MergeSpaceCellSelection& cells = temp_track->get_all_cells();
	auto it1 = find(cells.begin(),cells.end(),cell);
	if (it1 != cells.end()){
	  vertex->Add(temp_track);
	}
      }else{
	continue;
      }
    }
  }

}

void WireCell2dToy::ToyTracking::Associate(){
  WCTrackSelection saved_tracks;
  for (int i=0;i!=vertices.size();i++){
    for (int j=0;j!=vertices.at(i)->get_tracks().size();j++){
      saved_tracks.push_back(vertices.at(i)->get_tracks().at(j));
    }
  }


  // Now need to match tracks with vertices ... 
  for (int i= 0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    MergeSpaceCell *cell = vertex->get_msc();
    WCTrackSelection& tracks = vertex->get_tracks();
    for (int j=0;j!=saved_tracks.size();j++){
      WCTrack *temp_track = saved_tracks.at(j);
      auto it = find(tracks.begin(),tracks.end(),temp_track);
      if (it == tracks.end()){
	MergeSpaceCellSelection& cells = temp_track->get_all_cells();
	auto it1 = find(cells.begin(),cells.end(),cell);
	if (it1 != cells.end()){
	  vertex->Add(temp_track);
	}
      }else{
	continue;
      }
    }
  }
}


WireCell2dToy::ToyTracking::~ToyTracking(){
  for (int i=0;i!=tracks.size();i++){
    delete tracks.at(i);
  }
  for (int i=0;i!=vertices.size();i++){
    delete vertices.at(i);
  }
}
