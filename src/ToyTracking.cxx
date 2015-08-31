#include "WireCell2dToy/ToyTracking.h"

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
  //CheckVertices(toycrawler);
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
    //if (i==20)
    bool success = vertex->FindVertex();

    if (success){
      if (!ExamineVertex(vertex,toycrawler)){
    	bool success1 = vertex->FindVertex();
    	if (success1){
    	  ExamineVertex(vertex,toycrawler);
    	}
      }
    }

    std::cout << i << " Vertex " << vertex->get_ntracks() << " " << success << " " << vertex->Center().x/units::cm << " " << vertex->Center().y/units::cm << " " << vertex->Center().z/units::cm << " " << std::endl;
    // for (int j=0;j!=vertex->get_ntracks();j++){
    //   std::cout << vertex->get_tracks().at(j)->get_end_scells().at(0)->Get_Center().x/units::cm << " " 
    // 		<< vertex->get_tracks().at(j)->get_end_scells().at(1)->Get_Center().x/units::cm << " " 
    // 		<< vertex->get_tracks().at(j)->get_all_cells().front()->Get_Center().x/units::cm << " " 
    // 		<< vertex->get_tracks().at(j)->get_all_cells().back()->Get_Center().x/units::cm << " " 
    // 		<< std::endl;
    // }

  }
  CheckVertices(toycrawler);
  OrganizeTracks(); 


  
  
  //Now do fine tracking??? 

}



bool WireCell2dToy::ToyTracking::ExamineVertex(WCVertex* vertex, WireCell2dToy::ToyCrawler& toycrawler){
  Point vertex_location = vertex->Center();
  MergeSpaceCellSelection cells1;
  for (int i=0;i!=vertex->get_tracks().size();i++){
    WCTrack *track = vertex->get_tracks().at(i);
    for (int j=0;j!=track->get_all_cells().size();j++){
      MergeSpaceCell *cell1 = track->get_all_cells().at(j);
      cells1.push_back(cell1);
    }
  }
  
  MergeSpaceCell *vertex_cell = toycrawler.GetClosestMSC(vertex_location,cells1);
  
  // std::cout << vertex_cell << " " << cells1.size() << std::endl;

  //if (vertex_cell!=0)
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
      
      MergeSpaceCell *cell2 = cell1;
      while( (cell2->Get_Center().x - cell1->Get_Center().x) * (cell2->Get_Center().x - center->Get_Center().x) < 0 || cell2 == cell1){
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
      }
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



void WireCell2dToy::ToyTracking::CheckVertices(WireCell2dToy::ToyCrawler& toycrawler){
  
  MergeSpaceCellMap& mcells_map = toycrawler.Get_mcells_map();
  
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    MergeSpaceCell *center = vertex->get_msc();
    Point pcenter = vertex->Center();
    WCTrackSelection tracks = vertex->get_tracks();

    //deal the case where center is vertex ...
    // if (center->Get_Center().x == pcenter.x &&
    // 	center->Get_Center().y == pcenter.y &&
    // 	center->Get_Center().z == pcenter.z){
    for (int j=0;j!=tracks.size();j++){
      WCTrack *track = tracks.at(j);
      MergeSpaceCellSelection& cells = track->get_all_cells();
      auto it = find(cells.begin(),cells.end(),center);
      if (it == cells.end()){
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
      //std::cout << "Xin " << type << " " << cell->Get_Center().x/units::cm << " " << cell->Get_Center().y/units::cm << " " << cell->Get_Center().z/units::cm << std::endl;
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
