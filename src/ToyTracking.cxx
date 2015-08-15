#include "WireCell2dToy/ToyTracking.h"

using namespace WireCell;

WireCell2dToy::ToyTracking::ToyTracking(WireCell2dToy::ToyCrawler& toycrawler){
  CreateVertices(toycrawler);
  RemoveSame(); // get rid of duplicated ones
  

  // MergeVertices();   // merge things that are common ...
  // OrganizeTracks();  // improve the end points nearby and break things
  // Associate();    //associate the rest
   
  // MergeVertices(0); //merge things, but require the center must be at the end of common track
  // OrganizeTracks(); //improve the end points nearby and break things
  // Associate(); //associate the rest
  

  // MergeVertices(2);  // merge some vertices with distance cut only 
  // Crawl();
  // OrganizeTracks(); //associate the rest
  // Associate(); //associate the rest
  
  
  // MergeVertices(2);  // merge some vertices together, allow single track
  // BreakTracks();    //improve the end points and break things
  // OrganizeTracks(); //associate the rest


  // Associate();  //associate the rest .. 
  // CleanUpVertex(); //do some final clean up ... 
  

  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    std::cout << i << " " << vertex->get_ntracks() << " " << vertex->Center().x/units::cm << " " << vertex->Center().y/units::cm << " " << vertex->Center().z/units::cm << " " << std::endl;
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
  // first shif the vertex locations
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    vertex->OrganizeTracks();
  }
}

void WireCell2dToy::ToyTracking::BreakTracks(){
  // Now need to break the track?
  
  
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
	  if (vertex2->get_ntracks()==1) continue;
	if (vertex2 == vertex1) continue;
  	auto it2 = find(to_be_removed.begin(),to_be_removed.end(),vertex2);
  	if (it2 == to_be_removed.end()){
  	  //std::cout << i << " " << j << " " << vertex2->get_ntracks() << std::endl;
  	  if (vertex1->AddVertex(vertex2, flag)){
	    // std::cout << "remove2 " << vertex1->Center().x/units::cm << " " <<
  	    // vertex1->Center().y/units::cm << " " << vertex1->Center().z/units::cm << " " <<
  	    // vertex2->Center().x/units::cm << " " <<
  	    // vertex2->Center().y/units::cm << " " << vertex2->Center().z/units::cm << " " <<std::endl;
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
    WCTrack *track = new WCTrack(*toycrawler.Get_allMCT().at(i));
    tracks.push_back(track);
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


WireCell2dToy::ToyTracking::~ToyTracking(){
  for (int i=0;i!=tracks.size();i++){
    delete tracks.at(i);
  }
  for (int i=0;i!=vertices.size();i++){
    delete vertices.at(i);
  }
}
