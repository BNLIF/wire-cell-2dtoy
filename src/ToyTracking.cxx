#include "WireCell2dToy/ToyTracking.h"

using namespace WireCell;
WireCell2dToy::ToyTracking::ToyTracking(WireCell2dToy::ToyCrawler& toycrawler){
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

  // Now need to clean up vertices ... 
  for (int i= 0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    MergeSpaceCell *cell = vertex->get_msc();
    WCTrackSelection& tracks = vertex->get_tracks();

    for (int j=0;j!=saved_tracks.size();j++){
      WCTrack *temp_track = saved_tracks.at(j);
      auto it = find(tracks.begin(),tracks.end(),temp_track);
      if (it == tracks.end()){
	MergeClusterTrack& mct = temp_track->get_mct();
	MergeSpaceCellSelection& cells = mct.Get_allmcells();
	auto it1 = find(cells.begin(),cells.end(),cell);
	if (it1 != cells.end()){
	  vertex->Add(temp_track);
	}
      }else{
	continue;
      }
    }
    
    // std::cout << vertex->get_ntracks() << std::endl;
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
