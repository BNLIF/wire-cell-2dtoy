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

  MergeSpaceCellSelection vertex_candidate;
  
  //start to judge vertex ... 
  for (auto it = msc_wct_map.begin();it!=msc_wct_map.end();it++){
    MergeSpaceCell *cell = it->first;
    WCTrackSelection temp_tracks = it->second;
    for (int i=0;i!=temp_tracks.size();i++){
      WCTrack *temp_track = temp_tracks.at(i);
      int type = temp_track->TrackType(*cell);
      //std::cout << "Xin " << type << " " << cell->Get_Center().x/units::cm << " " << cell->Get_Center().y/units::cm << " " << cell->Get_Center().z/units::cm << std::endl;
      if (type==2){
      	WCVertex *vertex = new WCVertex(*cell);
      	vertices.push_back(vertex);
      	break;
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
