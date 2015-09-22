#include "WireCell2dToy/ToyCosmic.h"
#include "TVector3.h"

using namespace WireCell;

WireCell2dToy::ToyCosmic::ToyCosmic(WireCell2dToy::ToyTrackingSelection& trackings)
  : trackings(trackings)
{
  ToyTrackingSelection used_trackings;
  std::vector<ToyTrackingSelection> cosmic_candidates;
  
  while(used_trackings.size() != trackings.size()){
    ToyTracking *curr_tracking;
    for (int i=0;i!=trackings.size();i++){
      curr_tracking = trackings.at(i);
      auto it = find(used_trackings.begin(),used_trackings.end(),curr_tracking);
      if (it == used_trackings.end())
	break;
    }
    used_trackings.push_back(curr_tracking);
    
    //push it in anyway
    ToyTrackingSelection temp;
    temp.push_back(curr_tracking);
    cosmic_candidates.push_back(temp);

    //loop over all the trackings to add in
    int flag = 1;
    while(flag){
      flag = 0;
      for (int i=0;i!=trackings.size();i++){
	ToyTracking *curr_tracking = trackings.at(i);
	auto it = find(used_trackings.begin(),used_trackings.end(),curr_tracking);
	if (it != used_trackings.end()) continue;
	for (int j=0;j!=temp.size();j++){
	  if (IsConnected(temp.at(j),curr_tracking)){
	    temp.push_back(curr_tracking);
	    used_trackings.push_back(curr_tracking);
	    flag = 1;
	    break;
	  }
	}
	if (flag == 1) break;
      }
    } // end while loop
  }
  
  std::cout << trackings.size() << " " << cosmic_candidates.size() << std::endl;
}

bool WireCell2dToy::ToyCosmic::IsConnected(ToyTracking *tracking1, ToyTracking *tracking2){
  WCTrackSelection& tracking1_tracks = tracking1->get_good_tracks();
  WCTrackSelection& tracking2_tracks = tracking2->get_good_tracks();
  
  MergeSpaceCellSelection tracking1_mcells;
  MergeSpaceCellSelection tracking2_mcells;
  
  std::map<WCTrack*, std::vector<float>> tracks_angle_map;

  for (int i = 0;i!=tracking1_tracks.size();i++){
    WCTrack *track = tracking1_tracks.at(i);
    for (int j=0;j!=track->get_centerVP_cells().size();j++){
      tracking1_mcells.push_back(track->get_centerVP_cells().at(j));
    }
    
    tracks_angle_map[track] = track->get_direction();
    //std::vector<float> abc = track->get_direction();
    //std::cout << abc.at(0) << " " << abc.at(1) << " " << abc.at(2) << std::endl;
  }
  
  for (int i = 0;i!=tracking2_tracks.size();i++){
    WCTrack *track = tracking2_tracks.at(i);
    for (int j=0;j!=track->get_centerVP_cells().size();j++){
      tracking2_mcells.push_back(track->get_centerVP_cells().at(j));
    }
    tracks_angle_map[track] = track->get_direction();
    //std::vector<float> abc = track->get_direction();
    //std::cout << abc.at(0) << " " << abc.at(1) << " " << abc.at(2) << std::endl;
  }

  
  
  // std::cout << tracking1_tracks.size() << " " << tracking2_tracks.size() << " " << tracking1_mcells.size() << " " << tracking2_mcells.size() << std::endl;


  return false;
}

WireCell2dToy::ToyCosmic::~ToyCosmic(){
  
}
