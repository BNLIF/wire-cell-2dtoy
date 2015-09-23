#include "WireCell2dToy/ToyCosmic.h"
#include "TVector3.h"

using namespace WireCell;

WireCell2dToy::ToyCosmic::ToyCosmic(WireCell2dToy::ToyTrackingSelection& trackings)
  : trackings(trackings)
{
  ToyTrackingSelection used_trackings;
  cosmic_candidates.clear();
  
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
	for (int j=0;j!=cosmic_candidates.back().size();j++){
	  if (IsConnected(cosmic_candidates.back().at(j),curr_tracking)){
	    cosmic_candidates.back().push_back(curr_tracking);
	    used_trackings.push_back(curr_tracking);
	    flag = 1;
	    break;
	  }
	}
	if (flag == 1) break;
      }
    } // end while loop
  }
  
  //  std::cout << trackings.size() << " " << cosmic_candidates.size() << std::endl;

  for (int i = 0;i!=cosmic_candidates.size();i++){
    WCCosmic *cosmic = new WCCosmic(cosmic_candidates.at(i));
    cosmics.push_back(cosmic);
  }


  //further merge things
  int flag = 1;
  while(flag){
    flag = 0;
    WCCosmicSelection temp;
    for (int i=0;i!=cosmics.size();i++){
      WCCosmic *cosmic1 = cosmics.at(i);
      if (cosmic1->get_points().size()==0) continue;
      // if (cosmic1->get_points().size()==0){
      // 	flag = 1;
      // 	temp.push_back(cosmic1);
      // 	break;
      // }
      for (int j=0;j!=cosmics.size();j++){
	WCCosmic *cosmic2 = cosmics.at(j);
	if (cosmic1 == cosmic2) continue;
	if (cosmic2->get_points().size() == 0) continue;
	// if (cosmic2->get_points().size() == 0){
	//   flag = 1;
	//   temp.push_back(cosmic1);
	//   break;
	// }
	  
	// judge wheter merge
	float dis = cosmic1->cal_dist(cosmic2);
	float angle = cosmic1->cal_costh(cosmic2);
	
	if (dis < 10*units::cm && fabs(angle) > 0.94){
	  Point p1_f = cosmic1->get_points().front();
	  Point p1_b = cosmic1->get_points().back();
	  Point p2_f = cosmic2->get_points().front();
	  Point p2_b = cosmic2->get_points().back();
	  
	  int flag1 = 0;
	  if (sqrt(pow(p1_f.x-p2_f.x,2) +pow(p1_f.y-p2_f.y,2) + pow(p1_f.z-p2_f.z,2)) < 10*units::cm){
	    flag1 = 1;
	  }
	  if (flag1 == 0)
	    if (sqrt(pow(p1_f.x-p2_b.x,2) +pow(p1_f.y-p2_b.y,2) + pow(p1_f.z-p2_b.z,2)) < 10*units::cm){
	      flag1 = 1;
	    }
	  if(flag1 == 0)
	    if (sqrt(pow(p1_b.x-p2_f.x,2) +pow(p1_b.y-p2_f.y,2) + pow(p1_b.z-p2_f.z,2)) < 10*units::cm){
	      flag1 = 1;
	    }
	  if (flag1==0)
	    if (sqrt(pow(p1_b.x-p2_b.x,2) +pow(p1_b.y-p2_b.y,2) + pow(p1_b.z-p2_b.z,2)) < 10*units::cm){
	      flag1 = 1;
	    }

	  if (flag1==1){
	    cosmic1->Add(cosmic2);
	    temp.push_back(cosmic2);
	    flag = 1;
	    break;
	  }
	}
      }
      if (flag==1) break;
    }
    
    for (int i=0;i!=temp.size();i++){
      auto it = find(cosmics.begin(),cosmics.end(),temp.at(i));
      delete *it;
      cosmics.erase(it);
    }
  }
  

  // judge if anything is cosmics ... if not delete them ...  
  WCCosmicSelection temp;
  for (int i=0;i!=cosmics.size();i++){
    if (!cosmics.at(i)->IsCosmic())
      temp.push_back(cosmics.at(i));
  }
  for (int i=0;i!=temp.size();i++){
    auto it = find(cosmics.begin(),cosmics.end(),temp.at(i));
    delete *it;
    cosmics.erase(it);
  }
  
  // Now loop over all other toytrackings to see if any one is saved

}

bool WireCell2dToy::ToyCosmic::IsConnected(ToyTracking *tracking1, ToyTracking *tracking2){
  //return false;
  WCTrackSelection& tracking1_tracks = tracking1->get_good_tracks();
  WCTrackSelection& tracking2_tracks = tracking2->get_good_tracks();
  
  std::map<WCTrack*, std::vector<float>> tracks_angle_map;
  std::map<WCTrack*, std::vector<float>> tracks_pos_map;
  
  for (int i = 0;i!=tracking1_tracks.size();i++){
    WCTrack *track = tracking1_tracks.at(i);
    tracks_angle_map[track] = track->get_direction();
    tracks_pos_map[track] = track->get_position();
  }
  
  for (int i = 0;i!=tracking2_tracks.size();i++){
    WCTrack *track = tracking2_tracks.at(i);
    tracks_angle_map[track] = track->get_direction();
    tracks_pos_map[track] = track->get_position();
  }

  for (int i=0;i!=tracking1_tracks.size();i++){
    WCTrack *track1 = tracking1_tracks.at(i);
    TVector3 dir1(tracks_angle_map[track1].at(0),
		  tracks_angle_map[track1].at(1),
		  tracks_angle_map[track1].at(2));
    TVector3 pos1(tracks_pos_map[track1].at(0),
		  tracks_pos_map[track1].at(1),
		  tracks_pos_map[track1].at(2));
    
    // Point track1_p1(tracks_pos_map[track1].at(0),tracks_pos_map[track1].at(1),tracks_pos_map[track1].at(2));
    // Point track1_p2(tracks_pos_map[track1].at(0) + tracks_direction_map[track1].at(0),
    // 		    tracks_pos_map[track1].at(1) + tracks_direction_map[track1].at(1),
    // 		    tracks_pos_map[track1].at(2) + tracks_direction_map[track1].at(2));
    // Line l1(track1_p1, track1_p2);
    
    for (int j=0;j!=tracking2_tracks.size();j++){
      WCTrack *track2 = tracking2_tracks.at(j);
      TVector3 dir2(tracks_angle_map[track2].at(0),
		    tracks_angle_map[track2].at(1),
		    tracks_angle_map[track2].at(2));
      TVector3 pos2(tracks_pos_map[track2].at(0),
		    tracks_pos_map[track2].at(1),
		    tracks_pos_map[track2].at(2));

      TVector3 perp1 = dir1.Cross(dir2);
      TVector3 perp2 = pos1 - pos2;
      float abc = perp1.Dot(perp2);

      float dis=1e9;
      if (perp1.Mag() == 0){
	TVector3 perp3 = perp2.Cross(dir1);
	dis = perp3.Mag()/dir1.Mag();
      }else{
	dis = fabs(abc/perp1.Mag());
      }
      float angle = fabs(dir1.Dot(dir2)/dir1.Mag()/dir2.Mag());

      //std::cout << "qx: " << dis << " " << angle << std::endl;

      float dis_cut=3.0*units::cm;
      if (dis < 5*units::cm && angle >= 0.97){
       	dis_cut = 12*units::cm;
      }else if (dis < 10*units::cm && angle > 0.94){
	dis_cut = 6 * units::cm;
      }

      // MergeSpaceCell *mcell1_f = track1->get_centerVP_cells().front();
      // MergeSpaceCell *mcell1_b = track1->get_centerVP_cells().back();
      // MergeSpaceCell *mcell2_f = track2->get_centerVP_cells().front();
      // MergeSpaceCell *mcell2_b = track2->get_centerVP_cells().back();

      // if (IsConnected(mcell1_f,mcell2_f,dis_cut)) return true;
      // if (IsConnected(mcell1_f,mcell2_b,dis_cut)) return true;
      // if (IsConnected(mcell1_b,mcell2_f,dis_cut)) return true;
      // if (IsConnected(mcell1_b,mcell2_b,dis_cut)) return true;

      // check five cells
      for (int k1 = 0; k1!= track1->get_centerVP_cells().size();k1++){
	if (k1>=5 && k1 + 5 < track1->get_centerVP_cells().size()) continue;
	MergeSpaceCell *mcell1 = track1->get_centerVP_cells().at(k1);
	for (int k2 = 0; k2!= track2->get_centerVP_cells().size();k2++){
	  if (k2>=5 && k2 + 5 < track2->get_centerVP_cells().size()) continue;
	  MergeSpaceCell *mcell2 = track2->get_centerVP_cells().at(k2);
	  if (IsConnected(mcell1,mcell2,dis_cut)) return true;
	}
      }
    }
  }
  
  // std::cout << tracking1_tracks.size() << " " << tracking2_tracks.size() << " " << tracking1_mcells.size() << " " << tracking2_mcells.size() << std::endl;
  
  return false;
}

bool WireCell2dToy::ToyCosmic::IsConnected(MergeSpaceCell *mcell1, MergeSpaceCell *mcell2, float dis_cut){
  float dy1 = mcell1->get_dy();
  float dz1 = mcell1->get_dz();
  Point p1 = mcell1->Get_Center();

  float dy2 = mcell2->get_dy();
  float dz2 = mcell2->get_dz();
  Point p2 = mcell2->Get_Center();

  int flag_x = 0;
  int flag_y = 0; 
  int flag_z = 0;
  
  if (fabs(p1.x-p2.x)<dis_cut) flag_x = 1;
  if (flag_x == 0 ) return false;
  
  float min_y = 1e9;
  float max_y = -1e9;
  
  if (min_y > p1.y - dy1) min_y = p1.y - dy1;
  if (max_y < p1.y + dy1) max_y = p1.y + dy1;  
  if (min_y > p2.y - dy2) min_y = p2.y - dy2;
  if (max_y < p2.y + dy2) max_y = p2.y + dy2;

  if (max_y - min_y < dis_cut + dy1*2 + dy2*2)
    flag_y = 1;
  
  if (flag_y ==0 )return false;
  
  float min_z = 1e9;
  float max_z = -1e9;
  
  
  if (min_z > p1.z - dz1) min_z = p1.z - dz1;
  if (max_z < p1.z + dz1) max_z = p1.z + dz1;  
  if (min_z > p2.z - dz2) min_z = p2.z - dz2;
  if (max_z < p2.z + dz2) max_z = p2.z + dz2;

  if (max_z - min_z < dis_cut + dz1*2 + dz2*2)
    flag_z = 1;

  
  if (flag_x == 1 && flag_y == 1 && flag_z == 1) 
    return true;
  
  return false;
}


WireCell2dToy::ToyCosmic::~ToyCosmic(){
  
}
