#include "WCP2dToy/ToyTracking.h"
#include "WCP2dToy/ToyWalking.h"
#include "WCPData/Singleton.h"
#include "WCPData/TPCParams.h"

#include "TVector3.h"
#include "WCPData/Line.h"


using namespace WCP;

void WCP2dToy::ToyTracking::create_new_tracks_from_leftover(WCP2dToy::ToyCrawler& toycrawler){
  // find out the mcells not in the current tracks ...

  MergeSpaceCellSelection all_track_mscs;
  for (int i=0;i!=tracks.size();i++){
    auto it = find(good_tracks.begin(),good_tracks.end(),tracks.at(i));
    auto it1 = find(parallel_tracks.begin(),parallel_tracks.end(),tracks.at(i));
    if (it != good_tracks.end() || it1 != parallel_tracks.end()){
    //    if (tracks.at(i)->get_centerVP().size() > 0 ){
      all_track_mscs.insert(all_track_mscs.end(),tracks.at(i)->get_centerVP_cells().begin(),tracks.at(i)->get_centerVP_cells().end()) ; 
      //save associated cells as well
      all_track_mscs.insert(all_track_mscs.end(),tracks.at(i)->get_associated_cells().begin(),tracks.at(i)->get_associated_cells().end());
    }
  }
  
  MergeSpaceCellSelection leftover_mcells;
  
  MergeSpaceCellMap mcells_map = toycrawler.Get_mcells_map();
  for (auto it = mcells_map.begin(); it!= mcells_map.end();it ++ ){
    MergeSpaceCell *mscell = it->first;
    auto it1 = find(all_track_mscs.begin(),all_track_mscs.end(),mscell);
    if (it1 == all_track_mscs.end())
      leftover_mcells.push_back(mscell);
  }
  
  // std::cout << leftover_mcells.size() << std::endl;
 
  // go through everything and form clusters ...
  std::vector<MergeSpaceCellSelection> clusters;
  while (leftover_mcells.size() > 0){
    MergeSpaceCell *mcell = leftover_mcells.back();
    leftover_mcells.pop_back();
    
    MergeSpaceCellSelection cluster;
    
    cluster.push_back(mcell);

    MergeSpaceCellSelection to_be_searched;
    to_be_searched.push_back(mcell);
    
    while(to_be_searched.size()>0){
      MergeSpaceCell *mcell1 = to_be_searched.back();
      to_be_searched.pop_back();

      MergeSpaceCellSelection& mcells = mcells_map[mcell1];
      for (int i1 = 0; i1!=mcells.size();i1++){
	MergeSpaceCell *mcell2 = mcells.at(i1);
	auto it = find(leftover_mcells.begin(),leftover_mcells.end(),mcell2);
	auto it1 = find(cluster.begin(),cluster.end(),mcell2);
	if (it != leftover_mcells.end() ){
	  leftover_mcells.erase(it);
	  cluster.push_back(mcell2);
	  to_be_searched.push_back(mcell2);
	}
      }

    }
    clusters.push_back(cluster);
    
  }

  // std::cout << leftover_mcells.size() << " " << clusters.size() << std::endl;
  // int sum = 0;
  // for (int i=0;i!=clusters.size();i++){
  //   std::cout << i << " " << clusters.at(i).size() << std::endl;
  //   sum += clusters.at(i).size();
  // }
  // std::cout << sum << std::endl;
  // if cluster large than certain size, create tracks 
  for (int i=0;i!=clusters.size();i++){
    MergeSpaceCellSelection mcells = clusters.at(i);
    int sum = 0;
    Point center(0,0,0);
    for (int j=0;j!=mcells.size();j++){
      center.x += mcells.at(j)->Get_Center().x *mcells.at(j)->Get_all_spacecell().size();
      center.y += mcells.at(j)->Get_Center().y *mcells.at(j)->Get_all_spacecell().size();
      center.z += mcells.at(j)->Get_Center().z *mcells.at(j)->Get_all_spacecell().size();
      sum += mcells.at(j)->Get_all_spacecell().size();
    }
    if (mcells.size()>=3&& sum > 50){
      // form tracks ... 
      // find the center
      center.x /= sum;
      center.y /= sum;
      center.z /= sum;
      // find the furthest mcell
      float furthest_dist = 0;
      MergeSpaceCell *furthest_mcell = 0;
      SpaceCell *furthest_cell = 0;
      
      for (int j=0;j!=mcells.size();j++){
	MergeSpaceCell *mcell = mcells.at(j);
	for (int k=0;k!=mcell->Get_all_spacecell().size();k++){
	  SpaceCell *cell = mcell->Get_all_spacecell().at(k);
	  float dist = sqrt(pow(cell->x() - center.x,2)
			    +pow(cell->y() - center.y,2)
			    +pow(cell->z() - center.z,2));
	
	  if (dist > furthest_dist){
	    furthest_dist = dist;
	    furthest_mcell = mcell;
	    furthest_cell = cell;
	  }
	}
      }
      
      //      std::cout << furthest_mcell << " " << mcells.size() << std::endl;

      
      // create the vertex
      
      WCVertex *vertex = new WCVertex(*furthest_mcell);
      
      //      std::cout << "xin 1:" << mcells_map[furthest_mcell].size() << std::endl;
      
      Point p(furthest_cell->x(),furthest_cell->y(),furthest_cell->z());
      vertex->set_center(p);
      vertices.push_back(vertex);

      //examine all the cells to find out outlier 
      for (int j=0;j!=mcells.size();j++){
	MergeSpaceCell *mcell = mcells.at(j);
	new_mcells_map[mcell] = mcell;
      }

      // do parallel tracking ... 
      parallel_tracking(vertex,mcells,toycrawler);
      
    }
    //    std::cout << i << " " << sum << " " << mcells.size() << std::endl; 
  }
}


void WCP2dToy::ToyTracking::IterateMergeTracks(WCP::MergeSpaceCellMap& mcells_map){
  

  int flag = 1;
  std::cout << "Iterate Merge Tracks! " << std::endl;
  while(flag){
    int prev_ntracks = good_tracks.size()+parallel_tracks.size();
    // std::cout << prev_ntracks << " a " << good_tracks.size() << std::endl;
    MergeTracks(mcells_map);
    // std::cout << prev_ntracks << " b " << good_tracks.size() << std::endl;
    if (prev_ntracks == good_tracks.size()+parallel_tracks.size()) flag = 0;
  }

  flag = 1;
  while(flag){
    int prev_ntracks = good_tracks.size()+parallel_tracks.size();
    MergeTracks_no_shared_vertex(mcells_map);
    MergeTracks_no_shared_vertex(mcells_map,1);
    MergeTracks_no_shared_vertex(mcells_map,2);
    if (prev_ntracks == good_tracks.size()+parallel_tracks.size()) flag = 0;
  }

  
  // for (int i=0;i!=vertices.size();i++){
  //   if (mcells_map[vertices.at(i)->get_msc()].size()==0){
  //     std::cout << i << " " << std::endl;
  //   }
  // }



  for (int i=0;i!=good_tracks.size();i++){
    WCTrack *track = good_tracks.at(i);
    WCVertexSelection vertices = wct_wcv_map[track];
    std::cout << "R: " << vertices.size() << " " << vertices.at(0)->Center().x/units::cm << " " << vertices.at(0)->Center().y/units::cm << " " << vertices.at(0)->Center().z/units::cm << " "
  	      << vertices.at(1)->Center().x/units::cm << " "  << vertices.at(1)->Center().y/units::cm << " "  << vertices.at(1)->Center().z/units::cm << " " << 
      std::endl;
  }

  for (int i=0;i!=parallel_tracks.size();i++){
    WCTrack *track = parallel_tracks.at(i);
    WCVertexSelection vertices = wct_wcv_map[track];
    std::cout << "P: " << vertices.size() << " " << vertices.at(0)->Center().x/units::cm << " " << vertices.at(0)->Center().y/units::cm << " " << vertices.at(0)->Center().z/units::cm << " "
  	      << vertices.at(1)->Center().x/units::cm << " "  << vertices.at(1)->Center().y/units::cm << " "  << vertices.at(1)->Center().z/units::cm << " " << 
      std::endl;
  }

  
}

void WCP2dToy::ToyTracking::MergeTracks_no_shared_vertex(WCP::MergeSpaceCellMap& mcells_map, int type){
  std::cout << "Merge Tracks shared no vertex " << type << std::endl;
  
  WCTrackSelection all_tracks;
  for (int i=0;i!=good_tracks.size();i++){
    if (good_tracks.at(i)->get_centerVP_cells().size()>0)
      all_tracks.push_back(good_tracks.at(i));
  }
  for (int i=0;i!=parallel_tracks.size();i++){
    if (parallel_tracks.at(i)->get_centerVP_cells().size()>0)
      all_tracks.push_back(parallel_tracks.at(i));
  }
  //std::cout << all_tracks.size() << std::endl;

   ClusterTrackSelection ct_tracks;
  for (int i=0;i!=all_tracks.size();i++){
    WCTrack *track = all_tracks.at(i);
    MergeSpaceCellSelection& mcells = track->get_centerVP_cells();
    //std::cout << i << " " << mcells.size() << std::endl;
    ClusterTrack *ct = new ClusterTrack(mcells.at(0));
    for (int j=1;j< mcells.size();j++){
      ct->AddMSCell_anyway(mcells.at(j));
    }
    ct_tracks.push_back(ct);
  }
  // find the pair vertex that are close or if they share mcell?  
  // satisfy the angle cut
  std::vector<WCTrackSelection> merge_tracks_array;
  std::vector<WCVertexSelection> merge_vertices_array;
  

  WCTrackSelection used_tracks;
  for (int i=0;i!=all_tracks.size();i++){
    WCTrack *track1 = all_tracks.at(i);
    used_tracks.push_back(track1);
    WCVertexSelection vertices1 = wct_wcv_map[track1];
    for (int j=0;j!=all_tracks.size();j++){
      if (i==j) continue;
      WCTrack *track2 = all_tracks.at(j);
      if (find(used_tracks.begin(),used_tracks.end(),track2) !=used_tracks.end()) continue;
      WCVertexSelection vertices2 = wct_wcv_map[track2];

      for (int i1 = 0;i1!=vertices1.size();i1++){
	WCVertex *vertex1 = vertices1.at(i1);
	for (int j1 = 0; j1!=vertices2.size(); j1++){
	  WCVertex *vertex2 = vertices2.at(j1);
	  
	  // find out the other vertex
	  WCVertex *vertex1_other = 0;
	  for (int i2 = 0;i2!=vertices1.size();i2++){
	    vertex1_other = vertices1.at(i2);
	    if (vertex1_other != vertex1) break;
	    }
	  WCVertex *vertex2_other = 0;
	  for (int i3 = 0;i3!=vertices2.size();i3++){
	    vertex2_other = vertices2.at(i3);
	    if (vertex2_other != vertex2) break;
	  }
	  if (vertex1_other==0 || vertex2_other ==0) continue;
	  
	  float dis_1 = sqrt(pow(vertex1->Center().x - vertex2->Center().x,2)+
			      pow(vertex1->Center().y - vertex2->Center().y,2)+
			      pow(vertex1->Center().z - vertex2->Center().z,2));
	  float dis_2 = sqrt(pow(vertex1->Center().x - vertex2_other->Center().x,2)+
			      pow(vertex1->Center().y - vertex2_other->Center().y,2)+
			      pow(vertex1->Center().z - vertex2_other->Center().z,2));
	  float dis_3 = sqrt(pow(vertex1_other->Center().x - vertex2->Center().x,2)+
			      pow(vertex1_other->Center().y - vertex2->Center().y,2)+
			      pow(vertex1_other->Center().z - vertex2->Center().z,2));
	  float dis_4 = sqrt(pow(vertex1_other->Center().x - vertex2_other->Center().x,2)+
			      pow(vertex1_other->Center().y - vertex2_other->Center().y,2)+
			      pow(vertex1_other->Center().z - vertex2_other->Center().z,2));
	  TVector3 kk1(vertex1_other->Center().x-vertex1->Center().x,
		       vertex1_other->Center().y-vertex1->Center().y,
		       vertex1_other->Center().z-vertex1->Center().z);
	  TVector3 kk2(vertex2_other->Center().x-vertex2->Center().x,
		       vertex2_other->Center().y-vertex2->Center().y,
		       vertex2_other->Center().z-vertex2->Center().z);

	  //std::cout << dis << " " << vertex1->Center().x/units::cm << " " << vertex2->Center().x/units::cm << " " << track1->get_range() << " " << track2->get_range() << " " << vertex1_other->Center().x/units::cm << " " << vertex2_other->Center().x/units::cm << " " << track1 << " " << track2 << " " << std::endl;

	  if (dis_1 > dis_4 || kk1.Angle(kk2) <3.1415926/2.) continue;
	  if (vertex1 == vertex2) continue;

	  if (type == 0){
	    float dis = sqrt(pow(vertex1->Center().x-vertex2->Center().x,2) + 
			     pow(vertex1->Center().y-vertex2->Center().y,2) + 
			     pow(vertex1->Center().z-vertex2->Center().z,2)); 
	    if (dis > 5*units::cm) continue;
	    //std::cout << dis << " " << vertex1->Center().x/units::cm << " " << vertex2->Center().x/units::cm << " " << track1->get_range() << " " << track2->get_range() << " " << vertex1_other->Center().x/units::cm << " " << vertex2_other->Center().x/units::cm << " " << track1 << " " << track2 << std::endl;
	    
	    
	  }else if (type == 1 || type==2){
	    
	    
	    if (type == 1){
	      
	      float dis = fabs(vertex1->Center().x - vertex2->Center().x);
	      if (dis > 3*units::cm) continue;
	      if (track1->get_range()<3*units::cm || track2->get_range()<3*units::cm) continue;
	      //std::cout << dis << " " << vertex1->Center().x/units::cm << " " << vertex2->Center().x/units::cm << " " << track1->get_range() << " " << track2->get_range() << " " << vertex1_other->Center().x/units::cm << " " << vertex2_other->Center().x/units::cm << " " << track1 << " " << track2 << std::endl;
	    }
	    if (type==2){
	      float dis =  sqrt(pow(vertex1->Center().x-vertex2->Center().x,2) + 
			     pow(vertex1->Center().y-vertex2->Center().y,2) + 
			     pow(vertex1->Center().z-vertex2->Center().z,2)); 
	      
	      //std::cout << dis << " " << track1->get_range() << " " << track2->get_range() << std::endl;
	      if (dis > track1->get_range() || dis > track2->get_range()) continue;
	      
	      int flag_overlap = 0;
	      float min_dis = 1e9;
	      for (int i2 = 0;i2!=track1->get_centerVP_cells().size();i2++){
		MergeSpaceCell *mcell1 = track1->get_centerVP_cells().at(i2);
		for (int j2 = 0;j2!=track2->get_centerVP_cells().size();j2++){
		  MergeSpaceCell *mcell2 = track2->get_centerVP_cells().at(j2);
		  float dis10 = sqrt(pow(mcell1->Get_Center().x-mcell2->Get_Center().x,2)+
				     pow(mcell1->Get_Center().y-mcell2->Get_Center().y,2)+
				     pow(mcell1->Get_Center().z-mcell2->Get_Center().z,2));
		  if (dis10 < min_dis){
		    min_dis = dis10;
		  }
		}
	      }

	      for (int i2 = 0;i2!=track1->get_centerVP_cells().size();i2++){
		MergeSpaceCell *mcell1 = track1->get_centerVP_cells().at(i2);
		for (int j2 = 0;j2!=track2->get_centerVP_cells().size();j2++){
		  MergeSpaceCell *mcell2 = track2->get_centerVP_cells().at(j2);
		  if (mcell1 == mcell2) {
		    flag_overlap ++;
		    break;
		  }
		}
	      }
	      //std::cout << flag_overlap << " " << track1->get_range() << " " << track2->get_range() << " " << min_dis/units::cm << std::endl;
	      
	      if (flag_overlap<=1) {
		if (min_dis < 1*units::cm && track1->get_range() > 5*units::cm && track2->get_range()>5*units::cm){
		}else{
		  continue;
		}
	      }
		

	    }
	    
	    
	   
	    
	    TVector3 abc1(vertex1_other->Center().x - vertex1->Center().x,
			  vertex1_other->Center().y - vertex1->Center().y,
			  vertex1_other->Center().y - vertex1->Center().z);
	    TVector3 abc2(vertex2->Center().x - vertex1->Center().x,
			  vertex2->Center().y - vertex1->Center().y,
			  vertex2->Center().y - vertex1->Center().z);

	    TVector3 abc3(vertex2_other->Center().x - vertex2->Center().x,
			  vertex2_other->Center().y - vertex2->Center().y,
			  vertex2_other->Center().y - vertex2->Center().z);
	    TVector3 abc4(vertex1->Center().x - vertex2->Center().x,
			  vertex1->Center().y - vertex2->Center().y,
			  vertex1->Center().y - vertex2->Center().z);

	    // std::cout << vertex1->Center().x/units::cm << " " << vertex1_other->Center().x/units::cm << " " 
	    // 	      << vertex2->Center().x/units::cm << " " << vertex2_other->Center().x/units::cm << " " 
	    // 	      << abc1.Angle(abc2)/3.1415926*180. << " " << abc3.Angle(abc4)/3.1415926*180. << std::endl;
	    
	    if ((abc1.Angle(abc2)/3.1415926*180.<15&& abc3.Angle(abc4)/3.1415926*180. < 15.)||
		sqrt(pow(abc1.Angle(abc2)/3.1415926*180.,2)+pow(abc3.Angle(abc4)/3.1415926*180.,2))<20.){
	    }else{
	      continue;
	    }
	  }
	  
	  auto it1 = find(all_tracks.begin(),all_tracks.end(),track1);
	  auto it2 = find(all_tracks.begin(),all_tracks.end(),track2);

	  Point center1;
	  center1.x = vertex1->Center().x;
	  center1.y = vertex1->Center().y;
	  center1.z = vertex1->Center().z;
	  
	  Point center2;
	  center2.x = vertex2->Center().x;
	  center2.y = vertex2->Center().y;
	  center2.z = vertex2->Center().z;

	  ct_tracks.at(it1 - all_tracks.begin())->SC_Hough(center1);
	  ct_tracks.at(it2 - all_tracks.begin())->SC_Hough(center2);

	  double theta1 = ct_tracks.at(it1 - all_tracks.begin())->Get_Theta();
	  double phi1 = ct_tracks.at(it1 - all_tracks.begin())->Get_Phi();

	  double theta2 = ct_tracks.at(it2 - all_tracks.begin())->Get_Theta();
	  double phi2 = ct_tracks.at(it2 - all_tracks.begin())->Get_Phi();
	  
	  TVector3 abc1(sin(theta1)*cos(phi1),sin(theta1)*sin(phi1),cos(theta1));
	  TVector3 abc2(sin(theta2)*cos(phi2),sin(theta2)*sin(phi2),cos(theta2));
	  TVector3 abc3(vertex1_other->Center().x - vertex1->Center().x,
	  		vertex1_other->Center().y - vertex1->Center().y,
	  		vertex1_other->Center().z - vertex1->Center().z);
	  TVector3 abc4(vertex2_other->Center().x - vertex2->Center().x,
	  		vertex2_other->Center().y - vertex2->Center().y,
	  		vertex2_other->Center().z - vertex2->Center().z);

	  double theta_abc = abc1.Angle(abc2);
	  double theta_abc1 = abc3.Angle(abc4);
	  //std::cout << ct_tracks.at(it1 - all_tracks.begin())->Get_allmcells().size() << " " << ct_tracks.at(it2 - all_tracks.begin())->Get_allmcells().size() << " " << center1.x/units::cm << " " << center2.x/units::cm << std::endl;
	  
	  // std::cout << vertex1->Center().x/units::cm << " " << vertex1->Center().y/units::cm << " " << vertex1->Center().z/units::cm << " " << 
	  //   vertex1_other->Center().x/units::cm << " " << vertex1_other->Center().y/units::cm << " " << vertex1_other->Center().z/units::cm << std::endl;
	  // std::cout << vertex2->Center().x/units::cm << " " << vertex2->Center().y/units::cm << " " << vertex2->Center().z/units::cm << " " << 
	  //   vertex2_other->Center().x/units::cm << " " << vertex2_other->Center().y/units::cm << " " << vertex2_other->Center().z/units::cm << std::endl;
	  // for (int kk = 0 ; kk!=ct_tracks.at(it1 - all_tracks.begin())->Get_allmcells().size();kk++){
	  //   std::cout << ct_tracks.at(it1 - all_tracks.begin())->Get_allmcells().at(kk)->Get_Center().x/units::cm << std::endl;
	  // }
	  // std::cout << "haha " << std::endl;
	  // for (int kk = 0 ; kk!=ct_tracks.at(it2 - all_tracks.begin())->Get_allmcells().size();kk++){
	  //   std::cout << ct_tracks.at(it2 - all_tracks.begin())->Get_allmcells().at(kk)->Get_Center().x/units::cm << std::endl;
	  // }
	  

	  // std::cout << center1.x/units::cm << " " << center1.y/units::cm << " " << center1.z/units::cm << std::endl;
	  // std::cout << center2.x/units::cm << " " << center2.y/units::cm << " " << center2.z/units::cm << std::endl;
	  // std::cout << abc1.x() << " " << abc1.y() << " " << abc1.z() << std::endl;
	  // std::cout << abc2.x() << " " << abc2.y() << " " << abc2.z() << std::endl;
	  //std::cout << theta1 << " " << phi1 << " " << theta2 << " " << phi2 << std::endl;
	  //std::cout << theta_abc/3.1415926*180. << std::endl;
	  //	  if (fabs(theta1+theta2-3.1415926) < 20./180.*3.1415926 && 
	  //  fabs(fabs(phi1-phi2)-3.1415926)<20./180.*3.1415926 &&
	  //   sqrt(pow(theta1+theta2-3.1415926,2) + pow(fabs(phi1-phi2)-3.1415926,2)) < 28./180.*3.1415926){
	  
	  // std::cout << theta_abc/3.1415926*180. << " " << theta_abc1/3.1415926*180. << std::endl;

	  if (fabs(theta_abc/3.1415926*180-180.)>15){
	    ct_tracks.at(it1 - all_tracks.begin())->SC_Hough(center1,15*units::cm);
	    ct_tracks.at(it2 - all_tracks.begin())->SC_Hough(center2,15*units::cm);
	    theta1 = ct_tracks.at(it1 - all_tracks.begin())->Get_Theta();
	    phi1 = ct_tracks.at(it1 - all_tracks.begin())->Get_Phi();

	    theta2 = ct_tracks.at(it2 - all_tracks.begin())->Get_Theta();
	    phi2 = ct_tracks.at(it2 - all_tracks.begin())->Get_Phi();
	  
	    abc1.SetXYZ(sin(theta1)*cos(phi1),sin(theta1)*sin(phi1),cos(theta1));
	    abc2.SetXYZ(sin(theta2)*cos(phi2),sin(theta2)*sin(phi2),cos(theta2));
	    theta_abc = abc1.Angle(abc2);
	  }
	  //std::cout << theta_abc/3.1415926*180. << " " << theta_abc1/3.1415926*180. << std::endl;
	  
	  if ((fabs(theta_abc/3.1415926*180-180.)<15 || fabs(theta_abc1/3.1415926*180-180.)<15)||
	      (fabs(theta_abc/3.1415926*180-180.)<25 && fabs(theta_abc1/3.1415926*180-180.)<25)){
	    if (merge_tracks_array.size()==0){
	      WCTrackSelection merge_tracks;
	      merge_tracks.push_back(track1);
	      merge_tracks.push_back(track2);
	      WCVertexSelection merge_vertices;
	      merge_vertices.push_back(vertex1);
	      merge_vertices.push_back(vertex2);
	      
	      merge_tracks_array.push_back(merge_tracks);
	      merge_vertices_array.push_back(merge_vertices);
	    }else{
	      int flag = 0;
	      for (int k = 0;k!=merge_tracks_array.size();k++){
		auto itt1 = find(merge_tracks_array.at(k).begin(),merge_tracks_array.at(k).end(),track1);
		auto itt2 = find(merge_tracks_array.at(k).begin(),merge_tracks_array.at(k).end(),track2);
		if (itt1 != merge_tracks_array.at(k).end() || itt2 != merge_tracks_array.at(k).end()){
		  if (itt1 == merge_tracks_array.at(k).end())
		    merge_tracks_array.at(k).push_back(track1);
		  if (itt2 == merge_tracks_array.at(k).end())
		    merge_tracks_array.at(k).push_back(track2);
		  auto itt3 = find(merge_vertices_array.at(k).begin(),merge_vertices_array.at(k).end(),vertex1);
		  auto itt4 = find(merge_vertices_array.at(k).begin(),merge_vertices_array.at(k).end(),vertex2);
		  if (itt3 == merge_vertices_array.at(k).end())
		    merge_vertices_array.at(k).push_back(vertex1);
		  if (itt4 == merge_vertices_array.at(k).end())
		    merge_vertices_array.at(k).push_back(vertex2);
		  flag = 1;
		  break;
		}
	      }
	      if (flag == 0){
		WCTrackSelection merge_tracks;
		merge_tracks.push_back(track1);
		merge_tracks.push_back(track2);
		WCVertexSelection merge_vertices;
		merge_vertices.push_back(vertex1);
		merge_vertices.push_back(vertex2);
		
		merge_tracks_array.push_back(merge_tracks);
		merge_vertices_array.push_back(merge_vertices);
	      }
		// WCTrackSelection& merge_tracks = merge_tracks.at(k);
		// WCVertexSelection& merge_vertices = merge_vertices.at(k);
	    }

	    // std::cout << "abc " << track1->get_range() << " " << track2->get_range() << std::endl;
	  }
	}
      }
    }
  }

  
  for (int i=0;i!=ct_tracks.size();i++){
    delete ct_tracks.at(i);
  }
 
  // std::cout << merge_tracks_array.size() << std::endl;
  // for (int i=0;i!=merge_tracks_array.size();i++){
  //   for (int j=0;j!=merge_tracks_array.at(i).size();j++){
  //     std::cout << merge_tracks_array.at(i).at(j) << " ";
  //   }
  //   std::cout << std::endl;
  // }
  

  // merge and do fine tracking ... 
  //std::cout << merge_tracks_array.size() << std::endl;
  for (int i=0;i!=merge_tracks_array.size();i++){
    // add new track to all vertices
    WCVertexSelection good_vertices1;
    // figure out the good vertices
    for (int j=0;j!=merge_tracks_array.at(i).size();j++){
      WCTrack *track1 = merge_tracks_array.at(i).at(j);
      WCVertexSelection temp_vertices = wct_wcv_map[track1];
      for (int k=0;k!=temp_vertices.size();k++){
	WCVertex* temp_vertex = temp_vertices.at(k);
	auto it1 = find(merge_vertices_array.at(i).begin(),merge_vertices_array.at(i).end(),temp_vertex);
	if (it1 == merge_vertices_array.at(i).end()){
	  auto it2 = find(good_vertices1.begin(),good_vertices1.end(),temp_vertex);
	  if (it2 == good_vertices1.end()){
	    for (int k1=0;k1!=temp_vertex->get_tracks().size();k1++){
	      WCTrack *temp_track = temp_vertex->get_tracks().at(k1);
	      auto it3 = find(merge_tracks_array.at(i).begin(),merge_tracks_array.at(i).end(),temp_track);
	      if (it3 != merge_tracks_array.at(i).end()){
		good_vertices1.push_back(temp_vertex);
		break;
	      }
	    }
	  }
	}
      }
    }
    
    if (good_vertices1.size()<=1){
    }else{
      
      std::vector<int> track_no;
      for (int j=0;j!=good_vertices1.size();j++){
	WCVertex* temp_vertex = good_vertices1.at(j);
       
	for (int k=0;k!=temp_vertex->get_tracks().size();k++){
	  WCTrack *temp_track = temp_vertex->get_tracks().at(k);
	  auto it1 = find(merge_tracks_array.at(i).begin(),merge_tracks_array.at(i).end(),temp_track);
	  if (it1 != merge_tracks_array.at(i).end()){
	    track_no.push_back(j);
	    break;
	  }
	}
      }
    
      // need to pick up two vertices ...
      if (track_no.size()>2){
	// //      std::cout << track_no.size() << std::endl;
	// // form a vector for the cells to be merged ... 
	// MergeSpaceCellSelection to_be_merged;
	// for (int j=0;j!=merge_tracks_array.at(i).size();j++){
	//   for (int k=0;k!=merge_tracks_array.at(i).at(j)->get_centerVP_cells().size();k++){
	//     MergeSpaceCell *mcell = merge_tracks_array.at(i).at(j)->get_centerVP_cells().at(k);
	//     auto it = find(to_be_merged.begin(),to_be_merged.end(),mcell);
	//     if (it == to_be_merged.end())
	//       to_be_merged.push_back(mcell);
	//   }
	// }
	// Point pcenter(0,0,0);
	// int psum = 0;
	// for (int j=0;j!=to_be_merged.size();j++){
	//   pcenter.x += to_be_merged.at(j)->Get_Center().x * to_be_merged.at(j)->Get_all_spacecell().size();
	//   pcenter.y += to_be_merged.at(j)->Get_Center().y * to_be_merged.at(j)->Get_all_spacecell().size();
	//   pcenter.z += to_be_merged.at(j)->Get_Center().z * to_be_merged.at(j)->Get_all_spacecell().size();
	//   psum += to_be_merged.at(j)->Get_all_spacecell().size();
	// }
	// pcenter.x /= psum;
	// pcenter.y /= psum;
	// pcenter.z /= psum;
	
	// float max_dist = 0;
	// Point max_point;
	// for (int j=0;j!=to_be_merged.size();j++){
	//   Point point = to_be_merged.at(j)->Get_Center();
	//   float dist = sqrt(pow(point.x - pcenter.x,2) + 
	// 		    pow(point.y - pcenter.y,2) + 
	// 		    pow(point.z - pcenter.z,2));
	//   if (dist > max_dist){
	//     max_dist = dist;
	//     max_point = point;
	//   }
	// }
	// TVector3 vec1(max_point.x - pcenter.x, max_point.y - pcenter.y, max_point.z - pcenter.z);
	
	int vertex_i1=track_no.at(0);
	int vertex_i2=track_no.at(1);
	
	// TVector3 vec2(good_vertices1.at(vertex_i1)->Center().x-good_vertices1.at(vertex_i2)->Center().x,
	// 	      good_vertices1.at(vertex_i1)->Center().y-good_vertices1.at(vertex_i2)->Center().y,
	// 	      good_vertices1.at(vertex_i1)->Center().z-good_vertices1.at(vertex_i2)->Center().z);
	
	// float max_dis = fabs(vec1.Dot(vec2));

	float max_dis = sqrt(pow(good_vertices1.at(vertex_i1)->Center().x-good_vertices1.at(vertex_i2)->Center().x,2)+
			     pow(good_vertices1.at(vertex_i1)->Center().y-good_vertices1.at(vertex_i2)->Center().y,2)+
			     pow(good_vertices1.at(vertex_i1)->Center().z-good_vertices1.at(vertex_i2)->Center().z,2));
	for (int i1 = 0; i1!=track_no.size();i1++){
	  for (int i2 = 0; i2!=track_no.size();i2++){
	    if (i2==i1) continue;
	    // TVector3 vec3(good_vertices1.at(track_no.at(i1))->Center().x-good_vertices1.at(track_no.at(i2))->Center().x,
	    // 		  good_vertices1.at(track_no.at(i1))->Center().y-good_vertices1.at(track_no.at(i2))->Center().y,
	    // 		  good_vertices1.at(track_no.at(i1))->Center().z-good_vertices1.at(track_no.at(i2))->Center().z);
	    // float dis = fabs(vec1.Dot(vec3));
	    float dis = sqrt(pow(good_vertices1.at(track_no.at(i1))->Center().x-good_vertices1.at(track_no.at(i2))->Center().x,2)+
	    		     pow(good_vertices1.at(track_no.at(i1))->Center().y-good_vertices1.at(track_no.at(i2))->Center().y,2)+
	    		     pow(good_vertices1.at(track_no.at(i1))->Center().z-good_vertices1.at(track_no.at(i2))->Center().z,2));
	    if (dis > max_dis){
	      //std::cout << dis << " " << max_dis << std::endl;
	      max_dis = dis;
	      vertex_i1 = track_no.at(i1);
	      vertex_i2 = track_no.at(i2);
	    }
	  }
	}
	track_no.clear();
	track_no.push_back(vertex_i1);
	track_no.push_back(vertex_i2);
      }
      

    MergeSpaceCellSelection mcells;
    
    
    //create a new track
    MergeSpaceCell* start_cell = good_vertices1.at(track_no.at(0))->get_msc();
    Point start_point = good_vertices1.at(track_no.at(0))->Center();
    MergeSpaceCell* end_cell = good_vertices1.at(track_no.at(1))->get_msc();
    Point end_point = good_vertices1.at(track_no.at(1))->Center();
    WCP2dToy::ToyWalking walking(start_cell,start_point,end_cell,end_point,mcells_map,5000);
    mcells = walking.get_cells();
    // std::cout << start_point.x/units::cm << " " << start_point.y/units::cm << " " << start_point.z/units::cm << " " << std::endl;
    // std::cout << end_point.x/units::cm << " " << end_point.y/units::cm << " " << end_point.z/units::cm << " " << std::endl;
    //if (mcells_map.find(start_cell) == mcells_map.end()){
    // std::cout <<"xin:1 " << mcells_map[start_cell].size() << std::endl;
    //   // }
    //   // if (mcells_map.find(end_cell) == mcells_map.end()){
    // std::cout <<"xin:2 " << mcells_map[end_cell].size() << std::endl;
    //   // }
    // std::cout << "xin: " << mcells.size() << " " << walking.get_counter() << " " << walking.get_global_counter() << std::endl;

    if (mcells.size()==0){
      MergeSpaceCellSet msc_set;
      for (int j=0;j!=merge_tracks_array.at(i).size();j++){
	WCTrack *track = merge_tracks_array.at(i).at(j);
	for (int k=0;k!=track->get_centerVP_cells().size();k++){
	  msc_set.insert(track->get_centerVP_cells().at(k));
	}
      }
      for (auto it = msc_set.begin(); it!= msc_set.end(); it++){
	mcells.push_back(*it);
      }
    }
    //    std::cout << mcells.size() << std::endl;
    
    if (mcells.size()==0) continue;
    WCTrack *track = new WCTrack(mcells);
    tracks.push_back(track);

   
    //    std::cout << good_vertices1.size() << std::endl;
    
    WCVertexSelection good_vertices2;    
    for (int j=0;j!=track_no.size();j++){
      WCVertex* temp_vertex = good_vertices1.at(track_no.at(j));
      good_vertices2.push_back(temp_vertex);
      temp_vertex->Add(track);
      //std::cout << temp_vertex << " " << temp_vertex->get_tracks().size() << std::endl;
      for (int k=0;k!=temp_vertex->get_tracks().size();k++){
	WCTrack *temp_track = temp_vertex->get_tracks().at(k);
	auto it1 = find(merge_tracks_array.at(i).begin(),merge_tracks_array.at(i).end(),temp_track);
	if (it1 != merge_tracks_array.at(i).end()){
	  temp_vertex->set_ky(track,temp_vertex->get_ky(temp_track));
	  temp_vertex->set_kz(track,temp_vertex->get_kz(temp_track));
	  break;
	}
      }
      wcv_wct_map[temp_vertex].push_back(track);
    }
     wct_wcv_map[track] = good_vertices2;

      
      // // add new track to all maps
      // std::vector<int> track_no;
      // for (int j=0;j!=good_vertices1.size();j++){
      // 	WCVertex* temp_vertex = good_vertices1.at(j);
      // 	temp_vertex->Add(track);
      // 	// figure out the original ky and kz ... 
      // 	for (int k=0;k!=temp_vertex->get_tracks().size();k++){
      // 	  WCTrack *temp_track = temp_vertex->get_tracks().at(k);
      // 	  auto it1 = find(merge_tracks_array.at(i).begin(),merge_tracks_array.at(i).end(),temp_track);
      // 	  if (it1 != merge_tracks_array.at(i).end()){
      // 	    temp_vertex->set_ky(track,temp_vertex->get_ky(temp_track));
      // 	    temp_vertex->set_kz(track,temp_vertex->get_kz(temp_track));
      // 	    //std::cout << j<< " abc " << std::endl;
      // 	    track_no.push_back(j);
      // 	    break;
      // 	  }
      // 	}
      // 	wcv_wct_map[temp_vertex].push_back(track);
      // }
      // wct_wcv_map[track] = good_vertices1;
      
      // do fine tracking for the new track
      Point p1 = good_vertices1.at(track_no.at(0))->Center();
      Point p2 = good_vertices1.at(track_no.at(1))->Center();
      int np1 = good_vertices1.at(track_no.at(0))->get_ntracks();
      int np2 = good_vertices1.at(track_no.at(1))->get_ntracks();
      
      //std::cout << good_vertices1.size() << std::endl;
      double ky1, kz1, ky2, kz2;
      ky1 = good_vertices1.at(track_no.at(0))->get_ky(track);
      kz1 = good_vertices1.at(track_no.at(0))->get_kz(track);
      ky2 = good_vertices1.at(track_no.at(1))->get_ky(track);
      kz2 = good_vertices1.at(track_no.at(1))->get_kz(track);
      
      track->fine_tracking(np1,p1,ky1,kz1,np2,p2,ky2,kz2,1);
      good_tracks.push_back(track);

      
      //deal with the other cells not included ...  ??? (to be done)



      // delete old tracks
      for (int j=0;j!=merge_tracks_array.at(i).size();j++){
	WCTrack *temp_track = merge_tracks_array.at(i).at(j);
	temp_track->reset_fine_tracking();
	auto it = find(good_tracks.begin(),good_tracks.end(),temp_track);
	if (it != good_tracks.end())
	  good_tracks.erase(it);
	auto it1 = find(parallel_tracks.begin(),parallel_tracks.end(),temp_track);
	if (it1 != parallel_tracks.end())
	  parallel_tracks.erase(it1);
	bad_tracks.push_back(temp_track);
      }
      
    }
  }

  
  
 

  // delete old vertices if none tracks
  // save good vertices
  good_vertices.clear();
  for (int i=0;i!=good_tracks.size();i++){
    for (int j=0;j!=wct_wcv_map[good_tracks.at(i)].size();j++){
      WCVertex *temp_vertex = wct_wcv_map[good_tracks.at(i)].at(j);
      auto it = find(good_vertices.begin(),good_vertices.end(),temp_vertex);
      if (it == good_vertices.end())
  	good_vertices.push_back(temp_vertex);
    }
  }
  for (int i=0;i!=parallel_tracks.size();i++){
    for (int j=0;j!=wct_wcv_map[parallel_tracks.at(i)].size();j++){
      WCVertex *temp_vertex = wct_wcv_map[parallel_tracks.at(i)].at(j);
      auto it = find(good_vertices.begin(),good_vertices.end(),temp_vertex);
      if (it == good_vertices.end())
  	good_vertices.push_back(temp_vertex);
    }
  }


  bad_vertices.clear();
  // move the rest vertices to bad vertices
  for (int i=0;i!=vertices.size();i++){
    WCVertex *temp_vertex = vertices.at(i);
    auto it = find(good_vertices.begin(),good_vertices.end(),temp_vertex);
    if (it == good_vertices.end())
      bad_vertices.push_back(temp_vertex);
  }



}


void WCP2dToy::ToyTracking::MergeTracks(WCP::MergeSpaceCellMap& mcells_map){
  std::cout << "Merge Tracks! " << std::endl;

  WCTrackSelection all_tracks;
  
  // First merge tracks according to vertices and direction (2 or 3 tracks, and more tracks?)
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    for (int j=0;j!=vertex->get_tracks().size();j++){
      WCTrack *track = vertex->get_tracks().at(j);

      auto it1 = find(good_tracks.begin(),good_tracks.end(),track);
      auto it2 = find(parallel_tracks.begin(),parallel_tracks.end(),track);
      if (it1 == good_tracks.end() && it2 == parallel_tracks.end()) continue;

      auto it = find(all_tracks.begin(),all_tracks.end(),track);
      if (it == all_tracks.end() && track->get_centerVP_cells().size()>0)
	all_tracks.push_back(track);
      // std::cout << i << " " << j << " " << vertex->Center().x/units::cm <<
      // " " << vertex->Center().y/units::cm << " " << vertex->Center().z/units::cm << " " << 
      // vertex->get_ky(track) << " " << vertex->get_kz(track) << std::endl;
    }
  }

  // std::cout << all_tracks.size() << std::endl;
  
  ClusterTrackSelection ct_tracks;
  for (int i=0;i!=all_tracks.size();i++){
    WCTrack *track = all_tracks.at(i);
    MergeSpaceCellSelection& mcells = track->get_centerVP_cells();
    //std::cout << i << " " << mcells.size() << std::endl;
    ClusterTrack *ct = new ClusterTrack(mcells.at(0));
    for (int j=1;j< mcells.size();j++){
      ct->AddMSCell_anyway(mcells.at(j));
    }
    ct_tracks.push_back(ct);
  }
  
  std::map<WCVertex*, float> map_vertex_theta;
  std::map<WCVertex*, float> map_vertex_phi;
  std::map<WCVertex*, WCTrackSelection> map_vertex_tracks;

  //std::cout << "abc" << std::endl;

  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);

    // std::cout << i << " " << vertex->get_tracks().size() << std::endl;
    
    Point center;
    center.x = vertex->get_msc()->Get_Center().x;
    center.y = vertex->get_msc()->Get_Center().y;
    center.z = vertex->get_msc()->Get_Center().z;
    WCTrackSelection used_tracks;

    int ngood_tracks=0;
    for (int j=0;j!=vertex->get_tracks().size();j++){
      WCTrack *track1 = vertex->get_tracks().at(j);
      auto it1 = find(good_tracks.begin(),good_tracks.end(),track1);
      auto it2 = find(parallel_tracks.begin(),parallel_tracks.end(),track1);
      if (it1 != good_tracks.end() || it2!=parallel_tracks.end())
	ngood_tracks ++;
    }
    
    for (int j=0;j!=vertex->get_tracks().size();j++){
      WCTrack *track1 = vertex->get_tracks().at(j);
      auto itt1 = find(all_tracks.begin(),all_tracks.end(),track1);
      if (itt1 == all_tracks.end()) continue;
      
      //      if (track1->get_centerVP_cells().size()==0) continue;

      auto it3 = find(used_tracks.begin(),used_tracks.end(),track1);
      if (it3 != used_tracks.end()) continue;
      used_tracks.push_back(track1);

      auto it1 = find(all_tracks.begin(),all_tracks.end(),track1);
      ct_tracks.at(it1 - all_tracks.begin())->SC_Hough(center);
      double theta1 = ct_tracks.at(it1 - all_tracks.begin())->Get_Theta();
      double phi1 = ct_tracks.at(it1 - all_tracks.begin())->Get_Phi();

      WCVertexSelection t1_vertices = wct_wcv_map[track1];
      WCVertex *new_vertex1 = 0;
      for (int i1 = 0;i1!=t1_vertices.size();i1++){
	if (t1_vertices.at(i1) == vertex) continue;
	if (find(t1_vertices.at(i1)->get_tracks().begin(),
		 t1_vertices.at(i1)->get_tracks().end(),track1)==t1_vertices.at(i1)->get_tracks().end()) continue;
	new_vertex1 = t1_vertices.at(i1);
	break;
      }


      if (track1->get_range() < 7*units::cm){
	
	
	TVector3 temp_abc1(new_vertex1->Center().x - vertex->Center().x,
			   new_vertex1->Center().y - vertex->Center().y,
			   new_vertex1->Center().z - vertex->Center().z);
	theta1 = temp_abc1.Theta();
	phi1 = temp_abc1.Phi();
      }

      for (int k=0;k!=vertex->get_tracks().size();k++){
	WCTrack *track2 = vertex->get_tracks().at(k);

	//	std::cout << track1->get_range() << " " << track2->get_range() <<std::endl;
	
	auto itt2 = find(all_tracks.begin(),all_tracks.end(),track2);
	if (itt2 == all_tracks.end()) continue;
      //if (track2->get_centerVP_cells().size()==0) continue;
	auto it4 = find(used_tracks.begin(),used_tracks.end(),track2);
	if (it4 != used_tracks.end()) continue;
	
	 auto it2 = find(all_tracks.begin(),all_tracks.end(),track2);
	 ct_tracks.at(it2 - all_tracks.begin())->SC_Hough(center);
	 double theta2 = ct_tracks.at(it2 - all_tracks.begin())->Get_Theta();
	 double phi2 = ct_tracks.at(it2 - all_tracks.begin())->Get_Phi();
	 
	 WCVertexSelection t2_vertices = wct_wcv_map[track2];
	 WCVertex *new_vertex2 = 0;
	 for (int i1 = 0;i1!=t2_vertices.size();i1++){
	   if (t2_vertices.at(i1) == vertex) continue;
	   if (find(t2_vertices.at(i1)->get_tracks().begin(),
		    t2_vertices.at(i1)->get_tracks().end(),track2)==t2_vertices.at(i1)->get_tracks().end()) continue;
	   new_vertex2 = t2_vertices.at(i1);
	   break;
	 }


	 if (track2->get_range() < 7*units::cm){
	   
	   TVector3 temp_abc2(new_vertex2->Center().x - vertex->Center().x,
			      new_vertex2->Center().y - vertex->Center().y,
			      new_vertex2->Center().z - vertex->Center().z);
	   theta2 = temp_abc2.Theta();
	   phi2 = temp_abc2.Phi();
	 }	 
	 
	 TVector3 abc1(sin(theta1)*cos(phi1),sin(theta1)*sin(phi1),cos(theta1));
	 TVector3 abc2(sin(theta2)*cos(phi2),sin(theta2)*sin(phi2),cos(theta2));
	 double theta_abc = abc1.Angle(abc2);
	 //std::cout << theta_abc/3.1415926*180 << std::endl;
	 
	 TVector3 abc3(new_vertex1->Center().x - vertex->Center().x,
		       new_vertex1->Center().y - vertex->Center().y,
		       new_vertex1->Center().z - vertex->Center().z);
	 TVector3 abc4(new_vertex2->Center().x - vertex->Center().x,
		       new_vertex2->Center().y - vertex->Center().y,
		       new_vertex2->Center().z - vertex->Center().z);
	 double theta_abc1 = abc3.Angle(abc4);

	 
	 //	 std::cout << ngood_tracks << " " << theta_abc/3.1415926*180. << " " << theta_abc1/3.1415926*180. << " " << fabs(theta1+theta2-3.1415926)/3.1415926*180. << " " << (fabs(phi1-phi2)-3.1415926)/3.1415926*180. << " " << track1->get_range() << " " << track2->get_range() << " " << vertex->Center().x/units::cm << " " <<  std::endl;
	 
	 int flag_qx = 0;
	 if (ngood_tracks==2){
	   //if (sqrt(pow(theta1+theta2-3.1415926,2) + pow(fabs(phi1-phi2)-3.1415926,2)) < 28./180.*3.1415926){
	   if (fabs(theta_abc/3.1415926*180-180.)<25. || fabs(theta_abc1/3.1415926*180.-180.)<25.){
	     flag_qx = 1;
	   }else{
	     if (track1->get_range() > 7*units::cm){
	       ct_tracks.at(it1 - all_tracks.begin())->SC_Hough(center,15*units::cm);
	       theta1 = ct_tracks.at(it1 - all_tracks.begin())->Get_Theta();
	       phi1 = ct_tracks.at(it1 - all_tracks.begin())->Get_Phi();
	     }
	     if (track2->get_range()>7*units::cm){
	       ct_tracks.at(it2 - all_tracks.begin())->SC_Hough(center,15*units::cm);
	       theta2 = ct_tracks.at(it2 - all_tracks.begin())->Get_Theta();
	       phi2 = ct_tracks.at(it2 - all_tracks.begin())->Get_Phi();
	     }
	     abc1.SetXYZ(sin(theta1)*cos(phi1),sin(theta1)*sin(phi1),cos(theta1));
	     abc2.SetXYZ(sin(theta2)*cos(phi2),sin(theta2)*sin(phi2),cos(theta2));
	     theta_abc = abc1.Angle(abc2);

	     // if (sqrt(pow(theta1+theta2-3.1415926,2) + pow(fabs(phi1-phi2)-3.1415926,2)) < 28./180.*3.1415926)
	     if (fabs(theta_abc/3.1415926*180.-180.)<25.)
	       flag_qx = 1;
	     
	     

	   }
	 }else if (ngood_tracks>2){
	   //  if (fabs(theta1+theta2-3.1415926) < 10./180.*3.1415926 && 
	   //  fabs(fabs(phi1-phi2)-3.1415926)<10./180.*3.1415926 &&
	   //  sqrt(pow(theta1+theta2-3.1415926,2) + pow(fabs(phi1-phi2)-3.1415926,2)) < 14./180.*3.1415926)
	   if (fabs(theta_abc/3.1415926*180.-180)<15.|| fabs(theta_abc1/3.1415926*180.-180.)<15.)
	     flag_qx = 1;
	 }
	 //std::cout << ngood_tracks << " " << fabs(theta1+theta2-3.1415926)/3.1415926*180. << " " << (fabs(phi1-phi2)-3.1415926)/3.1415926*180. << " " << track1->get_range() << " " << track2->get_range() << " " << vertex->Center().x/units::cm << " " << flag_qx << std::endl;

	 if (flag_qx == 1){

	   if (map_vertex_theta.find(vertex)==map_vertex_theta.end()){
	     map_vertex_theta[vertex] = theta1+theta2-3.1415926;
	     map_vertex_phi[vertex] = fabs(phi1-phi2)-3.1415926;
	     WCTrackSelection temp_tracks;
	     temp_tracks.push_back(track1);
	     temp_tracks.push_back(track2);
	     map_vertex_tracks[vertex] = temp_tracks;
	   }else{
	     WCTrack *track3 = map_vertex_tracks[vertex].at(0);
	     WCTrack *track4 = map_vertex_tracks[vertex].at(1);
	     
	     // std::cout << (theta1+theta2-3.1415926)/3.1415926*180. << " " << (fabs(phi1-phi2)-3.1415926)/3.1415926*180.
	     // 	       << " " << map_vertex_theta[vertex]/3.1415926*180. << " " << map_vertex_phi[vertex]/3.1415926*180.
	     // 	       << " " << track1->get_range() << " " << track2->get_range() 
	     // 	       << " " << track3->get_range() << " " << track4->get_range() << std::endl;
	     
	     // save the smallest??? 
	     if (pow(theta1+theta2-3.1415926,2) + pow(fabs(phi1-phi2)-3.1415926,2) < 
	         pow(map_vertex_theta[vertex],2) + pow(map_vertex_phi[vertex],2)){
	       map_vertex_theta[vertex] = theta1+theta2-3.1415926;
	       map_vertex_phi[vertex] = fabs(phi1-phi2)-3.1415926;
	       map_vertex_tracks[vertex].clear();
	       map_vertex_tracks[vertex].push_back(track1);
	       map_vertex_tracks[vertex].push_back(track2);
	     }
	     
	   }
	 }
	 
	 // std::cout << i << " " << center.x/units::cm << " " << center.y/units::cm << " " << center.z/units::cm 
	 // 	   << " " << theta1+theta2-3.1415926 << " " << fabs(phi1-phi2)-3.1415926 << std::endl;

      }
    }
  }
  
  //std::cout << "abc:" << std::endl;

  for (int i=0;i!=ct_tracks.size();i++){
    delete ct_tracks.at(i);
  }

  // for (auto it = map_vertex_theta.begin(); it!= map_vertex_theta.end(); it++){
  //     WCVertex *vertex = it->first;
  //     WCTrack *track1 = map_vertex_tracks[vertex].at(0);
  //     WCTrack *track2 = map_vertex_tracks[vertex].at(1);
  //     std::cout << vertex << " " << track1 << " " << track2 << " " << wcv_wct_map[vertex].size() << std::endl;
  // }
  

  std::vector<WCVertexSelection> merge_vertices;
  std::vector<WCTrackSelection> merge_tracks;

  WCVertexSelection used_vertices;
  int prev_vertices = 0;
  int flag = 1;

  while(used_vertices.size() < map_vertex_theta.size()){
    prev_vertices = used_vertices.size();


    for (auto it = map_vertex_theta.begin(); it!= map_vertex_theta.end(); it++){
      WCVertex *vertex = it->first;
      WCTrack *track1 = map_vertex_tracks[vertex].at(0);
      WCTrack *track2 = map_vertex_tracks[vertex].at(1);
      
      auto it1 = find(used_vertices.begin(),used_vertices.end(),vertex);
      if (it1 != used_vertices.end()) continue;
      
      if (0<1){
	if (merge_vertices.size()==0 || flag == 1){
	  WCVertexSelection temp_vertices;
	  temp_vertices.push_back(vertex);
	  WCTrackSelection temp_tracks;
	  temp_tracks.push_back(track1);
	  temp_tracks.push_back(track2);
	  merge_vertices.push_back(temp_vertices);
	  merge_tracks.push_back(temp_tracks);
	  used_vertices.push_back(vertex);
	  flag = 0;
	}else{
	  // loop through exiting stuff
	  for (int i=0;i!=merge_vertices.size();i++){
	    WCVertexSelection& temp_vertices = merge_vertices.at(i);
	    WCTrackSelection& temp_tracks = merge_tracks.at(i);
	    auto it1 = find(temp_tracks.begin(),temp_tracks.end(),track1);
	    auto it2 = find(temp_tracks.begin(),temp_tracks.end(),track2);
	    if (it1 != temp_tracks.end() || it2 != temp_tracks.end()){
	      temp_vertices.push_back(vertex);
	      if (it1 == temp_tracks.end())
		temp_tracks.push_back(track1);
	      if (it2 == temp_tracks.end())
		temp_tracks.push_back(track2);
	      used_vertices.push_back(vertex);
	      break;
	    }
	    
	  }
	}
      }else{
	used_vertices.push_back(vertex);
      }
      
     
      // std::cout << it->first->Center().x/units::cm << " " << it->first->Center().y/units::cm 
      // 	      << " " << it->first->Center().z/units::cm << " " << 
      //   it->second/3.1415926*180. << " " << map_vertex_phi[it->first]/3.1415926*180. << std::endl;
    }
    
    if (prev_vertices == used_vertices.size()){
      flag = 1;
    }

  }
  
  // std::cout << merge_tracks.size() << " " << merge_tracks.at(0).size() << std::endl;

  
  // Now create a new track and remove the old tracks and vertices ... 
  for (int i=0;i!=merge_tracks.size();i++){
     // add new track to all vertices
    WCVertexSelection good_vertices1;
    // figure out the good vertices
    for (int j=0;j!=merge_tracks.at(i).size();j++){
      WCTrack *track1 = merge_tracks.at(i).at(j);
      WCVertexSelection temp_vertices = wct_wcv_map[track1];
      for (int k=0;k!=temp_vertices.size();k++){
	WCVertex* temp_vertex = temp_vertices.at(k);
	auto it1 = find(merge_vertices.at(i).begin(),merge_vertices.at(i).end(),temp_vertex);
	if (it1 == merge_vertices.at(i).end()){
	  auto it2 = find(good_vertices1.begin(),good_vertices1.end(),temp_vertex);
	  if (it2 == good_vertices1.end()){
	    for (int k1=0;k1!=temp_vertex->get_tracks().size();k1++){
	      WCTrack *temp_track = temp_vertex->get_tracks().at(k1);
	      auto it3 = find(merge_tracks.at(i).begin(),merge_tracks.at(i).end(),temp_track);
	      if (it3 != merge_tracks.at(i).end()){
		good_vertices1.push_back(temp_vertex);
		break;
	      }
	    }
	  }
	}
      }
    }
    
    
    std::vector<int> track_no;
    for (int j=0;j!=good_vertices1.size();j++){
      WCVertex* temp_vertex = good_vertices1.at(j);
      for (int k=0;k!=temp_vertex->get_tracks().size();k++){
	WCTrack *temp_track = temp_vertex->get_tracks().at(k);
	auto it1 = find(merge_tracks.at(i).begin(),merge_tracks.at(i).end(),temp_track);
	if (it1 != merge_tracks.at(i).end()){
	  track_no.push_back(j);
	  break;
	}
      }
    }
    
    // need to pick up two vertices ...
    if (track_no.size()>2){
      int vertex_i1=track_no.at(0);
      int vertex_i2=track_no.at(1);
      float max_dis = sqrt(pow(good_vertices1.at(vertex_i1)->Center().x-good_vertices1.at(vertex_i2)->Center().x,2)+
			   pow(good_vertices1.at(vertex_i1)->Center().y-good_vertices1.at(vertex_i2)->Center().y,2)+
			   pow(good_vertices1.at(vertex_i1)->Center().z-good_vertices1.at(vertex_i2)->Center().z,2));
      for (int i1 = 0; i1!=track_no.size();i1++){
	for (int i2 = 0; i2!=track_no.size();i2++){
	  if (i2==i1) continue;
	  float dis = sqrt(pow(good_vertices1.at(track_no.at(i1))->Center().x-good_vertices1.at(track_no.at(i2))->Center().x,2)+
			   pow(good_vertices1.at(track_no.at(i1))->Center().y-good_vertices1.at(track_no.at(i2))->Center().y,2)+
			   pow(good_vertices1.at(track_no.at(i1))->Center().z-good_vertices1.at(track_no.at(i2))->Center().z,2));
	  if (dis > max_dis){
	    max_dis = dis;
	    vertex_i1 = track_no.at(i1);
	    vertex_i2 = track_no.at(i2);
	  }
	}
      }
      track_no.clear();
      track_no.push_back(vertex_i1);
      track_no.push_back(vertex_i2);
    }
    
    

    //    std::cout << "abc: "<<  good_vertices1.size() << " " << track_no.size() << " " << merge_tracks.at(i).size() << " " << merge_vertices.at(i).size() << std::endl;
    // for (int i1 = 0;i1!=merge_vertices.at(i).size();i1++){
    //   std::cout << merge_vertices.at(i).at(i1) << std::endl;
    // }
    // for (int i1 =0;i1!=merge_tracks.at(i).size();i1++){
    //   std::cout << merge_tracks.at(i).at(i1) << " " << wct_wcv_map[merge_tracks.at(i).at(i1)].size() << " " << wct_wcv_map[merge_tracks.at(i).at(i1)].at(0) << " " << wct_wcv_map[merge_tracks.at(i).at(i1)].at(1) << std::endl;
    // }
    
    // creat a new track
    
    //
    MergeSpaceCellSelection mcells;
    
    MergeSpaceCell* start_cell = good_vertices1.at(track_no.at(0))->get_msc();
    Point start_point = good_vertices1.at(track_no.at(0))->Center();
    MergeSpaceCell* end_cell = good_vertices1.at(track_no.at(1))->get_msc();
    Point end_point = good_vertices1.at(track_no.at(1))->Center();
    WCP2dToy::ToyWalking walking(start_cell,start_point,end_cell,end_point,mcells_map,5000);
    mcells = walking.get_cells();
    // std::cout << "xin: " << mcells.size() << std::endl;
  
    //std::cout << mcells.size() << std::endl;
    if (mcells.size()==0){
      MergeSpaceCellSet msc_set;
      for (int j=0;j!=merge_tracks.at(i).size();j++){
	WCTrack *track = merge_tracks.at(i).at(j);
	for (int k=0;k!=track->get_centerVP_cells().size();k++){
	  msc_set.insert(track->get_centerVP_cells().at(k));
	}
      }
      for (auto it = msc_set.begin(); it!= msc_set.end(); it++){
	mcells.push_back(*it);
      }
    }

    WCTrack *track = new WCTrack(mcells);
    tracks.push_back(track);
    //std::cout << i << " " << merge_tracks.size() << " " << msc_set.size() << " " << mcells.size() << std::endl;


    WCVertexSelection good_vertices2;    
    for (int j=0;j!=track_no.size();j++){
      WCVertex* temp_vertex = good_vertices1.at(track_no.at(j));
      good_vertices2.push_back(temp_vertex);
      temp_vertex->Add(track);
      //std::cout << temp_vertex << " " << temp_vertex->get_tracks().size() << std::endl;
      for (int k=0;k!=temp_vertex->get_tracks().size();k++){
	WCTrack *temp_track = temp_vertex->get_tracks().at(k);
	auto it1 = find(merge_tracks.at(i).begin(),merge_tracks.at(i).end(),temp_track);
	if (it1 != merge_tracks.at(i).end()){
	  temp_vertex->set_ky(track,temp_vertex->get_ky(temp_track));
	  temp_vertex->set_kz(track,temp_vertex->get_kz(temp_track));
	  break;
	}
      }
      wcv_wct_map[temp_vertex].push_back(track);
    }
     wct_wcv_map[track] = good_vertices2;
     
    
    // //std::cout << good_vertices1.size() << std::endl;
    // // add new track to all maps
    
    // for (int j=0;j!=good_vertices1.size();j++){
    //   WCVertex* temp_vertex = good_vertices1.at(j);
    //   temp_vertex->Add(track);
    //   // figure out the original ky and kz ... 
    //   for (int k=0;k!=temp_vertex->get_tracks().size();k++){
    // 	WCTrack *temp_track = temp_vertex->get_tracks().at(k);
    // 	auto it1 = find(merge_tracks.at(i).begin(),merge_tracks.at(i).end(),temp_track);
    // 	if (it1 != merge_tracks.at(i).end()){
    // 	  temp_vertex->set_ky(track,temp_vertex->get_ky(temp_track));
    // 	  temp_vertex->set_kz(track,temp_vertex->get_kz(temp_track));
    // 	  //std::cout << j<< " abc " << std::endl;
    // 	  break;
    // 	}
    //   }
    //   wcv_wct_map[temp_vertex].push_back(track);
    // }
    // wct_wcv_map[track] = good_vertices1;

    

    // do fine tracking for the new track
    Point p1 = good_vertices1.at(track_no.at(0))->Center();
    Point p2 = good_vertices1.at(track_no.at(1))->Center();
    int np1 = good_vertices1.at(track_no.at(0))->get_ntracks();
    int np2 = good_vertices1.at(track_no.at(1))->get_ntracks();
    
    //std::cout << good_vertices1.size() << std::endl;
    double ky1, kz1, ky2, kz2;
    ky1 = good_vertices1.at(track_no.at(0))->get_ky(track);
    kz1 = good_vertices1.at(track_no.at(0))->get_kz(track);
    ky2 = good_vertices1.at(track_no.at(1))->get_ky(track);
    kz2 = good_vertices1.at(track_no.at(1))->get_kz(track);
    

    track->fine_tracking(np1,p1,ky1,kz1,np2,p2,ky2,kz2,1);
    good_tracks.push_back(track);
    // fine_tracking(1);
  }

  //  std::cout <<"abc:" << std::endl;

  for (int i=0;i!=merge_tracks.size();i++){
    // delete old tracks
    for (int j=0;j!=merge_tracks.at(i).size();j++){
      WCTrack *temp_track = merge_tracks.at(i).at(j);
      temp_track->reset_fine_tracking();
      auto it = find(good_tracks.begin(),good_tracks.end(),temp_track);
      if (it != good_tracks.end())
     	good_tracks.erase(it);
      auto it1 = find(parallel_tracks.begin(),parallel_tracks.end(),temp_track);
      if (it1 != parallel_tracks.end())
     	parallel_tracks.erase(it1);
      bad_tracks.push_back(temp_track);
    }
  }
  // delete old vertices if none tracks
  // save good vertices
  good_vertices.clear();
  for (int i=0;i!=good_tracks.size();i++){
    for (int j=0;j!=wct_wcv_map[good_tracks.at(i)].size();j++){
      WCVertex *temp_vertex = wct_wcv_map[good_tracks.at(i)].at(j);
      auto it = find(good_vertices.begin(),good_vertices.end(),temp_vertex);
      if (it == good_vertices.end())
	good_vertices.push_back(temp_vertex);
    }
  }
  for (int i=0;i!=parallel_tracks.size();i++){
    for (int j=0;j!=wct_wcv_map[parallel_tracks.at(i)].size();j++){
      WCVertex *temp_vertex = wct_wcv_map[parallel_tracks.at(i)].at(j);
      auto it = find(good_vertices.begin(),good_vertices.end(),temp_vertex);
      if (it == good_vertices.end())
	good_vertices.push_back(temp_vertex);
    }
  }


  bad_vertices.clear();
  // move the rest vertices to bad vertices
  for (int i=0;i!=vertices.size();i++){
    WCVertex *temp_vertex = vertices.at(i);
    auto it = find(good_vertices.begin(),good_vertices.end(),temp_vertex);
    if (it == good_vertices.end())
      bad_vertices.push_back(temp_vertex);
  }
  
  // need another special merging code ...
  

}

void WCP2dToy::ToyTracking::fill_maps(){
  if (tracks_angle_map.size()==0){
    for (int i = 0;i!=good_tracks.size();i++){
      WCTrack *track = good_tracks.at(i);
      tracks_angle_map[track] = track->get_direction();
      tracks_pos_map[track] = track->get_position();
    }
  }
}

bool WCP2dToy::ToyTracking::IsContained(){
  MergeSpaceCellSelection all_cells;
  for (int i=0;i!=tracks.size();i++){
    for (int j=0;j!=tracks.at(i)->get_all_cells().size();j++){
      auto it = find(all_cells.begin(), all_cells.end(), tracks.at(i)->get_all_cells().at(j));
      if (it == all_cells.end()) all_cells.push_back(tracks.at(i)->get_all_cells().at(j));
    }
  }
  
  for (int i=0;i!=all_cells.size();i++){
    if (all_cells.at(i)->Get_Center().x >= 1*units::cm && all_cells.at(i)->Get_Center().x <= (256 - 1) * units::cm
     && fabs(all_cells.at(i)->Get_Center().y) <= 233 * units::cm/2. - 10*units::cm
     && all_cells.at(i)->Get_Center().z >= 10 * units::cm && all_cells.at(i)->Get_Center().z <= 10.3692 *100 * units::cm - 10 * units::cm) 
      return true;
  }
  return false;
}

WCP2dToy::ToyTracking::ToyTracking(WCP2dToy::ToyCrawler& toycrawler, int tracking_type){

  
  CreateVertices(toycrawler);

  //std::cout << vertices.size() << " " << tracks.size() << std::endl;

  RemoveSame(); // get rid of duplicated ones
  
  //std::cout << "abc: " << std::endl;

  MergeVertices();   // merge things that are common ...
  OrganizeTracks();  // improve the end points nearby and break things
  Associate();    //associate the rest

  //std::cout << "abc: " << std::endl;

  MergeVertices(0); //merge things, but require the center must be at the end of common track
  //CheckVertices(toycrawler);
  OrganizeTracks(); //improve the end points nearby and break things
  Associate(); //associate the rest
  
  //std::cout << "abc: " << std::endl;

  MergeVertices(2);  // merge some vertices with distance cut only 
  CheckVertices(toycrawler);
  Crawl();
  OrganizeTracks(); //associate the rest
  Associate(); //associate the rest

  RemoveSameTrack();

  //std::cout << "abc: " << std::endl;
    
  MergeVertices(2);  // merge some vertices together, allow single track
  // std::cout << "abc: " << std::endl;
  CheckVertices(toycrawler);
  // std::cout << "abc: " << std::endl;
  BreakTracks();    //improve the end points and break things
  // std::cout << "abc: " << std::endl;
  OrganizeTracks(); //associate the rest
  //std::cout << "abc: " << std::endl;
  
  // MergeVertices(2);  // merge some vertices together, allow single track
  // CheckVertices(toycrawler);
  // OrganizeTracks(); //associate the rest

  RemoveSameTrack();

  Associate();  //associate the rest .. 
  CleanUpVertex(); //do some final clean up ...
  
  //std::cout << "abc: " << std::endl;

  if (tracking_type == 0 ){
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
      
      //  std::cout << i << " Vertex " << vertex->get_ntracks() << " " << success << " " << vertex->Center().x/units::cm << " " << vertex->Center().y/units::cm << " " << vertex->Center().z/units::cm << " " << vertex->get_msc()->Get_Center().x/units::cm << std::endl;
      // for (int j=0;j!=vertex->get_ntracks();j++){
      //   std::cout << vertex->get_tracks().at(j)->get_end_scells().at(0)->Get_Center().x/units::cm << " " 
      // 		<< vertex->get_tracks().at(j)->get_end_scells().at(1)->Get_Center().x/units::cm << " " 
      // 		<< vertex->get_tracks().at(j)->get_all_cells().front()->Get_Center().x/units::cm << " " 
      // 		<< vertex->get_tracks().at(j)->get_all_cells().back()->Get_Center().x/units::cm << " " 
      // 		<< std::endl;
      // }
    }
  }else if(tracking_type==1){
     for (int i=0;i!=vertices.size();i++){
      WCVertex *vertex = vertices.at(i);
      
      bool success = vertex->FindVertex();
      if (!success){
     	vertex->reset_center();
      }else{
	ExamineVertex(vertex,toycrawler);
      }
     }
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

  

  std::map<WCVertex*,int> vertex_track_map;
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    vertex_track_map[vertex] = vertex->get_ntracks();
  }

  //deal with wiggle tracks
  //std::cout << "Deal with Wiggle Tracks " << std::endl;
  deal_wiggle_tracks();
  CheckVertices(toycrawler); // redefine all the tracks ... 
  update_maps();
  
 
  //std::cout << wiggle_vertices.size() << std::endl;

  //  second round vertex fitting ... 
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    bool success = vertex->get_fit_success();
        
    std::cout << i << " Vertex " << vertex->get_ntracks() << " " << success << " " << vertex->Center().x/units::cm << " " << vertex->Center().y/units::cm << " " << vertex->Center().z/units::cm << " " << vertex->get_msc()->Get_Center().x/units::cm << std::endl;
    
    if (vertex_track_map[vertex] != vertex->get_ntracks()){
      vertex->FindVertex(1);
    //std::cout << i << " Vertex " << vertex << " " << vertex->get_msc() << " " << vertex->get_ntracks() << " " << success << " " << vertex->Center().x/units::cm << " " << vertex->Center().y/units::cm << " " << vertex->Center().z/units::cm << " " << vertex->get_msc()->Get_Center().x/units::cm << std::endl;
    }
    
    
    // for (int j=0;j!=vertex->get_ntracks();j++){
    //   std::cout << j << " " << vertex->get_tracks().at(j)->get_end_scells().at(0)->Get_Center().x/units::cm 
    // 		<< " " << vertex->get_tracks().at(j)->get_end_scells().at(1)->Get_Center().x/units::cm  << std::endl;
    // }
  }
  CheckVertices(toycrawler); // redefine all the tracks ... 

  

  if (vertices.size() > 0 ) {
    
    

    if (tracking_type == 0){
      //Now do fine tracking for existing tracks
      std::cout << "FineTracking " << std::endl; 
      update_maps();
      fine_tracking();
      cleanup_bad_tracks();
      update_maps1();

      

      // Now need to figure out how to judge whether this is a shower or track ... 
      // just from the number of tracks and the connectivities?
      bool shower_flag = IsThisShower(toycrawler);
      //shower_flag = false;
      
      // return;
      std::cout << "Shower? " << shower_flag << " Vertices " << vertices.size() << std::endl;
      // separate various issues ... 
      if (!shower_flag){
	//not a shower
      	std::cout << "Grown single track " << std::endl;
	if (grow_track_fill_gap(toycrawler)){
	  std::cout << "FineTracking Again" << std::endl; 
	  update_maps(1);
	  fine_tracking();
	}
	form_parallel_tiny_tracks(toycrawler);

	update_maps(1);
	fine_tracking(1);

	// create new tracks ... 
	create_new_tracks_from_leftover(toycrawler);
      }else{
      	//is a shower  
      	// std::cout << "Grown single track " << std::endl;
      	if (grow_track_fill_gap(toycrawler)){
      	  std::cout << "FineTracking Again" << std::endl; 
      	  update_maps(1);
      	  fine_tracking();
      	}
      	//std::cout << "Test Shower only" << std::endl; 
	
      	// Judge vertex with multiple tracks ...
      	if (track_shower_reco(toycrawler)){
      	  std::cout << "Track + Shower " << std::endl; 
      	  Cleanup_showers();
      	  std::cout << "Parallel Tracking " << std::endl; 
      	  // do the rest of fine tracking ... 
      	  form_parallel_tiny_tracks(toycrawler);
      	  update_maps(1);
      	  fine_tracking(1);
      	}else{
      	  std::cout << "Shower only" << std::endl; 
      	  // Judge vertex for single shower ...
      	  single_shower_reco(toycrawler);
      	  Cleanup_showers();
      	}
      }
    }else if (tracking_type == 1){
      //Now do fine tracking for existing tracks
      std::cout << "FineTracking " << std::endl; 
      update_maps();
      std::cout << "FineTracking 1" << std::endl; 
      fine_tracking();
      std::cout << "FineTracking 2" << std::endl; 
      update_maps1();      
      std::cout << "Grow tracks" << std::endl;
      //cosmic ray tuned
      if (grow_track_fill_gap(toycrawler)){
  	CheckVertices(toycrawler); 
       	std::cout << "FineTracking Again" << std::endl; 
       	update_maps(1);
       	fine_tracking();
      }
      
      //std::cout << "Parallel Tracking " << std::endl; 
      // form_parallel_tiny_tracks(toycrawler);
      // update_map(1);
      // fine_tracking(1);
    }
  }else{
    // if no vertices are found ... 
    //  do another parallel track finding without vertex... 
    //std::cout << "ABC: " << std::endl;
    form_parallel_tiny_tracks_wovertex(toycrawler);
    update_maps(1);
    fine_tracking(1);
    //std::cout << "ABC: " << std::endl;
  }



  if (tracking_type==1){
    std::cout << "Special tracking for cosmic " << std::endl; 
    // deal with the case where there are all black 
    if (good_tracks.size() == 0 ){
      std::cout << "Cosmic all " << std::endl;
      cosmic_finder_all(toycrawler);
      std::cout << "FineTracking Again All" << std::endl; 
      if (vertices.size() > 0){
	update_maps();
	fine_tracking();
      }
    }else{
      // deal with the case where there are some black cases
      std::cout << "Cosmic part " << std::endl;
      cosmic_finder_part(toycrawler);
      std::cout << "FineTracking Again Part" << std::endl; 
      if (vertices.size() > 0 ){
	update_maps();
	fine_tracking();
      }
    }
    
  }





}

void WCP2dToy::ToyTracking::cosmic_finder_part(WCP2dToy::ToyCrawler& toycrawler){
  MergeSpaceCellMap& mcells_map = toycrawler.Get_mcells_map();
  MergeSpaceCellSelection all_cells;  
  
  for (auto it =toycrawler.Get_mcells_map().begin();it!=toycrawler.Get_mcells_map().end();it++){
    MergeSpaceCell *mcell = it->first;
    int flag = 0;
    for (int i=0;i!=good_tracks.size();i++){
      WCTrack *track = good_tracks.at(i);
      auto it = find(track->get_centerVP_cells().begin(),
		     track->get_centerVP_cells().end(),
		     mcell);
      if (it != track->get_centerVP_cells().end()){
	flag = 1;
	break;
      }
      //save everything not inside the good tracks
    }
    if (flag == 0){
      all_cells.push_back(mcell);
    }
  }

  std::vector<MergeSpaceCellSelection> cluster_msc;

  MergeSpaceCellSelection used_cell;
  MergeSpaceCellSelection must_cell = all_cells;
  int start_num = 0;
  while(used_cell.size() != must_cell.size() && start_num < must_cell.size()){
    MergeSpaceCell* start_cell=0;
    for (int i=start_num; i<must_cell.size();i++){
      start_cell = must_cell.at(i);
      auto it = find(used_cell.begin(),used_cell.end(),start_cell);
      if (it == used_cell.end()){
	start_num = i+1;
	break;
      }
    }
    
    ToyNWalking toynwalking(start_cell,mcells_map,used_cell,must_cell);
    MergeSpaceCellSelection cells = toynwalking.get_cells();
    
    used_cell.insert(used_cell.end(),cells.begin(),cells.end());

    int sum = 0;
    for (int i=0;i!=cells.size();i++){
      sum += cells.at(i)->Get_all_spacecell().size();
    }
    if (sum >= 150 && cells.size() >3 || sum > 500 && cells.size() > 2 ||
	sum >=1000 && cells.size() > 1 || sum >=2000 && cells.size() ==1)
      cluster_msc.push_back(cells);

    // std::cout << used_cell.size() << " " << must_cell.size() << " " 
    // 	      << start_num << " " << cluster_msc.size() << std::endl; 
  }
  

  //std::cout << cluster_msc.size() << std::endl;
  // for (int i=0;i!=cluster_msc.size();i++){
  //     MergeSpaceCellSelection& mscs_1 = cluster_msc.at(i);
  //     std::cout << mscs_1.size() << std::endl;
  // }

  for (int i=0;i!=cluster_msc.size();i++){
    Point p(0,0,0);
    int sum = 0;
    all_cells.clear();
    all_cells = cluster_msc.at(i);
    for (int j = 0; j!= all_cells.size();j++){
      MergeSpaceCell *mcell = all_cells.at(j);
      p.x += mcell->Get_Center().x * mcell->Get_all_spacecell().size(); 
      p.y += mcell->Get_Center().y * mcell->Get_all_spacecell().size(); 
      p.z += mcell->Get_Center().z * mcell->Get_all_spacecell().size(); 
      sum += mcell->Get_all_spacecell().size(); 
    }
    p.x/=sum;
    p.y/=sum;
    p.z/=sum;

    ClusterTrack ct(all_cells.front());
    for (int i=1;i< all_cells.size();i++){
      ct.AddMSCell_anyway(all_cells.at(i));
    }
    ct.SC_Hough(p);
    float theta = ct.Get_Theta();
    float phi = ct.Get_Phi();
    
    Point p1(p.x + sin(theta)*cos(phi), p.y + sin(theta)*sin(phi), p.z + cos(theta));
    Line l1(p,p1);
    TVector3 dir = l1.vec();
    
    
    float max_dis = 0;
    float min_dis = 0;
    MergeSpaceCell *max_cell = all_cells.at(0);
    MergeSpaceCell *min_cell = all_cells.at(0);
    for (int i=0;i!=all_cells.size();i++){
      MergeSpaceCell *mcell = all_cells.at(i);
      TVector3 dir1(mcell->Get_Center().x-p.x,mcell->Get_Center().y-p.y,mcell->Get_Center().z-p.z);
      float dis = dir.Dot(dir1);
      float dis1 = l1.closest_dis(mcell->Get_Center());
      if (dis >0 ) dis = dis - dis1;
      if (dis <0 ) dis = dis + dis1;
      
      if (dis > max_dis){
	max_dis = dis;
	max_cell = mcell;
      }
      if (dis < min_dis){
	min_dis = dis;
	min_cell = mcell;
      }
    }
    
    Point max_point = max_cell->Get_Center();
    Point min_point = min_cell->Get_Center();
    for (int i=0;i!=max_cell->Get_all_spacecell().size();i++){
      SpaceCell *cell = max_cell->Get_all_spacecell().at(i);
      TVector3 dir1(cell->x()-p.x,cell->y()-p.y,cell->z()-p.z);
      float dis = dir.Dot(dir1);
      Point p2(cell->x(),cell->y(),cell->z());
      float dis1 = l1.closest_dis(p2);
      if (dis >0 ) dis = dis - dis1;
      if (dis <0 ) dis = dis + dis1;
      if (dis > max_dis){
	max_dis = dis;
	max_point = p2;
      }
    }
    
    for (int i=0;i!=min_cell->Get_all_spacecell().size();i++){
      SpaceCell *cell = min_cell->Get_all_spacecell().at(i);
      TVector3 dir1(cell->x()-p.x,cell->y()-p.y,cell->z()-p.z);
      float dis = dir.Dot(dir1);
      Point p2(cell->x(),cell->y(),cell->z());
      float dis1 = l1.closest_dis(p2);
      if (dis >0 ) dis = dis - dis1;
      if (dis <0 ) dis = dis + dis1;
      if (dis < min_dis){
	min_dis = dis;
	min_point = p2;
      }
    }
    
    
    MergeSpaceCellSelection track_mcells;
    WCP2dToy::ToyWalking walking(max_cell,max_point,min_cell,min_point,mcells_map);
    track_mcells = walking.get_cells();
    
    //std::cout << track_mcells.size() << " " << walking.get_counter() << " " << walking.get_global_counter() << std::endl;
    if (track_mcells.size() > 0){
      double ky, kz;
      if (max_point.x == p.x){
	ky = 0;
	kz = 0;
      }else{
	ky = (max_point.y-p.y)/(max_point.x-p.x);
	kz = (max_point.z-p.z)/(max_point.x-p.x);
      }
      
      
      WCTrack *track = new WCTrack(track_mcells);
      tracks.push_back(track);
      
      WCVertex *vertex1 = new WCVertex(*max_cell);
      vertex1->set_center(max_point);
      vertex1->Add(track);
      vertex1->set_ky(track,ky);
      vertex1->set_kz(track,kz);
      vertices.push_back(vertex1);
      
      
      if (min_point.x == p.x){
	ky = 0;
	kz = 0;
      }else{
	ky = (min_point.y-p.y)/(min_point.x-p.x);
	kz = (min_point.z-p.z)/(min_point.x-p.x);
      }
      
      WCVertex *vertex2 = new WCVertex(*min_cell);
      vertex2->set_center(min_point);
      vertex2->Add(track);
      vertex2->set_ky(track,ky);
      vertex2->set_kz(track,kz);
      vertices.push_back(vertex2);


      //std::cout << tracks.size() << " " << vertices.size() << std::endl;
    }
  }

}


void WCP2dToy::ToyTracking::cosmic_finder_all(WCP2dToy::ToyCrawler& toycrawler){
  MergeSpaceCellMap& mcells_map = toycrawler.Get_mcells_map();
  MergeSpaceCellSelection all_cells;  
  Point p(0,0,0);
  int sum=0;
  for (auto it =toycrawler.Get_mcells_map().begin();it!=toycrawler.Get_mcells_map().end();it++){
    MergeSpaceCell *mcell = it->first;
    all_cells.push_back(mcell);
    p.x += mcell->Get_Center().x * mcell->Get_all_spacecell().size(); 
    p.y += mcell->Get_Center().y * mcell->Get_all_spacecell().size(); 
    p.z += mcell->Get_Center().z * mcell->Get_all_spacecell().size(); 
    sum += mcell->Get_all_spacecell().size(); 
  }

  if (sum >0){
    p.x/=sum;
    p.y/=sum;
    p.z/=sum;
    
    ClusterTrack ct(all_cells.front());
    for (int i=1;i< all_cells.size();i++){
      ct.AddMSCell_anyway(all_cells.at(i));
    }
    ct.SC_Hough(p);
    float theta = ct.Get_Theta();
    float phi = ct.Get_Phi();
    
    Point p1(p.x + sin(theta)*cos(phi), p.y + sin(theta)*sin(phi), p.z + cos(theta));
    Line l1(p,p1);
    TVector3 dir = l1.vec();
    
    
    float max_dis = 0;
    float min_dis = 0;
    MergeSpaceCell *max_cell = all_cells.at(0);
    MergeSpaceCell *min_cell = all_cells.at(0);
    for (int i=0;i!=all_cells.size();i++){
      MergeSpaceCell *mcell = all_cells.at(i);
      TVector3 dir1(mcell->Get_Center().x-p.x,mcell->Get_Center().y-p.y,mcell->Get_Center().z-p.z);
      float dis = dir.Dot(dir1);
      float dis1 = l1.closest_dis(mcell->Get_Center());
      if (dis >0 ) dis = dis - dis1;
      if (dis <0 ) dis = dis + dis1;
      
      if (dis > max_dis){
	max_dis = dis;
	max_cell = mcell;
      }
      if (dis < min_dis){
	min_dis = dis;
	min_cell = mcell;
      }
    }
    
    Point max_point = max_cell->Get_Center();
    Point min_point = min_cell->Get_Center();
    for (int i=0;i!=max_cell->Get_all_spacecell().size();i++){
      SpaceCell *cell = max_cell->Get_all_spacecell().at(i);
      TVector3 dir1(cell->x()-p.x,cell->y()-p.y,cell->z()-p.z);
      float dis = dir.Dot(dir1);
      Point p2(cell->x(),cell->y(),cell->z());
      float dis1 = l1.closest_dis(p2);
      if (dis >0 ) dis = dis - dis1;
      if (dis <0 ) dis = dis + dis1;
      if (dis > max_dis){
	max_dis = dis;
	max_point = p2;
      }
    }
    
    for (int i=0;i!=min_cell->Get_all_spacecell().size();i++){
      SpaceCell *cell = min_cell->Get_all_spacecell().at(i);
      TVector3 dir1(cell->x()-p.x,cell->y()-p.y,cell->z()-p.z);
      float dis = dir.Dot(dir1);
      Point p2(cell->x(),cell->y(),cell->z());
      float dis1 = l1.closest_dis(p2);
      if (dis >0 ) dis = dis - dis1;
      if (dis <0 ) dis = dis + dis1;
      if (dis < min_dis){
	min_dis = dis;
	min_point = p2;
      }
    }
    
    
    MergeSpaceCellSelection track_mcells;
    WCP2dToy::ToyWalking walking(max_cell,max_point,min_cell,min_point,mcells_map);
    track_mcells = walking.get_cells();
    
    //std::cout << track_mcells.size() << " " << walking.get_counter() << " " << walking.get_global_counter() << std::endl;
    if (track_mcells.size() > 0){
      double ky, kz;
      if (max_point.x == p.x){
	ky = 0;
	kz = 0;
      }else{
	ky = (max_point.y-p.y)/(max_point.x-p.x);
	kz = (max_point.z-p.z)/(max_point.x-p.x);
      }
      
      
      WCTrack *track = new WCTrack(track_mcells);
      tracks.push_back(track);
      
      WCVertex *vertex1 = new WCVertex(*max_cell);
      vertex1->set_center(max_point);
      vertex1->Add(track);
      vertex1->set_ky(track,ky);
      vertex1->set_kz(track,kz);
      vertices.push_back(vertex1);
      
      
      if (min_point.x == p.x){
	ky = 0;
	kz = 0;
      }else{
	ky = (min_point.y-p.y)/(min_point.x-p.x);
	kz = (min_point.z-p.z)/(min_point.x-p.x);
      }
      
      WCVertex *vertex2 = new WCVertex(*min_cell);
      vertex2->set_center(min_point);
      vertex2->Add(track);
      vertex2->set_ky(track,ky);
      vertex2->set_kz(track,kz);
      vertices.push_back(vertex2);
      
      
      // std::cout << all_cells.size() << " " << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << " " << 
      //   theta << " " << phi << " " << max_dis << " " << min_dis << " " << track_mcells.size() << std::endl;
    }
  }
}

void WCP2dToy::ToyTracking::RemoveSameTrack(){
  WCTrackSelection to_be_removed;
  //WCTrackSelection parent_track;
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    for (int j=0;j!=vertex->get_tracks().size();j++){
      WCTrack *track1 = vertex->get_tracks().at(j);
      auto it1 = find(to_be_removed.begin(),to_be_removed.end(),track1);
      if (it1 != to_be_removed.end()) continue;
      
      for (int k=0;k!=vertex->get_tracks().size();k++){
	WCTrack *track2 = vertex->get_tracks().at(k);
	if (track2 == track1) continue;
	auto it2 = find(to_be_removed.begin(),to_be_removed.end(),track2);
	if (it2 != to_be_removed.end()) continue;
	
	if (track1->Inside(track2)){
	  to_be_removed.push_back(track1);
	  // parent_track.push_back(track2);
	}

      }
    }
  }

  WCVertexSelection vertex_to_removed;
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    for (int j= 0; j!=to_be_removed.size();j++){
      auto it =find(vertex->get_tracks().begin(),
		    vertex->get_tracks().end(),
		    to_be_removed.at(j));
      if (it != vertex->get_tracks().end())
	vertex->get_tracks().erase(it);
    }
    
    if (vertex->get_ntracks()==0){
      vertex_to_removed.push_back(vertex);
    }else{
      int flag = 1;
      for (int j=0;j!=vertex->get_tracks().size();j++){
	if (vertex->get_msc() == vertex->get_tracks().at(j)->get_end_scells().front() || 
	    vertex->get_msc() == vertex->get_tracks().at(j)->get_end_scells().back()){
	  flag = 0;
	}
      }
      if (flag == 1)
      vertex_to_removed.push_back(vertex);
    }
  }
  //  std::cout << to_be_removed.size() << " " << vertex_to_removed.size() << std::endl;

  for (int i=0;i!=vertex_to_removed.size();i++){
    WCVertex *vertex = vertex_to_removed.at(i);
    
      auto it = find(vertices.begin(),vertices.end(),vertex);
      vertices.erase(it);
      delete vertex;
  }

}




void WCP2dToy::ToyTracking::single_shower_reco(WCP2dToy::ToyCrawler& toycrawler){
  // MergeSpaceCellSelection all_cells;  
  
  // for (auto it =toycrawler.Get_mcells_map().begin();it!=toycrawler.Get_mcells_map().end();it++){
  //   all_cells.push_back(it->first);
  // }
  //  std::cout << all_cells.size() << std::endl;
  
  MergeSpaceCellMap& mcells_map = toycrawler.Get_mcells_map();

  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex =vertices.at(i);
    
    // find vertex containing on track
    int ntrack = 0;
    WCTrack *track = 0;
    for (int j=0;j!=vertex->get_ntracks();j++){
      auto it1 = find(good_tracks.begin(),good_tracks.end(),vertex->get_tracks().at(j));
      auto it2 = find(bad_tracks.begin(),bad_tracks.end(),vertex->get_tracks().at(j));
      if (it1 != good_tracks.end() || it2 != bad_tracks.end()){
	ntrack ++;
	track = vertex->get_tracks().at(j);
      }
    }
    
    // form shower
    if (ntrack == 1){
      MergeSpaceCellSelection empty_cells;
      WCShower *shower = new WCShower(vertex,track,empty_cells,mcells_map);
    
      // std::cout << vertex->Center().x/units::cm << " " << 
      // 	vertex->Center().y/units::cm << " " <<
      // 	vertex->Center().z/units::cm << std::endl;
    // judge if agrees with shower
      if (shower->IsShower(empty_cells)){
	showers.push_back(shower);
      }else{
	delete shower;
      }
    }
  }
  // pick up the furthest one 
  
  float max_dis = 0;
  WCShower *max_shower = 0;
  for (int i=0;i!=showers.size();i++){
    WCShower *shower = showers.at(i);
    if (shower->distance_to_center() > max_dis){
      max_dis = shower->distance_to_center();
      max_shower = shower;
    }
  }

  int flag_qx = 1;
  while(flag_qx){
    flag_qx = 0;
    for (int i=0;i!=showers.size();i++){
      WCShower *shower = showers.at(i);
      if (shower != max_shower){
  	delete shower;
  	showers.erase(showers.begin()+i);
  	flag_qx = 1;
  	break;
      }
    }
  }
  

  for (int i=0;i!=showers.size();i++){
   WCShower *shower = showers.at(i);
   std::cout << "Showers "<< i << " " << " " <<   shower->distance_to_center() << " " << 
      shower->get_all_cells().size() << " " 
  	      << shower->get_vertex()->Center().x/units::cm << " " 
  	      << shower->get_vertex()->Center().y/units::cm << " " 
  	      << shower->get_vertex()->Center().z/units::cm << " " 
  	      << std::endl;
  }

}



void WCP2dToy::ToyTracking::Cleanup_showers(){
  WCVertexSelection not_remove_vertex;
  WCTrackSelection not_remove_track;
  for (int i = 0;i!=showers.size();i++){
    not_remove_vertex.push_back(showers.at(i)->get_vertex());
    not_remove_track.push_back(showers.at(i)->get_track());
  }

  WCVertexSelection to_be_remove_vertex;
  WCTrackSelection to_be_remove_track;
  for (int i=0;i!=good_tracks.size();i++){
    auto it = find(not_remove_track.begin(),not_remove_track.end(),good_tracks.at(i));
    if (it == not_remove_track.end()){
      int flag_qx = 0;
      for (int j=0;j!=showers.size();j++){
	if (showers.at(j)->Contain(good_tracks.at(i))){
	  flag_qx = 1;
	  break;
	}
      }
      if (flag_qx == 1){
	to_be_remove_track.push_back(good_tracks.at(i));
	// for (int j=0;j!=vertices.size();j++){
	//   WCVertex *vertex = vertices.at(j);
	//   WCTrackSelection& vtracks = vertex->get_tracks();
	//   auto it = find(vtracks.begin(),vtracks.end(),good_tracks.at(i));
	//   if (it != vtracks.end()){
	//     vtracks.erase(it);
	//   }
	//   if (vtracks.size()==0){
	//     auto it1 = find(to_be_remove_vertex.begin(),to_be_remove_vertex.end(),vertices.at(j));
	//     if (it1 == to_be_remove_vertex.end())
	//       to_be_remove_vertex.push_back(vertices.at(j));
	//   }
	// }
      }
    }
  }

  for (int i=0;i!=to_be_remove_track.size();i++){
    auto it = find(good_tracks.begin(),good_tracks.end(),to_be_remove_track.at(i));
    bad_tracks.push_back(to_be_remove_track.at(i));
    //to_be_remove_track.at(i)->reset_fine_tracking();
    good_tracks.erase(it);
  }

  // to_be_remove_track.clear();
  
  
  // for (int i=0;i!=bad_tracks.size();i++){
  //   auto it = find(not_remove_track.begin(),not_remove_track.end(),bad_tracks.at(i));
  //   if (it == not_remove_track.end()){
  //     int flag_qx = 0;
  //     for (int j=0;j!=showers.size();j++){
  // 	if (showers.at(j)->Contain(bad_tracks.at(i))){
  // 	  flag_qx = 1;
  // 	  break;
  // 	}
  //     }
  //     if (flag_qx == 1){
  // 	to_be_remove_track.push_back(bad_tracks.at(i));
  // 	for (int j=0;j!=vertices.size();j++){
  // 	  WCVertex *vertex = vertices.at(j);
  // 	  WCTrackSelection& vtracks = vertex->get_tracks();
  // 	  auto it = find(vtracks.begin(),vtracks.end(),bad_tracks.at(i));
  // 	  if (it != vtracks.end()){
  // 	    vtracks.erase(it);
  // 	  }
  // 	  if (vtracks.size()==0){
  // 	    auto it1 = find(to_be_remove_vertex.begin(),to_be_remove_vertex.end(),vertices.at(j));
  // 	    if (it1 == to_be_remove_vertex.end())
  // 	      to_be_remove_vertex.push_back(vertices.at(j));
  // 	  }
  // 	}
  //     }
  //   }
  // }

  // for (int i=0;i!=to_be_remove_track.size();i++){
  //   auto it = find(bad_tracks.begin(),bad_tracks.end(),to_be_remove_track.at(i));
  //   to_be_remove_track.at(i)->reset_fine_tracking();
  //   bad_tracks.erase(it);
  // }

 

  // for (int i=0;i!=to_be_remove_vertex.size();i++){
  //   WCVertex *vertex = to_be_remove_vertex.at(i);
  //   auto it = find(vertices.begin(),vertices.end(),vertex);
  //   if (it != vertices.end()){
  //     vertices.erase(it);
  //     delete vertex;
  //   }
  // }
}



bool  WCP2dToy::ToyTracking::track_shower_reco(WCP2dToy::ToyCrawler& toycrawler){
  bool result = false;

  // find vertex with distance cut
  std::vector<WCVertexSelection> possible_vertex;
  std::vector<WCTrackSelection> possible_track;

  MergeSpaceCellSelection good_cells;
  for (int i= 0;i!=good_tracks.size();i++){
    WCTrack *track = good_tracks.at(i);
    for (int j=0;j!=track->get_centerVP_cells().size();j++){
      good_cells.push_back(track->get_centerVP_cells().at(j));
    }
  }

  WCVertexSelection used_vertices;

  //put bad vertices in first
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    int flag = 0;
    for (int j=0;j!=vertex->get_tracks().size();j++){
      auto it = find(bad_tracks.begin(),bad_tracks.end(),vertex->get_tracks().at(j));
      if (it == bad_tracks.end()){
  	flag = 1;
      }
    }
    if (flag == 0 )
      used_vertices.push_back(vertex);
  }

  


  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    auto it = find(used_vertices.begin(),used_vertices.end(),vertex);

    if (it == used_vertices.end()){
      used_vertices.push_back(vertex);
      // do something about the distance ... 
      WCVertexSelection temp;
      temp.push_back(vertex);
      for (int j=0;j!=vertices.size();j++){
	WCVertex* vertex1 = vertices.at(j);
	auto it1 =find(used_vertices.begin(),used_vertices.end(),vertex1);
	if (it1 == used_vertices.end() && vertex1 != vertex){

	  //Need to look through every thing inside the temp
	  
	  for (int k=0;k!=temp.size();k++){
	    WCVertex *vertex2 = temp.at(k);
	    float dis = sqrt(pow(vertex2->Center().y-vertex1->Center().y,2)
			     +pow(vertex2->Center().z-vertex1->Center().z,2));
	    float dis1 = fabs(vertex2->Center().x-vertex1->Center().x);

	    MergeSpaceCell *vertex_mcell = vertex2->get_msc();
	    MergeSpaceCell *vertex1_mcell = vertex1->get_msc();
	    
	    int flag = 0;
	    if (dis < Singleton<TPCParams>::Instance().get_pitch()*3.3 && dis1 < Singleton<TPCParams>::Instance().get_ts_width()*1.1 && vertex2 != vertex1){ 
	      flag = 1;
	    }else if (vertex2->get_msc() == vertex1->get_msc() && vertex2 != vertex1){
	      flag = 1;
	    }else if (fabs(vertex2->get_msc()->Get_Center().x - vertex1->get_msc()->Get_Center().x) < Singleton<TPCParams>::Instance().get_ts_width()*1.1 
		      && vertex_mcell->Overlap(*vertex1_mcell)
		      && vertex2 != vertex1){
	      flag = 1;
	    }

	    // std::cout << vertex2->get_msc()->Get_Center().x/units::cm << " " << 
	    //   vertex1->get_msc()->Get_Center().x/units::cm << " " << dis << " " 
	    // 	      << flag << std::endl;
	    // std::cout << vertex1 << " " << vertex2 << " " << vertex_mcell->Overlap(*vertex1_mcell) << " " << flag << std::endl;
	    
	    if (flag == 1){
	      used_vertices.push_back(vertex1);
	      temp.push_back(vertex1);
	      break;
	    }
	  }
	}
      }

      // now save all the vertices that are close into temp 


      WCTrackSelection temp_tracks;
      std::vector<TVector3*> temp_dirs;
      
      //judge if satisfy requirement
      //     if (temp.size()>1 || temp.at(0)->get_ntracks()>1){
      for (int j = 0; j!= temp.size();j++){
	for (int k=0;k!=temp.at(j)->get_ntracks();k++){
	  WCTrack* track =  temp.at(j)->get_tracks().at(k);
	  auto it2 = find(good_tracks.begin(),good_tracks.end(),track);
	  auto it3 = find(bad_tracks.begin(),bad_tracks.end(),track);
	  if (it2 != good_tracks.end() || it3 != bad_tracks.end()){
	    temp_tracks.push_back(track);
	    TVector3 *vec = new TVector3(1,temp.at(j)->get_ky(track),temp.at(j)->get_kz(track));
	    temp_dirs.push_back(vec);
	  }
	}
      }

      

      

      // // look at the bad tracks ...
      // for (int j=0;j!=bad_tracks.size();j++){
      // 	int flag_qx = 0;
      // 	for (int k=0;k!=temp.size();k++){
      // 	  float dis1 = sqrt(pow(temp.at(k)->Center().x - bad_tracks.at(j)->get_centerVP_cells().front()->Get_Center().x,2)
      // 			    + pow(temp.at(k)->Center().y - bad_tracks.at(j)->get_centerVP_cells().front()->Get_Center().y,2)
      // 			    + pow(temp.at(k)->Center().z - bad_tracks.at(j)->get_centerVP_cells().front()->Get_Center().z,2));
      // 	  if (dis1 < 1 * units::cm){
      // 	    flag_qx = 1;
      // 	    break;
      // 	  }

      // 	  float dis2 = sqrt(pow(temp.at(k)->Center().x - bad_tracks.at(j)->get_centerVP_cells().back()->Get_Center().x,2)
      // 			    + pow(temp.at(k)->Center().y - bad_tracks.at(j)->get_centerVP_cells().back()->Get_Center().y,2)
      // 			    + pow(temp.at(k)->Center().z - bad_tracks.at(j)->get_centerVP_cells().back()->Get_Center().z,2));
      // 	  if (dis2 < 1*units::cm){
      // 	    flag_qx = 1;
      // 	    break;
      // 	  }
      // 	}
	
      // 	if (flag_qx==1){
      // 	  temp_tracks.push_back(bad_tracks.at(j));
      // 	  TVector3 *vec = new TVector3(0,0,0);
      // 	  temp_dirs.push_back(vec);
      // 	}
      // }

      
      // 


      if (temp_tracks.size() == 2){
	TVector3 *vec1 = temp_dirs.at(0);
	TVector3 *vec2 = temp_dirs.at(1);
	
	double a = (*vec1) * (*vec2);
	
	if (vec1->Mag()!=0 & vec2->Mag()!=0)
	  a = a / vec1->Mag() / vec2->Mag();
       
	// std::cout << temp.at(0)->Center().x/units::cm << " " 
	// 	  << temp.at(0)->Center().y/units::cm << " " 
	// 	  << temp.at(0)->Center().z/units::cm << " " << a << std::endl;
	if (a > -0.9 && a < 0.9){
	  possible_vertex.push_back(temp);
	  possible_track.push_back(temp_tracks);
	}
      }else if (temp_tracks.size()>2){
	possible_vertex.push_back(temp);
	possible_track.push_back(temp_tracks);
      }
      
      for (int j=0;j!=temp_dirs.size();j++){
	delete temp_dirs.at(j);
      }
      
	// std::cout << vertex->Center().x/units::cm << " " 
	// 	  << vertex->Center().y/units::cm << " " 
	// 	  << vertex->Center().z/units::cm << " " << " " 
	// 	  << temp.size() << " " << temp_tracks.size() << std::endl;
	  
	//}
    }
  }


  // Now prepare to remove things ... 
  
  
  for (int i=0;i!=possible_vertex.size();i++){
    WCVertexSelection temp_vertex = possible_vertex.at(i); 
    WCTrackSelection temp_track = possible_track.at(i);
    
    for (int j=0;j!=temp_track.size();j++){
      WCTrack *spec_track = temp_track.at(j);
      
      // find the special vertex
      WCVertex *spec_vertex = 0;
      for (int k=0;k!=temp_vertex.size();k++){
	auto it = find(temp_vertex.at(k)->get_tracks().begin(),temp_vertex.at(k)->get_tracks().end(), spec_track);
	if (it != temp_vertex.at(k)->get_tracks().end()){
	  spec_vertex = temp_vertex.at(k);
	  break;
	}
      }
      if (spec_vertex == 0){
	float min_dis = 1e9;
	for (int k=0;k!=temp_vertex.size();k++){
	  float dis1 = sqrt(pow(temp_vertex.at(k)->Center().x - spec_track->get_centerVP_cells().front()->Get_Center().x,2)+
			    pow(temp_vertex.at(k)->Center().y - spec_track->get_centerVP_cells().front()->Get_Center().y,2)+
			    pow(temp_vertex.at(k)->Center().z - spec_track->get_centerVP_cells().front()->Get_Center().z,2));

	  if (dis1 < min_dis){
	    min_dis = dis1;
	    spec_vertex = temp_vertex.at(k);
	  }

	  float dis2 = sqrt(pow(temp_vertex.at(k)->Center().x - spec_track->get_centerVP_cells().back()->Get_Center().x,2)+
			    pow(temp_vertex.at(k)->Center().y - spec_track->get_centerVP_cells().back()->Get_Center().y,2)+
			    pow(temp_vertex.at(k)->Center().z - spec_track->get_centerVP_cells().back()->Get_Center().z,2));

	  if (dis2 < min_dis){
	    min_dis = dis2;
	    spec_vertex = temp_vertex.at(k);
	  }
	}
      }


      // form excluded cells ...
      MergeSpaceCellSelection exclude_cells;
      WCTrackSelection exclude_tracks;
      //save all exclude tracks ... 
      for (int k=0;k!=temp_track.size();k++){
	WCTrack *curr_track = temp_track.at(k);
	if (curr_track != spec_track){ 
	  exclude_tracks.push_back(curr_track);
	}
      }
      
      int flag_qx = 1;
      while(flag_qx){
      	flag_qx = 0;
      	for (int k=0;k!=good_tracks.size();k++){
      	  WCTrack *curr_track = good_tracks.at(k);
      	  if (curr_track != spec_track){
      	    auto it = find(exclude_tracks.begin(),exclude_tracks.end(),
      			   curr_track);
      	    if (it == exclude_tracks.end()){
      	      for (int k1=0;k1!= exclude_tracks.size();k1++){
      		if (curr_track->IsConnected(exclude_tracks.at(k1))){
      		  exclude_tracks.push_back(curr_track);
      		  break;
      		}
      	      }
      	    }
      	  }
      	}
      }

      for (int k=0;k!=exclude_tracks.size();k++){
       	for (int k1 = 0; k1!=exclude_tracks.at(k)->get_centerVP_cells().size();k1++){
	  auto it = find(exclude_cells.begin(),exclude_cells.end(),exclude_tracks.at(k)->get_centerVP_cells().at(k1));
	  if (it == exclude_cells.end())
	    exclude_cells.push_back(exclude_tracks.at(k)->get_centerVP_cells().at(k1));
	}
      }

      WCShower *shower = new WCShower(spec_vertex,spec_track,exclude_cells,toycrawler.Get_mcells_map());

      bool flag_shower = shower->IsShower(good_cells);

      

      // std::cout << "Xin: " << spec_vertex->Center().x/units::cm << " "
      // 		<< spec_vertex->Center().y/units::cm << " " 
      // 		<< spec_vertex->Center().z/units::cm << " " 
      // 		<< flag_shower << " " << good_cells.size() << " " << shower->get_all_cells().size()  << " " <<  std::endl;
      
      // flag_shower = false;

      if (flag_shower&& shower->IsShower_CheckVertex()){
	showers.push_back(shower);
	// spec_track, spec_vertex, exclude_cells ... 
      }else{
	delete shower;
      }
      


      // }
      //Judge ... 

    } // end loop of temp_tracks ... 

    // Add in a shower without any tracks ... 
    for (int j=0;j!=temp_vertex.size();j++){
      WCVertex *spec_vertex = temp_vertex.at(j);
      int ngood_track = 0;
      WCTrackSelection exclude_tracks;
      for (int k=0;k!=spec_vertex->get_tracks().size();k++){
	WCTrack *spec_track = spec_vertex->get_tracks().at(k);
	auto itq = find(good_tracks.begin(),good_tracks.end(),spec_track);
	if (itq != good_tracks.end()){
	  exclude_tracks.push_back(spec_track);
	}
      }
      // Need to contain more than one good tracks
      if (exclude_tracks.size() <2) continue; // not working 
      
      // work out exclude cells 
      MergeSpaceCellSelection exclude_cells;
      int flag_qx = 1;
      while(flag_qx){
      	flag_qx = 0;
      	for (int k=0;k!=good_tracks.size();k++){
      	  WCTrack *curr_track = good_tracks.at(k);
	  auto it = find(exclude_tracks.begin(),exclude_tracks.end(),
			 curr_track);
	  if (it == exclude_tracks.end()){
	    for (int k1=0;k1!= exclude_tracks.size();k1++){
	      if (curr_track->IsConnected(exclude_tracks.at(k1))){
		exclude_tracks.push_back(curr_track);
		break;
	      }
	    }
      	  }
      	}
      }
      
      for (int k=0;k!=exclude_tracks.size();k++){
       	for (int k1 = 0; k1!=exclude_tracks.at(k)->get_centerVP_cells().size();k1++){
	  auto it = find(exclude_cells.begin(),exclude_cells.end(),exclude_tracks.at(k)->get_centerVP_cells().at(k1));
	  if (it == exclude_cells.end())
	    exclude_cells.push_back(exclude_tracks.at(k)->get_centerVP_cells().at(k1));
	}
      }
      // at least one adajcent cell not inside good tracks 
      MergeSpaceCellSelection& ad_mcells = toycrawler.Get_mcells_map()[spec_vertex->get_msc()];
      int flag_cont = 0;
      for (int k=0;k!=ad_mcells.size();k++){
	auto it = find(exclude_cells.begin(),exclude_cells.end(),ad_mcells.at(k));
	if (it == exclude_cells.end()){
	  flag_cont = 1;
	  break;
	}
      }
      // Set track to zero, and need to make sure the other code are compatible
      if (flag_cont ==1){
	WCShower *shower = new WCShower(spec_vertex,0,exclude_cells,toycrawler.Get_mcells_map());
	bool flag_shower = shower->IsShower(good_cells);
	
	if (flag_shower && shower->IsShower_CheckVertex()){
	  showers.push_back(shower);
	  // spec_track, spec_vertex, exclude_cells ... 
	}else{
	  delete shower;
	}
      }
    }
  }

  // find the maximum shower, and judge if other shower belong to it?
  int max_cells = 0;
  WCShower *max_shower = 0;
  for (int i=0;i!=showers.size();i++){
    WCShower *shower = showers.at(i);
    if (shower->get_all_cells().size() > max_cells){
      max_cells = shower->get_all_cells().size();
      max_shower = shower;
    }
  }
  // go through the rest of shower to judge if they belong to max_shower
  int flag_qx = 1;
  while(flag_qx){
    flag_qx = 0;
    for (int i=0;i!=showers.size();i++){
      WCShower *shower = showers.at(i);
      if (shower != max_shower){
  	if (shower->IsContained(max_shower)){
  	  delete shower;
  	  showers.erase(showers.begin()+i);
  	  flag_qx = 1;
  	  break;
  	}
      }
    }
  }

  // find the second largest shower
  if (showers.size() >2){
    int second_max_cells = 0;
    WCShower *second_max_shower = 0;
    
    for (int i=0;i!=showers.size();i++){
      WCShower *shower = showers.at(i);
      if (shower->get_all_cells().size() > second_max_cells && shower != max_shower){
  	second_max_cells = shower->get_all_cells().size();
  	second_max_shower = shower;
      }
    }

    flag_qx = 1;
    while(flag_qx){
      flag_qx = 0;
      for (int i=0;i!=showers.size();i++){
  	WCShower *shower = showers.at(i);
  	if (shower != max_shower && shower != second_max_shower){
  	  if (shower->IsContained(second_max_shower)){
  	    delete shower;
  	    showers.erase(showers.begin()+i);
  	    flag_qx = 1;
  	    break;
  	  }
  	}
      }
    }
  }

  //Now for the remaining showers, judge if any of two showers are overlapping, 
  //if so, delete the one with larger rms ... 
  flag_qx = 1;
  while(flag_qx){
    flag_qx = 0;
    for (int i=0;i!=showers.size();i++){
      WCShower *shower1 = showers.at(i);
      for (int j=0;j!=showers.size();j++){
  	WCShower *shower2 = showers.at(j);
  	if (shower1!=shower2){
  	  if (shower1->Overlap(shower2)){
  	    shower1->SC_Hough(shower1->get_vertex()->Center());
  	    shower2->SC_Hough(shower2->get_vertex()->Center());
  	    float rms1 = shower1->SC_proj_Hough(shower1->get_vertex()->Center());
  	    float rms2 = shower2->SC_proj_Hough(shower2->get_vertex()->Center());
  	    //std::cout << rms1 << " " << rms2 << " " << shower1->get_all_cells().size() << " " << shower2->get_all_cells().size() << std::endl;
  	    if (rms1 < rms2){
  	      delete shower2;
  	      showers.erase(showers.begin()+j);
  	    }else{
  	      delete shower1;
  	      showers.erase(showers.begin()+i);
  	    }
  	    flag_qx = 1;
  	    break;
  	  }
  	}
      }
      if (flag_qx == 1) break;
    }
  }

  WCShowerSelection bad_showers;

  for (int i=0;i!=showers.size();i++){
    WCShower *shower = showers.at(i);
    int nvertex_outside = 0;
    for (int k=0;k!=vertices.size();k++){
      auto it = find(shower->get_all_cells().begin(),
		     shower->get_all_cells().end(),
		     vertices.at(k)->get_msc());
      if (it == shower->get_all_cells().end())
	nvertex_outside ++;
    }

    if (nvertex_outside >=8){
      bad_showers.push_back(shower);
    }
  }

  for (int i=0;i!=bad_showers.size();i++){
    WCShower *shower = bad_showers.at(i);
    auto it = find(showers.begin(),showers.end(),shower);
    showers.erase(it);
  }


  //  std::cout << showers.size() << std::endl;
  for (int i=0;i!=showers.size();i++){
    WCShower *shower = showers.at(i);
    std::cout << "Showers: "<< i << " " << " " <<  shower->IsShower(good_cells) << " " << 
      shower->get_all_cells().size() << " " 
  	      << shower->get_vertex()->Center().x/units::cm << " " 
  	      << shower->get_vertex()->Center().y/units::cm << " " 
  	      << shower->get_vertex()->Center().z/units::cm << " " 
  	      << std::endl;
  }
  
  
  if (showers.size() >0)
    result = true;

  return result;
}


void WCP2dToy::ToyTracking::form_parallel_tiny_tracks_wovertex(WCP2dToy::ToyCrawler& toycrawler){
  //form all the cells
  MergeSpaceCellMap& mcells_map = toycrawler.Get_mcells_map();
  MergeSpaceCellSelection mcells;
  for (auto it = mcells_map.begin(); it != mcells_map.end(); it++){
      mcells.push_back(it->first);
  }
  
  //examine all the cells to find out outlier 
  for (int i=0;i!=mcells.size();i++){
    MergeSpaceCell *mcell = mcells.at(i);
    new_mcells.push_back(mcell);
    new_mcells_map[mcell] = mcell;
  }

  std::vector<MergeSpaceCellSelection> cluster_msc;
  // Need to cluster these by whether they are connected ...
  for (int i=0;i!=new_mcells.size();i++){
    MergeSpaceCell *mcell = new_mcells.at(i);
    //std::cout << mcell->Get_Center().x/units::cm << std::endl;
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
   	  if (fabs(mcell1->Get_Center().x - mcell->Get_Center().x) < mcell1->thickness() + 0.2*units::mm ){
	    if (mcell1->Overlap(*mcell)){
	      mscs.push_back(mcell);
	      flag = 1;
	      break;
	    }
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
	    
  	    if (fabs(mcell1->Get_Center().x - mcell2->Get_Center().x) < mcell1->thickness() + 0.2*units::mm){ 
  	      if (mcell1->Overlap(*mcell2)){

  		cluster_msc.at(i).insert(cluster_msc.at(i).end(),cluster_msc.at(j).begin(),cluster_msc.at(j).end());
  		cluster_msc.erase(cluster_msc.begin() + j);
		
  		//std::cout << flag << std::endl;
  		flag = 1;
  		break;
  	      }
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

  //std::cout << "Xin: " << cluster_msc.size() << " " << new_mcells.size() << std::endl;

  for (int i=0;i!=cluster_msc.size();i++){
    // Create a vertex which is furthest away from the center ... 
    // find the center
    MergeSpaceCellSelection& mscs_1 = cluster_msc.at(i);

    float sum = 0;
    Point center;
    center.x = 0;
    center.y = 0;
    center.z = 0;

    

    for (int j=0;j!=mscs_1.size();j++){
      sum += mscs_1.at(j)->Get_all_spacecell().size();
      center.x += mscs_1.at(j)->Get_Center().x*mscs_1.at(j)->Get_all_spacecell().size();
      center.y += mscs_1.at(j)->Get_Center().y*mscs_1.at(j)->Get_all_spacecell().size();
      center.z += mscs_1.at(j)->Get_Center().z*mscs_1.at(j)->Get_all_spacecell().size();
    }
    center.x /= sum;
    center.y /= sum;
    center.z /= sum;

    // find the furthest mcell, and cell
    float furthest_dist = 0;
    MergeSpaceCell *furthest_mcell = 0;
    SpaceCell *furthest_cell = 0;

    //std::cout << i << " " << " " << mscs_1.size() << std::endl;
    

    for (int j=0;j!=mscs_1.size();j++){
      MergeSpaceCell *mcell = mscs_1.at(j);
      for (int k=0;k!=mscs_1.at(j)->Get_all_spacecell().size();k++){
	SpaceCell *cell = mscs_1.at(j)->Get_all_spacecell().at(k);
	
	float dist = sqrt(pow(cell->x() - center.x,2)
			  +pow(cell->y() - center.y,2)
			  +pow(cell->z() - center.z,2));
	
	if (dist > furthest_dist){
	  furthest_dist = dist;
	  furthest_mcell = mcell;
	  furthest_cell = cell;
	}
      }
    }
    // create the vertex
    WCVertex *vertex = new WCVertex(*furthest_mcell);
    Point p(furthest_cell->x(),furthest_cell->y(),furthest_cell->z());
    vertex->set_center(p);
    vertices.push_back(vertex);

    // do parallel tracking ... 
    parallel_tracking(vertex,cluster_msc.at(i),toycrawler);
  }

}


void WCP2dToy::ToyTracking::form_parallel_tiny_tracks(WCP2dToy::ToyCrawler& toycrawler){
  
  //form all the cells
  MergeSpaceCellMap& mcells_map = toycrawler.Get_mcells_map();
  MergeSpaceCellSelection mcells;
  for (auto it = mcells_map.begin(); it != mcells_map.end(); it++){
    int flag = 0;
    for (int j=0;j!=showers.size();j++){
      auto it1 = find(showers.at(j)->get_all_cells().begin(),
		      showers.at(j)->get_all_cells().end(),
		      it->first);
      if (it1 != showers.at(j)->get_all_cells().end()){
	flag = 1;
	break;
      }
    }
    if (flag == 0)
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
      nmcell->set_id(mcell->get_id());
      //nmcell->set_mcell(mcell->get_mcell());
      int flag2 = 0;

      for (int j=0;j!=mcell->Get_all_spacecell().size();j++){
	SpaceCell *cell = mcell->Get_all_spacecell().at(j);
	int flag1 = 0;

	for (int k=0;k!=track_no.size();k++){
	  double dist = good_tracks.at(track_no.at(k))->dist_proj(mcell,cell);
	  double dist1 = good_tracks.at(track_no.at(k))->dist(mcell,cell);
	  if (dist < 2 * Singleton<TPCParams>::Instance().get_pitch()
	      && dist1 < 20 * units::mm){
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
   	  if (fabs(mcell1->Get_Center().x - mcell->Get_Center().x) < mcell1->thickness() + 0.2*units::mm ){
	    if (mcell1->Overlap(*mcell)){
	      mscs.push_back(mcell);
	      flag = 1;
	      break;
	    }
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
	    
  	    if (fabs(mcell1->Get_Center().x - mcell2->Get_Center().x) < mcell1->thickness() + 0.2*units::mm){ 
	      if (mcell1->Overlap(*mcell2)){

		cluster_msc.at(i).insert(cluster_msc.at(i).end(),cluster_msc.at(j).begin(),cluster_msc.at(j).end());
		cluster_msc.erase(cluster_msc.begin() + j);
		
		//std::cout << flag << std::endl;
		flag = 1;
		break;
	      }
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
    WCVertex *vertex = vertices.at(0);
    int vertex_ntrack = 0;
    WCVertex *temp_vertex;
    int flag = 0;

    for (int k=0;k!=cluster_msc.at(i).size();k++){
      MergeSpaceCell *mcell1 = cluster_msc.at(i).at(k);
      MergeSpaceCell *mcell2 = new_mcells_map[mcell1];
      for (int j=0;j!=vertices.size();j++){
	temp_vertex = vertices.at(j);
	
	// vertex must contain a good track to start with 
	int flag1 = 0;
	for (int i1 = 0;i1!=temp_vertex->get_ntracks();i1++){
	  auto it = find(good_tracks.begin(),good_tracks.end(),temp_vertex->get_tracks().at(i1));
	  if (it != good_tracks.end()){
	    flag1 ++;
	  }
	}
       
	if (flag1 >0){
	  MergeSpaceCell *vertex_cell = temp_vertex->get_msc();
	  if (fabs(vertex_cell->Get_Center().x - mcell2->Get_Center().x) < 2*vertex_cell->thickness() + 0.5*units::mm && vertex_cell->Overlap(*mcell2)){

	    if (flag1 > vertex_ntrack){
	      flag = 1;
	      vertex = temp_vertex;
	      vertex_ntrack = flag1;
	    }

	  }
	}
      }
      if (flag == 1 && vertex_ntrack>=2)
	break;
    }
    
    //std::cout << "Xin: " << flag << " " << vertex->Center().x/units::cm << " "  << vertex->get_ntracks()<< std::endl;

    if (flag == 1){
      // if yes, do parallel track finder ... 
      //MergeSpaceCellSelection mcells = cluster_msc.at(i);
      parallel_tracking(vertex,cluster_msc.at(i),toycrawler);
      //std::cout << flag << " " << vertex->Center().x/units::cm << " " << std::endl;
    }else{
        // If not, do short track, just save the nmcell etc
      WCTrack *track = new WCTrack(cluster_msc.at(i));
      short_tracks.push_back(track);
    }
  }

  // std::cout << short_tracks.size() << " " << parallel_tracks.size() << std::endl;
  //  

}

void WCP2dToy::ToyTracking::parallel_tracking(WCVertex *vertex, MergeSpaceCellSelection &mcells, WCP2dToy::ToyCrawler& toycrawler){
 
  // re-organize mcells to do a merge, merge anything that are close ... 
  MergeSpaceCellMap& mcells_map_old = toycrawler.Get_mcells_map();  
  
  MergeSpaceCellMap mcells_map; // new map ... 
  WCP::MergeSpaceCellMap1 new_mcells_map1; 

  WCTrack *track_1st=0;
  WCTrack *track_2nd=0;

  bool first_track = false;
  bool second_track = false;
  float first_track_direction;
  float second_track_direction;

  MergeSpaceCell *vertex_cell = vertex->get_msc();

  int flag_qx = 0;
  MergeSpaceCellSelection mcells_old; // old cells ... 
  for (int i=0;i!=mcells.size();i++){
    auto it = find(mcells_old.begin(),mcells_old.end(),new_mcells_map[mcells.at(i)]);
    if (it == mcells_old.end()){
      mcells_old.push_back(new_mcells_map[mcells.at(i)]);
      // if (new_mcells_map[mcells.at(i)] == vertex_cell && new_mcells_map[mcells.at(i)]!=mcells.at(i)) flag_qx = 1;
      if (new_mcells_map[mcells.at(i)] == vertex_cell) flag_qx = 1;
    }
  }
  
  // force to not generate new cells ... 
  // this option (flag_qx == 1)  is in conflict with the track merging ... 
  flag_qx = 0;

  
  //  std::cout << mcells.size() << " " << mcells_old.size() << " " << flag_qx << std::endl;
  // mcells = mcells_old;
  // mcells_map = mcells_map_old;

  if (flag_qx == 0){
    mcells_map = mcells_map_old;
  }else{
    // Create a new set of merged cells including vertex and anything between, 
    std::vector<MergeSpaceCellSelection> mcells_merge;
    MergeSpaceCellSelection used_mcells_merge;
    while (used_mcells_merge.size() < mcells_old.size()){
      if (mcells_merge.size() == 0 ){
	MergeSpaceCellSelection temp;
	temp.push_back(mcells_old.at(0));
	used_mcells_merge.push_back(mcells_old.at(0));
	mcells_merge.push_back(temp);
      }
      
      int i = mcells_merge.size()-1;
      for (int j=0;j!=mcells_merge.at(i).size();j++){
	MergeSpaceCell *mcell1 = mcells_merge.at(i).at(j);
	for (int k=1;k<mcells_old.size();k++){
	  MergeSpaceCell *mcell = mcells_old.at(k);
	  auto it = find(used_mcells_merge.begin(),used_mcells_merge.end(),mcell);
	  if (it != used_mcells_merge.end()) continue;
	  
	  
	  if (fabs(mcell->Get_Center().x - mcell1->Get_Center().x)<0.1*units::mm){
	    
	    if (mcell->Overlap(*mcell1,1.5)){
	      mcells_merge.at(i).push_back(mcell);
	      used_mcells_merge.push_back(mcell);
	    }
	  }
	}
      }
      
      for (int k=1; k<mcells_old.size();k++){
	MergeSpaceCell *mcell = mcells_old.at(k);
	auto it = find(used_mcells_merge.begin(),used_mcells_merge.end(),mcell);
	if (it != used_mcells_merge.end()) continue;
	MergeSpaceCellSelection temp;
	temp.push_back(mcells_old.at(k));
	used_mcells_merge.push_back(mcells_old.at(k));
	mcells_merge.push_back(temp);
	break;
      }
    }
    //    std::cout << mcells_old.size() << " " << mcells_merge.size() << std::endl;
    //    form a new set of mergecells
    //update vertex_cell ... 
    MergeSpaceCellSelection new_merge_cells;
    for (int i = 0; i !=mcells_merge.size();i++){
      MergeSpaceCell *new_merge_cell = new MergeSpaceCell();
      for (int j=0;j!=mcells_merge.at(i).size();j++){
	// std::cout << i << " " << j << " " << mcells_merge.at(i).at(j)->Get_Center().x/units::cm
	// 	  << " " << mcells_merge.at(i).at(j)->Get_Center().y/units::cm
	// 	  << " " << mcells_merge.at(i).at(j)->Get_Center().z/units::cm << std::endl;

	new_mcells_map1[mcells_merge.at(i).at(j)] = new_merge_cell;
	for (int k=0;k!=mcells_merge.at(i).at(j)->Get_all_spacecell().size();k++){
	  new_merge_cell->AddSpaceCell(mcells_merge.at(i).at(j)->Get_all_spacecell().at(k));
	}
      }
      new_merge_cells.push_back(new_merge_cell);
    }
    //std::cout << new_merge_cells.size() << std::endl;
    // form mcells_map
    for (int i=0;i!=new_merge_cells.size();i++){
      MergeSpaceCell *mcell1 = new_merge_cells.at(i);
      MergeSpaceCellSelection temp;
      for (int j=0;j!=new_merge_cells.size();j++){
	MergeSpaceCell *mcell2 = new_merge_cells.at(j);
	if (fabs(mcell1->Get_Center().x-mcell2->Get_Center().x) < mcell1->thickness() + 0.5*units::mm && 
	    fabs(mcell1->Get_Center().x-mcell2->Get_Center().x) > 0.5*units::mm){
	  if (mcell1->Overlap(*mcell2))
	    temp.push_back(mcell2);
	}
      }
      mcells_map[mcell1] = temp;
    }
    
    //std::cout << vertex_cell << std::endl;
    //vertex_cell = new_mcells_map1[vertex_cell];
    //std::cout << vertex_cell <<std::endl;
  }
  
  

  // find the furthest point and mcell, 
  Point p = vertex->Center();
  MergeSpaceCellSelection used_mcells; 
  
  double max_dis1 = 0;
  
  MergeSpaceCell *max_mcell;
  
  for (int i=0;i!=mcells.size();i++){
    MergeSpaceCell *mcell = mcells.at(i);
    
    Point center = mcell->Get_Center();
    mcell->CalMinMax();
    double dy = mcell->get_dy();
    double dz = mcell->get_dz();
    
    double max_dis = 0;
    double dis[5];
    dis[0] = sqrt(pow(p.x-center.x,2)+pow(p.y-center.y,2)+pow(p.z-center.z,2));
    dis[1] = sqrt(pow(p.x-center.x,2)+pow(p.y-center.y-dy,2)+pow(p.z-center.z-dz,2));
    dis[2] = sqrt(pow(p.x-center.x,2)+pow(p.y-center.y+dy,2)+pow(p.z-center.z-dz,2));
    dis[3] = sqrt(pow(p.x-center.x,2)+pow(p.y-center.y-dy,2)+pow(p.z-center.z+dz,2));
    dis[4] = sqrt(pow(p.x-center.x,2)+pow(p.y-center.y+dy,2)+pow(p.z-center.z+dz,2));
    
    for (int j=0;j!=5;j++){
      if (max_dis < dis[j]) max_dis = dis[j];
    }
    
    // std::cout << max_dis << " " << center.x/units::cm << " " << center.y/units::cm << " " << center.z/units::cm << " " << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << " " << dy/units::cm << " " << dz/units::cm << std::endl;
    
    if (max_dis > max_dis1){
      max_dis1 = max_dis;
      max_mcell = mcell;
    }
  }
  
  max_mcell = new_mcells_map[max_mcell];   //original cell

  max_dis1 = 0;
  Point max_point;
  //find the furtherst point 
  for (int i=0;i!=max_mcell->Get_all_spacecell().size();i++){
    SpaceCell *cell = max_mcell->Get_all_spacecell().at(i);
    double dis = sqrt(pow(p.x-cell->x(),2)+pow(p.y-cell->y(),2)+pow(p.z-cell->z(),2));
    if (dis > max_dis1){
      max_dis1 = dis;
      max_point.x = cell->x();
      max_point.y = cell->y();
      max_point.z = cell->z();
    }
  }




  // if (flag_qx == 0){
  // }else{ // new cellls ... 
  //   max_mcell = new_mcells_map1[new_mcells_map[max_mcell]]; 
  // }
  // std::cout << max_point.x/units::cm << " " << max_point.y/units::cm << " " << max_point.z/units::cm << std::endl;
  // std::cout << max_mcell->Get_Center().x << " " << max_mcell->Get_Center().y << std::endl;
  // walk back to the mcell containing the vertex ...  Progressive tracking ...  
  
  
  MergeSpaceCellSelection track_mcells;
  double dist;
  
  if (flag_qx == 0){
    // max_cell --> vertex_cell, find the shortest path 
    WCP2dToy::ToyWalking walking(max_mcell,max_point,vertex_cell,vertex->Center(),mcells_map);
    track_mcells = walking.get_cells();
    dist = walking.get_dist();
  }else{ // map over the new cells ... 
    WCP2dToy::ToyWalking walking(new_mcells_map1[max_mcell],max_point,new_mcells_map1[vertex_cell],vertex->Center(),mcells_map);
    dist = walking.get_dist();
    // How to deal with these??? 
    track_mcells = walking.get_cells();
  }
  
  if (dist < 1e9){
    double ky, kz;
    if (max_point.x == p.x){
      ky = 0;
      kz = 0;
    }else{
      ky = (max_point.y-p.y)/(max_point.x-p.x);
      kz = (max_point.z-p.z)/(max_point.x-p.x);
    }
    
    // judge the angle ... 
    // find the track that contain max_cell
    // judge angle between these two tracks
    float min_angle = 1e9;
    for (int qx = 0;qx!=vertex->get_tracks().size();qx++){
      WCTrack *track1 = vertex->get_tracks().at(qx);
      //auto it = find(track1->get_centerVP_cells().begin(),track1->get_centerVP_cells().end(),max_mcell);
      //if (it != track1->get_centerVP_cells().end()){
      TVector3 vec1(max_point.x-p.x,max_point.y-p.y,max_point.z-p.z);
      TVector3 vec2;
      
      for (int qx1 = 0; qx1!= wct_wcv_map[track1].size();qx1++){
  	if (wct_wcv_map[track1].at(qx1) != vertex ){
  	  vec2.SetXYZ(wct_wcv_map[track1].at(qx1)->Center().x - vertex->Center().x,
  		      wct_wcv_map[track1].at(qx1)->Center().y - vertex->Center().y,
  		      wct_wcv_map[track1].at(qx1)->Center().z - vertex->Center().z);
  	  if (vec1.Angle(vec2)/3.1415926*180. < min_angle){
  	    min_angle = vec1.Angle(vec2)/3.1415926*180.;
  	  }
  	  //	  std::cout << mcells.size() << "Angle: " << vec1.Angle(vec2)/3.1415926*180. << std::endl;
  	  break;
  	}
      }
           
    }

    //  std::cout << mcells.size() << " " << min_angle << std::endl;

    if (mcells.size()==1 && min_angle < 25 || min_angle < 5){
    }else{
      
      WCTrack *track = 0;
      if (flag_qx ==0){
	track = new WCTrack(track_mcells);
      }else{
	// How to deal with this ??? 
	track = new WCTrack(track_mcells);
      }
      
      // for (int i1 = 0; i1!= track_mcells.size();i1++){
      // 	std::cout << track_mcells.at(i1)->Get_Center().x/units::cm << " " 
      // 		  << track_mcells.at(i1)->Get_Center().y/units::cm << " " 
      // 		  << track_mcells.at(i1)->Get_Center().z/units::cm << " " 
      // 		  << std::endl;
      // }
      // MergeSpaceCellSelection temp_track_mcells;
      // for (int i1 = 0; i1!= track_mcells.size();i1++){
      // 	if (new_mcells_map.find(track_mcells.at(i1))==new_mcells_map.end()){
      // 	  temp_track_mcells.push_back(track_mcells.at(i1));
      // 	}else{
      // 	  temp_track_mcells.push_back(new_mcells_map[track_mcells.at(i1)]);
      // 	}
      // }
      // WCTrack *track = new WCTrack(temp_track_mcells);

      //put everything into used cells ... (new cells or old cells ...)
      used_mcells.insert(used_mcells.end(),track_mcells.begin(),track_mcells.end());
      first_track = true;
      track_1st = track;
      if (track_mcells.size() >=2){
	first_track_direction = track_mcells.at(track_mcells.size()-2)->Get_Center().x - 
					track_mcells.at(track_mcells.size()-1)->Get_Center().x;
      }else{
	first_track_direction = 0;
      }

      vertex->Add(track);
      vertex->set_ky(track,ky);
      vertex->set_kz(track,kz);
      tracks.push_back(track);
      parallel_tracks.push_back(track);
      
      // std::cout << vertex->Center().x/units::cm << " " << vertex->Center().y/units::cm << " " << vertex->Center().z/units::cm << " " << std::endl;
      // std::cout << max_point.x/units::cm << " " << max_point.y/units::cm << " " << max_point.z/units::cm << std::endl;

      WCVertex *vertex1 = new WCVertex(*max_mcell);
      vertex1->set_center(max_point);
      vertex1->Add(track);
      vertex1->set_ky(track,ky);
      vertex1->set_kz(track,kz);
      vertices.push_back(vertex1);

      WCVertexSelection qvertex;
      qvertex.push_back(vertex);
      qvertex.push_back(vertex1);
      wct_wcv_map[track] = qvertex; 
      
      fine_tracking(1);
      // //check existing cells ... 
      //       for (int i1 = 0; i1!= track_mcells.size();i1++){
      // 	MergeSpaceCell *mcell = track_mcells.at(i1);
      // 	for (int j1 = 0; j1!= track_mcells.at(i1)->Get_all_spacecell().size();j1++){
      // 	  SpaceCell *cell = track_mcells.at(i1)->Get_all_spacecell().at(j1);
      // 	  double dist = track->dist_proj(mcell,cell);
      // 	  if (dist > 4.5*units::mm)
      // 	    std::cout << track->dist_proj(mcell,cell)/units::mm << " " << mcell->Get_Center().x/units::cm  << " " << mcell->Get_Center().y/units::cm << " " << mcell->Get_Center().z/units::cm <<
      // std::endl;
      // 	}
      //       }

    }
  }
  // std::cout << track_mcells.size() << std::endl;
  // form vertex ... 

 


  
  if (first_track){
    // what about the duplicated tracks ??? 
    // mcells_map
    // used_mcells
    MergeSpaceCellSelection mcells_save;
    for (int i = 0;i!=mcells.size();i++){
      MergeSpaceCell *mcell = new_mcells_map[mcells.at(i)];
      if (flag_qx == 1){
  	mcell = new_mcells_map1[mcell];
      }
      auto it = find(used_mcells.begin(),used_mcells.end(),mcell);
      int flag = 0;
      if (it == used_mcells.end()){
  	for (int j = 0;j!=mcells_map[mcell].size();j++){
  	  auto it1 = find(used_mcells.begin(),used_mcells.end(),mcells_map[mcell].at(j));
	  
  	  int flag1 = 0;
  	  for (int k=0;k!=mcells.size();k++){
  	    MergeSpaceCell *mcell1 = new_mcells_map[mcells.at(k)];
  	    if (flag_qx == 1)
  	      mcell1 = new_mcells_map1[mcell1];
  	    if ( mcells_map[mcell].at(j) == mcell1  ){
  	      flag1 = 1;
  	      break;
  	    }
  	  }
  	  if ( ((it1 == used_mcells.end() && flag1 == 1) || mcells_map[mcell].at(j) == vertex_cell) 
  	       && (mcell->Get_Center().x-vertex_cell->Get_Center().x) * first_track_direction < 0 ){
  	    flag = 1;
  	    break;
  	  }
  	}
      }
      if (flag == 1){
  	mcells_save.push_back(mcells.at(i));
      }
    }





    

    // std::cout << mcells.size() << " " << mcells_save.size() << std::endl;
    // for (int i=0;i!=mcells_save.size();i++){
    //   std::cout << mcells_save.at(i)->Get_Center().x/units::cm << " " 
    // 	      << mcells_save.at(i)->Get_Center().y/units::cm << " " 
    // 	      << mcells_save.at(i)->Get_Center().z/units::cm << " " <<std::endl;
    // }
    
    
    // try to find the second track only ...
    if (mcells_save.size() > 0){
      max_dis1 = 0;
      for (int i=0;i!=mcells_save.size();i++){
  	MergeSpaceCell *mcell = mcells_save.at(i);
	
  	Point center = mcell->Get_Center();
  	mcell->CalMinMax();
  	double dy = mcell->get_dy();
  	double dz = mcell->get_dz();
	
  	double max_dis = 0;
  	double dis[5];
  	dis[0] = sqrt(pow(p.x-center.x,2)+pow(p.y-center.y,2)+pow(p.z-center.z,2));
  	dis[1] = sqrt(pow(p.x-center.x,2)+pow(p.y-center.y-dy,2)+pow(p.z-center.z-dz,2));
  	dis[2] = sqrt(pow(p.x-center.x,2)+pow(p.y-center.y+dy,2)+pow(p.z-center.z-dz,2));
  	dis[3] = sqrt(pow(p.x-center.x,2)+pow(p.y-center.y-dy,2)+pow(p.z-center.z+dz,2));
  	dis[4] = sqrt(pow(p.x-center.x,2)+pow(p.y-center.y+dy,2)+pow(p.z-center.z+dz,2));
	
  	for (int j=0;j!=5;j++){
  	  if (max_dis < dis[j]) max_dis = dis[j];
  	}
	
  	// std::cout << max_dis << " " << center.x/units::cm << " " << center.y/units::cm << " " << center.z/units::cm << " " << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << " " << dy/units::cm << " " << dz/units::cm << std::endl;
	
  	if (max_dis > max_dis1){
  	  max_dis1 = max_dis;
  	  max_mcell = mcell;
  	}
      }
      
      if (flag_qx == 0){
  	max_mcell = new_mcells_map[max_mcell];
      }else{
  	max_mcell = new_mcells_map1[new_mcells_map[max_mcell]];
      }
      max_dis1 = 0;
      
      //find the furtherst point 
      for (int i=0;i!=max_mcell->Get_all_spacecell().size();i++){
  	SpaceCell *cell = max_mcell->Get_all_spacecell().at(i);
  	double dis = sqrt(pow(p.x-cell->x(),2)+pow(p.y-cell->y(),2)+pow(p.z-cell->z(),2));
  	if (dis > max_dis1){
  	  max_dis1 = dis;
  	  max_point.x = cell->x();
  	  max_point.y = cell->y();
  	  max_point.z = cell->z();
  	}
      }
      
      //std::cout << max_point.x << " " << max_point.y << " " <<max_point.z << std::endl;
      
      
      // max_cell --> vertex_cell, find the shortest path 
      WCP2dToy::ToyWalking walking1(max_mcell,max_point,vertex_cell,vertex->Center(),mcells_map);
      track_mcells = walking1.get_cells();
      dist = walking1.get_dist();
      
      //std::cout << dist << " " << track_mcells.size() << std::endl;
      
      if (dist < 1e9){
  	double ky, kz;
  	if (max_point.x == p.x){
  	  ky = 0;
  	  kz = 0;
  	}else{
  	  ky = (max_point.y-p.y)/(max_point.x-p.x);
  	  kz = (max_point.z-p.z)/(max_point.x-p.x);
  	}
	
  	// judge the angle ... 
  	// find the track that contain max_cell
  	// judge angle between these two tracks
  	float min_angle = 1e9;
  	for (int qx = 0;qx!=vertex->get_tracks().size();qx++){
  	  WCTrack *track1 = vertex->get_tracks().at(qx);
  	  //auto it = find(track1->get_centerVP_cells().begin(),track1->get_centerVP_cells().end(),max_mcell);
  	  //if (it != track1->get_centerVP_cells().end()){
  	  TVector3 vec1(max_point.x-p.x,max_point.y-p.y,max_point.z-p.z);
  	  TVector3 vec2;
	  
  	  for (int qx1 = 0; qx1!= wct_wcv_map[track1].size();qx1++){
  	    if (wct_wcv_map[track1].at(qx1) != vertex ){
  	      vec2.SetXYZ(wct_wcv_map[track1].at(qx1)->Center().x - vertex->Center().x,
  			  wct_wcv_map[track1].at(qx1)->Center().y - vertex->Center().y,
  			  wct_wcv_map[track1].at(qx1)->Center().z - vertex->Center().z);
  	      if (vec1.Angle(vec2)/3.1415926*180. < min_angle){
  		min_angle = vec1.Angle(vec2)/3.1415926*180.;
  	      }
  	      //	  std::cout << mcells.size() << "Angle: " << vec1.Angle(vec2)/3.1415926*180. << std::endl;
  	      break;
  	    }
  	  }
	  
  	}
	
  	if (mcells.size()==1 && min_angle < 25 || min_angle < 5){
  	}else{
  	  WCTrack *track = new WCTrack(track_mcells);
  	  second_track = true;
  	  track_2nd = track;

  	  if (track_mcells.size() >=2){
  	    second_track_direction = track_mcells.at(track_mcells.size()-2)->Get_Center().x - 
  					    track_mcells.at(track_mcells.size()-1)->Get_Center().x;
  	  }else{
  	    second_track_direction = 0;
  	  }

  	  //put everything into used cells ...
  	  used_mcells.insert(used_mcells.end(),track_mcells.begin(),track_mcells.end());
	  
  	  vertex->Add(track);
  	  vertex->set_ky(track,ky);
  	  vertex->set_kz(track,kz);
  	  tracks.push_back(track);
  	  parallel_tracks.push_back(track);
	  
  	  WCVertex *vertex1 = new WCVertex(*max_mcell);
  	  vertex1->set_center(max_point);
  	  vertex1->Add(track);
  	  vertex1->set_ky(track,ky);
  	  vertex1->set_kz(track,kz);
  	  vertices.push_back(vertex1);

  	  WCVertexSelection qvertex;
  	  qvertex.push_back(vertex);
  	  qvertex.push_back(vertex1);
  	  wct_wcv_map[track] = qvertex;

  	  fine_tracking(1);
  	}
      }
    }
  }



  //deal with the remaining cells ... 
  
  if (first_track){
    MergeSpaceCellSelection mcells_save;
    for (int i = 0;i!=mcells.size();i++){
      MergeSpaceCell *mcell = new_mcells_map[mcells.at(i)];
      if (flag_qx == 1){
  	mcell = new_mcells_map1[mcell];
      }
      auto it = find(used_mcells.begin(),used_mcells.end(),mcell);
      int flag = 0;
      if (it == used_mcells.end()){
  	for (int j = 0;j!=mcells_map[mcell].size();j++){
  	  auto it1 = find(used_mcells.begin(),used_mcells.end(),mcells_map[mcell].at(j));
	  
  	  int flag1 = 0;
  	  for (int k=0;k!=mcells.size();k++){
  	    MergeSpaceCell *mcell1 = new_mcells_map[mcells.at(k)];
  	    if (flag_qx == 1)
  	      mcell1 = new_mcells_map1[mcell1];
  	    if ( mcells_map[mcell].at(j) == mcell1  ){
  	      flag1 = 1;
  	      break;
  	    }
  	  }
  	  if ( (it1 == used_mcells.end() && flag1 == 1)  ){
  	    flag = 1;
  	    break;
  	  }
  	}
      }
      
      // good and try to look at this cell 
      if (flag == 1){
  	int status_flag = 0;

  	if (first_track){
  	  // deal with the parallel track
  	  WCTrack *track =track_1st;
  	  // find out the closest merged cell
  	  double closest_dist = 1e9;
  	  MergeSpaceCell *closest_mcell = 0;
  	  for (int i1 = 0;i1!=track->get_all_cells().size();i1++){
  	    double dist = sqrt(pow(mcells.at(i)->Get_Center().y - track->get_all_cells().at(i1)->Get_Center().y,2)
  			       + pow(mcells.at(i)->Get_Center().z - track->get_all_cells().at(i1)->Get_Center().z,2));
  	    if (dist < closest_dist && 
  		fabs(mcells.at(i)->Get_Center().x - track->get_all_cells().at(i1)->Get_Center().x) < 
  		Singleton<TPCParams>::Instance().get_ts_width()*1.1
  		){
  	      closest_dist = dist;
  	      closest_mcell = track->get_all_cells().at(i1);
  	    }
  	  }
	  
  	  if (closest_mcell!=0){
  	    // check distance
  	    int flag1 = 0;
  	    for (int i1 = 0; i1!= mcells.at(i)->Get_all_spacecell().size();i1++){
  	      SpaceCell *cell = mcells.at(i)->Get_all_spacecell().at(i1);
  	      double dist = track->dist_proj(closest_mcell,cell);
  	      if (dist > 1.5 * Singleton<TPCParams>::Instance().get_ts_width()
  		  && (mcell->Get_Center().x-vertex_cell->Get_Center().x) * first_track_direction > 0
  		  || dist > 0.75 * Singleton<TPCParams>::Instance().get_ts_width()
  		  && (mcell->Get_Center().x-vertex_cell->Get_Center().x) * first_track_direction < 0
  		  ){
  		//std::cout << dist/units::mm << " " << mcells.at(i)->Get_Center().x/units::cm  << " " << mcells.at(i)->Get_Center().y/units::cm << " " << mcells.at(i)->Get_Center().z/units::cm << std::endl;
  		flag1 = 1;
  		break;
  	      }
  	      // 
  	    }
  	    if (flag1 == 0){
  	      track->get_associated_cells().push_back(mcells.at(i));
  	      status_flag = 1;
  	    }
  	  }
  	}
	
  	if (status_flag == 0 && second_track){
  	  // deal with the parallel track
  	  WCTrack *track =track_2nd;
  	  // find out the closest merged cell
  	  double closest_dist = 1e9;
  	  MergeSpaceCell *closest_mcell = 0;
  	  for (int i1 = 0;i1!=track->get_all_cells().size();i1++){
  	    double dist = sqrt(pow(mcells.at(i)->Get_Center().y - track->get_all_cells().at(i1)->Get_Center().y,2)
  			       + pow(mcells.at(i)->Get_Center().z - track->get_all_cells().at(i1)->Get_Center().z,2));
  	    if (dist < closest_dist && 
  		fabs(mcells.at(i)->Get_Center().x - track->get_all_cells().at(i1)->Get_Center().x) < 
  		Singleton<TPCParams>::Instance().get_ts_width()*1.1
  		){
  	      closest_dist = dist;
  	      closest_mcell = track->get_all_cells().at(i1);
  	    }
  	  }
	  
  	  if(closest_mcell !=0){
  	    // check distance
  	    int flag1 = 0;
  	    for (int i1 = 0; i1!= mcells.at(i)->Get_all_spacecell().size();i1++){
  	      SpaceCell *cell = mcells.at(i)->Get_all_spacecell().at(i1);
  	      double dist = track->dist_proj(closest_mcell,cell);
  	      if (dist > 1.5 * Singleton<TPCParams>::Instance().get_ts_width()
  		  && (mcell->Get_Center().x-vertex_cell->Get_Center().x) * second_track_direction > 0
  		  || dist > 0.75 * Singleton<TPCParams>::Instance().get_ts_width()
  		  && (mcell->Get_Center().x-vertex_cell->Get_Center().x) * second_track_direction < 0
  		  ){
  		//std::cout << dist/units::mm << " " << mcells.at(i)->Get_Center().x/units::cm  << " " << mcells.at(i)->Get_Center().y/units::cm << " " << mcells.at(i)->Get_Center().z/units::cm << std::endl;
  		flag1 = 1;
  		break;
  	      }
  	      // 
  	    }
  	    if (flag1 == 0){
  	      track->get_associated_cells().push_back(mcells.at(i));
  	      status_flag = 2;
  	    }
  	  }
  	}

  	if (status_flag == 0){
  	  mcells_save.push_back(mcells.at(i));
  	}
      }
    }

    // do short tracks ... //should we connect these and help??? 
    for (int i=0;i!=mcells_save.size();i++){
      MergeSpaceCellSelection temp;
      temp.push_back(mcells_save.at(i));
      WCTrack *track = new WCTrack(temp);
      short_tracks.push_back(track);
    //   std::cout << mcells_save.at(i)->Get_Center().x/units::cm << " " 
    // 		<< mcells_save.at(i)->Get_Center().y/units::cm << " " 
    // 		<< mcells_save.at(i)->Get_Center().z/units::cm << " " <<std::endl;
    }

  }
	  



}


bool WCP2dToy::ToyTracking::grow_track_fill_gap(WCP2dToy::ToyCrawler& toycrawler){
  bool result = false;
  MergeSpaceCellMap& mcells_map = toycrawler.Get_mcells_map();
  
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    if (vertex->get_ntracks() == 1){
      WCTrack *track = vertex->get_tracks().at(0);
      
      //only extend good track
      // auto itqx = find(good_tracks.begin(),good_tracks.end(),track);
      // if (itqx == good_tracks.end()) continue;
      // 

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
	
	MergeSpaceCell *prev_cell1 = 0;
	for (int j=0;j!=cells.size();j++){
	  MergeSpaceCell *cell = cells.at(j);
	  auto it = find(track->get_all_cells().begin(), track->get_all_cells().end(), cell);
	  if (it != track->get_all_cells().end()){
	    prev_cell1 = cell;
	    break;
	  }
	}

	//std::cout << track->get_all_cells().size() << std::endl;

	for (int j=0;j!=cells.size();j++){
	  MergeSpaceCell *cell = cells.at(j);
	  auto it = find(track->get_all_cells().begin(), track->get_all_cells().end(), cell);
	  if (it == track->get_all_cells().end() && fabs(cell->Get_Center().x - final_cell->Get_Center().x)>0.1*units::mm && (cell->Get_Center().x - final_cell->Get_Center().x)*(final_cell->Get_Center().x - prev_cell1->Get_Center().x)>0){
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
	    // if (itqx != bad_tracks.end()){
	    //   track->get_centerVP_cells().insert(track->get_all_cells().begin(),saved_cells.at(qx));
	    // }
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

  CheckVertices(toycrawler); 


  WCVertexSelection temp_vertices;
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex1 = vertices.at(i);
    auto it = find(temp_vertices.begin(),temp_vertices.end(),vertex1);
    if (it != temp_vertices.end()) continue;
    for (int j=0;j!=vertices.size();j++){
      WCVertex *vertex2 = vertices.at(j);
      auto it1 = find(temp_vertices.begin(),temp_vertices.end(),vertex2);
      if (it1 != temp_vertices.end()) continue;
      if (vertex1 != vertex2 && vertex1->get_msc() == vertex2->get_msc()){
  	if (vertex1->get_ntracks() > vertex2->get_ntracks()){
  	  // vertex1->AddVertex(vertex2);
  	  for (int k=0;k!=vertex2->get_tracks().size();k++){
  	    WCTrack *track = vertex2->get_tracks().at(k);
  	    auto it2 = find(vertex1->get_tracks().begin(),
  			    vertex1->get_tracks().end(),
  			    track);
  	    if (it2 == vertex1->get_tracks().end()){
  	      vertex1->get_tracks().push_back(track);
  	      vertex1->set_ky(track,vertex2->get_ky(track));
  	      vertex1->set_kz(track,vertex2->get_kz(track));
  	    }
  	    //std::cout << vertex1->get_ntracks() << " " << vertex2->get_ntracks() << std::endl;
  	    vertex1->FindVertex();
  	    //std::cout << vertex1->get_ntracks() << " " << vertex2->get_ntracks() << " " << vertex1->get_fit_success() << std::endl;
  	  }
  	  temp_vertices.push_back(vertex2);
  	}else{
  	  for (int k=0;k!=vertex1->get_tracks().size();k++){
  	    WCTrack *track = vertex1->get_tracks().at(k);
  	    auto it2 = find(vertex2->get_tracks().begin(),
  			    vertex2->get_tracks().end(),
  			    track);
  	    if (it2 == vertex2->get_tracks().end()){
  	      vertex2->get_tracks().push_back(track);
  	      vertex2->set_ky(track,vertex1->get_ky(track));
  	      vertex2->set_kz(track,vertex1->get_kz(track));
  	    }
  	    //std::cout << vertex2->get_ntracks() << " " << vertex1->get_ntracks() << std::endl;
  	    vertex2->FindVertex();
  	    //std::cout << vertex2->get_ntracks() << " " << vertex1->get_ntracks() << " " << vertex2->get_fit_success() << std::endl;
  	  }
  	  // vertex2->AddVertex(vertex1);
  	  temp_vertices.push_back(vertex1);
  	}
      }
      //std::cout << vertices.at(i) << " " << vertices.at(i)->get_msc() << std::endl;
    }
    //std::cout << i << " " << vertices.size() << std::endl;
  }
  
  for (int i=0;i!=temp_vertices.size();i++){
    auto it = find(vertices.begin(),vertices.end(),temp_vertices.at(i));
    if (it != vertices.end()){
      delete temp_vertices.at(i);
      vertices.erase(it);
    }
  }

  return result;
}


void WCP2dToy::ToyTracking::deal_wiggle_tracks(){
  
  //wiggle_vertices.clear();
  
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
	  if (dis < Singleton<TPCParams>::Instance().get_pitch() * 6.7
	      && fabs(end_p1.x-p.x) < Singleton<TPCParams>::Instance().get_ts_width() * 2.1 ){
	    vertex_cell1 = track->get_end_scells().at(0);
	    vertex_cell2 = track->get_end_scells().at(1);
	    flag = 1;
	    break;
	  }
	  
	  dis = sqrt(pow(end_p2.x - p.x,2) + pow(end_p2.y - p.y,2) + pow(end_p2.z - p.z,2));
	  //std::cout << dis/units::cm << std::endl;
	  if (dis < Singleton<TPCParams>::Instance().get_pitch() * 6.7
	      && fabs(end_p2.x-p.x) < Singleton<TPCParams>::Instance().get_ts_width() * 2.1){
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

	  //wiggle_vertices.push_back(vertex);
	  
	  WCVertex *vertex1 = new WCVertex(*vertex_cell2);
	  vertex1->get_tracks().push_back(track);
	  vertices.push_back(vertex1);
	  vertex1->FindVertex();
	  //wiggle_vertices.push_back(vertex);
	}
	//std::cout << vertices.size() << std::endl;
	
	//      std::cout << track->get_end_scells().at(0)->Get_Center().x/units::cm << " " << 
	//	track->get_end_scells().at(1)->Get_Center().x/units::cm << std::endl;
      }
    }
  }

  
}

bool WCP2dToy::ToyTracking::IsThisShower(WCP2dToy::ToyCrawler& toycrawler){
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
      // std::cout << "a: " << time << std::endl;
    }else{
      time_mcells_set2.insert(time);
      //std::cout << "b: " << time << std::endl;
    }
  }
  
  //1.5 max x, min x, max/min y, max/minz
  float min_x=1e9, max_x=-1e9;
  float min_y=1e9, max_y=-1e9;
  float min_z=1e9, max_z=-1e9;
  for (auto it = mcells_map.begin(); it!= mcells_map.end();it ++ ){
    MergeSpaceCell *mscell = it->first;
    if (min_x > mscell->Get_Center().x) min_x = mscell->Get_Center().x;
    if (max_x < mscell->Get_Center().x) max_x = mscell->Get_Center().x;
    if (min_y > mscell->Get_Center().y-mscell->get_dy()) min_y = mscell->Get_Center().y-mscell->get_dy();
    if (max_y < mscell->Get_Center().y+mscell->get_dy()) max_y = mscell->Get_Center().y+mscell->get_dy();
    if (min_z > mscell->Get_Center().z-mscell->get_dz()) min_z = mscell->Get_Center().z-mscell->get_dz();
    if (max_z < mscell->Get_Center().z+mscell->get_dz()) max_z = mscell->Get_Center().z+mscell->get_dz();
  }
  
  
  
  // std::cout << fabs(max_x-min_x) << " " << sqrt(pow(max_y-min_y,2)+pow(max_z-min_z,2)) << std::endl;
  
  
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


  std::cout << "Number of Tracks " << ntracks << " "  << time_mcells_set.size() << " " << track_cluster.size() << " " << time_mcells_set1.size() << " " << time_mcells_set2.size() << " " << fabs(max_x-min_x)/sqrt(pow(max_y-min_y,2)+pow(max_z-min_z,2)) << std::endl; 
  
  
  
  
  //need better tuning .... 
  // no tracks 
  if (ntracks == 0  ){
    if (time_mcells_set.size() >= 4){
      return true;
    }else{
      return false;
    }
  }else{
    //unacounted ones spam more than 10
    if (time_mcells_set.size() >=10){
      // special deal with the parallel tracks
      if (fabs(max_x-min_x)/sqrt(pow(max_y-min_y,2)+pow(max_z-min_z,2))>0.1){
	
	// find out the structure of the remaing stuff...
	MergeSpaceCellSelection remaining_mcells;
	for (auto it = mcells_map.begin(); it!= mcells_map.end();it ++ ){
	  MergeSpaceCell *mscell = it->first;
	  auto it1 = find(all_track_mscs.begin(),all_track_mscs.end(),mscell);
	  if (it1 == all_track_mscs.end()) 
	    remaining_mcells.push_back(mscell);
	}
	// find center
	Point p1(0,0,0);
	int sum = 0;
	for (int i =0;i!=remaining_mcells.size();i++){
	  MergeSpaceCell *mcell = remaining_mcells.at(i);
	  p1.x += mcell->Get_Center().x * mcell->Get_all_spacecell().size();
	  p1.y += mcell->Get_Center().y * mcell->Get_all_spacecell().size();
	  p1.z += mcell->Get_Center().z * mcell->Get_all_spacecell().size();
	  sum += mcell->Get_all_spacecell().size();
	}
	p1.x /= sum;
	p1.y /= sum;
	p1.z /= sum;
	
	// find the furthest one, define the axis
	float max_dis = 0;
	Point max_point;
	for (int i=0;i!=remaining_mcells.size();i++){
	  MergeSpaceCell *mcell = remaining_mcells.at(i);
	  float dis = sqrt(pow(mcell->Get_Center().x - p1.x,2) + 
			   pow(mcell->Get_Center().y - p1.y,2) + 
			   pow(mcell->Get_Center().z - p1.z,2));
	  if (dis > max_dis){
	    max_dis = dis;
	    max_point = mcell->Get_Center();
	  }
	}
	Line l1(max_point,p1);
	TVector3 dir = l1.vec();
	float transverse_dis=0;
	float max_long_dis=-1e9;
	float min_long_dis=1e9;
	// find the maximum distance along it and transverse to it
	for (int i=0; i!=remaining_mcells.size();i++){
	  Point p2 = remaining_mcells.at(i)->Get_Center();
	  if (l1.closest_dis(p2) > transverse_dis)
	    transverse_dis = l1.closest_dis(p2);
	  TVector3 v1(p2.x - p1.x, p2.y-p1.y,p2.z-p1.z);
	  if (v1.Dot(dir)/dir.Mag() > max_long_dis)
	    max_long_dis = v1.Dot(dir)/dir.Mag() ;
	  if (v1.Dot(dir)/dir.Mag() < min_long_dis)
	    min_long_dis = v1.Dot(dir)/dir.Mag() ;
	}
	
	// form ratio and judge .. 
	//std::cout << transverse_dis << " " << transverse_dis * 2/(max_long_dis - min_long_dis) << std::endl;
	
	//std::cout << remaining_mcells.size() << std::endl;
	if (transverse_dis * 2/(max_long_dis - min_long_dis) >0.1)
	  return true;
      } 
    }else{
      if (track_cluster.size() >=3 && time_mcells_set.size() >=4){
	if (fabs(max_x-min_x)/sqrt(pow(max_y-min_y,2)+pow(max_z-min_z,2))>0.1){
	  // find out the structure of the remaing stuff...
	  MergeSpaceCellSelection remaining_mcells;
	  for (auto it = mcells_map.begin(); it!= mcells_map.end();it ++ ){
	    MergeSpaceCell *mscell = it->first;
	    auto it1 = find(all_track_mscs.begin(),all_track_mscs.end(),mscell);
	    if (it1 == all_track_mscs.end()) 
	      remaining_mcells.push_back(mscell);
	  }
	  // find center
	  Point p1(0,0,0);
	  int sum = 0;
	  for (int i =0;i!=remaining_mcells.size();i++){
	    MergeSpaceCell *mcell = remaining_mcells.at(i);
	    p1.x += mcell->Get_Center().x * mcell->Get_all_spacecell().size();
	    p1.y += mcell->Get_Center().y * mcell->Get_all_spacecell().size();
	    p1.z += mcell->Get_Center().z * mcell->Get_all_spacecell().size();
	    sum += mcell->Get_all_spacecell().size();
	  }
	  p1.x /= sum;
	  p1.y /= sum;
	  p1.z /= sum;
	  
	  // find the furthest one, define the axis
	  float max_dis = 0;
	  Point max_point;
	  for (int i=0;i!=remaining_mcells.size();i++){
	    MergeSpaceCell *mcell = remaining_mcells.at(i);
	    float dis = sqrt(pow(mcell->Get_Center().x - p1.x,2) + 
			     pow(mcell->Get_Center().y - p1.y,2) + 
			     pow(mcell->Get_Center().z - p1.z,2));
	    if (dis > max_dis){
	      max_dis = dis;
	      max_point = mcell->Get_Center();
	    }
	  }
	  Line l1(max_point,p1);
	  TVector3 dir = l1.vec();
	  float transverse_dis=0;
	  float max_long_dis=-1e9;
	  float min_long_dis=1e9;
	  // find the maximum distance along it and transverse to it
	  for (int i=0; i!=remaining_mcells.size();i++){
	    Point p2 = remaining_mcells.at(i)->Get_Center();
	    if (l1.closest_dis(p2) > transverse_dis)
	      transverse_dis = l1.closest_dis(p2);
	    TVector3 v1(p2.x - p1.x, p2.y-p1.y,p2.z-p1.z);
	    if (v1.Dot(dir)/dir.Mag() > max_long_dis)
	      max_long_dis = v1.Dot(dir)/dir.Mag() ;
	    if (v1.Dot(dir)/dir.Mag() < min_long_dis)
	      min_long_dis = v1.Dot(dir)/dir.Mag() ;
	  }
	  
	  // form ratio and judge .. 
	  //std::cout << transverse_dis << " " << transverse_dis * 2/(max_long_dis - min_long_dis) << std::endl;
	  
	  //std::cout << remaining_mcells.size() << std::endl;
	  if (transverse_dis * 2/(max_long_dis - min_long_dis) >0.1)
	    return true;
	  
	}
      }else if (track_cluster.size()==2 && time_mcells_set.size() >=10){
	if (fabs(max_x-min_x)/sqrt(pow(max_y-min_y,2)+pow(max_z-min_z,2))>0.1){
	  // find out the structure of the remaing stuff...
	  MergeSpaceCellSelection remaining_mcells;
	  for (auto it = mcells_map.begin(); it!= mcells_map.end();it ++ ){
	    MergeSpaceCell *mscell = it->first;
	    auto it1 = find(all_track_mscs.begin(),all_track_mscs.end(),mscell);
	    if (it1 == all_track_mscs.end()) 
	      remaining_mcells.push_back(mscell);
	  }
	  // find center
	  Point p1(0,0,0);
	  int sum = 0;
	  for (int i =0;i!=remaining_mcells.size();i++){
	    MergeSpaceCell *mcell = remaining_mcells.at(i);
	    p1.x += mcell->Get_Center().x * mcell->Get_all_spacecell().size();
	    p1.y += mcell->Get_Center().y * mcell->Get_all_spacecell().size();
	    p1.z += mcell->Get_Center().z * mcell->Get_all_spacecell().size();
	    sum += mcell->Get_all_spacecell().size();
	  }
	  p1.x /= sum;
	  p1.y /= sum;
	  p1.z /= sum;
	  
	  // find the furthest one, define the axis
	  float max_dis = 0;
	  Point max_point;
	  for (int i=0;i!=remaining_mcells.size();i++){
	    MergeSpaceCell *mcell = remaining_mcells.at(i);
	    float dis = sqrt(pow(mcell->Get_Center().x - p1.x,2) + 
			     pow(mcell->Get_Center().y - p1.y,2) + 
			     pow(mcell->Get_Center().z - p1.z,2));
	    if (dis > max_dis){
	      max_dis = dis;
	      max_point = mcell->Get_Center();
	    }
	  }
	  Line l1(max_point,p1);
	  TVector3 dir = l1.vec();
	  float transverse_dis=0;
	  float max_long_dis=-1e9;
	  float min_long_dis=1e9;
	  // find the maximum distance along it and transverse to it
	  for (int i=0; i!=remaining_mcells.size();i++){
	    Point p2 = remaining_mcells.at(i)->Get_Center();
	    if (l1.closest_dis(p2) > transverse_dis)
	      transverse_dis = l1.closest_dis(p2);
	    TVector3 v1(p2.x - p1.x, p2.y-p1.y,p2.z-p1.z);
	    if (v1.Dot(dir)/dir.Mag() > max_long_dis)
	      max_long_dis = v1.Dot(dir)/dir.Mag() ;
	    if (v1.Dot(dir)/dir.Mag() < min_long_dis)
	      min_long_dis = v1.Dot(dir)/dir.Mag() ;
	  }
	  
	  // form ratio and judge .. 
	  //std::cout << transverse_dis << " " << transverse_dis * 2/(max_long_dis - min_long_dis) << std::endl;
	  
	  //std::cout << remaining_mcells.size() << std::endl;
	  if (transverse_dis * 2/(max_long_dis - min_long_dis) >0.1)
	    return true;
	  
	}
      }else{
	if (time_mcells_set1.size() > time_mcells_set2.size() * 0.75 && time_mcells_set1.size() > 5){
	  
	  //special deal with the parallel tracks 
	  if (fabs(max_x-min_x)/sqrt(pow(max_y-min_y,2)+pow(max_z-min_z,2))>0.1){
	    // find out the structure of the remaing stuff...
	    MergeSpaceCellSelection remaining_mcells;
	    for (auto it = mcells_map.begin(); it!= mcells_map.end();it ++ ){
	      MergeSpaceCell *mscell = it->first;
	      auto it1 = find(all_track_mscs.begin(),all_track_mscs.end(),mscell);
	      if (it1 == all_track_mscs.end()) 
		remaining_mcells.push_back(mscell);
	    }
	    // find center
	    Point p1(0,0,0);
	    int sum = 0;
	    for (int i =0;i!=remaining_mcells.size();i++){
	      MergeSpaceCell *mcell = remaining_mcells.at(i);
	      p1.x += mcell->Get_Center().x * mcell->Get_all_spacecell().size();
	      p1.y += mcell->Get_Center().y * mcell->Get_all_spacecell().size();
	      p1.z += mcell->Get_Center().z * mcell->Get_all_spacecell().size();
	      sum += mcell->Get_all_spacecell().size();
	    }
	    p1.x /= sum;
	    p1.y /= sum;
	    p1.z /= sum;
	    
	    // find the furthest one, define the axis
	    float max_dis = 0;
	    Point max_point;
	    for (int i=0;i!=remaining_mcells.size();i++){
	      MergeSpaceCell *mcell = remaining_mcells.at(i);
	      float dis = sqrt(pow(mcell->Get_Center().x - p1.x,2) + 
			       pow(mcell->Get_Center().y - p1.y,2) + 
			       pow(mcell->Get_Center().z - p1.z,2));
	      if (dis > max_dis){
		max_dis = dis;
		max_point = mcell->Get_Center();
	      }
	    }
	    Line l1(max_point,p1);
	    TVector3 dir = l1.vec();
	    float transverse_dis=0;
	    float max_long_dis=-1e9;
	    float min_long_dis=1e9;
	    // find the maximum distance along it and transverse to it
	    for (int i=0; i!=remaining_mcells.size();i++){
	      Point p2 = remaining_mcells.at(i)->Get_Center();
	      if (l1.closest_dis(p2) > transverse_dis)
		transverse_dis = l1.closest_dis(p2);
	      TVector3 v1(p2.x - p1.x, p2.y-p1.y,p2.z-p1.z);
	      if (v1.Dot(dir)/dir.Mag() > max_long_dis)
		max_long_dis = v1.Dot(dir)/dir.Mag() ;
	      if (v1.Dot(dir)/dir.Mag() < min_long_dis)
		min_long_dis = v1.Dot(dir)/dir.Mag() ;
	    }
	    
	    // form ratio and judge .. 
	    //std::cout << transverse_dis << " " << transverse_dis * 2/(max_long_dis - min_long_dis) << std::endl;
	    
	    //std::cout << remaining_mcells.size() << std::endl;
	    if (transverse_dis * 2/(max_long_dis - min_long_dis) >0.1)
	      return true;
	    
	  }
	}
      }
    }
  }
  
  

  return false;
  


}


void WCP2dToy::ToyTracking::cleanup_bad_tracks(){
  WCVertexSelection to_be_removed;

  for (int i=0;i!=good_tracks.size();i++){
    if (good_tracks.at(i)->IsBadTrack()){      
      bad_tracks.push_back(good_tracks.at(i));
      //bad_tracks.at(i)->reset_fine_tracking();
      
      // for (int j=0;j!=vertices.size();j++){
      //   WCVertex *vertex = vertices.at(j);
      // 	WCTrackSelection& vtracks = vertex->get_tracks();
      // 	auto it = find(vtracks.begin(),vtracks.end(),good_tracks.at(i));
      
      // if (it != vtracks.end()){
      //   vtracks.erase(it);
      // }
      // if (vtracks.size()==0){
      //   //eliminate this vertex
      //   auto it1 = find(to_be_removed.begin(),to_be_removed.end(),vertices.at(j));
      //   if (it1 == to_be_removed.end())
      //     to_be_removed.push_back(vertices.at(j));
      // }
      // }
    }  
  }

  for (int i=0;i!=bad_tracks.size();i++){
    auto it = find(good_tracks.begin(),good_tracks.end(),bad_tracks.at(i));
    good_tracks.erase(it);
  }
  
  //deal with vertices ... 
  


  // for (int i=0;i!=vertices.size();i++){
  //   auto it1 = find(to_be_removed.begin(),to_be_removed.end(),vertices.at(i));
  //   if (vertices.at(i)->get_ntracks()==0 && it1==to_be_removed.end())
  //     to_be_removed.push_back(vertices.at(i));
  // }

  // for (int i=0;i!=to_be_removed.size();i++){
  //   WCVertex *vertex = to_be_removed.at(i);
  //   auto it = find(vertices.begin(),vertices.end(),vertex);
  //   if (it != vertices.end()){
  //     vertices.erase(it);
  //     delete vertex;
  //   }
  // }

}

void WCP2dToy::ToyTracking::update_maps1(){
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


void WCP2dToy::ToyTracking::update_maps(int flag){
  wct_wcv_map.clear();
  wcv_wct_map.clear();

  for (int i=0; i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    WCTrackSelection tracks = vertex->get_tracks();
    wcv_wct_map[vertex] = tracks;
    
    for (int j=0;j!=tracks.size();j++){
      WCTrack *track = tracks.at(j);

      // IF FLAG == 1, REQUIRE THE TRACK TO BE IN GOOD TRACK LIST
      if (flag==1){
	auto it = find(good_tracks.begin(),good_tracks.end(),track);
	auto it1 = find(parallel_tracks.begin(),parallel_tracks.end(),track);
	if (it == good_tracks.end() && it1 == parallel_tracks.end()) continue;
      }
      //

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

void WCP2dToy::ToyTracking::fine_tracking(int flag){
  for (int i=0;i!=tracks.size();i++){
    //std::cout << i << " " << tracks.size() << std::endl;
    WCTrack *track = tracks.at(i);
    if (wct_wcv_map.find(track)!=wct_wcv_map.end()){
      // std::cout << wct_wcv_map[track].size() << std::endl;
      if (wct_wcv_map[track].size()==2){
	//std::cout << i << "abc1 " << std::endl;
	WCVertex *vertex1 = wct_wcv_map[track].at(0);
	//std::cout << i << "abc2 " << std::endl;
	WCVertex *vertex2 = wct_wcv_map[track].at(1);
	//	std::cout << i << "abc3 " << std::endl;
	if (vertex1 == vertex2) continue;
	Point p1 = vertex1->Center();
	Point p2 = vertex2->Center();

	// std::cout << p1.x/units::cm << " " << p1.y/units::cm << " " <<p1.z/units::cm 
	// 	  << " " << vertex1->get_msc()->Get_Center().x/units::cm << " " << vertex1->get_msc()->Get_Center().y/units::cm
	// 	  << " " << vertex1->get_msc()->Get_Center().z/units::cm 
	// 	  << " " << p2.x/units::cm << " " << p2.y/units::cm << " " <<p2.z/units::cm 
	//           << " " << vertex2->get_msc()->Get_Center().x/units::cm << " " << vertex2->get_msc()->Get_Center().y/units::cm
	// 	  << " " << vertex2->get_msc()->Get_Center().z/units::cm 
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
	if (flag == 0 )
	  track->reset_fine_tracking();

	//	std::cout << np1 << " " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << ky1 << " " << ky2 << " " << np2 << " " << 
	//	 p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm <<
	// " " << kz1 << " " << kz2 << std::endl;

	track->fine_tracking(np1,p1,ky1,kz1,np2,p2,ky2,kz2,flag);
	if (track->get_centerVP().size() > 0){
	  auto itt = find(good_tracks.begin(),good_tracks.end(),track);
	  auto itt1 = find(parallel_tracks.begin(),parallel_tracks.end(),track);
	  if (itt == good_tracks.end() && itt1 == parallel_tracks.end())
	    good_tracks.push_back(track);
	}
	//std::cout << i << "abc5 " << std::endl;
      }
    }
  }
}

bool WCP2dToy::ToyTracking::ExamineVertex(WCVertex* vertex, WCP2dToy::ToyCrawler& toycrawler){
  Point vertex_location = vertex->Center();
  
  if (fabs(vertex_location.x - vertex->get_msc()->Get_Center().x) < Singleton<TPCParams>::Instance().get_ts_width() * 5 ){
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
	//std::cout << cell1 << " " << cell2 << std::endl;

	while( (cell2->Get_Center().x - cell1->Get_Center().x) * (cell2->Get_Center().x - center->Get_Center().x) < 0 || cell2 == cell1){
	  
	  
	  //	std::cout << cell2->Get_Center().x/units::cm << " " << cell2 << std::endl;
	  
	  // start to crawl ... 
	  MergeSpaceCellSelection cells2 = mcells_map[cell2];
	  
	  MergeSpaceCell *min_cell3=0;
	  double min = 1e9;
	  
	  for (int k=0;k!=cells2.size();k++){
	    MergeSpaceCell *cell3 = cells2.at(k);
	    if (dis * (cell3->Get_Center().x - cell2->Get_Center().x) < 0)
	      continue;
	    
	    double dy = cell3->get_dy();
	    double dz = cell3->get_dz();
	    
	    if (dy == 0) dy = Singleton<TPCParams>::Instance().get_pitch()/2.; 
	    if (dz == 0) dz = Singleton<TPCParams>::Instance().get_pitch()/2.; 
	    
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



void WCP2dToy::ToyTracking::CheckVertices(WCP2dToy::ToyCrawler& toycrawler){
  
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
	    
	    if (dy == 0) dy = Singleton<TPCParams>::Instance().get_pitch()/2.;
	    if (dz == 0) dz = Singleton<TPCParams>::Instance().get_pitch()/2.;
	    
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

  
  // Check if any of the two vertices are the same, if so merge them
  int flag_qx = 1;
  while(flag_qx){
    flag_qx = 0;
    for (int i=0;i!=vertices.size();i++){
      WCVertex *vertex1 = vertices.at(i);
      for (int j=0;j!=vertices.size();j++){
  	WCVertex *vertex2 = vertices.at(j);
  	if (vertex1!=vertex2){
  	  if (vertex1->get_ntracks() > vertex2->get_ntracks()){
  	    if (vertex1->MergeVertex(vertex2)){
	      //std::cout << vertex2->get_ntracks() << " " << vertex1->get_ntracks() << std::endl;
  	      flag_qx = 1;
  	      vertices.erase(vertices.begin()+j);
  	      break;
  	    }
  	  }else{
  	    if (vertex2->MergeVertex(vertex1)){
	      //std::cout << vertex2->get_ntracks() << " " << vertex1->get_ntracks() << std::endl;
  	      flag_qx = 1;
  	      vertices.erase(vertices.begin()+i);
  	      break;
  	    }
  	  }
  	}
      }
      if (flag_qx==1)
  	break;
    }
  }


}


void WCP2dToy::ToyTracking::CleanUpVertex(){
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


void WCP2dToy::ToyTracking::OrganizeTracks(int flag){
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


void WCP2dToy::ToyTracking::Crawl(){
  // first shift the vertex locations
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);

    //if (i==1){
    WCVertex *vertex1 = vertex->OrganizeTracks();
    
    if (vertex1 !=0){
      vertices.push_back(vertex1);
    }
    //std::cout << vertex->get_ntracks() << " " << vertex->Center().x/units::cm << " " << vertex->Center().y/units::cm << " " << vertex->Center().z/units::cm << std::endl;
    //}
  }
}

void WCP2dToy::ToyTracking::BreakTracks(){
  // Now need to break the track?
  //std::cout << "abc" << std::endl;

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

  //std::cout << "abc: " << std::endl;
  // second break the track if there is a direction change ... 
  break_tracks.clear();
  WCTrackSelection finished_tracks;
  WCVertexSelection NewVertices;
  
  for (int i=0;i!=vertices.size();i++){
    WCVertex *vertex = vertices.at(i);
    //std::cout << "abc1: " << std::endl;
    break_tracks = vertex->BreakTracksAngle(finished_tracks);
    //std::cout << "abc2: " << std::endl;
    if (break_tracks.size()>0){

      for (int j=1;j<break_tracks.size();j++){
      	tracks.push_back(break_tracks.at(j));
      }
      // Form new vertices
      for (int j=0;j<break_tracks.size()-1;j++){
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
  //std::cout << "abc: " << std::endl;

  for (int i=0;i!=NewVertices.size();i++){
    vertices.push_back(NewVertices.at(i));
  }


}


void WCP2dToy::ToyTracking::MergeVertices(int flag){
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


void WCP2dToy::ToyTracking::RemoveSame(){
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
  	  std::cout << "remove " << vertex1->Center().x/units::cm << " " <<
  	    vertex1->Center().y/units::cm << " " << vertex1->Center().z/units::cm << " " <<
  	    vertex2->Center().x/units::cm << " " <<
  	    vertex2->Center().y/units::cm << " " << vertex2->Center().z/units::cm << " " <<std::endl;
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



void WCP2dToy::ToyTracking::CreateVertices(ToyCrawler& toycrawler){
  std::map<MergeSpaceCell*,WCTrackSelection> msc_wct_map;

  //  std::cout << toycrawler.Get_allMCT().size() << std::endl;

   // fill in tracks ... 
  for (int i=0;i!=toycrawler.Get_allMCT().size();i++){
    toycrawler.Get_allMCT().at(i)->Organize();
    
    // int flag = 0;
    // for (int j=0;j!=tracks.size();j++){
    // }
    // if (flag == 1) continue;

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

void WCP2dToy::ToyTracking::Associate(){
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


WCP2dToy::ToyTracking::~ToyTracking(){
  for (int i=0;i!=tracks.size();i++){
    delete tracks.at(i);
  }
  for (int i=0;i!=vertices.size();i++){
    delete vertices.at(i);
  }
}
