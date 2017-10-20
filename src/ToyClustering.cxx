#include "WireCell2dToy/ToyClustering.h"

using namespace WireCell;

void WireCell2dToy::Clustering_jump_gap_cosmics(WireCell::PR3DClusterSelection& live_clusters, double dis){
  for (size_t i=0;i!=live_clusters.size();i++){
    PR3DCluster* cluster_1 = live_clusters.at(i);
    for (size_t j=i+1;j<live_clusters.size();j++){
      PR3DCluster* cluster_2 = live_clusters.at(j);
      std::pair<SlimMergeGeomCell*,SlimMergeGeomCell*> mcells = WireCell2dToy::Get_Closest_MCells_Clusters(cluster_1,cluster_2,dis);
      
    }
  }
}

std::pair<WireCell::SlimMergeGeomCell*, WireCell::SlimMergeGeomCell*> WireCell2dToy::Get_Closest_MCells_Clusters(WireCell::PR3DCluster *cluster1, WireCell::PR3DCluster *cluster2, double dis){
  cluster1->Create_point_cloud();
  cluster2->Create_point_cloud();

  
}


void WireCell2dToy::Clustering_live_dead(WireCell::PR3DClusterSelection& live_clusters, WireCell::PR3DClusterSelection& dead_clusters){

   // dead to live clusters mapping ... 
   std::map<PR3DCluster*,std::vector<PR3DCluster*>> dead_live_cluster_mapping;
   std::map<PR3DCluster*,std::vector<SMGCSelection>> dead_live_mcells_mapping;
   
   // form map between live and dead clusters ... 
   for (size_t i = 0; i!=live_clusters.size(); i++){
     for (size_t j = 0; j!= dead_clusters.size(); j++){
       SMGCSelection mcells = (live_clusters.at(i))->Is_Connected(dead_clusters.at(j),2);
       int live_cluster_id = live_clusters.at(i)->get_cluster_id();
       int dead_cluster_id = dead_clusters.at(j)->get_cluster_id();
       // if (live_cluster_id==5 || live_cluster_id == 21 || live_cluster_id == 76){
       // 	 if (dead_cluster_id==8 || dead_cluster_id==100)
       // 	   std::cout << live_cluster_id << " " << dead_cluster_id << " " << mcells.size() << std::endl;
       // }
       if ( mcells.size()>0 ){
	 if (dead_live_cluster_mapping.find(dead_clusters.at(j))==dead_live_cluster_mapping.end() ){
	   std::vector<PR3DCluster*> temp_clusters;
	   temp_clusters.push_back(live_clusters.at(i));
	   dead_live_cluster_mapping[dead_clusters.at(j)] = temp_clusters;
	   std::vector<SMGCSelection> temp_mcells;
	   temp_mcells.push_back(mcells);
	   dead_live_mcells_mapping[dead_clusters.at(j)] = temp_mcells;
	 }else{
	   dead_live_cluster_mapping[dead_clusters.at(j)].push_back(live_clusters.at(i));
	   dead_live_mcells_mapping[dead_clusters.at(j)].push_back(mcells);
	 }
       }
	   
       // if (mcells.size()>0)
       // 	 std::cout << mcells.size() << " " <<
       // 	   live_clusters.at(i)->get_cluster_id() << " "
       // 		   << live_clusters.at(i)->get_num_mcells() << " "
       // 		   << live_clusters.at(i)->get_num_time_slices() << " " 
       // 		   << dead_clusters.at(j)->get_cluster_id() << " "
       // 		   << dead_clusters.at(j)->get_num_mcells() << " "
       // 		   << dead_clusters.at(j)->get_num_time_slices() << std::endl;
     }
   }
  
  std::set<std::pair<PR3DCluster*, PR3DCluster*>> tested_pairs;
  std::set<std::pair<PR3DCluster*, PR3DCluster*>> to_be_merged_pairs;

  
  
  for (auto it = dead_live_cluster_mapping.begin(); it!= dead_live_cluster_mapping.end(); it++){
    PR3DCluster* the_dead_cluster = (*it).first;
    std::vector<PR3DCluster*> connected_live_clusters = (*it).second;
    std::vector<SMGCSelection> connected_live_mcells = dead_live_mcells_mapping[the_dead_cluster];
    if (connected_live_clusters.size()>1){
      //      std::cout << the_dead_cluster->get_cluster_id() << " " << connected_live_clusters.size() << std::endl;

      // std::cout << "Candidates: ";
      // for (size_t i=0; i!= connected_live_clusters.size(); i++){
      //   PR3DCluster* cluster_1 = connected_live_clusters.at(i);
      // 	std::cout << cluster_1->get_cluster_id() << " ";
      // }
      // std::cout << std::endl;
      
      for (size_t i=0; i!= connected_live_clusters.size(); i++){
        PR3DCluster* cluster_1 = connected_live_clusters.at(i);
	//	std::cout << cluster_1->get_cluster_id() << " ";
	SMGCSelection mcells_1 = connected_live_mcells.at(i);
	for (size_t j=i+1;j<connected_live_clusters.size(); j++){
	  PR3DCluster* cluster_2 = connected_live_clusters.at(j);
	  SMGCSelection mcells_2 = connected_live_mcells.at(j);

	  if (tested_pairs.find(std::make_pair(cluster_1,cluster_2))==tested_pairs.end()){
	    cluster_1->Create_point_cloud();
	    cluster_2->Create_point_cloud();
	    //std::cout << cluster_1->get_cluster_id() << " " << cluster_2->get_cluster_id() << std::endl;
	    tested_pairs.insert(std::make_pair(cluster_1,cluster_2));
	    tested_pairs.insert(std::make_pair(cluster_2,cluster_1));
	    // starting the test ... 

	    bool flag_merge = false;
	    for (auto it1 = mcells_1.begin(); it1!= mcells_1.end(); it1++){
	      SlimMergeGeomCell* mcell_1 = (*it1);
	      Point mcell1_center = mcell_1->center();
	      // double theta, phi;
	      // std::pair<double,double> angles_1 = cluster_1->HoughTrans(mcell1_center,30*units::cm);
	      // theta = 3.1415926 - angles_1.first;
	      // phi = 3.1415926 + angles_1.second;
	      TVector3 dir1 = cluster_1->calc_dir(mcell1_center,mcell1_center,30*units::cm);
	      for (auto it2 = mcells_2.begin(); it2!=mcells_2.end(); it2++){
		SlimMergeGeomCell* mcell_2 = (*it2);
		Point mcell2_center = mcell_2->center();
		TVector3 dir2 = cluster_2->calc_dir(mcell1_center,mcell2_center,30*units::cm);
		double angle_diff = (3.1415926-dir1.Angle(dir2))/3.1415926*180.;
		double dis = sqrt(pow(mcell1_center.x-mcell2_center.x,2)
				  +pow(mcell1_center.y-mcell2_center.y,2)
				  +pow(mcell1_center.z-mcell2_center.z,2));
		// if (angle_diff < 30)
		//std::cout << "Xin: " << angle_diff << " " << dis/units::cm << std::endl;
		if (dis <= 3*units::cm && angle_diff <= 60 ||
		      dis <= 10*units::cm && angle_diff <=20 ||
		    angle_diff<10  ||
		    mcells_1.size()==1 && mcells_2.size()==1 && dis < 15*units::cm){
		  flag_merge = true;
		  break;
		}
		// std::cout << mcell1_center.x/units::cm << " " << mcell1_center.y/units::cm << " " << mcell1_center.z/units::cm << " " << mcell2_center.x/units::cm << " " << mcell2_center.y/units::cm << " " << mcell2_center.z/units::cm << std::endl;
		// if (cluster_1->get_cluster_id()==5 && cluster_2->get_cluster_id()==21)
		//   std::cout << sqrt(pow(mcell2_center.x - mcell1_center.x,2)+pow(mcell2_center.y - mcell1_center.y,2)+pow(mcell2_center.z - mcell1_center.z,2)) << " " << (mcell2_center.x - mcell1_center.x)/sqrt(pow(mcell2_center.x - mcell1_center.x,2)+pow(mcell2_center.y - mcell1_center.y,2)+pow(mcell2_center.z - mcell1_center.z,2)) << " " << (mcell2_center.y - mcell1_center.y)/sqrt(pow(mcell2_center.x - mcell1_center.x,2)+pow(mcell2_center.y - mcell1_center.y,2)+pow(mcell2_center.z - mcell1_center.z,2)) << " " << (mcell2_center.z - mcell1_center.z)/sqrt(pow(mcell2_center.x - mcell1_center.x,2)+pow(mcell2_center.y - mcell1_center.y,2)+pow(mcell2_center.z - mcell1_center.z,2))  << " " << sin(theta)*cos(phi) << " " << sin(theta)*sin(phi) << " " << cos(theta) << std::endl;
		// if (WireCell2dToy::IsCrossing(mcell1_center,theta,phi,mcell_2,gds)){
		//   flag_merge = true;
		//   break;
		// }
	      }
	      if (flag_merge) break;
	    }

	    if (!flag_merge){
	      for (auto it1 = mcells_2.begin(); it1!= mcells_2.end(); it1++){
		SlimMergeGeomCell* mcell_1 = (*it1);
		Point mcell1_center = mcell_1->center();
		//double theta, phi;
		TVector3 dir1 = cluster_1->calc_dir(mcell1_center,mcell1_center,30*units::cm);
		for (auto it2 = mcells_1.begin(); it2!=mcells_1.end(); it2++){
		  SlimMergeGeomCell* mcell_2 = (*it2);
		  Point mcell2_center = mcell_2->center();
		  TVector3 dir2 = cluster_2->calc_dir(mcell1_center,mcell2_center,30*units::cm);
		  double angle_diff = (3.1415926-dir1.Angle(dir2))/3.1415926*180.;
		  double dis = sqrt(pow(mcell1_center.x-mcell2_center.x,2)
				    +pow(mcell1_center.y-mcell2_center.y,2)
				    +pow(mcell1_center.z-mcell2_center.z,2));
		  if (dis <= 3*units::cm && angle_diff <= 60 ||
		      dis <= 10*units::cm && angle_diff <=20 ||
		      angle_diff<10  ||
		      mcells_1.size()==1 && mcells_2.size()==1 && dis < 15*units::cm ){
		    flag_merge = true;
		    break;
		  }
		  //		  if (angle_diff < 30)
		  // std::cout << "Xin: " << angle_diff << " " << dis/units::cm << std::endl;
		  
		}
		if (flag_merge) break;
	      }
	    }
	    
	    
	    if (flag_merge)
	      to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	  }
	}
      }
      //  std::cout << std::endl;
    }
  }
  
  std::vector<std::set<PR3DCluster*>> merge_clusters;
  for (auto it = to_be_merged_pairs.begin(); it!=to_be_merged_pairs.end(); it++){
    PR3DCluster *cluster1 = (*it).first;
    PR3DCluster *cluster2 = (*it).second;
    //std::cout << cluster1 << " " << cluster2 << " " << cluster1->get_cluster_id() << " " << cluster2->get_cluster_id() << std::endl;
    
    bool flag_new = true;
    std::vector<std::set<PR3DCluster*>> temp_set;
    for (auto it1 = merge_clusters.begin(); it1!=merge_clusters.end(); it1++){
      std::set<PR3DCluster*>& clusters = (*it1);
      if (clusters.find(cluster1)!=clusters.end() ||
	  clusters.find(cluster2)!=clusters.end()){
	clusters.insert(cluster1);
	clusters.insert(cluster2);
	flag_new = false;
	temp_set.push_back(clusters);
	//break;
      }
    }
    if (flag_new){
      std::set<PR3DCluster*> clusters;
      clusters.insert(cluster1);
      clusters.insert(cluster2);
      merge_clusters.push_back(clusters);
    }
    if (temp_set.size()>1){
      // merge them further ...
      for (size_t i=1;i!=temp_set.size();i++){
	for (auto it1 = temp_set.at(i).begin(); it1!= temp_set.at(i).end(); it1++){
	  temp_set.at(0).insert(*it1);
	}
	merge_clusters.erase(find(merge_clusters.begin(),merge_clusters.end(),temp_set.at(i)));
      }
    }
  }

  // merge clusters into new clusters, delete old clusters 
  for (auto it = merge_clusters.begin(); it!=merge_clusters.end();it++){
    std::set<PR3DCluster*>& clusters = (*it);
    PR3DCluster *ncluster = new PR3DCluster((*clusters.begin())->get_cluster_id());
    live_clusters.push_back(ncluster);
    for (auto it1 = clusters.begin(); it1!=clusters.end(); it1++){
      PR3DCluster *ocluster = *(it1);
      // std::cout << ocluster->get_cluster_id() << " ";
      SMGCSelection& mcells = ocluster->get_mcells();
      for (auto it2 = mcells.begin(); it2!=mcells.end(); it2++){
  	SlimMergeGeomCell *mcell = (*it2);
  	//std::cout << ocluster->get_cluster_id() << " " << mcell << std::endl;
  	int time_slice = mcell->GetTimeSlice();
  	ncluster->AddCell(mcell,time_slice);
      }
      live_clusters.erase(find(live_clusters.begin(), live_clusters.end(), ocluster));
      delete ocluster;
    }
    //std::cout << std::endl;
  }
  
  //std::cout << merge_clusters.size() << std::endl;
  
  //std::cout << to_be_merged_pairs.size() << std::endl;
  // 
  //   std::cout << "Pair: " << (*it).first->get_cluster_id() << " " << (*it).second->get_cluster_id() << std::endl;
  // }


  
  
}


// bool WireCell2dToy::IsCrossing(WireCell::Point& p, double theta, double phi, WireCell::SlimMergeGeomCell *mcell, WireCellSst::GeomDataSource& gds){
//   // dx = sintheta*cosphi*dl
//   // dy = sintheta*sinphi*dl
//   // dz = costheta*dl
//   Point pc = mcell->center();
//   Point p_test;
//   if (sin(theta)*cos(phi)!=0){
//     double dl = (pc.x - p.x)/(sin(theta)*cos(phi));
//     if (dl<0) return false;
    
//     p_test.x = pc.x;
//     p_test.y = p.y + sin(theta)*sin(phi)*dl;
//     p_test.z = p.z + cos(theta)*dl;

//     // std::cout << pc.x/units::cm << " " << p.x/units::cm << " " << dl/units::cm << " " << sqrt(pow(p_test.x - pc.x,2)+pow(p_test.y - pc.y,2)+pow(p_test.z - pc.z,2))/units::cm << std::endl;
    
//   }else{
//     p_test.x = pc.x;
//     double dl = (pc.y-p.y)*sin(theta)*sin(phi) + (pc.z-p.z)*cos(theta);
//     if (dl<0) return false;
//     p_test.y = p.y + sin(theta)*sin(phi)*dl;
//     p_test.z = p.z + cos(theta)*dl;
//   }
//   const GeomWire *uwire = gds.closest(p_test,WirePlaneType_t(0));
//   const GeomWire *vwire = gds.closest(p_test,WirePlaneType_t(1));
//   const GeomWire *wwire = gds.closest(p_test,WirePlaneType_t(2));

//   GeomWireSelection uwires = mcell->get_uwires();
//   GeomWireSelection vwires = mcell->get_vwires();
//   GeomWireSelection wwires = mcell->get_wwires();
  
//   // std::cout << uwires.front()->index() << " " << uwire->index() << " " << uwires.back()->index() << " "
//   // 	    << vwires.front()->index() << " " << vwire->index() << " " << vwires.back()->index() << " "
//   // 	    << wwires.front()->index() << " " << wwire->index() << " " << wwires.back()->index() << " " << std::endl;
    
  
  

//   if (uwire->index() >= uwires.front()->index() && uwire->index() <= uwires.back()->index() &&
//       vwire->index() >= vwires.front()->index() && vwire->index() <= vwires.back()->index() &&
//       wwire->index() >= wwires.front()->index() && wwire->index() <= wwires.back()->index()){
//     return true;
//   }else{
//     return false;
//   }
  
// }
