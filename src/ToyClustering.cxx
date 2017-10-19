#include "WireCell2dToy/ToyClustering.h"

using namespace WireCell;

void WireCell2dToy::Clustering_live_dead(WireCell::PR3DClusterSelection& live_clusters, WireCell::PR3DClusterSelection& dead_clusters, std::map<WireCell::PR3DCluster*,std::vector<WireCell::PR3DCluster*>>& dead_live_cluster_mapping, std::map<WireCell::PR3DCluster*,std::vector<WireCell::SMGCSelection>>& dead_live_mcells_mapping, WireCellSst::GeomDataSource& gds){
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

  //std::cout << to_be_merged_pairs.size() << std::endl;
  // for (auto it = to_be_merged_pairs.begin(); it!=to_be_merged_pairs.end(); it++){
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
