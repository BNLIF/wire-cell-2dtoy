#include "WireCell2dToy/ImprovePR3DCluster.h"

using namespace WireCell;

WireCell::PR3DCluster* WireCell2dToy::Improve_PR3DCluster(WireCell::PR3DCluster* cluster, ToyCTPointCloud& ct_point_cloud,WireCellSst::GeomDataSource& gds){

  std::map<int,std::set<int>> u_time_chs; // time chs
  std::map<int,std::set<int>> v_time_chs; // time chs
  std::map<int,std::set<int>> w_time_chs; // time chs
  std::map<std::pair<int,int>,double> time_ch_charge_map;
  std::map<std::pair<int,int>,double> time_ch_charge_err_map;

  // fill all map according to existing mcells
  SMGCSelection& old_mcells = cluster->get_mcells();
  for (auto it=old_mcells.begin(); it!=old_mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    GeomWireSelection& uwires = mcell->get_uwires();
    GeomWireSelection& vwires = mcell->get_vwires();
    GeomWireSelection& wwires = mcell->get_wwires();
    int time_slice = mcell->GetTimeSlice();

    if (u_time_chs.find(time_slice)==u_time_chs.end()){
      std::set<int> uchs;
      u_time_chs[time_slice] = uchs;
      std::set<int> vchs;
      v_time_chs[time_slice] = vchs;
      std::set<int> wchs;
      w_time_chs[time_slice] = wchs;
    }
    

    for (auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
      const GeomWire *wire = (*it1);
      double charge = mcell->Get_Wire_Charge(wire);
      double charge_err = mcell->Get_Wire_Charge_Err(wire);
      u_time_chs[time_slice].insert(wire->channel());
      time_ch_charge_map[std::make_pair(time_slice,wire->channel())] = charge;
      time_ch_charge_err_map[std::make_pair(time_slice,wire->channel())] = charge_err;
      // std::cout << time_slice << " " << wire->channel() << " " << charge << " " << charge_err << std::endl;
    }

    for (auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
      const GeomWire *wire = (*it1);
      double charge = mcell->Get_Wire_Charge(wire);
      double charge_err = mcell->Get_Wire_Charge_Err(wire);
      v_time_chs[time_slice].insert(wire->channel());
      time_ch_charge_map[std::make_pair(time_slice,wire->channel())] = charge;
      time_ch_charge_err_map[std::make_pair(time_slice,wire->channel())] = charge_err;
      // std::cout << time_slice << " " << wire->channel() << " " << charge << " " << charge_err << std::endl;
    }

    for (auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
      const GeomWire *wire = (*it1);
      double charge = mcell->Get_Wire_Charge(wire);
      double charge_err = mcell->Get_Wire_Charge_Err(wire);
      w_time_chs[time_slice].insert(wire->channel());
      time_ch_charge_map[std::make_pair(time_slice,wire->channel())] = charge;
      time_ch_charge_err_map[std::make_pair(time_slice,wire->channel())] = charge_err;
      // std::cout << time_slice << " " << wire->channel() << " " << charge << " " << charge_err << std::endl;
    }
  }
  //std::cout << u_time_chs.size() << " " << v_time_chs.size() << " " << w_time_chs.size() << " " << time_ch_charge_map.size() << std::endl;

  {
    // add in missing pieces based on trajectory points
    std::list<WCPointCloud<double>::WCPoint>& wcps = cluster->get_path_wcps();
    std::vector<Point> path_pts;
    double low_dis_limit = 0.3*units::cm;
    for (auto it = wcps.begin(); it!=wcps.end(); it++){
      Point p((*it).x,(*it).y,(*it).z);
      if (path_pts.size()==0){
	path_pts.push_back(p);
      }else{
	double dis = sqrt(pow(p.x-path_pts.back().x,2) +
			  pow(p.y-path_pts.back().y,2) +
			  pow(p.z-path_pts.back().z,2) );
	if (dis < low_dis_limit ){
	  path_pts.push_back(p);
	}else{
	  int ncount = int(dis/low_dis_limit)+1;
	  
	  for (int i=0; i != ncount; i++){
	    Point p1;
	    p1.x = path_pts.back().x + (p.x - path_pts.back().x) * (i+1)/ncount;
	    p1.y = path_pts.back().y + (p.y - path_pts.back().y) * (i+1)/ncount;
	    p1.z = path_pts.back().z + (p.z - path_pts.back().z) * (i+1)/ncount;
	    path_pts.push_back(p1);
	  }
	}
      }
    }

    for (auto it=path_pts.begin(); it!=path_pts.end(); it++){
      std::vector<int> results = ct_point_cloud.convert_3Dpoint_time_ch((*it));
      // +- 2 time_slices;
      // +- 2 wire channels ...
    }
    std::cout << path_pts.size() << std::endl;
    
  }
  
  
  // recreate the merged wires

  // recreate the merge cells

  // examine the newly create merged cells

  
  // create a new cluster ...
  
  
  
  //Point p(150*units::cm, 35*units::cm, 532*units::cm);
  //std::vector<int> results = ct_point_cloud.convert_3Dpoint_time_ch(p);
  //std::cout << results.at(0) << " " << results.at(1) << " " << results.at(2) << " " << results.at(3) << std::endl;
  
  return 0;
}
