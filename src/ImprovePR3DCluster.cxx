#include "WCP2dToy/ImprovePR3DCluster.h"
#include "WCP2dToy/LowmemTiling.h"


using namespace WCP;

WCP::PR3DCluster* WCP2dToy::Improve_PR3DCluster(WCP::PR3DCluster* cluster, ToyCTPointCloud& ct_point_cloud,WCPSst::GeomDataSource& gds){

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

  // for (auto it = time_ch_charge_map.begin(); it!= time_ch_charge_map.end(); it++){
  //   std::cout << it->first.first << " " << it->first.second << std::endl;
  // }
  
  
  //  std::cout << u_time_chs.size() << " Xin1: " << v_time_chs.size() << " " << w_time_chs.size() << " " << time_ch_charge_map.size() << std::endl;

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

    std::pair<int,int> uch_limits = ct_point_cloud.get_uch_limits();
    std::pair<int,int> vch_limits = ct_point_cloud.get_vch_limits();
    std::pair<int,int> wch_limits = ct_point_cloud.get_wch_limits();

    std::vector<bool> path_pts_flag;
    for (auto it=path_pts.begin(); it!=path_pts.end(); it++){
      std::vector<int> results = ct_point_cloud.convert_3Dpoint_time_ch((*it));

      int nu = 0;
      int nv = 0;
      int nw = 0;

      if (time_ch_charge_map.find(std::make_pair(results.at(0),results.at(1)))!=time_ch_charge_map.end()) nu ++;
      if (time_ch_charge_map.find(std::make_pair(results.at(0),results.at(1)-1))!=time_ch_charge_map.end()) nu +=2;
      if (time_ch_charge_map.find(std::make_pair(results.at(0),results.at(1)+1))!=time_ch_charge_map.end()) nu ++;

      if (time_ch_charge_map.find(std::make_pair(results.at(0),results.at(2)))!=time_ch_charge_map.end()) nv ++;
      if (time_ch_charge_map.find(std::make_pair(results.at(0),results.at(2)-1))!=time_ch_charge_map.end()) nv+=2;
      if (time_ch_charge_map.find(std::make_pair(results.at(0),results.at(2)+1))!=time_ch_charge_map.end()) nv++;
      
      if (time_ch_charge_map.find(std::make_pair(results.at(0),results.at(3)))!=time_ch_charge_map.end()) nw++;
      if (time_ch_charge_map.find(std::make_pair(results.at(0),results.at(3)-1))!=time_ch_charge_map.end()) nw += 2;
      if (time_ch_charge_map.find(std::make_pair(results.at(0),results.at(3)+1))!=time_ch_charge_map.end()) nw ++;

      if (nu>0&&nv>0&&nw>0&&nu+nv+nw>=6){
	path_pts_flag.push_back(true);
      }else{
	path_pts_flag.push_back(false);
      }
      
      //  std::cout << "Path: " << (*it).x/units::cm << " " << (*it).y/units::cm << " " << (*it).z/units::cm << " " << path_pts_flag.back() << std::endl;
    }

    // std::cout << path_pts.size() << " " << path_pts_flag.size() << " " << cluster->get_cluster_id() << std::endl;
    
    for (size_t i=0;i!=path_pts.size();i++){
      std::vector<int> results = ct_point_cloud.convert_3Dpoint_time_ch(path_pts.at(i));

      if (i==0){
	if (path_pts_flag.at(i) && path_pts_flag.at(i+1)) continue;
      }else if (i+1==path_pts.size()){
	if (path_pts_flag.at(i) && path_pts_flag.at(i-1)) continue;
      }else{
	if (path_pts_flag.at(i-1) && path_pts_flag.at(i) && path_pts_flag.at(i+1)) continue;
      }
      
    
      
      for (int time_slice = results.at(0)-3; time_slice<=results.at(0)+3; time_slice ++){
	if (time_slice <0) continue;
	
	if (u_time_chs.find(time_slice)==u_time_chs.end()){
	  std::set<int> uchs;
	  u_time_chs[time_slice] = uchs;
	  std::set<int> vchs;
	  v_time_chs[time_slice] = vchs;
	  std::set<int> wchs;
	  w_time_chs[time_slice] = wchs;
	}
	for (int ch = results.at(1)-3; ch<=results.at(1)+3; ch++){
	  if (ch < uch_limits.first || ch > uch_limits.second ||
	      fabs(ch - results.at(1)) + fabs(time_slice -results.at(0))>3)
	    continue;
	  u_time_chs[time_slice].insert(ch);
	  if (time_ch_charge_map.find(std::make_pair(time_slice,ch))==time_ch_charge_map.end()){
	    time_ch_charge_map[std::make_pair(time_slice,ch)]=0;
	    time_ch_charge_err_map[std::make_pair(time_slice,ch)]=0;
	    //std::cout << time_slice << " " << ch << std::endl;
	  }
	}
	for (int ch = results.at(2)-3; ch<=results.at(2)+3; ch++){
	  if (ch < vch_limits.first || ch > vch_limits.second ||
	      fabs(ch - results.at(2)) + fabs(time_slice -results.at(0))>3)
	    continue;
	  v_time_chs[time_slice].insert(ch);
	  if (time_ch_charge_map.find(std::make_pair(time_slice,ch))==time_ch_charge_map.end()){
	    time_ch_charge_map[std::make_pair(time_slice,ch)]=0;
	    time_ch_charge_err_map[std::make_pair(time_slice,ch)]=0;
	    //std::cout << time_slice << " " << ch << std::endl;
	  }
	}
	for (int ch = results.at(3)-3; ch<=results.at(3)+3; ch++){
	  if (ch < wch_limits.first || ch > wch_limits.second ||
	      fabs(ch - results.at(3)) + fabs(time_slice -results.at(0))>3)
	    continue;
	  w_time_chs[time_slice].insert(ch);
	  if (time_ch_charge_map.find(std::make_pair(time_slice,ch))==time_ch_charge_map.end()){
	    time_ch_charge_map[std::make_pair(time_slice,ch)]=0;
	    time_ch_charge_err_map[std::make_pair(time_slice,ch)]=0;
	    //std::cout << time_slice << " " << ch << std::endl;
	  }
	}	
      }
    }
    //std::cout << path_pts.size() << std::endl;
    //   std::cout << u_time_chs.size() << " Xin2: " << v_time_chs.size() << " " << w_time_chs.size() << " " << time_ch_charge_map.size() << std::endl;
  }

  
  WCP2dToy::WCPHolder *WCholder = new WCPHolder();
  for (auto it = u_time_chs.begin(); it!= u_time_chs.end(); it++){
    int time_slice = it->first;
    WCP2dToy::LowmemTiling tiling(time_slice,gds,*WCholder);
    // recreate the merged wires
    // recreate the merge cells
    tiling.init_good_cells(u_time_chs, v_time_chs, w_time_chs);  
  }
  
  
  // examine the newly create merged cells
  std::map<int,SMGCSelection> old_time_mcells_map;
  for (auto it = old_mcells.begin(); it!=old_mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    int time_slice = mcell->GetTimeSlice();
    if (old_time_mcells_map.find(time_slice)==old_time_mcells_map.end()){
      SMGCSelection mcells;
      mcells.push_back(mcell);
      old_time_mcells_map[time_slice] = mcells;
    }else{
      old_time_mcells_map[time_slice].push_back(mcell);
    }
  }
  
  std::map<int,SMGCSelection> new_time_mcells_map;
  SMGCSet new_mcells_set;
  
  GeomCellSelection& temp_cells = WCholder->get_cells();
  for (auto it = temp_cells.begin(); it!=temp_cells.end(); it++){
    SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);
    int time_slice = mcell->GetTimeSlice();
    bool flag_good = false;
    
    if (!flag_good){
      // -1 time slice
      if (old_time_mcells_map.find(time_slice-1) != old_time_mcells_map.end()){
	for (auto it1 = old_time_mcells_map[time_slice-1].begin(); it1!=old_time_mcells_map[time_slice-1].end(); it1++){
	  if (mcell->Overlap_fast((*it1))){
	    flag_good = true;
	    break;
	  }
	}
      }
    }
    if (!flag_good){
      // same time_slice
      if (old_time_mcells_map.find(time_slice) != old_time_mcells_map.end()){
	for (auto it1 = old_time_mcells_map[time_slice].begin(); it1!=old_time_mcells_map[time_slice].end(); it1++){
	  if (mcell->Overlap_fast((*it1))){
	    flag_good = true;
	    break;
	  }
	}
      }
    }
    if (!flag_good){
      // +1 time_slice
      if (old_time_mcells_map.find(time_slice+1) != old_time_mcells_map.end()){
	for (auto it1 = old_time_mcells_map[time_slice+1].begin(); it1!=old_time_mcells_map[time_slice+1].end(); it1++){
	  if (mcell->Overlap_fast((*it1))){
	    flag_good = true;
	    break;
	  }
	}
      }
    }
    

    if (flag_good){
      if (new_time_mcells_map.find(time_slice)==new_time_mcells_map.end()){
	SMGCSelection mcells;
	mcells.push_back(mcell);
	new_time_mcells_map[time_slice] = mcells;
      }else{
	new_time_mcells_map[time_slice].push_back(mcell);
      }
    }else{
      new_mcells_set.insert(mcell);
    }
    
  }
  
  int prev_num_mcells = new_mcells_set.size()+1;
  while(new_mcells_set.size() != prev_num_mcells && new_mcells_set.size() > 0){
    prev_num_mcells = new_mcells_set.size();

    SMGCSelection temp_mcells;
    
    //start to work on it ...
    for (auto it= new_mcells_set.begin(); it!=new_mcells_set.end(); it++){
      SlimMergeGeomCell *mcell = (*it);
      int time_slice = mcell->GetTimeSlice();
      bool flag_good = false;

      if (!flag_good){
	// -1 time slice
	if (new_time_mcells_map.find(time_slice-1) != new_time_mcells_map.end()){
	  for (auto it1 = new_time_mcells_map[time_slice-1].begin(); it1!=new_time_mcells_map[time_slice-1].end(); it1++){
	    if (mcell->Overlap_fast((*it1))){
	      flag_good = true;
	      break;
	    }
	  }
	}
      }
      if (!flag_good){
	// same time_slice
	if (new_time_mcells_map.find(time_slice) != new_time_mcells_map.end()){
	  for (auto it1 = new_time_mcells_map[time_slice].begin(); it1!=new_time_mcells_map[time_slice].end(); it1++){
	    if (mcell->Overlap_fast((*it1))){
	      flag_good = true;
	      break;
	    }
	  }
	}
      }
      if (!flag_good){
	// +1 time_slice
	if (new_time_mcells_map.find(time_slice+1) != new_time_mcells_map.end()){
	  for (auto it1 = new_time_mcells_map[time_slice+1].begin(); it1!=new_time_mcells_map[time_slice+1].end(); it1++){
	    if (mcell->Overlap_fast((*it1))){
	      flag_good = true;
	      break;
	    }
	  }
	}
      }

      if (flag_good){
	temp_mcells.push_back(mcell);
      }
    }

    for (auto it = temp_mcells.begin(); it!=temp_mcells.end(); it++){
      SlimMergeGeomCell *mcell = (*it);
      int time_slice = mcell->GetTimeSlice();

      if (new_time_mcells_map.find(time_slice)==new_time_mcells_map.end()){
	SMGCSelection mcells;
	mcells.push_back(mcell);
	new_time_mcells_map[time_slice] = mcells;
      }else{
	new_time_mcells_map[time_slice].push_back(mcell);
      }
      new_mcells_set.erase(mcell);
    }
    //  std::cout << new_mcells_set.size() << std::endl;
  }

  for (auto it=new_mcells_set.begin(); it!=new_mcells_set.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    delete mcell;
  }
  new_mcells_set.clear();

  // create a new cluster ...
  PR3DCluster *new_cluster = new PR3DCluster(cluster->get_cluster_id());
  for (auto it = new_time_mcells_map.begin(); it!= new_time_mcells_map.end(); it++){
    int time_slice = it->first;
    SMGCSelection& temp_mcells = it->second;
    for (auto it1 = temp_mcells.begin(); it1!=temp_mcells.end(); it1++){
      SlimMergeGeomCell *mcell = (*it1);
      new_cluster->AddCell(mcell,time_slice);
    }
  }
  
  //std::cout << cluster->get_cluster_id() << " " << old_mcells.size() << " " << u_time_chs.size() << " " << WCholder->get_ncell() << " " << WCholder->get_nwire() << " " << new_mcells_set.size() << std::endl;
  
  
  //Point p(150*units::cm, 35*units::cm, 532*units::cm);
  //std::vector<int> results = ct_point_cloud.convert_3Dpoint_time_ch(p);
  //std::cout << results.at(0) << " " << results.at(1) << " " << results.at(2) << " " << results.at(3) << std::endl;
  
  return new_cluster;
}
