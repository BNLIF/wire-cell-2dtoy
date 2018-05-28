#include "WireCell2dToy/ToyFiducial.h"

using namespace WireCell;

int pnpoly(std::vector<double>& vertx, std::vector<double>& verty, double testx, double testy)
{
  int i, j, c = 0;
  for (i = 0, j = int(vertx.size())-1; i < int(vertx.size()); j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
      c = !c;
  }
  return c;
}


WireCell2dToy::ToyFiducial::ToyFiducial(int dead_region_ch_ext, double offset_t, double offset_u, double offset_v, double offset_w, double slope_t, double slope_u, double slope_v, double slope_w, double angle_u, double angle_v, double angle_w, double boundary_dis_cut, double top, double bottom, double upstream, double downstream, double anode, double cathode)
  : dead_region_ch_ext(dead_region_ch_ext)
  , offset_t(offset_t)
  , offset_u(offset_u)
  , offset_v(offset_v)
  , offset_w(offset_w)
  , slope_t(slope_t)
  , slope_u(slope_u)
  , slope_v(slope_v)
  , slope_w(slope_w)
  , angle_u(angle_u)
  , angle_v(angle_v)
  , angle_w(angle_w)
  , m_top(top)
  , m_bottom(bottom)
  , m_upstream(upstream)
  , m_downstream(downstream)
  , m_anode(anode)
  , m_cathode(cathode)
{
  m_sc_bottom_1_y=-116*units::cm;
  m_sc_bottom_1_x=80*units::cm;

  m_sc_bottom_2_y=-99*units::cm;
  m_sc_bottom_2_x=256*units::cm;

  m_sc_top_1_y = 118*units::cm;
  m_sc_top_1_x = 100*units::cm;

  m_sc_top_2_y = 103*units::cm;
  m_sc_top_2_x = 256*units::cm;

  m_sc_upstream_1_z = 0*units::cm;
  m_sc_upstream_1_x = 120*units::cm;

  m_sc_upstream_2_z = 11*units::cm;
  m_sc_upstream_2_x = 256*units::cm;

  m_sc_downstream_1_z=1037*units::cm;
  m_sc_downstream_1_x=120*units::cm;

  m_sc_downstream_2_z=1026*units::cm;
  m_sc_downstream_2_x=256*units::cm;

  //
  boundary_xy_x.clear(); boundary_xy_y.clear();
  boundary_xy_x.push_back(m_anode + boundary_dis_cut); boundary_xy_y.push_back(m_bottom + boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_bottom_1_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_bottom_1_y+boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_bottom_2_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_bottom_2_y+boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_top_2_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_top_2_y-boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_top_1_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_top_1_y-boundary_dis_cut);
  boundary_xy_x.push_back(m_anode + boundary_dis_cut); boundary_xy_y.push_back(m_top - boundary_dis_cut);
  // boundary_xy_x.push_back(m_anode + boundary_dis_cut); boundary_xy_y.push_back(m_bottom + boundary_dis_cut);
  
  // for (size_t i=0;i!=boundary_xy_x.size();i++){
  //   std::cout << boundary_xy_x.at(i)/units::cm << " XY " << boundary_xy_y.at(i)/units::cm << std::endl;
  // }
  
  boundary_xz_x.clear(); boundary_xz_z.clear();
  boundary_xz_x.push_back(m_anode + boundary_dis_cut); boundary_xz_z.push_back(m_upstream + boundary_dis_cut+1*units::cm);
  boundary_xz_x.push_back(m_sc_upstream_1_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_upstream_1_z + boundary_dis_cut+1*units::cm);
  boundary_xz_x.push_back(m_sc_upstream_2_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_upstream_2_z + boundary_dis_cut+1*units::cm);
  boundary_xz_x.push_back(m_sc_downstream_2_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_downstream_2_z - boundary_dis_cut-1*units::cm);
  boundary_xz_x.push_back(m_sc_downstream_1_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_downstream_1_z - boundary_dis_cut-1*units::cm);
  boundary_xz_x.push_back(m_anode + boundary_dis_cut); boundary_xz_z.push_back(m_downstream - boundary_dis_cut-1*units::cm);
  // boundary_xz_x.push_back(m_anode + boundary_dis_cut); boundary_xz_z.push_back(m_upstream + boundary_dis_cut+2*units::cm);

  // for (size_t i=0;i!=boundary_xz_x.size();i++){
  //   std::cout << boundary_xz_x.at(i)/units::cm << " XZ " << boundary_xz_z.at(i)/units::cm << std::endl;
  // }
}

bool WireCell2dToy::ToyFiducial::check_tgm(WireCell::FlashTPCBundle *bundle, double offset_x, WireCell::ToyCTPointCloud& ct_point_cloud){

  PR3DCluster *main_cluster = bundle->get_main_cluster();
  Opflash *flash = bundle->get_flash();

  std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = main_cluster->get_extreme_wcps();

  // int max_group = 0;
  // int max_count = out_vec_wcps.at(max_group).size();

  // for (size_t i=1; i!=out_vec_wcps.size();i++){
  //   if (out_vec_wcps.at(i).size() > max_count){
  //     max_group = i;
  //     max_count = out_vec_wcps.at(max_group).size();
  //   }
  // }

  TVector3 drift_dir(1,0,0);
  // hard coded for U and V plane ... 
  TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
  TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
  TVector3 W_dir(0,1,0);

  
  //  std::cout << "Flash: " << flash->get_flash_id() << std::endl;

  // take a look at the first point ...
  for (size_t i=0;i!=out_vec_wcps.size();i++){
    bool flag_p1_inside = true;
    int p1_index = -1;
    for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
      Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
      flag_p1_inside = flag_p1_inside && inside_fiducial_volume(p1,offset_x);
      if (!flag_p1_inside){
	p1_index = j;
	break;
      }
    }
    
    
    // loop through the remaining groups and check ...
    for (size_t k=i+1;k!=out_vec_wcps.size();k++){
      bool flag_p2_inside = true;
      int p2_index = -1;
      for(size_t j=0;j!=out_vec_wcps.at(k).size();j++){
	Point p2(out_vec_wcps.at(k).at(j).x,out_vec_wcps.at(k).at(j).y,out_vec_wcps.at(k).at(j).z);
	flag_p2_inside = flag_p2_inside && inside_fiducial_volume(p2,offset_x);
	if (!flag_p2_inside){
	  p2_index = j;
	  break;
	}
      }
      

      // if (main_cluster->get_cluster_id()==13){
       
      // }
      //  std::cout << main_cluster->get_cluster_id() << " " << i << " " << k << " " << p1_index << " " << p2_index << " " << flag_p1_inside << " " << flag_p2_inside << std::endl;
	
      if ((!flag_p1_inside) && (!flag_p2_inside)){
	// if not a neutrino candidate ... to be worked out ...
	 // Point p1(out_vec_wcps.at(i).at(p1_index).x,out_vec_wcps.at(i).at(p1_index).y,out_vec_wcps.at(i).at(p1_index).z);
	 // Point p2(out_vec_wcps.at(k).at(p2_index).x,out_vec_wcps.at(k).at(p2_index).y,out_vec_wcps.at(k).at(p2_index).z);
	 // std::cout << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << inside_fiducial_volume(p1,offset_x) << " A " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << inside_fiducial_volume(p2,offset_x) << " " << offset_x/units::cm << std::endl;


	//	std::cout << main_cluster->get_cluster_id() << " " << (out_vec_wcps.at(i).at(p1_index).x-offset_x)/units::cm << " " << out_vec_wcps.at(i).at(p1_index).y/units::cm << " " << out_vec_wcps.at(i).at(p1_index).z/units::cm << " " ;
	//std::cout << (out_vec_wcps.at(k).at(p2_index).x-offset_x)/units::cm << " " << out_vec_wcps.at(k).at(p2_index).y/units::cm << " " << out_vec_wcps.at(k).at(p2_index).z/units::cm << " " <<  flag_p1_inside << " " << flag_p2_inside << " " << out_vec_wcps.size() << std::endl;

	//std::cout << p1_index << " " << p2_index << std::endl;
	
	// check two points in between
	bool flag_check = false;
	for (int kk=0;kk!=3;kk++){
	  Point p3(out_vec_wcps.at(i).at(p1_index).x+ (out_vec_wcps.at(k).at(p2_index).x - out_vec_wcps.at(i).at(p1_index).x)/4.*(kk+1),
		   out_vec_wcps.at(i).at(p1_index).y+ (out_vec_wcps.at(k).at(p2_index).y - out_vec_wcps.at(i).at(p1_index).y)/4.*(kk+1),
		   out_vec_wcps.at(i).at(p1_index).z+ (out_vec_wcps.at(k).at(p2_index).z - out_vec_wcps.at(i).at(p1_index).z)/4.*(kk+1));
	  flag_check = flag_check || inside_fiducial_volume(p3,offset_x);
	}
	//std::cout << flag_check << " " << out_vec_wcps.size() << std::endl;
	
	if (flag_check){
	  if (flash->get_type()==2){
	    if (!check_neutrino_candidate(main_cluster,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x))
	      return true;
	  }else{
	    return true; // through going muon ...
	  }
	}else{
	  if (out_vec_wcps.size()==2)
	    return true;
	}
      }else{
	Vector main_dir = main_cluster->get_PCA_axis(0);
	TVector3 dir_main(main_dir.x,main_dir.y,main_dir.z);
	TVector3 dir_test(out_vec_wcps.at(i).at(0).x-out_vec_wcps.at(k).at(0).x,
			  out_vec_wcps.at(i).at(0).y-out_vec_wcps.at(k).at(0).y,
			  out_vec_wcps.at(i).at(0).z-out_vec_wcps.at(k).at(0).z);

	//	std::cout << main_cluster->get_cluster_id() << " " << fabs((3.1415926/2.-dir_test.Angle(dir_main))/3.1415926*180.) << std::endl;

	
	
	if (fabs((3.1415926/2.-dir_test.Angle(dir_main))/3.1415926*180.)>75 || i==0 && k==1){
	  // check dead region ...
	  bool flag_p1_inside_p = flag_p1_inside;
	  if (flag_p1_inside_p){
	    Point p1(out_vec_wcps.at(i).at(0).x,out_vec_wcps.at(i).at(0).y,out_vec_wcps.at(i).at(0).z);
	    TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
	    dir *= (-1);
	    
	    //	    std::cout << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << dir.X() << " " << dir.Y() << dir.Z() << std::endl;

	    // check U and V and W
	    TVector3 dir_1(0,dir.Y(),dir.Z());
	    double angle1 = dir_1.Angle(U_dir);
	    TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
	    double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;

	    double angle2 = dir_1.Angle(V_dir);
	    TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
	    double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;

	    double angle3 = dir_1.Angle(W_dir);
	    TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
	    double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;

	    //   std::cout << "A: " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << angle1_1 << " " << angle2_1 << " " << angle3_1 << std::endl;

	    if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
	      flag_p1_inside_p = flag_p1_inside_p && check_signal_processing(p1,dir,ct_point_cloud,1*units::cm,offset_x);
	    }
	    
	    if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 )
	      flag_p1_inside_p= flag_p1_inside_p && check_dead_volume(p1,dir,1*units::cm,offset_x);
	  }
	  
	  bool flag_p2_inside_p = flag_p2_inside;
	  if (flag_p2_inside_p){
	    Point p2(out_vec_wcps.at(k).at(0).x,out_vec_wcps.at(k).at(0).y,out_vec_wcps.at(k).at(0).z);
	    TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
	    dir *= (-1);
	    
	    //	    std::cout << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << dir.X() << " " << dir.Y() << dir.Z() << std::endl;
	    TVector3 dir_1(0,dir.Y(),dir.Z());
	    double angle1 = dir_1.Angle(U_dir);
	    TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
	    double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;

	    double angle2 = dir_1.Angle(V_dir);
	    TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
	    double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;

	    double angle3 = dir_1.Angle(W_dir);
	    TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
	    double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;

	    //	    std::cout << "B: " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " <<  angle1_1 << " " << angle2_1 << " " << angle3_1 << std::endl;

	    if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
	      flag_p2_inside_p = flag_p2_inside_p && check_signal_processing(p2,dir,ct_point_cloud,1*units::cm,offset_x);
	    }
	    
	    if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 )
	      flag_p2_inside_p=flag_p2_inside_p && check_dead_volume(p2,dir,1*units::cm,offset_x);
	  }
	  
	  if ((!flag_p1_inside_p) && (!flag_p2_inside_p)){
	    if (flash->get_type()==2){
	      if (!check_neutrino_candidate(main_cluster,out_vec_wcps.at(i).at(0),out_vec_wcps.at(k).at(0),offset_x))
		return true;
	    }else{
	      return true;
	    }
	  }
	  // check signal processing ...

	// {
	// 	if (flag_p1_inside)
	// 	  ;
	
	// 	if (flag_p2_inside)
	// 	  ;
	
	// }

	}
	
	
	
      }
      
    }
  }
  

  return false;

  // also check against the dead channel ...  
  // // check the fiducial volume ...
  // std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = main_cluster->get_main_axis_wcps();
  // Point p1(wcps.first.x,wcps.first.y,wcps.first.z);
  // Point p2(wcps.second.x,wcps.second.y,wcps.second.z);
  // //double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
  // bool flag_inside_p1 = inside_fiducial_volume(p1,offset_x);
  // bool flag_inside_p2 = inside_fiducial_volume(p2,offset_x);
  // //std::cout << main_cluster->get_cluster_id() << " " << (p1.x-offset_x)/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << (p2.x-offset_x)/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << fid->inside_fiducial_volume(p1,offset_x) << " " << fid->inside_fiducial_volume(p2,offset_x) << std::endl;

  
  
  // // check the dead region ... 
  // if (flag_inside_p1){
  //   // define a local direction ...
  //   TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
  //   dir *= (-1);
  //   flag_inside_p1=check_dead_volume(p1,dir,1*units::cm,offset_x);
  // }
  // if (flag_inside_p2){
  //   // define a  local direction ...
  //   TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
  //   dir *= (-1);
  //   flag_inside_p2=check_dead_volume(p2,dir,1*units::cm,offset_x);
  // }
  // return (!flag_inside_p1)&&(!flag_inside_p2);


  // bool flag_2nd = true;
       // {
	 
	 
       // 	 if ((!flag_inside_p1)&&(!flag_inside_p2)){
       // 	   event_type |= 1UL << 3; // through going muon ... 
       // 	   flag_2nd = false;
       // 	 }
	 
       // 	 if (flag_2nd && ((!flag_inside_p1)|| (!flag_inside_p2) )){
       // 	   // check the fiducial volume ...
       // 	   std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = main_cluster->get_extreme_wcps();
	   
       // 	   Point p1(wcps.first.x,wcps.first.y,wcps.first.z);
       // 	   Point p2(wcps.second.x,wcps.second.y,wcps.second.z);
	   
       // 	   flag_inside_p1 = fid->inside_fiducial_volume(p1,offset_x);
       // 	   flag_inside_p2 = fid->inside_fiducial_volume(p2,offset_x);

       // 	   std::cout << main_cluster->get_cluster_id() << " " << (p1.x-offset_x)/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << (p2.x-offset_x)/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << fid->inside_fiducial_volume(p1,offset_x) << " " << fid->inside_fiducial_volume(p2,offset_x) << std::endl;
	   
       // 	   // check the dead region ...
       // 	   if (flag_inside_p1){
       // 	     // define a local direction ...
       // 	     TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
       // 	     dir *= (-1);
       // 	     flag_inside_p1=fid->check_dead_volume(p1,dir,1*units::cm,offset_x);
       // 	   }
       // 	   if (flag_inside_p2){
       // 	     // define a  local direction ...
       // 	     TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
       // 	     dir *= (-1);
       // 	     flag_inside_p2=fid->check_dead_volume(p2,dir,1*units::cm,offset_x);
       // 	   }
	   
       // 	   if ((!flag_inside_p1)&&(!flag_inside_p2)){
       // 	     event_type |= 1UL << 3; // through going muon ... 
       // 	   }
       // 	 }

  
  return false;
}

bool WireCell2dToy::ToyFiducial::check_neutrino_candidate(WireCell::PR3DCluster *main_cluster,WCPointCloud<double>::WCPoint& wcp1 ,WCPointCloud<double>::WCPoint& wcp2, double offset_x){
  main_cluster->Create_graph();
  main_cluster->dijkstra_shortest_paths(wcp1);
  main_cluster->cal_shortest_path(wcp2);

  std::list<WCPointCloud<double>::WCPoint>& path_wcps = main_cluster->get_path_wcps();

  PointVector path_wcps_vec;  
  // if (fine_tracking_path.size()==0){
  double low_dis_limit = 0.5*units::cm;
  for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
    if (path_wcps_vec.size()==0){
      Point p((*it).x,(*it).y,(*it).z);
      path_wcps_vec.push_back(p);
    }else{
      double dis = sqrt(pow((*it).x - path_wcps_vec.back().x,2)
			+pow((*it).y - path_wcps_vec.back().y,2)
			+pow((*it).z - path_wcps_vec.back().z,2));
      if (dis > low_dis_limit){
	Point p((*it).x,(*it).y,(*it).z);
	path_wcps_vec.push_back(p);
      }
    }
  }    
  // }else{
  //   path_wcps_vec = fine_tracking_path;
  // }
  
  
  // if (cluster_id == 13){
  //   std::cout << wcp1.x/units::cm << " " << wcp1.y/units::cm << " " << wcp1.z/units::cm << " " << wcp2.x/units::cm << " " << wcp2.y/units::cm << " " << wcp2.z/units::cm << std::endl;
  int count = 0;
  double max_angle=0;
  Point max_point(0,0,0);
  TVector3 drift_dir(1,0,0);
  for (size_t i=5;i+5<path_wcps_vec.size();i++){
    TVector3 dir1(path_wcps_vec.at(i).x - path_wcps_vec.at(i-5).x,
		  path_wcps_vec.at(i).y - path_wcps_vec.at(i-5).y,
		  path_wcps_vec.at(i).z - path_wcps_vec.at(i-5).z);
    TVector3 dir2(path_wcps_vec.at(i).x - path_wcps_vec.at(i+5).x,
		  path_wcps_vec.at(i).y - path_wcps_vec.at(i+5).y,
		  path_wcps_vec.at(i).z - path_wcps_vec.at(i+5).z);
    
    TVector3 dir3, dir4, dir5, dir6;
    {
      PointVector pts;
      double temp_x = 0;
      double temp_y = 0;
      double temp_z = 0;
      double temp_count = 0;
      for (size_t j=1;j!=15;j++){
       	if (i>=j){
	  Point pt(path_wcps_vec.at(i-j).x,path_wcps_vec.at(i-j).y,path_wcps_vec.at(i-j).z);
	  
	  if (j<=12&&j>2){
	    temp_x += pt.x;
	    temp_y += pt.y;
	    temp_z += pt.z;
	    temp_count ++;
	  }
	  pts.push_back(pt);
	}
	Point pt(path_wcps_vec.at(i).x,path_wcps_vec.at(i).y,path_wcps_vec.at(i).z);
	dir3 = main_cluster->calc_PCA_dir(pt,pts);
	dir5.SetXYZ(temp_x/temp_count - path_wcps_vec.at(i).x,
		    temp_y/temp_count - path_wcps_vec.at(i).y,
		    temp_z/temp_count - path_wcps_vec.at(i).z);
	if (dir3.Angle(dir1)>3.1415926/2.)
	  dir3 *= -1;
      }
    }
    {
      PointVector pts;
      double temp_x = 0;
      double temp_y = 0;
      double temp_z = 0;
      double temp_count = 0;
      for (size_t j=1;j!=15;j++){
	if (i+j<path_wcps_vec.size()){
	  Point pt(path_wcps_vec.at(i+j).x,path_wcps_vec.at(i+j).y,path_wcps_vec.at(i+j).z);
	  if (j<=12&&j>2){
	    temp_x += pt.x;
	    temp_y += pt.y;
	    temp_z += pt.z;
	    temp_count ++;
	  }
	  pts.push_back(pt);
	}
      }
      Point pt(path_wcps_vec.at(i).x,path_wcps_vec.at(i).y,path_wcps_vec.at(i).z);
      dir4 = main_cluster->calc_PCA_dir(pt,pts);
      dir6.SetXYZ(temp_x/temp_count - path_wcps_vec.at(i).x,
      		  temp_y/temp_count - path_wcps_vec.at(i).y,
      		  temp_z/temp_count - path_wcps_vec.at(i).z);
      if (dir4.Angle(dir2)>3.1415926/2.)
	dir4 *= -1;
    }

    int cut1 = 0;
    if ((3.1415926 - dir1.Angle(dir2))/3.1415926*180.>30) cut1++;
    if ((3.1415926 - dir3.Angle(dir4))/3.1415926*180.>30) cut1++;
    if ((3.1415926 - dir5.Angle(dir6))/3.1415926*180.>30) cut1++;
    int cut2 = 0;
    if (fabs(3.1415926/2.-drift_dir.Angle(dir1-dir2))/3.1415926*180. > 5) cut2++;
    if (fabs(3.1415926/2.-drift_dir.Angle(dir3-dir4))/3.1415926*180. > 5) cut2++;
    if (fabs(3.1415926/2.-drift_dir.Angle(dir5-dir6))/3.1415926*180. > 5) cut2++;
    
    // std::cout << i << " " << path_wcps_vec.at(i).x/units::cm << " " << path_wcps_vec.at(i).y/units::cm << " " << path_wcps_vec.at(i).z/units::cm << " " << (3.1415926 - dir1.Angle(dir2))/3.1415926*180. << " " << (3.1415926 - dir3.Angle(dir4))/3.1415926*180. << " " << (3.1415926 - dir5.Angle(dir6))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(dir1-dir2))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(dir3-dir4))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(dir5-dir6))/3.1415926*180. << " " << cut1 << " " << cut2 << std::endl;

   
    
    
    if (cut1>=3 && cut2>=2){
      if ((3.1415926 - dir3.Angle(dir4))/3.1415926*180. > max_angle){
	max_angle = (3.1415926 - dir3.Angle(dir4))/3.1415926*180.;
	max_point = path_wcps_vec.at(i);
      }
      
      count ++;
      if (count >=3){
	TVector3 temp1(path_wcps_vec.at(i).x-wcp1.x,
		       path_wcps_vec.at(i).y-wcp1.y,
		       path_wcps_vec.at(i).z-wcp1.z);
	TVector3 temp2(path_wcps_vec.at(i).x-wcp2.x,
		       path_wcps_vec.at(i).y-wcp2.y,
		       path_wcps_vec.at(i).z-wcp2.z);

	//	std::cout << "A: " << (3.1415926-temp1.Angle(temp2))/3.1415926*180. << " " << temp1.Mag()/units::cm << " " << temp2.Mag()/units::cm << std::endl;


	if ((3.1415926-temp1.Angle(temp2))/3.1415926*180. >30 ||
	    (3.1415926-temp1.Angle(temp2))/3.1415926*180. >25 && temp1.Mag()>15*units::cm && temp2.Mag()>15*units::cm){

	  if ((!inside_fiducial_volume(max_point,offset_x)) || // must be in fiducial
	      inside_dead_region(max_point)&&(3.1415926-temp1.Angle(temp2))/3.1415926*180<45 // not in dead_volume
	      ){ // should not too close to anode 
	  }else{
	    return true;
	  }
	}
	
      }
    }else{
      count = 0 ;
      max_angle = 0;
      max_point.x = 0;
      max_point.y = 0;
      max_point.z = 0;
    }
  }
  
  return false;
}


WireCell2dToy::ToyFiducial::~ToyFiducial(){
}

bool WireCell2dToy::ToyFiducial::check_signal_processing(WireCell::Point& p, TVector3& dir, WireCell::ToyCTPointCloud& ct_point_cloud, double step, double offset_x){

  if (dir.Mag()==0){
    return true;
  }else{
    dir *= 1./dir.Mag();
    Point temp_p = p;

    int num_points = 0;
    int num_points_dead = 0;

    //  std::cerr << temp_p.x/units::cm << " " << temp_p.y/units::cm << " " << temp_p.z/units::cm << std::endl;
    
    while(inside_fiducial_volume(temp_p,offset_x)){
      num_points ++;
      //if (inside_dead_region(temp_p))
      //	num_points_dead ++;

      //      std::cerr << temp_p.x/units::cm << " " << temp_p.y/units::cm << " " << temp_p.z/units::cm << " ";
      
      WireCell::CTPointCloud<double> cloud_u = ct_point_cloud.get_closest_points(temp_p,1.2*units::cm,0);
      WireCell::CTPointCloud<double> cloud_v = ct_point_cloud.get_closest_points(temp_p,1.2*units::cm,1);
      WireCell::CTPointCloud<double> cloud_w = ct_point_cloud.get_closest_points(temp_p,1.2*units::cm,2);
      
      //      std::cerr << cloud_u.pts.size() << " " << cloud_v.pts.size() << " " << cloud_w.pts.size() << std::endl;

      if (cloud_u.pts.size()>0 || cloud_v.pts.size()>0 || cloud_w.pts.size() > 0 || inside_dead_region(temp_p))
      	num_points_dead++;
      
      if (num_points - num_points_dead >=5) return true;
	
      temp_p.x += dir.X() * step;
      temp_p.y += dir.Y() * step;
      temp_p.z += dir.Z() * step;
    }

    //std::cout << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << " " << num_points << " " << num_points_dead << std::endl;
    
    if (num_points_dead > 0.8*num_points){
    	return false;
    }else{
    	return true;
    }
  }
  
  return true;
}

bool WireCell2dToy::ToyFiducial::check_dead_volume(WireCell::Point& p, TVector3& dir, double step, double offset_x){
  if (!inside_fiducial_volume(p,offset_x)){
    return false;
  }else{
    if (dir.Mag()==0){
      return true;
    }else{
      dir *= 1./dir.Mag();
      Point temp_p = p;
      int num_points = 0;
      int num_points_dead = 0;
      while(inside_fiducial_volume(temp_p,offset_x)){

	num_points ++;
	if (inside_dead_region(temp_p))
	  num_points_dead ++;

	if (num_points - num_points_dead >=4) return true;
	
	temp_p.x += dir.X() * step;
	temp_p.y += dir.Y() * step;
	temp_p.z += dir.Z() * step;
      }

      // std::cout << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << " " << num_points << " " << num_points_dead << std::endl;
      
      if (num_points_dead > 0.9*num_points){
	return false;
      }else{
	return true;
      }
      
    }
  }
}


bool WireCell2dToy::ToyFiducial::inside_fiducial_volume(WireCell::Point& p, double offset_x){

  int c1 = pnpoly(boundary_xy_x, boundary_xy_y, p.x-offset_x, p.y);
  int c2 = pnpoly(boundary_xz_x, boundary_xz_z, p.x-offset_x, p.z);

  //  std::cout << (p.x-offset_x)/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << std::endl;
  //std::cout << c1 << " " << c2 << std::endl;
  
  if (c1 && c2){
    return true;
  }else{
    return false;
  }
}

bool WireCell2dToy::ToyFiducial::inside_dead_region(WireCell::Point& p){
  // convert the position into U, V, W, and T number ...
  int time_slice = p.x * slope_t + offset_t;
  double pos_u = cos(angle_u) * p.z - sin(angle_u) *p.y;
  double pos_v = cos(angle_v) * p.z - sin(angle_v) *p.y;
  double pos_w = cos(angle_w) * p.z - sin(angle_w) *p.y;
  int ch_u = pos_u * slope_u + offset_u;
  int ch_v = pos_v * slope_v + offset_v + 2400;
  int ch_w = pos_w * slope_w + offset_w + 4800;

  //std::cout << ch_u << " " << ch_v << " " << ch_w << " " << time_slice << std::endl;
  //  std::cout << slope_w << " " << offset_w << " " << pos_w << std::endl;
  
  if (time_slice <0 || time_slice >=2398) return false;
  if (ch_u <0 || ch_u>=2400)  return false;
  if (ch_v <2400 || ch_v>=4800)  return false;
  if (ch_w <4800 || ch_w>=8256)  return false;

  std::set<SlimMergeGeomCell*> dead_u_mcells;
  std::set<SlimMergeGeomCell*> dead_v_mcells;
  std::set<SlimMergeGeomCell*> dead_w_mcells;

  if (ch_mcell_set_map.find(ch_u)!=ch_mcell_set_map.end())
    dead_u_mcells = ch_mcell_set_map[ch_u];
  if (ch_mcell_set_map.find(ch_v)!=ch_mcell_set_map.end())
    dead_v_mcells = ch_mcell_set_map[ch_v];
  if (ch_mcell_set_map.find(ch_w)!=ch_mcell_set_map.end())
    dead_w_mcells = ch_mcell_set_map[ch_w];
  
  // std::cout << ch_u << " " << ch_v << " " << ch578d-_w << " " << dead_u_mcells.size() << " " << dead_v_mcells.size() << " " << dead_w_mcells.size() << std::endl;
  
  // find the dead region given the U, V, and W number
  std::set<SlimMergeGeomCell*> results;
  for (auto it = dead_u_mcells.begin(); it!=dead_u_mcells.end(); it++){
    if (dead_v_mcells.find(*it)!=dead_v_mcells.end()){
      // compare UV sets
      results.insert(*it);
    }else if (dead_w_mcells.find(*it)!=dead_w_mcells.end()){
      // compare UW sets
      results.insert(*it);
    }
  }
  // compare VW sets ...
  for (auto it = dead_v_mcells.begin(); it!=dead_v_mcells.end(); it++){
    if (dead_w_mcells.find(*it)!=dead_w_mcells.end()){
      // compare UW sets
      results.insert(*it);
    }
  }
  
  
  // Check the T number for the remaining T ... 
  for (auto it = results.begin(); it!=results.end(); it++){
    if (mcell_time_map[*it].first <= time_slice &&
	time_slice <= mcell_time_map[*it].second)
      return true;
  }
  
  return false;
}


void WireCell2dToy::ToyFiducial::AddDeadRegion(WireCell::SlimMergeGeomCell* mcell, std::vector<int>& time_slices){

  mcells.push_back(mcell);
  int start_time = time_slices.front() - dead_region_ch_ext ;
  int end_time = time_slices.back() + dead_region_ch_ext;
  mcell_time_map[mcell] = std::make_pair(start_time, end_time);

  GeomWireSelection& uwires = mcell->get_uwires();
  GeomWireSelection& vwires = mcell->get_vwires();
  GeomWireSelection& wwires = mcell->get_wwires();

  // std::cout << uwires.size() << " " << vwires.size() << " " << wwires.size() << " " << start_time << " " << end_time << std::endl;
  
  std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();

  if (find(bad_planes.begin(), bad_planes.end(), WirePlaneType_t(0))!=bad_planes.end()){
    int start_ch = uwires.front()->channel() - dead_region_ch_ext;
    if (start_ch <0) start_ch = 0;
    int end_ch = uwires.back()->channel() + dead_region_ch_ext;
    if (end_ch>=2400) end_ch = 2399;
    for (int i = start_ch; i<=end_ch;i++){
      if (ch_mcell_set_map.find(i)==ch_mcell_set_map.end()){
	std::set<SlimMergeGeomCell*> mcells_set;
	mcells_set.insert(mcell);
	ch_mcell_set_map[i] = mcells_set;
      }else{
	ch_mcell_set_map[i].insert(mcell);
      }
    }
  }
  if (find(bad_planes.begin(), bad_planes.end(), WirePlaneType_t(1))!=bad_planes.end()){
    int start_ch = vwires.front()->channel() - dead_region_ch_ext;
    if (start_ch <2400) start_ch = 2400;
    int end_ch = vwires.back()->channel() + dead_region_ch_ext;
    if (end_ch>=4800) end_ch = 4799;
    for (int i = start_ch; i<=end_ch;i++){
      if (ch_mcell_set_map.find(i)==ch_mcell_set_map.end()){
	std::set<SlimMergeGeomCell*> mcells_set;
	mcells_set.insert(mcell);
	ch_mcell_set_map[i] = mcells_set;
      }else{
	ch_mcell_set_map[i].insert(mcell);
      }
    }
  }
  if (find(bad_planes.begin(), bad_planes.end(), WirePlaneType_t(2))!=bad_planes.end()){
    int start_ch = wwires.front()->channel() - dead_region_ch_ext;
    if (start_ch <4800) start_ch = 4800;
    int end_ch = wwires.back()->channel() + dead_region_ch_ext;
    if (end_ch>=8256) end_ch = 8255;
    for (int i = start_ch; i<=end_ch;i++){
      if (ch_mcell_set_map.find(i)==ch_mcell_set_map.end()){
	std::set<SlimMergeGeomCell*> mcells_set;
	mcells_set.insert(mcell);
	ch_mcell_set_map[i] = mcells_set;
      }else{
	ch_mcell_set_map[i].insert(mcell);
      }
    }
  }
  
  
}
