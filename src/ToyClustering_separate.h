#include <boost/graph/connected_components.hpp>

bool sortbysec(const std::pair<PR3DCluster*,double> &a,
	       const std::pair<PR3DCluster*,double> &b){
  return (a.second > b.second);
}



bool WireCell2dToy::JudgeSeparateDec_1(WireCell::PR3DCluster* cluster, TVector3& drift_dir){
  // get the main axis
  cluster->Calc_PCA();
  TVector3 dir1(cluster->get_PCA_axis(0).x,cluster->get_PCA_axis(0).y,cluster->get_PCA_axis(0).z);
  TVector3 dir2(cluster->get_PCA_axis(1).x,cluster->get_PCA_axis(1).y,cluster->get_PCA_axis(1).z);
  TVector3 dir3(cluster->get_PCA_axis(2).x,cluster->get_PCA_axis(2).y,cluster->get_PCA_axis(2).z);
  
  double angle1 = fabs(dir2.Angle(drift_dir)-3.1415926/2.)/3.1415926*180.;
  double angle2 = fabs(dir3.Angle(drift_dir)-3.1415926/2.)/3.1415926*180.;
  double ratio1 = cluster->get_PCA_value(1)/cluster->get_PCA_value(0) ;
  double ratio2 = cluster->get_PCA_value(2)/cluster->get_PCA_value(0) ;

  //std::cout << " " << ratio1 << " " << ratio2 << " " << angle1 << " " << angle2 << std::endl;
  
  if (ratio1 > pow(10,exp(1.38115-1.19312*pow(angle1,1./3.))-2.2) ||
      ratio2 > pow(10,exp(1.38115-1.19312*pow(angle2,1./3.))-2.2) || ratio1 > 0.75)
    return true;
  return false;
}

bool WireCell2dToy::JudgeSeparateDec_2(WireCell::PR3DCluster* cluster, TVector3& drift_dir, std::vector<WCPointCloud<double>::WCPoint>& boundary_points, std::vector<WCPointCloud<double>::WCPoint>& independent_points){
  cluster->Create_point_cloud();
  ToyPointCloud* cloud = cluster->get_point_cloud();
  boundary_points = cloud->get_hull();
  std::vector<WCPointCloud<double>::WCPoint> hy_points;
  std::vector<WCPointCloud<double>::WCPoint> ly_points;
  std::vector<WCPointCloud<double>::WCPoint> hz_points;
  std::vector<WCPointCloud<double>::WCPoint> lz_points;
  std::vector<WCPointCloud<double>::WCPoint> hx_points;
  std::vector<WCPointCloud<double>::WCPoint> lx_points;

  std::set<int> independent_surfaces;
  	
  for (size_t j=0;j!=boundary_points.size();j++){
    if (j==0){
      hy_points.push_back(boundary_points.at(j));
      ly_points.push_back(boundary_points.at(j));
      hz_points.push_back(boundary_points.at(j));
      lz_points.push_back(boundary_points.at(j));
      hx_points.push_back(boundary_points.at(j));
      lx_points.push_back(boundary_points.at(j));
    }else{
      Point test_p(boundary_points.at(j).x,boundary_points.at(j).y,boundary_points.at(j).z);
      if (cluster->get_num_points(test_p,15*units::cm) > 100){
	if (boundary_points.at(j).y > hy_points.at(0).y) hy_points.at(0) = boundary_points.at(j);
	if (boundary_points.at(j).y < ly_points.at(0).y) ly_points.at(0) = boundary_points.at(j);
	if (boundary_points.at(j).x > hx_points.at(0).x) hx_points.at(0) = boundary_points.at(j);
	if (boundary_points.at(j).x < lx_points.at(0).x) lx_points.at(0) = boundary_points.at(j);
	if (boundary_points.at(j).z > hz_points.at(0).z) hz_points.at(0) = boundary_points.at(j);
	if (boundary_points.at(j).z < lz_points.at(0).z) lz_points.at(0) = boundary_points.at(j);
      }
    }
  }
  
  
  
  
  if (hy_points.at(0).y > 101.5*units::cm){
    for (size_t j=0;j!=boundary_points.size();j++){
      if (boundary_points.at(j).y > 101.5*units::cm){
	bool flag_save = true;
	for (size_t k=0;k!=hy_points.size();k++){
	  double dis = sqrt(pow(hy_points.at(k).x-boundary_points.at(j).x,2)+pow(hy_points.at(k).y-boundary_points.at(j).y,2)+pow(hy_points.at(k).z-boundary_points.at(j).z,2));
	  if (dis <25*units::cm){
	    if (boundary_points.at(j).y > hy_points.at(k).y)
	      hy_points.at(k) = boundary_points.at(j);
	    flag_save = false;
	  }
	}
	if(flag_save)
	  hy_points.push_back(boundary_points.at(j));
      }
    }
  }
  
  if (ly_points.at(0).y < -99.5*units::cm){
    for (size_t j=0;j!=boundary_points.size();j++){
      if (boundary_points.at(j).y < -99.5*units::cm){
	bool flag_save = true;
	for (size_t k=0;k!=ly_points.size();k++){
	  double dis = sqrt(pow(ly_points.at(k).x-boundary_points.at(j).x,2)+pow(ly_points.at(k).y-boundary_points.at(j).y,2)+pow(ly_points.at(k).z-boundary_points.at(j).z,2));
	  if (dis <25*units::cm){
	    if (boundary_points.at(j).y < ly_points.at(k).y)
	      ly_points.at(k) = boundary_points.at(j);
	    flag_save = false;
	  }
	}
	if(flag_save)
	  ly_points.push_back(boundary_points.at(j));
      }
    }
  }
  if (hz_points.at(0).z > 1022*units::cm){
    for (size_t j=0;j!=boundary_points.size();j++){
      if (boundary_points.at(j).z > 1022*units::cm){
	bool flag_save = true;
	for (size_t k=0;k!=hz_points.size();k++){
	  double dis = sqrt(pow(hz_points.at(k).x-boundary_points.at(j).x,2)+pow(hz_points.at(k).y-boundary_points.at(j).y,2)+pow(hz_points.at(k).z-boundary_points.at(j).z,2));
	  if (dis <25*units::cm){
	    if (boundary_points.at(j).z > hz_points.at(k).z)
	      hz_points.at(k) = boundary_points.at(j);
	    flag_save = false;
	  }
	}
	if(flag_save)
	  hz_points.push_back(boundary_points.at(j));
      }
    }
  }
  if (lz_points.at(0).z <15*units::cm){
    for (size_t j=0;j!=boundary_points.size();j++){
      if (boundary_points.at(j).z < 15*units::cm){
	bool flag_save = true;
	for (size_t k=0;k!=lz_points.size();k++){
	  double dis = sqrt(pow(lz_points.at(k).x-boundary_points.at(j).x,2)+pow(lz_points.at(k).y-boundary_points.at(j).y,2)+pow(lz_points.at(k).z-boundary_points.at(j).z,2));
	  if (dis <25*units::cm){
	    if (boundary_points.at(j).z < lz_points.at(k).z)
	      lz_points.at(k) = boundary_points.at(j);
	    flag_save = false;
	  }
	}
	if(flag_save)
	  lz_points.push_back(boundary_points.at(j));
      }
    }
  }
  
  int num_outside_points = 0;
  int num_outx_points = 0;
  
  for (size_t j=0;j!=hy_points.size();j++){
    if (hy_points.at(j).x >=0*units::cm && hy_points.at(j).x <=256*units::cm &&
	hy_points.at(j).y >=-99.5*units::cm && hy_points.at(j).y <=101.5*units::cm &&
	hy_points.at(j).z >= 15*units::cm && hy_points.at(j).z <= 1022*units::cm)
      continue;
    
    bool flag_save = true;
    for (size_t k=0;k!=independent_points.size();k++){
      double dis = sqrt(pow(hy_points.at(j).x - independent_points.at(k).x,2)+pow(hy_points.at(j).y - independent_points.at(k).y,2)+pow(hy_points.at(j).z - independent_points.at(k).z,2));
      if (dis < 30*units::cm)
	flag_save = false;
    }
    if (flag_save){
      independent_points.push_back(hy_points.at(j));
      if (hy_points.at(j).y > 104*units::cm){
	independent_surfaces.insert(0);
      }else if (hy_points.at(j).y <-99.5*units::cm){
	independent_surfaces.insert(1);
      }else if (hy_points.at(j).z  > 1025*units::cm){
	independent_surfaces.insert(2);
      }else if (hy_points.at(j).z < 12*units::cm){
	independent_surfaces.insert(3);
      }else if (hy_points.at(j).x > 257*units::cm){
	independent_surfaces.insert(4);
      }else if (hy_points.at(j).x < -1*units::cm){
	independent_surfaces.insert(5);
      }

      
      if (hy_points.at(j).y > 104*units::cm || hy_points.at(j).y <-99.5*units::cm ||
	  hy_points.at(j).z < 12*units::cm || hy_points.at(j).z > 1025*units::cm ||
	  hy_points.at(j).x < -1*units::cm || hy_points.at(j).x > 257*units::cm )
	num_outside_points ++;
      if (hy_points.at(j).x < -1*units::cm || hy_points.at(j).x > 257*units::cm )
	num_outx_points++;
    }
    
  }
  for (size_t j=0;j!=ly_points.size();j++){
    if (ly_points.at(j).x >=0*units::cm && ly_points.at(j).x <=256*units::cm &&
	ly_points.at(j).y >=-99.5*units::cm && ly_points.at(j).y <=101.5*units::cm &&
	ly_points.at(j).z >= 15*units::cm && ly_points.at(j).z <= 1022*units::cm)
      continue;

    
    bool flag_save = true;
    for (size_t k=0;k!=independent_points.size();k++){
      double dis = sqrt(pow(ly_points.at(j).x - independent_points.at(k).x,2)+pow(ly_points.at(j).y - independent_points.at(k).y,2)+pow(ly_points.at(j).z - independent_points.at(k).z,2));
      if (dis < 30*units::cm)
	flag_save = false;
    }
    if (flag_save){
      independent_points.push_back(ly_points.at(j));
      

      if (ly_points.at(j).y <-99.5*units::cm){
	independent_surfaces.insert(1);
      }else if (ly_points.at(j).y > 104*units::cm){
	independent_surfaces.insert(0);
      }else if (ly_points.at(j).z  > 1025*units::cm){
	independent_surfaces.insert(2);
      }else if (ly_points.at(j).z < 12*units::cm){
	independent_surfaces.insert(3);
      }else if (ly_points.at(j).x > 257*units::cm){
	independent_surfaces.insert(4);
      }else if (ly_points.at(j).x < -1*units::cm){
	independent_surfaces.insert(5);
      }

      
      if (ly_points.at(j).y > 104*units::cm || ly_points.at(j).y <-99.5*units::cm ||
	  ly_points.at(j).z < 12*units::cm || ly_points.at(j).z > 1025*units::cm ||
	  ly_points.at(j).x < -1*units::cm || ly_points.at(j).x > 257*units::cm )
	num_outside_points ++;
      if (ly_points.at(j).x < -1*units::cm || ly_points.at(j).x > 257*units::cm )
	num_outx_points++;
    }
  }
  for (size_t j=0;j!=hz_points.size();j++){

    if (hz_points.at(j).x >=0*units::cm && hz_points.at(j).x <=256*units::cm &&
	hz_points.at(j).y >=-99.5*units::cm && hz_points.at(j).y <=101.5*units::cm &&
	hz_points.at(j).z >= 15*units::cm && hz_points.at(j).z <= 1022*units::cm)
      continue;
    
    bool flag_save = true;
    for (size_t k=0;k!=independent_points.size();k++){
      double dis = sqrt(pow(hz_points.at(j).x - independent_points.at(k).x,2)+pow(hz_points.at(j).y - independent_points.at(k).y,2)+pow(hz_points.at(j).z - independent_points.at(k).z,2));
      if (dis < 30*units::cm)
	flag_save = false;
    }
    if (flag_save){
      independent_points.push_back(hz_points.at(j));
            
      if (hz_points.at(j).z  > 1025*units::cm){
	independent_surfaces.insert(2);
      }else if (hz_points.at(j).z < 12*units::cm){
	independent_surfaces.insert(3);
      }else if (hz_points.at(j).y > 104*units::cm){
	independent_surfaces.insert(0);
      }else if (hz_points.at(j).y <-99.5*units::cm){
	independent_surfaces.insert(1);
      }else if (hz_points.at(j).x > 257*units::cm){
	independent_surfaces.insert(4);
      }else if (hz_points.at(j).x < -1*units::cm){
	independent_surfaces.insert(5);
      }

      
      if (hz_points.at(j).y > 104*units::cm || hz_points.at(j).y <-99.5*units::cm ||
	  hz_points.at(j).z < 12*units::cm || hz_points.at(j).z > 1025*units::cm ||
	  hz_points.at(j).x < -1*units::cm || hz_points.at(j).x > 257*units::cm )
	num_outside_points ++;
      if (hz_points.at(j).x < -1*units::cm || hz_points.at(j).x > 257*units::cm )
	num_outx_points ++;
    }
  }
  for (size_t j=0;j!=lz_points.size();j++){

    if (lz_points.at(j).x >=0*units::cm && lz_points.at(j).x <=256*units::cm &&
	lz_points.at(j).y >=-99.5*units::cm && lz_points.at(j).y <=101.5*units::cm &&
	lz_points.at(j).z >= 15*units::cm && lz_points.at(j).z <= 1022*units::cm)
      continue;
    
    bool flag_save = true;
    for (size_t k=0;k!=independent_points.size();k++){
      double dis = sqrt(pow(lz_points.at(j).x - independent_points.at(k).x,2)+pow(lz_points.at(j).y - independent_points.at(k).y,2)+pow(lz_points.at(j).z - independent_points.at(k).z,2));
      if (dis < 30*units::cm)
	flag_save = false;
    }
    if (flag_save){
      independent_points.push_back(lz_points.at(j));

      if (lz_points.at(j).z < 12*units::cm){
	independent_surfaces.insert(3);
      }else if (lz_points.at(j).z  > 1025*units::cm){
	independent_surfaces.insert(2);
      }else if (lz_points.at(j).y > 104*units::cm){
	independent_surfaces.insert(0);
      }else if (lz_points.at(j).y <-99.5*units::cm){
	independent_surfaces.insert(1);
      }else if (lz_points.at(j).x > 257*units::cm){
	independent_surfaces.insert(4);
      }else if (lz_points.at(j).x < -1*units::cm){
	independent_surfaces.insert(5);
      }


      
      if (lz_points.at(j).y > 104*units::cm || lz_points.at(j).y <-99.5*units::cm ||
	  lz_points.at(j).z < 12*units::cm || lz_points.at(j).z > 1025*units::cm ||
	  lz_points.at(j).x < -1*units::cm || lz_points.at(j).x > 257*units::cm )
	num_outside_points ++;
      if (lz_points.at(j).x < -1*units::cm || lz_points.at(j).x > 257*units::cm )
	num_outx_points++;
    }
  }
  for (size_t j=0;j!=hx_points.size();j++){

    if (hx_points.at(j).x >=0*units::cm && hx_points.at(j).x <=256*units::cm &&
	hx_points.at(j).y >=-99.5*units::cm && hx_points.at(j).y <=101.5*units::cm &&
	hx_points.at(j).z >= 15*units::cm && hx_points.at(j).z <= 1022*units::cm)
      continue;
    
    bool flag_save = true;
    for (size_t k=0;k!=independent_points.size();k++){
      double dis = sqrt(pow(hx_points.at(j).x - independent_points.at(k).x,2)+pow(hx_points.at(j).y - independent_points.at(k).y,2)+pow(hx_points.at(j).z - independent_points.at(k).z,2));
      if (dis < 30*units::cm)
	flag_save = false;
    }
    if (flag_save){
      independent_points.push_back(hx_points.at(j));
      
      if (hx_points.at(j).y > 104*units::cm || hx_points.at(j).y <-99.5*units::cm ||
	  hx_points.at(j).z < 12*units::cm || hx_points.at(j).z > 1025*units::cm ||
	  hx_points.at(j).x < -1*units::cm || hx_points.at(j).x > 257*units::cm )
	num_outside_points ++;
      if (hx_points.at(j).x < -1*units::cm || hx_points.at(j).x > 257*units::cm ){
	
	num_outx_points++;
      }

      if (lx_points.at(j).x > 257*units::cm){
	independent_surfaces.insert(4);
      }else if (lx_points.at(j).x < -1*units::cm){
	independent_surfaces.insert(5);
      }else if (lx_points.at(j).y > 104*units::cm){
	independent_surfaces.insert(0);
      }else if (lx_points.at(j).y <-99.5*units::cm){
	independent_surfaces.insert(1);
      }else if (lx_points.at(j).z  > 1025*units::cm){
	independent_surfaces.insert(2);
      }else if (lx_points.at(j).z < 12*units::cm){
	independent_surfaces.insert(3);
      }
      
    }
  }
  for (size_t j=0;j!=lx_points.size();j++){

    if (lx_points.at(j).x >=0*units::cm && lx_points.at(j).x <=256*units::cm &&
	lx_points.at(j).y >=-99.5*units::cm && lx_points.at(j).y <=101.5*units::cm &&
	lx_points.at(j).z >= 15*units::cm && lx_points.at(j).z <= 1022*units::cm)
      continue;
    
    bool flag_save = true;
    for (size_t k=0;k!=independent_points.size();k++){
      double dis = sqrt(pow(lx_points.at(j).x - independent_points.at(k).x,2)+pow(lx_points.at(j).y - independent_points.at(k).y,2)+pow(lx_points.at(j).z - independent_points.at(k).z,2));
      if (dis < 30*units::cm)
	flag_save = false;
    }
    if (flag_save){
      independent_points.push_back(lx_points.at(j));
  
      if (lx_points.at(j).y > 104*units::cm || lx_points.at(j).y <-99.5*units::cm ||
	  lx_points.at(j).z < 12*units::cm || lx_points.at(j).z > 1025*units::cm ||
	  lx_points.at(j).x < -1*units::cm || lx_points.at(j).x > 257*units::cm )
	num_outside_points ++;
      if (lx_points.at(j).x < -1*units::cm || lx_points.at(j).x > 257*units::cm ){
	num_outx_points++;
	
      }

      if (lx_points.at(j).x < -1*units::cm){
	independent_surfaces.insert(5);
      }else if (lx_points.at(j).x > 257*units::cm){
	independent_surfaces.insert(4);
      }else if (lx_points.at(j).y > 104*units::cm){
	independent_surfaces.insert(0);
      }else if (lx_points.at(j).y <-99.5*units::cm){
	independent_surfaces.insert(1);
      }else if (lx_points.at(j).z  > 1025*units::cm){
	independent_surfaces.insert(2);
      }else if (lx_points.at(j).z < 12*units::cm){
	independent_surfaces.insert(3);
      }

      
    }
  }
  
  int num_far_points = 0;
  
  
  if (independent_points.size()==2&&independent_surfaces.size()>1){
    TVector3 dir_1(independent_points.at(1).x - independent_points.at(0).x, independent_points.at(1).y - independent_points.at(0).y, independent_points.at(1).z - independent_points.at(0).z);
    dir_1.SetMag(1);
    for (size_t j=0;j!=boundary_points.size();j++){
      TVector3 dir_2(boundary_points.at(j).x - independent_points.at(0).x, boundary_points.at(j).y - independent_points.at(0).y, boundary_points.at(j).z - independent_points.at(0).z);
      double angle_12 = dir_1.Angle(dir_2);
      TVector3 dir_3 = dir_2 - dir_1 * dir_2.Mag() * cos(angle_12);
      double angle_3 = dir_3.Angle(drift_dir);
      //std::cout << dir_3.Mag()/units::cm << " " << fabs(angle_3-3.1415926/2.)/3.1415926*180. << " " << fabs(dir_3.X()/units::cm) << std::endl;
      if (fabs(angle_3-3.1415926/2.)/3.1415926*180.<7.5){
	if (fabs(dir_3.X()/units::cm)>14*units::cm)
	  num_far_points ++;
	if (fabs(dir_1.Angle(drift_dir)-3.1415926/2.)/3.1415926*180. > 15){
	  if (dir_3.Mag() > 20*units::cm)
	    num_far_points ++;
	}
      }else{
	if (dir_3.Mag() > 20*units::cm)
	  num_far_points ++;
      }
    }

   
    // find the middle points and close distance ...
    Point middle_point((independent_points.at(1).x + independent_points.at(0).x)/2.,
		       (independent_points.at(1).y + independent_points.at(0).y)/2.,
		       (independent_points.at(1).z + independent_points.at(0).z)/2.);
    double middle_dis = cloud->get_closest_dis(middle_point);
    if (middle_dis > 25*units::cm){
      num_far_points = 0;
    }
    
    
  }

  
  
  
  // std::cout <<  cluster->get_cluster_id() << " " << hy_points.size() << " " << ly_points.size() << " " << hz_points.size() << " " << lz_points.size() <<  " " << hx_points.size() << " " << lz_points.size() << " " << num_outside_points << " " << num_outx_points << " " << independent_points.size() << " " << num_far_points << " " << independent_surfaces.size() << std::endl;
  /* for (auto it=independent_surfaces.begin(); it!=independent_surfaces.end(); it++){ */
  /*   std::cout << (*it) << std::endl; */
  /* } */

  /* for (auto it = independent_points.begin(); it!=independent_points.end(); it++){ */
  /*   std::cout << (*it).x/units::cm << " " << (*it).y/units::cm << " " << (*it).z/units::cm << std::endl; */
  /* } */
  
  
  if ((num_outside_points > 1 && independent_surfaces.size()>1
       || num_outx_points>0) &&
      (independent_points.size()>2 ||
       independent_points.size()==2 && num_far_points > 0))
    return true;
  return false;
}

std::vector<WireCell::PR3DCluster*> WireCell2dToy::Separate_1(WireCell::PR3DCluster* cluster,std::vector<WireCell::WCPointCloud<double>::WCPoint>& boundary_points, std::vector<WireCell::WCPointCloud<double>::WCPoint>& independent_points){

  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();

  TVector3 dir_drift(1,0,0);
  TVector3 dir_cosmic(0,1,0);
  TVector3 dir_beam(0,0,1);
  
  ToyPointCloud *temp_cloud = new ToyPointCloud(angle_u,angle_v,angle_w);
  
  ToyPointCloud* cloud = cluster->get_point_cloud();

  Point cluster_center(cluster->get_center().x, cluster->get_center().y, cluster->get_center().z);

  //std::cout << cluster->get_PCA_value(0) << " " << cluster->get_PCA_value(1) << " " << cluster->get_PCA_value(2) << " " << cluster->get_PCA_axis(0) << " " << cluster->get_PCA_axis(1) << std::endl;
  
  TVector3 main_dir,second_dir;
  main_dir.SetXYZ(cluster->get_PCA_axis(0).x,cluster->get_PCA_axis(0).y,cluster->get_PCA_axis(0).z);
  second_dir.SetXYZ(cluster->get_PCA_axis(1).x,cluster->get_PCA_axis(1).y,cluster->get_PCA_axis(1).z);

  // special case, if one of the cosmic is very close to the beam direction
  if ( cluster->get_PCA_value(1) > 0.08 * cluster->get_PCA_value(0) &&
       fabs(main_dir.Angle(dir_beam)-3.1415926/2.) > 75/180.*3.1415926 &&
       fabs(second_dir.Angle(dir_cosmic)-3.1415926/2.) > 60/180.*3.1415926){
    main_dir = second_dir;
  }
  //  std::cout << main_dir.Angle(dir_beam)/3.1415926*180. << " " << second_dir.Angle(dir_cosmic)/3.1415926*180. << " " << independent_points.size() << " " << std::endl;
  
  main_dir.SetMag(1);
  if (main_dir.Y()>0) main_dir *= -1;
    //  std::cout << cluster->get_PCA_value(0) << " " << cluster->get_PCA_value(1) << " " << cluster->get_PCA_value(2) << std::endl;
  
  double min_dis = 1e9;
  int min_index = 0;
  for (size_t j=0; j!= independent_points.size(); j++){
    TVector3 dir(independent_points.at(j).x - cluster_center.x, independent_points.at(j).y - cluster_center.y, independent_points.at(j).z -cluster_center.z);
    Point temp_p(independent_points.at(j).x, independent_points.at(j).y, independent_points.at(j).z);
    double dis = dir.Dot(main_dir);
    //std::cout << j << " " << dis << " " << dir.Mag() << " " << sqrt(dir.Mag()*dir.Mag() - dis*dis) << std::endl;
    if (dis < min_dis && cluster->get_num_points(temp_p,15*units::cm) > 100){
      min_dis = dis;
      min_index = j;
    }
  }

  /* std::cout << independent_points.at(min_index).x/units::cm << " " << independent_points.at(min_index).y/units::cm << "  " << independent_points.at(min_index).z/units::cm << " " << min_index << " " << min_dis/units::cm << " " << independent_points.size() << std::endl; */
  /* for (size_t k=0;k!=independent_points.size();k++){ */
  /*   std::cout << independent_points.at(k).x/units::cm << " " << independent_points.at(k).y/units::cm << "  " << independent_points.at(k).z/units::cm << std::endl; */
  /* } */
  
   
  WCPointCloud<double>::WCPoint start_wcpoint = independent_points.at(min_index);
  Point start_point(start_wcpoint.x, start_wcpoint.y, start_wcpoint.z);
  TVector3 dir = cluster->VHoughTrans(start_point,100*units::cm);
  dir.SetMag(1);

  TVector3 drift_dir(1,0,0);

  // std::cout  << " " << start_wcpoint.x/units::cm << " " << start_wcpoint.y/units::cm << " " << start_wcpoint.z/units::cm << " " <<  " " << dir.X() << " " << dir.Y() << " " << dir.Z() << std::endl;

  
  TVector3 inv_dir = dir * (-1);
  start_wcpoint = cluster->get_furthest_wcpoint(start_wcpoint,inv_dir,1*units::cm,0);
  WCPointCloud<double>::WCPoint end_wcpoint = cluster->get_furthest_wcpoint(start_wcpoint,dir);

  TVector3 test_dir(end_wcpoint.x - start_wcpoint.x, end_wcpoint.y - start_wcpoint.y, end_wcpoint.z - start_wcpoint.z);
  
  //std::cout  << " XQ1 " << start_wcpoint.x/units::cm << " " << start_wcpoint.y/units::cm << " " << start_wcpoint.z/units::cm << " " << end_wcpoint.x/units::cm << " " << end_wcpoint.y/units::cm << " " << end_wcpoint.z/units::cm << " " << dir.X() << " " << dir.Y() << " " << dir.Z() << std::endl;
  
  if (fabs(test_dir.Angle(drift_dir)-3.1415926/2.)<2.0*3.1415926/180.){
    cluster->adjust_wcpoints_parallel(start_wcpoint,end_wcpoint);
    // std::cout << "Parallel Case! " << " " << fabs(test_dir.Angle(drift_dir)-3.1415926/2.)/3.1415926*180. << std::endl;
  }
  
  //std::cout  << " XQ " << start_wcpoint.x/units::cm << " " << start_wcpoint.y/units::cm << " " << start_wcpoint.z/units::cm << " " << end_wcpoint.x/units::cm << " " << end_wcpoint.y/units::cm << " " << end_wcpoint.z/units::cm << " " << dir.X() << " " << dir.Y() << " " << dir.Z() << std::endl;
  
  
  cluster->dijkstra_shortest_paths(start_wcpoint);
  cluster->cal_shortest_path(end_wcpoint);
  
  std::list<WCPointCloud<double>::WCPoint>& path_wcps = cluster->get_path_wcps();
  std::vector<bool> flag_u_pts, flag_v_pts, flag_w_pts;
  std::vector<bool> flag1_u_pts, flag1_v_pts, flag1_w_pts;
  flag_u_pts.resize(cloud->get_cloud().pts.size(),false);
  flag_v_pts.resize(cloud->get_cloud().pts.size(),false);
  flag_w_pts.resize(cloud->get_cloud().pts.size(),false);

  flag1_u_pts.resize(cloud->get_cloud().pts.size(),false);
  flag1_v_pts.resize(cloud->get_cloud().pts.size(),false);
  flag1_w_pts.resize(cloud->get_cloud().pts.size(),false);
  
  PointVector pts;
  
  WCPointCloud<double>::WCPoint prev_wcp = path_wcps.front();
  for (auto it = path_wcps.begin(); it!=path_wcps.end();it++){

    // std::cout << "a: " << (*it).x/units::cm << " " << (*it).y/units::cm << " " << (*it).z/units::cm << std::endl;
    
    
    double dis = sqrt(pow((*it).x - prev_wcp.x,2) + pow((*it).y - prev_wcp.y,2) + pow((*it).z - prev_wcp.z,2));
    if (dis <=1.5*units::cm){
      Point current_pt((*it).x,(*it).y,(*it).z);
      pts.push_back(current_pt);
    }else{
      int num_points = int(dis/(1.5*units::cm))+1;
      double dis_seg = dis/num_points;
      for (int k=0;k!=num_points;k++){
	Point current_pt(prev_wcp.x + (k+1.)/num_points*((*it).x - prev_wcp.x),
			 prev_wcp.y + (k+1.)/num_points*((*it).y - prev_wcp.y),
			 prev_wcp.z + (k+1.)/num_points*((*it).z - prev_wcp.z));
	pts.push_back(current_pt);
	//	std::cout << "b: " << current_pt.x/units::cm << " " << current_pt.y/units::cm << " " << current_pt.z/units::cm << std::endl;
      }
    }
    prev_wcp = (*it);
  }
  temp_cloud->AddPoints(pts);
  temp_cloud->build_kdtree_index();
  
  for (size_t j=0;j!=flag_u_pts.size();j++){
    Point test_p;
    test_p.x = cloud->get_cloud().pts[j].x;
    test_p.y = cloud->get_cloud().pts[j].y;
    test_p.z = cloud->get_cloud().pts[j].z;
    std::pair<int,double> temp_results = temp_cloud->get_closest_2d_dis(test_p,0);
    double dis = temp_results.second;
    if (dis <= 1.5*units::cm){
      flag_u_pts.at(cloud->get_cloud().pts[j].index) = true;
    }
    if (dis <= 5.0*units::cm){
      flag1_u_pts.at(cloud->get_cloud().pts[j].index) = true;
    }
    temp_results = temp_cloud->get_closest_2d_dis(test_p,1);
    dis = temp_results.second;
    if (dis <= 1.5*units::cm){
      flag_v_pts.at(cloud->get_cloud().pts[j].index) = true;
    }
    if (dis <= 5.0*units::cm){
      flag1_v_pts.at(cloud->get_cloud().pts[j].index) = true;
    }
    temp_results = temp_cloud->get_closest_2d_dis(test_p,2);
    dis = temp_results.second;
    if (dis <= 1.5*units::cm){
      flag_w_pts.at(cloud->get_cloud().pts[j].index) = true;
    }
    if (dis <= 5.0*units::cm){
      flag1_w_pts.at(cloud->get_cloud().pts[j].index) = true;
    }
    /* if (fabs(test_p.y-76.7*units::cm)<2*units::cm && fabs(test_p.z-328.9*units::cm)<2*units::cm && fabs(test_p.x - 132.8*units::cm) < 2*units::cm){ */
    /*   std::pair<int,double> temp0 = temp_cloud->get_closest_2d_dis(test_p,0); */
    /*   std::pair<int,double> temp1 = temp_cloud->get_closest_2d_dis(test_p,1); */
    /*   std::pair<int,double> temp2 = temp_cloud->get_closest_2d_dis(test_p,2); */
    /*   std::cout << test_p.x/units::cm << " " << test_p.y/units::cm << " " << test_p.z/units::cm << " " << flag_u_pts.at(cloud->get_cloud().pts[j].index) << " " << flag_v_pts.at(cloud->get_cloud().pts[j].index) << " " << flag_w_pts.at(cloud->get_cloud().pts[j].index) << " " << flag1_u_pts.at(cloud->get_cloud().pts[j].index) << " " << flag1_v_pts.at(cloud->get_cloud().pts[j].index) << " " << flag1_w_pts.at(cloud->get_cloud().pts[j].index) << " " << temp0.second/units::cm << " " << pts[temp0.first].x/units::cm << " " << pts[temp0.first].y/units::cm << " " << pts[temp0.first].z/units::cm << " " << cloud->get_cloud().pts[j].mcell->GetTimeSlice() << std::endl; */
    /* } */
  }

  // special treatment of first and last point
  {
    std::vector<size_t> indices = cloud->get_closest_2d_index(pts.front(),3.0*units::cm,0);
    for (size_t k=0;k!=indices.size();k++){
      flag_u_pts.at(indices.at(k)) = true;
    }
    indices = cloud->get_closest_2d_index(pts.front(),3.0*units::cm,1);
    for (size_t k=0;k!=indices.size();k++){
      flag_v_pts.at(indices.at(k)) = true;
    }
    indices = cloud->get_closest_2d_index(pts.front(),3.0*units::cm,2);
    for (size_t k=0;k!=indices.size();k++){
      flag_w_pts.at(indices.at(k)) = true;
    }
    indices = cloud->get_closest_2d_index(pts.back(),3.0*units::cm,0);
    for (size_t k=0;k!=indices.size();k++){
      flag_u_pts.at(indices.at(k)) = true;
    }
    indices = cloud->get_closest_2d_index(pts.back(),3.0*units::cm,1);
    for (size_t k=0;k!=indices.size();k++){
      flag_v_pts.at(indices.at(k)) = true;
    }
    indices = cloud->get_closest_2d_index(pts.back(),3.0*units::cm,2);
    for (size_t k=0;k!=indices.size();k++){
      flag_w_pts.at(indices.at(k)) = true;
    }
    
  }
  
  
  SMGCSelection& mcells = cluster->get_mcells();
  std::map<SlimMergeGeomCell*, int> mcell_np_map;
  for (auto it=mcells.begin(); it!=mcells.end(); it++){
    mcell_np_map[*it] = 0;
  }
  for (size_t j=0;j!=flag_u_pts.size();j++){
    if (flag_u_pts.at(j) && flag_v_pts.at(j) && flag1_w_pts.at(j) ||
	flag_u_pts.at(j) && flag_w_pts.at(j) && flag1_v_pts.at(j) ||
	flag_w_pts.at(j) && flag_v_pts.at(j) && flag1_u_pts.at(j)){
      mcell_np_map[cloud->get_cloud().pts.at(j).mcell] ++;
    }
  }
  

  PR3DCluster *cluster1 = new PR3DCluster(1);
  PR3DCluster *cluster2 = new PR3DCluster(2);

  for (auto it=mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = *it;
    /* if (fabs(mcell->get_sampling_points().front().x-132.12*units::cm)<0.1*units::mm) */
    /*   std::cout << mcell_np_map[mcell] << " a " << mcell->get_sampling_points().size() << */
    /* 	" " << mcell->get_uwires().size() << " " <<  mcell->get_vwires().size() << " " <<  mcell->get_wwires().size() << std::endl; */
    
    if (mcell_np_map[mcell] > 0.5 * mcell->get_sampling_points().size() ||
	mcell_np_map[mcell] > 0.25 * mcell->get_sampling_points().size() && 
	mcell->get_uwires().size() + mcell->get_vwires().size() +  mcell->get_wwires().size() < 25){
      cluster1->AddCell(mcell,mcell->GetTimeSlice());
    }else{
      cluster2->AddCell(mcell,mcell->GetTimeSlice());
    }
  }

  std::vector<WireCell::PR3DCluster*> final_clusters;
  final_clusters.push_back(cluster1);
  
  std::vector<WireCell::PR3DCluster*> other_clusters = Separate_2(cluster2, 5*units::cm);
  delete cluster2;
  
  std::vector<WireCell::PR3DCluster*> saved_clusters;
  std::vector<WireCell::PR3DCluster*> to_be_merged_clusters;
  for (size_t i=0;i!=other_clusters.size();i++){

    // How to write???
    bool flag_save = false;
    std::vector<int> range_v1 = other_clusters.at(i)->get_uvwt_range();
    double length_1 = sqrt(2./3. * (pow(pitch_u*range_v1.at(0),2) + pow(pitch_v*range_v1.at(1),2) + pow(pitch_w*range_v1.at(2),2)) + pow(time_slice_width*range_v1.at(3),2));
    other_clusters.at(i)->Create_point_cloud();
    ToyPointCloud* temp_cloud1 = other_clusters.at(i)->get_point_cloud();
    std::tuple<int,int,double> temp_dis = temp_cloud1->get_closest_points(temp_cloud);
    
    if (length_1 < 30*units::cm && std::get<2>(temp_dis)<5*units::cm){
      
      int temp_total_points = other_clusters.at(i)->get_num_points();
      int temp_close_points = 0;
      for (size_t j=0;j!=temp_cloud1->get_num_points();j++){
	Point test_point(temp_cloud1->get_cloud().pts.at(j).x,temp_cloud1->get_cloud().pts.at(j).y,temp_cloud1->get_cloud().pts.at(j).z);
	if (temp_cloud->get_closest_dis(test_point) < 10*units::cm){
	  temp_close_points ++;
	}
      }
      //std::cout << temp_close_points << " " << temp_total_points << std::endl;
      if (temp_close_points > 0.8 * temp_total_points){
	saved_clusters.push_back(other_clusters.at(i));
	flag_save = true;
      }
      //std::vector<size_t> indices = cloud->get_closest_2d_index(pts.front(),2.4*units::cm,0);
      
    }
    
    
    if (!flag_save)
      to_be_merged_clusters.push_back(other_clusters.at(i));
  }

  
  cluster2 = new PR3DCluster(2);
  for (size_t i=0;i!=to_be_merged_clusters.size();i++){
    SMGCSelection& temp_mcells = to_be_merged_clusters.at(i)->get_mcells();
     for (auto it=temp_mcells.begin(); it!=temp_mcells.end(); it++){
       SlimMergeGeomCell *mcell = *it;
       cluster2->AddCell(mcell,mcell->GetTimeSlice());
     }
     delete to_be_merged_clusters.at(i);
  }
  
  final_clusters.push_back(cluster2);
  for (size_t i=0;i!=saved_clusters.size();i++){
    final_clusters.push_back(saved_clusters.at(i));
  }

  
  delete temp_cloud;


  /* int num_mcells = 0; */
  /* for (size_t i=0;i!=final_clusters.size();i++){ */
  /*   num_mcells += final_clusters.at(i)->get_num_mcells(); */
  /* } */
  /* std::cout << cluster->get_num_mcells() << " " << num_mcells << std::endl; */
  
  
  return final_clusters;
  
}
  


void WireCell2dToy::Clustering_separate(WireCell::PR3DClusterSelection& live_clusters, std::map<WireCell::PR3DCluster*,double>& cluster_length_map){
  TVector3 drift_dir(1,0,0);

  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();

  std::vector<PR3DCluster*> new_clusters;
  std::vector<PR3DCluster*> del_clusters;
  
  for (size_t i=0;i!=live_clusters.size();i++){
    PR3DCluster *cluster = live_clusters.at(i);
    
    
    if (cluster_length_map[cluster]> 100*units::cm){
      std::vector<WCPointCloud<double>::WCPoint> boundary_points;
      std::vector<WCPointCloud<double>::WCPoint> independent_points;
      //std::cout << cluster->get_cluster_id() << " " << WireCell2dToy::JudgeSeparateDec_1(cluster,drift_dir)  << " " << std::endl;
      if (WireCell2dToy::JudgeSeparateDec_2(cluster,drift_dir,boundary_points,independent_points)){
	if (WireCell2dToy::JudgeSeparateDec_1(cluster,drift_dir)){
	  std::cout << "Separate cluster " << cluster->get_cluster_id() << std::endl;
	  
	  
	  std::vector<PR3DCluster*> sep_clusters = WireCell2dToy::Separate_1(cluster,boundary_points,independent_points);
	  PR3DCluster* cluster1 = sep_clusters.at(0);
	  new_clusters.push_back(cluster1);
	  del_clusters.push_back(cluster);
	  for (size_t k=2;k<sep_clusters.size();k++){
	    new_clusters.push_back(sep_clusters.at(k));
	  }
	  
	  
	  std::vector<PR3DCluster*> temp_del_clusters;
	  PR3DCluster* cluster2 = sep_clusters.at(1);
	  std::vector<int> range_v1 = cluster2->get_uvwt_range();
	  double length_1 = sqrt(2./3. * (pow(pitch_u*range_v1.at(0),2) + pow(pitch_v*range_v1.at(1),2) + pow(pitch_w*range_v1.at(2),2)) + pow(time_slice_width*range_v1.at(3),2));
	  
	  PR3DCluster* final_sep_cluster = cluster2;
	  
	  
	  if (length_1 > 100*units::cm){
	    boundary_points.clear();
	    independent_points.clear();
	    if (WireCell2dToy::JudgeSeparateDec_1(cluster2,drift_dir) &&
		WireCell2dToy::JudgeSeparateDec_2(cluster2,drift_dir,boundary_points,independent_points)){
	      
	      // std::cout << "Separate 2nd level" << std::endl;
	      
	      std::vector<PR3DCluster*>  sep_clusters = WireCell2dToy::Separate_1(cluster2,boundary_points,independent_points);
	      PR3DCluster* cluster3 = sep_clusters.at(0);
	      new_clusters.push_back(cluster3);
	      temp_del_clusters.push_back(cluster2);
	      for (size_t k=2;k<sep_clusters.size();k++){
		new_clusters.push_back(sep_clusters.at(k));
	      }
	      
	      PR3DCluster* cluster4 = sep_clusters.at(1);
	      final_sep_cluster = cluster4;
	      range_v1 = cluster4->get_uvwt_range();
	      length_1 = sqrt(2./3. * (pow(pitch_u*range_v1.at(0),2) + pow(pitch_v*range_v1.at(1),2) + pow(pitch_w*range_v1.at(2),2)) + pow(time_slice_width*range_v1.at(3),2));
	      
	      if (length_1 > 100*units::cm){
		boundary_points.clear();
		independent_points.clear();
		if (WireCell2dToy::JudgeSeparateDec_1(cluster4,drift_dir) &&
		    WireCell2dToy::JudgeSeparateDec_2(cluster4,drift_dir,boundary_points,independent_points)){
		  //	std::cout << "Separate 3rd level" << std::endl;
		  
		  std::vector<PR3DCluster*>  sep_clusters = WireCell2dToy::Separate_1(cluster4,boundary_points,independent_points);
		  PR3DCluster* cluster5 = sep_clusters.at(0);
		  new_clusters.push_back(cluster5);
		  temp_del_clusters.push_back(cluster4);
		  for (size_t k=2;k<sep_clusters.size();k++){
		    new_clusters.push_back(sep_clusters.at(k));
		  }
		  PR3DCluster* cluster6 = sep_clusters.at(1);
		  final_sep_cluster = cluster6;
		}
	      }
	    }
	  }
	  
	  // temporary
	  range_v1 = final_sep_cluster->get_uvwt_range();
	  length_1 = sqrt(2./3. * (pow(pitch_u*range_v1.at(0),2) + pow(pitch_v*range_v1.at(1),2) + pow(pitch_w*range_v1.at(2),2)) + pow(time_slice_width*range_v1.at(3),2));
	  
	  //	std::cout << length_1/units::cm << std::endl;
	  
	  if (length_1 > 60*units::cm){
	    boundary_points.clear();
	    independent_points.clear();
	    WireCell2dToy::JudgeSeparateDec_1(final_sep_cluster,drift_dir);
	    WireCell2dToy::JudgeSeparateDec_2(final_sep_cluster,drift_dir,boundary_points,independent_points);
	    if (independent_points.size() > 0){
	      
	      // std::cout << "Separate final one" << std::endl;
	      
	      std::vector<PR3DCluster*> sep_clusters = WireCell2dToy::Separate_1(final_sep_cluster,boundary_points,independent_points);
	      PR3DCluster* cluster5 = sep_clusters.at(0);
	      new_clusters.push_back(cluster5);
	      temp_del_clusters.push_back(final_sep_cluster);
	      for (size_t k=2;k<sep_clusters.size();k++){
		new_clusters.push_back(sep_clusters.at(k));
	      }
	      
	      PR3DCluster* cluster6 = sep_clusters.at(1);
	      final_sep_cluster = cluster6;
	    }	    
	  }
	  
	  std::vector<PR3DCluster*> final_sep_clusters = Separate_2(final_sep_cluster);
	  for (auto it = final_sep_clusters.begin(); it!=final_sep_clusters.end(); it++){
	    new_clusters.push_back(*it);
	  }
	  temp_del_clusters.push_back(final_sep_cluster);
	  
	  
	  /* int num_mcells = 0; */
	  /* for (size_t i=0;i!=new_clusters.size();i++){ */
	  /*   num_mcells += new_clusters.at(i)->get_num_mcells(); */
	  /* } */
	  //std::cout << cluster->get_num_mcells() << " " << num_mcells << std::endl; 
	  
	  
	  
	  for (auto it = temp_del_clusters.begin(); it!= temp_del_clusters.end(); it++){
	    delete *it;
	  }
	}else{
	  std::cout << "Stripping Cluster " <<   cluster->get_cluster_id() << std::endl;

	  std::vector<PR3DCluster*> sep_clusters = WireCell2dToy::Separate_1(cluster,boundary_points,independent_points);
	  
	  PR3DCluster* cluster1 = sep_clusters.at(0);
	  new_clusters.push_back(cluster1);
	  del_clusters.push_back(cluster);
	  for (size_t k=2;k<sep_clusters.size();k++){
	    new_clusters.push_back(sep_clusters.at(k));
	  }
	  	  
	  std::vector<PR3DCluster*> temp_del_clusters;
	  PR3DCluster* cluster2 = sep_clusters.at(1);
	  std::vector<int> range_v1 = cluster2->get_uvwt_range();
	  double length_1 = sqrt(2./3. * (pow(pitch_u*range_v1.at(0),2) + pow(pitch_v*range_v1.at(1),2) + pow(pitch_w*range_v1.at(2),2)) + pow(time_slice_width*range_v1.at(3),2));
	  PR3DCluster* final_sep_cluster = cluster2;
	  
	  std::vector<PR3DCluster*> final_sep_clusters = Separate_2(final_sep_cluster);
	  for (auto it = final_sep_clusters.begin(); it!=final_sep_clusters.end(); it++){
	    new_clusters.push_back(*it);
	  }
	  temp_del_clusters.push_back(final_sep_cluster);
	  
	  for (auto it = temp_del_clusters.begin(); it!= temp_del_clusters.end(); it++){
	    delete *it;
	  }
	  
	}
      } // else ... 
    }
  }
  
  for (auto it=new_clusters.begin(); it!=new_clusters.end(); it++){
    PR3DCluster *ncluster = (*it);
    std::vector<int> range_v1 = ncluster->get_uvwt_range();
    double length_1 = sqrt(2./3. * (pow(pitch_u*range_v1.at(0),2) + pow(pitch_v*range_v1.at(1),2) + pow(pitch_w*range_v1.at(2),2)) + pow(time_slice_width*range_v1.at(3),2));
    cluster_length_map[ncluster] = length_1;
    live_clusters.push_back(ncluster);
  }
  for (auto it=del_clusters.begin(); it!=del_clusters.end(); it++){
    PR3DCluster *ocluster = (*it);
    cluster_length_map.erase(ocluster);
    live_clusters.erase(find(live_clusters.begin(), live_clusters.end(), ocluster));
    delete ocluster;
  }

 
  
   //std::cout << live_clusters.size() << " " << cluster_length_map.size() << std::endl;
   
}


std::vector<WireCell::PR3DCluster*> WireCell2dToy::Separate_2(WireCell::PR3DCluster* cluster, double dis_cut){
  
  std::map<int,SMGCSet>& time_cells_set_map = cluster->get_time_cells_set_map();
  SMGCSelection mcells = cluster->get_mcells();
  
  
  std::vector<int> time_slices;
  for (auto it1 = time_cells_set_map.begin(); it1!=time_cells_set_map.end(); it1++){
    time_slices.push_back((*it1).first);
  }
  
  std::vector<std::pair<SlimMergeGeomCell*,SlimMergeGeomCell*>> connected_mcells;
  for (size_t i=0; i!= time_slices.size(); i++){
    SMGCSet& mcells_set = time_cells_set_map[time_slices.at(i)];
    
    // create graph for points in mcell inside the same time slice
    if (mcells_set.size()>=2){
      for (auto it2 = mcells_set.begin(); it2!=mcells_set.end();it2++){
    	SlimMergeGeomCell *mcell1 = *it2;
    	auto it2p = it2;
    	if (it2p!=mcells_set.end()){
    	  it2p++;
    	  for (auto it3 = it2p; it3!=mcells_set.end(); it3++){
    	    SlimMergeGeomCell *mcell2 = *(it3);
	    if (mcell1->Overlap_fast(mcell2,5))
		connected_mcells.push_back(std::make_pair(mcell1,mcell2));
    	  }
    	}
      }
    }
    // create graph for points between connected mcells in adjacent time slices + 1, if not, + 2
    std::vector<SMGCSet> vec_mcells_set;
    if (i+1 < time_slices.size()){
      if (time_slices.at(i+1)-time_slices.at(i)==1){
    	vec_mcells_set.push_back(time_cells_set_map[time_slices.at(i+1)]);
    	if (i+2 < time_slices.size())
    	  if (time_slices.at(i+2)-time_slices.at(i)==2)
    	    vec_mcells_set.push_back(time_cells_set_map[time_slices.at(i+2)]);
      }else if (time_slices.at(i+1) - time_slices.at(i)==2){
    	vec_mcells_set.push_back(time_cells_set_map[time_slices.at(i+1)]);
      }
    }
    bool flag = false;
    for (size_t j=0; j!=vec_mcells_set.size(); j++){
      if (flag) break;
      SMGCSet& next_mcells_set = vec_mcells_set.at(j);
      for (auto it1 = mcells_set.begin(); it1!= mcells_set.end(); it1++){
    	SlimMergeGeomCell *mcell1 = (*it1);
    	for (auto it2 = next_mcells_set.begin(); it2!=next_mcells_set.end(); it2++){
    	  SlimMergeGeomCell *mcell2 = (*it2);
    	  if (mcell1->Overlap_fast(mcell2,2)){
    	    flag = true;
    	    connected_mcells.push_back(std::make_pair(mcell1,mcell2));
    	  }
    	}
      }
    }
  }

  // form ...
  
  const int N = mcells.size();
  MCUGraph *graph = new MCUGraph(N);
  
  std::map<SlimMergeGeomCell*, int> mcell_index_map;
  for (size_t i=0;i!=mcells.size();i++){
    SlimMergeGeomCell *curr_mcell = mcells.at(i);
    mcell_index_map[curr_mcell] = i;
    
    auto v = vertex(i, *graph); // retrieve vertex descriptor
    (*graph)[v].index = i;
  }
  
  for (auto it=connected_mcells.begin(); it!=connected_mcells.end(); it++){
    int index1 = mcell_index_map[it->first];
    int index2 = mcell_index_map[it->second];
    auto edge = add_edge(index1,index2,*graph);
    if (edge.second){
      (*graph)[edge.first].dist = 1;
    }
  }

  {
    std::vector<int> component(num_vertices(*graph));
    const int num = connected_components(*graph,&component[0]);

    if (num > 1){
      std::vector<ToyPointCloud*> pt_clouds;
      std::vector<std::vector<int>> vec_vec(num);
      for (int j=0;j!=num;j++){
	ToyPointCloud *pt_cloud = new ToyPointCloud();
	pt_clouds.push_back(pt_cloud);
      }
      std::vector<int>::size_type i;
      for (i=0;i!=component.size(); ++i){
	vec_vec.at(component[i]).push_back(i);
	SlimMergeGeomCell *mcell = mcells.at(i);
	pt_clouds.at(component[i])->AddPoints(mcell->get_sampling_points());
      }
      for (int j=0;j!=num;j++){
	pt_clouds.at(j)->build_kdtree_index();
      }
      
      for (int j=0;j!=num;j++){
	for (int k=j+1;k!=num;k++){
	  std::tuple<int,int,double> temp_results = pt_clouds.at(j)->get_closest_points(pt_clouds.at(k));
	  if (std::get<2>(temp_results)<dis_cut){
	    int index1 = vec_vec[j].front();
	    int index2 = vec_vec[k].front();
	    auto edge = add_edge(index1,index2,*graph);
	    if (edge.second){
	      (*graph)[edge.first].dist = 1;
	    }
	  }
	}
      }
      
      
      for (int j=0;j!=num;j++){
	delete pt_clouds.at(j);
      }
    }

    //std::cout << num << std::endl;
    
  }


  std::vector<WireCell::PR3DCluster*> final_clusters;
  {
    std::vector<int> component(num_vertices(*graph));
    const int num = connected_components(*graph,&component[0]);
    final_clusters.resize(num);
    for (size_t i=0;i!=num;i++){
      final_clusters.at(i) = new PR3DCluster(i);
    }
    
    std::vector<int>::size_type i;
    for (i=0;i!=component.size(); ++i){
      SlimMergeGeomCell *mcell = mcells.at(i);
      final_clusters[component[i]]->AddCell(mcell,mcell->GetTimeSlice());
    }
  }


  
  delete graph;
  
 
  /* int num_mcells = 0; */
  /* for (size_t i=0;i!=final_clusters.size();i++){ */
  /*   num_mcells += final_clusters.at(i)->get_num_mcells(); */
  /* } */
  /* std::cout << cluster->get_num_mcells() << " " << num_mcells << std::endl; */


  
  return final_clusters;
}
