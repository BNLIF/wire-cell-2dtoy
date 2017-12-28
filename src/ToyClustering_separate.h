void WireCell2dToy::Clustering_separate(WireCell::PR3DClusterSelection& live_clusters, std::map<WireCell::PR3DCluster*,double>& cluster_length_map){
  TVector3 drift_dir(1,0,0);

  for (size_t i=0;i!=live_clusters.size();i++){
    PR3DCluster *cluster = live_clusters.at(i);
    bool flag_sep = false;

    
    if (cluster_length_map[cluster]> 100*units::cm){
      // get the main axis
      cluster->Calc_PCA();

      
      
      TVector3 dir1(cluster->get_PCA_axis(0).x,cluster->get_PCA_axis(0).y,cluster->get_PCA_axis(0).z);
      TVector3 dir2(cluster->get_PCA_axis(1).x,cluster->get_PCA_axis(1).y,cluster->get_PCA_axis(1).z);
      TVector3 dir3(cluster->get_PCA_axis(2).x,cluster->get_PCA_axis(2).y,cluster->get_PCA_axis(2).z);

      double angle1 = fabs(dir2.Angle(drift_dir)-3.1415926/2.)/3.1415926*180.;
      double angle2 = fabs(dir3.Angle(drift_dir)-3.1415926/2.)/3.1415926*180.;
      double ratio1 = cluster->get_PCA_value(1)/cluster->get_PCA_value(0) ;
      double ratio2 = cluster->get_PCA_value(2)/cluster->get_PCA_value(0) ;
      if (ratio1 > pow(10,exp(1.38115-1.19312*pow(angle1,1./3.))-2.2) ||
	  ratio2 > pow(10,exp(1.38115-1.19312*pow(angle2,1./3.))-2.2))
	flag_sep = true;
      
    
      
      if (flag_sep) {
	cluster->Create_point_cloud();
	ToyPointCloud* cloud = cluster->get_point_cloud();
	std::vector<WCPointCloud<double>::WCPoint> boundary_points = cloud->get_hull();
	  std::vector<WCPointCloud<double>::WCPoint> hy_points;
	  std::vector<WCPointCloud<double>::WCPoint> ly_points;
	  std::vector<WCPointCloud<double>::WCPoint> hz_points;
	  std::vector<WCPointCloud<double>::WCPoint> lz_points;
	  std::vector<WCPointCloud<double>::WCPoint> hx_points;
	  std::vector<WCPointCloud<double>::WCPoint> lx_points;

	  for (size_t j=0;j!=boundary_points.size();j++){
	    if (j==0){
	      hy_points.push_back(boundary_points.at(j));
	      ly_points.push_back(boundary_points.at(j));
	      hz_points.push_back(boundary_points.at(j));
	      lz_points.push_back(boundary_points.at(j));
	      hx_points.push_back(boundary_points.at(j));
	      lx_points.push_back(boundary_points.at(j));
	    }else{
	      if (boundary_points.at(j).y > hy_points.at(0).y) hy_points.at(0) = boundary_points.at(j);
	      if (boundary_points.at(j).y < ly_points.at(0).y) ly_points.at(0) = boundary_points.at(j);
	      if (boundary_points.at(j).x > hx_points.at(0).x) hx_points.at(0) = boundary_points.at(j);
	      if (boundary_points.at(j).x < lx_points.at(0).x) lx_points.at(0) = boundary_points.at(j);
	      if (boundary_points.at(j).z > hz_points.at(0).z) hz_points.at(0) = boundary_points.at(j);
	      if (boundary_points.at(j).z < lz_points.at(0).z) lz_points.at(0) = boundary_points.at(j);
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
	     
	  if (ly_points.at(0).y < -101.5*units::cm){
	    for (size_t j=0;j!=boundary_points.size();j++){
	      if (boundary_points.at(j).y < -101.5*units::cm){
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
	   std::vector<WCPointCloud<double>::WCPoint> independent_points;
	   for (size_t j=0;j!=hy_points.size();j++){
	     bool flag_save = true;
	     for (size_t k=0;k!=independent_points.size();k++){
	       double dis = sqrt(pow(hy_points.at(j).x - independent_points.at(k).x,2)+pow(hy_points.at(j).y - independent_points.at(k).y,2)+pow(hy_points.at(j).z - independent_points.at(k).z,2));
	       if (dis < 30*units::cm)
		 flag_save = false;
	     }
	     if (flag_save){
	       independent_points.push_back(hy_points.at(j));
	       if (hy_points.at(j).y > 106.5*units::cm)
		 num_outside_points ++;
	     }
	   }
	   for (size_t j=0;j!=ly_points.size();j++){
	     bool flag_save = true;
	     for (size_t k=0;k!=independent_points.size();k++){
	       double dis = sqrt(pow(ly_points.at(j).x - independent_points.at(k).x,2)+pow(ly_points.at(j).y - independent_points.at(k).y,2)+pow(ly_points.at(j).z - independent_points.at(k).z,2));
	       if (dis < 30*units::cm)
		 flag_save = false;
	     }
	     if (flag_save){
	       independent_points.push_back(ly_points.at(j));
	       if (ly_points.at(j).y < -106.5*units::cm)
		 num_outside_points ++;
	     }
	   }
	   for (size_t j=0;j!=hz_points.size();j++){
	     bool flag_save = true;
	     for (size_t k=0;k!=independent_points.size();k++){
	       double dis = sqrt(pow(hz_points.at(j).x - independent_points.at(k).x,2)+pow(hz_points.at(j).y - independent_points.at(k).y,2)+pow(hz_points.at(j).z - independent_points.at(k).z,2));
	       if (dis < 30*units::cm)
		 flag_save = false;
	     }
	     if (flag_save){
	       independent_points.push_back(hz_points.at(j));
	       if (hz_points.at(j).z > 1027*units::cm)
		 num_outside_points ++;
	     }
	   }
	   for (size_t j=0;j!=lz_points.size();j++){
	     bool flag_save = true;
	     for (size_t k=0;k!=independent_points.size();k++){
	       double dis = sqrt(pow(lz_points.at(j).x - independent_points.at(k).x,2)+pow(lz_points.at(j).y - independent_points.at(k).y,2)+pow(lz_points.at(j).z - independent_points.at(k).z,2));
	       if (dis < 30*units::cm)
		 flag_save = false;
	     }
	     if (flag_save){
	       independent_points.push_back(lz_points.at(j));
	       if (lz_points.at(j).z < 10*units::cm)
		 num_outside_points ++;
	     }
	   }
	   for (size_t j=0;j!=hx_points.size();j++){
	     bool flag_save = true;
	     for (size_t k=0;k!=independent_points.size();k++){
	       double dis = sqrt(pow(hx_points.at(j).x - independent_points.at(k).x,2)+pow(hx_points.at(j).y - independent_points.at(k).y,2)+pow(hx_points.at(j).z - independent_points.at(k).z,2));
	       if (dis < 30*units::cm)
		 flag_save = false;
	     }
	     if (flag_save){
	       independent_points.push_back(hx_points.at(j));
	       if (hx_points.at(j).x > 257*units::cm)
		 num_outside_points ++;
	     }
	   }
	   for (size_t j=0;j!=lx_points.size();j++){
	     bool flag_save = true;
	     for (size_t k=0;k!=independent_points.size();k++){
	       double dis = sqrt(pow(lx_points.at(j).x - independent_points.at(k).x,2)+pow(lx_points.at(j).y - independent_points.at(k).y,2)+pow(lx_points.at(j).z - independent_points.at(k).z,2));
	       if (dis < 30*units::cm)
		 flag_save = false;
	     }
	     if (flag_save){
	       independent_points.push_back(lx_points.at(j));
	       if (lz_points.at(j).z < -1*units::cm)
		 num_outside_points ++;
	     }
	   }
	  
	   int num_far_points = 0;

	   
	   if (independent_points.size()==2){
	     TVector3 dir_1(independent_points.at(1).x - independent_points.at(0).x, independent_points.at(1).y - independent_points.at(0).y, independent_points.at(1).z - independent_points.at(0).z);
	     dir_1.SetMag(1);
	     for (size_t j=0;j!=boundary_points.size();j++){
		TVector3 dir_2(boundary_points.at(j).x - independent_points.at(0).x, boundary_points.at(j).y - independent_points.at(0).y, boundary_points.at(j).z - independent_points.at(0).z);
		double angle_12 = dir_1.Angle(dir_2);
		TVector3 dir_3 = dir_2 - dir_1 * dir_2.Mag() * cos(angle_12);
		double angle_3 = dir_3.Angle(drift_dir);
		if (fabs(angle_3-3.1415926/2.)/3.1415926*180.<7.5){
		  if (fabs(dir_3.X()/units::cm)>25*units::cm)
		    num_far_points ++;
		}else{
		  if (dir_3.Mag() > 40*units::cm)
		    num_far_points ++;
		}
	     }
	   }
	  
	   std::cout <<  cluster->get_cluster_id() << " " << hy_points.size() << " " << ly_points.size() << " " << hz_points.size() << " " << lz_points.size() <<  " " << hx_points.size() << " " << lz_points.size() << " " << num_outside_points << " " << independent_points.size() << " " << num_far_points << std::endl;

	   if (num_outside_points > 0 && (independent_points.size()>2 ||
					  independent_points.size()==2 && num_far_points > 0)){
	     std::cout << "Separate cluster " << cluster->get_cluster_id() << std::endl;
	     // separate ... 
	   }
	   /* std::cout << << " Xin " << cluster_length_map[cluster]/units::cm << " " << boundary_points.size() << " " << cluster->get_PCA_value(0) << " " << cluster->get_PCA_value(1)/cluster->get_PCA_value(0) << " " << cluster->get_PCA_value(2)/cluster->get_PCA_value(0) << " " << dir1.Angle(drift_dir) << " " << dir2.Angle(drift_dir) << " " << dir3.Angle(drift_dir) << std::endl;   */
      }
    }
  }

  
  
}
