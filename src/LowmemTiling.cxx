#include "WireCell2dToy/LowmemTiling.h"

using namespace WireCell;

WireCell2dToy::LowmemTiling::LowmemTiling(int time_slice, int nrebin, WireCell::GeomDataSource& gds,WireCell2dToy::WireCellHolder& holder)
  : gds(gds)
  , nrebin(nrebin)
  , time_slice(time_slice)
  , holder(holder)
{
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  
   
}

// bool WireCell2dToy::LowmemTiling::check_crossing(float dis_u, float u_pitch, float dis_v, float v_pitch,
// 						 float min_w, float max_w){
//   std::vector<Vector> psave(5);
  
//   for (int i=0;i!=5;i++){
//     bool flag;
//     if (i==0){
//       flag = gds.crossing_point(dis_u-u_pitch/2.,dis_v-v_pitch/2.,kUwire,kVwire, psave[0]);
//     }else if (i==1){
//       flag = gds.crossing_point(dis_u+u_pitch/2.,dis_v+v_pitch/2.,kUwire,kVwire, psave[0]);
//     }else if (i==2){
//       flag = gds.crossing_point(dis_u,dis_v,kUwire,kVwire, psave[0]);
//     }else if (i==3){
//       flag = gds.crossing_point(dis_u+u_pitch/2.,dis_v-v_pitch/2.,kUwire,kVwire, psave[0]);
//     }else if (i==4){
//       flag = gds.crossing_point(dis_u-u_pitch/2.,dis_v+v_pitch/2.,kUwire,kVwire, psave[0]);
//     }
 
//     if (flag){
//       if (psave[0].z > max_w && psave[0].z < max_w)
// 	return true;
//     }
//   }
  
//   return false;
// }

bool WireCell2dToy::LowmemTiling::check_crossing(const WireCell::GeomWire* wire1, const WireCell::GeomWire* wire2, float pitch1, float pitch2, WirePlaneType_t plane, float min, float max, float tolerance){
  
  float dis1 = gds.wire_dist(*wire1);
  float dis2 = gds.wire_dist(*wire2);  
  
  bool flag;
  float dis;
  std::vector<Vector> puv_save(5);
  
  flag = gds.crossing_point(dis1-pitch1/4.,dis2-pitch2/4.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1-pitch1/4.,dis2+pitch2/4.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1+pitch1/4.,dis2-pitch2/4.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1+pitch1/4.,dis2+pitch2/4.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1-pitch1/2.,dis2-pitch2/2.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1-pitch1/2.,dis2+pitch2/2.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1+pitch1/2.,dis2-pitch2/2.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1+pitch1/2.,dis2+pitch2/2.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1,dis2,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  
  return false;
}


void WireCell2dToy::LowmemTiling::test_crossing(PointVector& pcell, float dis_u, float dis_v, float bmin_w, float bmax_w, float u_pitch, float v_pitch, const GeomWire* uwire_1, const GeomWire* uwire_2, const GeomWire *vwire_1, const GeomWire *vwire_2, bool flag_test_u, bool flag_test_v){
  
  std::vector<Vector> puv_save(9);
  std::vector<Vector> pcross(9);
  float dis_u1, dis_v1;
  
  bool flag_qx = false;
  // min u and min v 
  bool flag[10];
  flag[0] = gds.crossing_point(dis_u-u_pitch/2.,dis_v-v_pitch/2.,kUwire,kVwire, puv_save[0]);
  flag[1] = gds.crossing_point(dis_u+u_pitch/2.,dis_v-v_pitch/2.,kUwire,kVwire, puv_save[1]);
  flag[2] = gds.crossing_point(dis_u-u_pitch/2.,dis_v+v_pitch/2.,kUwire,kVwire, puv_save[2]);
  flag[3] = gds.crossing_point(dis_u+u_pitch/2.,dis_v+v_pitch/2.,kUwire,kVwire, puv_save[3]);
  flag[4] = gds.crossing_point(dis_u,dis_v,kUwire,kVwire, puv_save[4]);
  // flag[5] = gds.crossing_point(dis_u-u_pitch/4.,dis_v-v_pitch/4.,kUwire,kVwire, puv_save[5]);
  // flag[6] = gds.crossing_point(dis_u+u_pitch/4.,dis_v-v_pitch/4.,kUwire,kVwire, puv_save[6]);
  // flag[7] = gds.crossing_point(dis_u-u_pitch/4.,dis_v+v_pitch/4.,kUwire,kVwire, puv_save[7]);
  // flag[8] = gds.crossing_point(dis_u+u_pitch/4.,dis_v+v_pitch/4.,kUwire,kVwire, puv_save[8]);
  
  float dis[9];
  for (int i=0;i!=5;i++){
    dis[i] = gds.wire_dist(puv_save[i],kYwire);
    if (flag[i] && dis[i] > bmin_w && dis[i] < bmax_w){
      pcell.push_back(puv_save[i]);
      flag_qx = true;
    }
  }
 

  
  if (!flag_qx){
    if (flag_test_u){
      // look at the U wire ...
      if (gds.crossing_point(bmin_w,dis_v,kYwire,kVwire,pcross[0])){
	// quick calculation
	const GeomWire *uwire_4 = gds.closest(pcross[0],WirePlaneType_t(0));
	for (int k=0;k!=5;k++){
	  int index = uwire_4->index()+(-2+k);
	  if (index < uwire_1->index() || 
	      index > uwire_2->index()) continue;
	  const GeomWire *uwire_3 = gds.by_planeindex(WirePlaneType_t(0), index);
	  dis_u1 = gds.wire_dist(*uwire_3);
	  flag[0] = gds.crossing_point(dis_u1-u_pitch/2.,dis_v-v_pitch/2.,kUwire,kVwire, puv_save[0]);
	  flag[1] = gds.crossing_point(dis_u1+u_pitch/2.,dis_v-v_pitch/2.,kUwire,kVwire, puv_save[1]);
	  flag[2] = gds.crossing_point(dis_u1-u_pitch/2.,dis_v+v_pitch/2.,kUwire,kVwire, puv_save[2]);
	  flag[3] = gds.crossing_point(dis_u1+u_pitch/2.,dis_v+v_pitch/2.,kUwire,kVwire, puv_save[3]);
	  flag[4] = gds.crossing_point(dis_u1,dis_v,kUwire,kVwire, puv_save[4]);
	  // flag[5] = gds.crossing_point(dis_u1-u_pitch/4.,dis_v-v_pitch/4.,kUwire,kVwire, puv_save[5]);
	  // flag[6] = gds.crossing_point(dis_u1+u_pitch/4.,dis_v-v_pitch/4.,kUwire,kVwire, puv_save[6]);
	  // flag[7] = gds.crossing_point(dis_u1-u_pitch/4.,dis_v+v_pitch/4.,kUwire,kVwire, puv_save[7]);
	  // flag[8] = gds.crossing_point(dis_u1+u_pitch/4.,dis_v+v_pitch/4.,kUwire,kVwire, puv_save[8]);
	  //flag_qx = false;
	  for (int i=0;i!=5;i++){
	    dis[i] = gds.wire_dist(puv_save[i],kYwire);
	    if (flag[i] && dis[i] > bmin_w && dis[i] < bmax_w){
	      pcell.push_back(puv_save[i]);
	      //  flag_qx = true;
	    }
	  }
	  
	  // if (flag_qx){
	  //   break;
	  // }
	}
      }
      
      if (gds.crossing_point(bmax_w,dis_v,kYwire,kVwire,pcross[0])){
	// quick calculation
	const GeomWire *uwire_4 = gds.closest(pcross[0],WirePlaneType_t(0));
	for (int k=0;k!=5;k++){
	  int index = uwire_4->index()+(-2+k);
	  if (index < uwire_1->index() || 
	      index > uwire_2->index()) continue;
	  const GeomWire *uwire_3 = gds.by_planeindex(WirePlaneType_t(0), index);
	  dis_u1 = gds.wire_dist(*uwire_3);
	  flag[0] = gds.crossing_point(dis_u1-u_pitch/2.,dis_v-v_pitch/2.,kUwire,kVwire, puv_save[0]);
	  flag[1] = gds.crossing_point(dis_u1+u_pitch/2.,dis_v-v_pitch/2.,kUwire,kVwire, puv_save[1]);
	  flag[2] = gds.crossing_point(dis_u1-u_pitch/2.,dis_v+v_pitch/2.,kUwire,kVwire, puv_save[2]);
	  flag[3] = gds.crossing_point(dis_u1+u_pitch/2.,dis_v+v_pitch/2.,kUwire,kVwire, puv_save[3]);
	  flag[4] = gds.crossing_point(dis_u1,dis_v,kUwire,kVwire, puv_save[4]);
	  // flag[5] = gds.crossing_point(dis_u1-u_pitch/4.,dis_v-v_pitch/4.,kUwire,kVwire, puv_save[5]);
	  // flag[6] = gds.crossing_point(dis_u1+u_pitch/4.,dis_v-v_pitch/4.,kUwire,kVwire, puv_save[6]);
	  // flag[7] = gds.crossing_point(dis_u1-u_pitch/4.,dis_v+v_pitch/4.,kUwire,kVwire, puv_save[7]);
	  // flag[8] = gds.crossing_point(dis_u1+u_pitch/4.,dis_v+v_pitch/4.,kUwire,kVwire, puv_save[8]);
	  //  flag_qx = false;
	  for (int i=0;i!=5;i++){
	    dis[i] = gds.wire_dist(puv_save[i],kYwire);
	    if (flag[i] && dis[i] > bmin_w && dis[i] < bmax_w){
	      pcell.push_back(puv_save[i]);
	      // flag_qx = true;
	    }
	  }

	  // if (flag_qx){
	  //   break;
	  // }
	}
      }
    }

    
    if (flag_test_v){
      // look at the V wire
      if (gds.crossing_point(bmin_w,dis_u,kYwire,kUwire,pcross[0])){
	const GeomWire *vwire_4 = gds.closest(pcross[0],WirePlaneType_t(1));
	for (int k=0;k!=5;k++){
	  int index = vwire_4->index()+(-2+k);
	  if (index < vwire_1->index() || 
	      index > vwire_2->index()) continue;
	  const GeomWire *vwire_3 = gds.by_planeindex(WirePlaneType_t(1), index);
	  dis_v1 = gds.wire_dist(*vwire_3);
	  flag[0] = gds.crossing_point(dis_u-u_pitch/2.,dis_v1-v_pitch/2.,kUwire,kVwire, puv_save[0]);
	  flag[1] = gds.crossing_point(dis_u+u_pitch/2.,dis_v1-v_pitch/2.,kUwire,kVwire, puv_save[1]);
	  flag[2] = gds.crossing_point(dis_u-u_pitch/2.,dis_v1+v_pitch/2.,kUwire,kVwire, puv_save[2]);
	  flag[3] = gds.crossing_point(dis_u+u_pitch/2.,dis_v1+v_pitch/2.,kUwire,kVwire, puv_save[3]);
	  flag[4] = gds.crossing_point(dis_u,dis_v1,kUwire,kVwire, puv_save[4]);
	  // flag[5] = gds.crossing_point(dis_u-u_pitch/4.,dis_v1-v_pitch/4.,kUwire,kVwire, puv_save[5]);
	  // flag[6] = gds.crossing_point(dis_u+u_pitch/4.,dis_v1-v_pitch/4.,kUwire,kVwire, puv_save[6]);
	  // flag[7] = gds.crossing_point(dis_u-u_pitch/4.,dis_v1+v_pitch/4.,kUwire,kVwire, puv_save[7]);
	  // flag[8] = gds.crossing_point(dis_u+u_pitch/4.,dis_v1+v_pitch/4.,kUwire,kVwire, puv_save[8]);
	  //flag_qx = false;
	  for (int i=0;i!=5;i++){
	    dis[i] = gds.wire_dist(puv_save[i],kYwire);
	    if (flag[i] && dis[i] > bmin_w && dis[i] < bmax_w ){
	      pcell.push_back(puv_save[i]);
	      //  flag_qx = true;
	    }
	  }
	  
	  // if (flag_qx){
	  //   break;
	  // }
	}
      }
    
      if (gds.crossing_point(bmax_w,dis_u,kYwire,kUwire,pcross[0])){
	const GeomWire *vwire_4 = gds.closest(pcross[0],WirePlaneType_t(1));
	for (int k=0;k!=5;k++){
	  int index = vwire_4->index()+(-2+k);
	  if (index < vwire_1->index() || 
	      index > vwire_2->index()) continue;
	  const GeomWire *vwire_3 = gds.by_planeindex(WirePlaneType_t(1), index);
	  dis_v1 = gds.wire_dist(*vwire_3);
	  flag[0] = gds.crossing_point(dis_u-u_pitch/2.,dis_v1-v_pitch/2.,kUwire,kVwire, puv_save[0]);
	  flag[1] = gds.crossing_point(dis_u+u_pitch/2.,dis_v1-v_pitch/2.,kUwire,kVwire, puv_save[1]);
	  flag[2] = gds.crossing_point(dis_u-u_pitch/2.,dis_v1+v_pitch/2.,kUwire,kVwire, puv_save[2]);
	  flag[3] = gds.crossing_point(dis_u+u_pitch/2.,dis_v1+v_pitch/2.,kUwire,kVwire, puv_save[3]);
	  flag[4] = gds.crossing_point(dis_u,dis_v1,kUwire,kVwire, puv_save[4]);
	  // flag[5] = gds.crossing_point(dis_u-u_pitch/4.,dis_v1-v_pitch/4.,kUwire,kVwire, puv_save[5]);
	  // flag[6] = gds.crossing_point(dis_u+u_pitch/4.,dis_v1-v_pitch/4.,kUwire,kVwire, puv_save[6]);
	  // flag[7] = gds.crossing_point(dis_u-u_pitch/4.,dis_v1+v_pitch/4.,kUwire,kVwire, puv_save[7]);
	  // flag[8] = gds.crossing_point(dis_u+u_pitch/4.,dis_v1+v_pitch/4.,kUwire,kVwire, puv_save[8]);
	  //flag_qx = false;
	  for (int i=0;i!=5;i++){
	    dis[i] = gds.wire_dist(puv_save[i],kYwire);
	    if (flag[i] && dis[i] > bmin_w && dis[i] < bmax_w){
	      pcell.push_back(puv_save[i]);
	      //flag_qx = true;
	    }
	  }
	  

	  // if (flag_qx){
	  //   break;
	  // }
	}
      }
    }
  }
  

}


WireCell::SlimMergeGeomCell* WireCell2dToy::LowmemTiling::create_slim_merge_cell(WireCell::MergeGeomWire *uwire, WireCell::MergeGeomWire *vwire, WireCell::MergeGeomWire *wwire){
  float u_pitch = gds.pitch(kUwire);
  float v_pitch = gds.pitch(kVwire);
  float w_pitch = gds.pitch(kYwire);
  
  float tolerance = 0.1 * units::mm / 2.;
  
  // find U plane wires
  float dis_u[3];
  const GeomWire *uwire_1 = uwire->get_allwire().front();
  const GeomWire *uwire_2 = uwire->get_allwire().back();
  dis_u[0] = gds.wire_dist(*uwire_1);
  dis_u[1] = gds.wire_dist(*uwire_2);
  // float bmin_u = dis_u[0] - u_pitch/2. - tolerance;
  // float bmax_u = dis_u[1] + u_pitch/2. + tolerance;

  // find V plane wires
  float dis_v[3];
  const GeomWire *vwire_1 = vwire->get_allwire().front();
  const GeomWire *vwire_2 = vwire->get_allwire().back();
  dis_v[0] = gds.wire_dist(*vwire_1);
  dis_v[1] = gds.wire_dist(*vwire_2);
  // float bmin_v = dis_v[0] - v_pitch/2. - tolerance;
  // float bmax_v = dis_v[1] + v_pitch/2. + tolerance;

  // define W plane range
  const GeomWire *wwire_1 = wwire->get_allwire().front();
  const GeomWire *wwire_2 = wwire->get_allwire().back();
  float bmin_w = gds.wire_dist(*wwire_1) - w_pitch/2. - tolerance;
  float bmax_w = gds.wire_dist(*wwire_2) + w_pitch/2. + tolerance;

  
  

  PointVector pcell;
  test_crossing(pcell, dis_u[0], dis_v[0], bmin_w, bmax_w, u_pitch, v_pitch, uwire_1, uwire_2, vwire_1, vwire_2,true,true);
  test_crossing(pcell, dis_u[0], dis_v[1], bmin_w, bmax_w, u_pitch, v_pitch, uwire_1, uwire_2, vwire_1, vwire_2,true,true);
  test_crossing(pcell, dis_u[1], dis_v[0], bmin_w, bmax_w, u_pitch, v_pitch, uwire_1, uwire_2, vwire_1, vwire_2,true,true);
  test_crossing(pcell, dis_u[1], dis_v[1], bmin_w, bmax_w, u_pitch, v_pitch, uwire_1, uwire_2, vwire_1, vwire_2,true,true);
  
  
  
  // find the max and min
  if (pcell.size() >=1){
    //Creat a cell and then get all the wires in ... 
    double u_max = -1e9, u_min = 1e9;
    double v_max = -1e9, v_min = 1e9;
    double w_max = -1e9, w_min = 1e9;
    // tolerance = 0.;
    for (int k=0;k!=pcell.size();k++){
      double udist = gds.wire_dist(pcell.at(k),WirePlaneType_t(0));
      if (udist + tolerance> u_max) u_max = udist + tolerance;
      if (udist - tolerance< u_min) u_min = udist - tolerance;
      double vdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(1));
      if (vdist + tolerance> v_max) v_max = vdist + tolerance;
      if (vdist - tolerance< v_min) v_min = vdist + tolerance;
      double wdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(2));
      if (wdist + tolerance> w_max) w_max = wdist + tolerance;
      if (wdist - tolerance< w_min) w_min = wdist - tolerance;
    }
    
    int ident = holder.get_ncell();
   
    

    const GeomWire* uwire_min = gds.closest(u_min,WirePlaneType_t(0));
    if (uwire_min->index() < uwire_1->index()){
      uwire_min = uwire_1;
    }else if (uwire_min->index() > uwire_2->index()){
      uwire_min = uwire_2;
    }
    const GeomWire* uwire_max = gds.closest(u_max,WirePlaneType_t(0));
    if (uwire_max->index() < uwire_1->index()){
      uwire_max = uwire_1;
    }else if (uwire_max->index() > uwire_2->index()){
      uwire_max = uwire_2;
    }
    const GeomWire* vwire_min = gds.closest(v_min,WirePlaneType_t(1));
    if (vwire_min->index() < vwire_1->index()){
      vwire_min = vwire_1;
    }else if (vwire_min->index() > vwire_2->index()){
      vwire_min = vwire_2;
    }
    const GeomWire* vwire_max = gds.closest(v_max,WirePlaneType_t(1));
    if (vwire_max->index() < vwire_1->index()){
      vwire_max = vwire_1;
    }else if (vwire_max->index() > vwire_2->index()){
      vwire_max = vwire_2;
    }
    const GeomWire* wwire_min = gds.closest(w_min,WirePlaneType_t(2));//.bounds(w_min,WirePlaneType_t(2)).second;
    if (wwire_min->index() < wwire_1->index()){
      wwire_min = wwire_1;
    }else if (wwire_min->index() > wwire_2->index()){
      wwire_min = wwire_2;
    }
    const GeomWire* wwire_max = gds.closest(w_max,WirePlaneType_t(2));//.bounds(w_max,WirePlaneType_t(2)).first;
    if (wwire_max->index() < wwire_1->index()){
      wwire_max = wwire_1;
    }else if (wwire_max->index() > wwire_2->index()){
      wwire_max = wwire_2;
    }
    u_min = gds.wire_dist(*uwire_min)-u_pitch/2.;
    u_max = gds.wire_dist(*uwire_max)+u_pitch/2.;
    v_min = gds.wire_dist(*vwire_min)-v_pitch/2.;
    v_max = gds.wire_dist(*vwire_max)+v_pitch/2.;
    w_min = gds.wire_dist(*wwire_min)-w_pitch/2.;
    w_max = gds.wire_dist(*wwire_max)+w_pitch/2.;
    
    // std::cout << uwire_min->index() << " " << uwire_max->index() << " " 
    // 	      << vwire_min->index() << " " << vwire_max->index() << " " 
    // 	      << wwire_min->index() << " " << wwire_max->index() << " " 
    // 	      << vwire_1->index() << " " << vwire_2->index() << " " 
    // 	      << w_min/units::m << " " << w_max / units::m << " " 
    // 	      << std::endl;
    
    
    GeomWireSelection ugroup,vgroup,wgroup;
    // find U group
    if (uwire_min!=0)
      ugroup.push_back(uwire_min);
    if (uwire_max!=uwire_min&&uwire_max!=0){
      ugroup.push_back(uwire_max);
      const GeomWire* uwire_p1 = gds.by_planeindex(WirePlaneType_t(0),uwire_min->index()+1);
      if (uwire_p1!=uwire_max&&uwire_p1!=0){
	ugroup.push_back(uwire_p1);
	const GeomWire* uwire_n1 = gds.by_planeindex(WirePlaneType_t(0),uwire_max->index()-1);
	if (uwire_n1!=uwire_p1&&uwire_n1!=0){
	  ugroup.push_back(uwire_n1);
	}
      }
    }
    // find V group
    if (vwire_min!=0)
      vgroup.push_back(vwire_min);
    if (vwire_max!=vwire_min&&vwire_max!=0){
      vgroup.push_back(vwire_max);
      const GeomWire* vwire_p1 = gds.by_planeindex(WirePlaneType_t(1),vwire_min->index()+1);
      if (vwire_p1!=vwire_max&&vwire_p1!=0){
	vgroup.push_back(vwire_p1);
	const GeomWire* vwire_n1 = gds.by_planeindex(WirePlaneType_t(1),vwire_max->index()-1);
	if (vwire_n1!=vwire_p1&&vwire_n1!=0){
	  vgroup.push_back(vwire_n1);
	}
      }
    }
    // find W group
    if (wwire_min!=0)
      wgroup.push_back(wwire_min);
    if (wwire_max!=wwire_min&&wwire_max!=0){
      wgroup.push_back(wwire_max);
      const GeomWire* wwire_p1 = gds.by_planeindex(WirePlaneType_t(2),wwire_min->index()+1);
      if (wwire_p1!=wwire_max&&wwire_p1!=0){
	wgroup.push_back(wwire_p1);
	const GeomWire* wwire_n1 = gds.by_planeindex(WirePlaneType_t(2),wwire_max->index()-1);
	if (wwire_n1!=wwire_p1&&wwire_n1!=0){
	  wgroup.push_back(wwire_n1);
	}
      }
    }
    // initialize map
    std::map<const GeomWire*, int> wiremap;
    for (auto it = ugroup.begin();it!=ugroup.end();it++){
      wiremap[*it]=0;
    }
    for (auto it = vgroup.begin();it!=vgroup.end();it++){
      wiremap[*it]=0;
    }
    for (auto it = wgroup.begin();it!=wgroup.end();it++){
      wiremap[*it]=0;
    }
    
    // std::cout << ugroup.size() << " " << vgroup.size() << " " << wgroup.size() << " " 
    // 	      << uwire_max->index()-uwire_min->index()+1 << " " 
    // 	      << vwire_max->index()-vwire_min->index()+1 << " " 
    // 	      << wwire_max->index()-wwire_min->index()+1 << " " 
    // 	      << std::endl;
    
    // check crossing 
    for (int i=0;i!=ugroup.size();i++){
      const GeomWire *uwire = ugroup.at(i);
      for (int j=0;j!=vgroup.size();j++){
	const GeomWire *vwire = vgroup.at(j);
	if (check_crossing(uwire, vwire, u_pitch, v_pitch, WirePlaneType_t(2), w_min, w_max, tolerance)){
	  wiremap[uwire]++;
	  wiremap[vwire]++;
	}
      }
    }
    for (int i=0;i!=wgroup.size();i++){
      const GeomWire *wwire = wgroup.at(i);
      for (int j=0;j!=vgroup.size();j++){
	const GeomWire *vwire = vgroup.at(j);
	if (check_crossing(wwire, vwire, w_pitch, v_pitch, WirePlaneType_t(0), u_min, u_max, tolerance)){
	  wiremap[wwire]++;
	  wiremap[vwire]++;
	}
      }
    }
    for (int i=0;i!=ugroup.size();i++){
      const GeomWire *uwire = ugroup.at(i);
      for (int j=0;j!=wgroup.size();j++){
	const GeomWire *wwire = wgroup.at(j);
	if (check_crossing(uwire, wwire, u_pitch, w_pitch, WirePlaneType_t(1), v_min, v_max, tolerance)){
	  wiremap[uwire]++;
	  wiremap[wwire]++;
	}
      }
    }

    
    // 	
    // 	for (int k=0;k!=wgroup.size();k++){
    // 	  const GeomWire *wwire = wgroup.at(k);
    // 	  if (check_crossing(uwire,vwire,wwire,u_pitch,v_pitch,w_pitch,tolerance)){
    // 	    wiremap[uwire]++;
    // 	    wiremap[vwire]++;
    // 	    wiremap[wwire]++;
    // 	  }
    // 	}
    //   }
    //   //      std::cout << wiremap[uwire] << " " << u_pitch << " " << tolerance << std::endl;
    // }

    int flag = 1;
    if (wiremap[uwire_min]==0 ){
      const GeomWire* wire_temp = gds.by_planeindex(WirePlaneType_t(0),uwire_min->index()+1);
      if (wiremap.find(wire_temp)!=wiremap.end() && wiremap[wire_temp]>0){ 
	uwire_min = wire_temp;
      }else{
	flag = 0;
      }
    }
    if (flag==1){
      if (wiremap[vwire_min]==0 ){
	const GeomWire* wire_temp = gds.by_planeindex(WirePlaneType_t(1),vwire_min->index()+1);
	if (wiremap.find(wire_temp)!=wiremap.end() && wiremap[wire_temp]>0){ 
	  vwire_min = wire_temp;
	}else{
	  flag = 0;
	}
      }
    }
    if (flag==1){
      if (wiremap[wwire_min]==0 ){
	const GeomWire* wire_temp = gds.by_planeindex(WirePlaneType_t(2),wwire_min->index()+1);
	if (wiremap.find(wire_temp)!=wiremap.end() && wiremap[wire_temp]>0){ 
	  wwire_min = wire_temp;
	}else{
	  flag = 0;
	}
      }
    }
    if (flag==1){
      if (wiremap[uwire_max]==0 ){
	const GeomWire* wire_temp = gds.by_planeindex(WirePlaneType_t(0),uwire_max->index()-1);
	if (wiremap.find(wire_temp)!=wiremap.end() && wiremap[wire_temp]>0){ 
	  uwire_max = wire_temp;
	}else{
	  flag = 0;
	}
      }
    }
    if (flag==1){
      if (wiremap[vwire_max]==0 ){
	const GeomWire* wire_temp = gds.by_planeindex(WirePlaneType_t(1),vwire_max->index()-1);
	if (wiremap.find(wire_temp)!=wiremap.end() && wiremap[wire_temp]>0){ 
	  vwire_max = wire_temp;
	}else{
	  flag = 0;
	}
      }
    }
    if (flag==1){
      if (wiremap[wwire_max]==0 ){
	const GeomWire* wire_temp = gds.by_planeindex(WirePlaneType_t(2),wwire_max->index()-1);
	if (wiremap.find(wire_temp)!=wiremap.end() && wiremap[wire_temp]>0){ 
	  wwire_max = wire_temp;
	}else{
	  flag = 0;
	}
      }
    }

   
   
    if (flag==1){
      SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident); 
      holder.AddCell(mcell);
      // Inser U    
      for (int k=uwire_min->index();k!=uwire_max->index()+1;k++){
	const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),k);
	mcell->AddWire(uwire,WirePlaneType_t(0));
      }
      //Insert V
      for (int k=vwire_min->index();k!=vwire_max->index()+1;k++){
	const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),k);
	mcell->AddWire(vwire,WirePlaneType_t(1));
      }
      // Insert W
      for (int k=wwire_min->index();k!=wwire_max->index()+1;k++){
	const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),k);
	mcell->AddWire(wwire,WirePlaneType_t(2));
      }
      
      return mcell;
    }
  }
  
  return 0;
}


void WireCell2dToy::LowmemTiling::init_good_cells(const WireCell::Slice& slice, std::vector<float>& uplane_rms, std::vector<float>& vplane_rms, std::vector<float>& wplane_rms){
  // form good wires group
  form_fired_merge_wires(slice);
  
  // create three good wire cells & two good wire + one bad wire cells
  // U/V/W = 1/1/1
  for (int i=0;i!=fired_wire_u.size();i++){
    MergeGeomWire *uwire = (MergeGeomWire *)fired_wire_u.at(i);
    for (int j=0;j!=fired_wire_v.size();j++){
      MergeGeomWire *vwire = (MergeGeomWire *)fired_wire_v.at(j);
      for (int k=0;k!=fired_wire_w.size();k++){
	MergeGeomWire *wwire = (MergeGeomWire *)fired_wire_w.at(k);
	SlimMergeGeomCell *mcell = create_slim_merge_cell(uwire,vwire,wwire);
	
	if (mcell !=0) three_good_wire_cells.push_back(mcell);
	
      }
    }
  }
  
  // //U/V/W = 1/1/0
  // for (int i=0;i!=fired_wire_u.size();i++){
  //   MergeGeomWire *uwire = (MergeGeomWire *)fired_wire_u.at(i);
  //   for (int j=0;j!=fired_wire_v.size();j++){
  //     MergeGeomWire *vwire = (MergeGeomWire *)fired_wire_v.at(j);
  //     for (int k=0;k!=bad_wire_w.size();k++){
  // 	MergeGeomWire *wwire = (MergeGeomWire *)bad_wire_w.at(k);
	
  // 	SlimMergeGeomCell *mcell = create_slim_merge_cell(uwire,vwire,wwire);
  // 	if (mcell !=0) two_good_wire_cells.push_back(mcell);
  //     }
  //   }
  // }
  
  // //U/V/W = 1/0/1
  // for (int i=0;i!=fired_wire_u.size();i++){
  //   MergeGeomWire *uwire = (MergeGeomWire *)fired_wire_u.at(i);
  //   for (int j=0;j!=bad_wire_v.size();j++){
  //     MergeGeomWire *vwire = (MergeGeomWire *)bad_wire_v.at(j);
  //     for (int k=0;k!=fired_wire_w.size();k++){
  // 	MergeGeomWire *wwire = (MergeGeomWire *)fired_wire_w.at(k);
	
  // 	SlimMergeGeomCell *mcell = create_slim_merge_cell(uwire,vwire,wwire);
  // 	if (mcell !=0) two_good_wire_cells.push_back(mcell);
  //     }
  //   }
  // }
  
  // //U/V/W = 0/1/1
  // for (int i=0;i!=bad_wire_u.size();i++){
  //   MergeGeomWire *uwire = (MergeGeomWire *)bad_wire_u.at(i);
  //   for (int j=0;j!=fired_wire_v.size();j++){
  //     MergeGeomWire *vwire = (MergeGeomWire *)fired_wire_v.at(j);
  //     for (int k=0;k!=fired_wire_w.size();k++){
  // 	MergeGeomWire *wwire = (MergeGeomWire *)fired_wire_w.at(k);
	
  // 	SlimMergeGeomCell *mcell = create_slim_merge_cell(uwire,vwire,wwire);
  // 	if (mcell !=0) two_good_wire_cells.push_back(mcell);
  //     }
  //   }
  // }
  

}




void WireCell2dToy::LowmemTiling::check_bad_cells(WireCell2dToy::LowmemTiling* tiling,WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map){
  
  int prev_time_slice = tiling->get_time_slice();

  int flag = 0;
  for (auto it = uplane_map.begin(); it!= uplane_map.end(); it++){
    int start_time_slice = it->second.first/ nrebin;
    int end_time_slice = it->second.second/ nrebin;
    if (prev_time_slice < start_time_slice && time_slice >= start_time_slice ||
	prev_time_slice <= end_time_slice && time_slice > end_time_slice )
      {
	flag = 1;
	break;
      }
    //    std::cout << start_time_slice << " " << end_time_slice << std::endl;
  }
  if (flag==0){
    for (auto it = vplane_map.begin(); it!= vplane_map.end(); it++){
      int start_time_slice = it->second.first/ nrebin;
      int end_time_slice = it->second.second/ nrebin;
      if (prev_time_slice < start_time_slice && time_slice >= start_time_slice ||
	prev_time_slice <= end_time_slice && time_slice > end_time_slice
	  ){
	flag = 1;
	break;
      }
      //    std::cout << start_time_slice << " " << end_time_slice << std::endl;
    }
  }
  if (flag==0){
    for (auto it = wplane_map.begin(); it!= wplane_map.end(); it++){
      int start_time_slice = it->second.first/ nrebin;
      int end_time_slice = it->second.second/ nrebin;
      if (prev_time_slice < start_time_slice && time_slice >= start_time_slice ||
	prev_time_slice <= end_time_slice && time_slice > end_time_slice
	  ){
	flag = 1;
	break;
      }
      //    std::cout << start_time_slice << " " << end_time_slice << std::endl;
    }
  }

  

  if (flag==1){
    std::cout << "Regenerate bad cells at time slice " << time_slice << std::endl;
    init_bad_cells(uplane_map,vplane_map,wplane_map);
  }else{
    //copy wires 
    for (int i=0;i!=tiling->get_bad_wire_u().size();i++){
      MergeGeomWire *mwire1 = (MergeGeomWire*)tiling->get_bad_wire_u().at(i);
      bad_wire_u.push_back(mwire1);
    }
    for (int i=0;i!=tiling->get_bad_wire_v().size();i++){
      MergeGeomWire *mwire1 = (MergeGeomWire*)tiling->get_bad_wire_v().at(i);
      bad_wire_v.push_back(mwire1);
    }
    for (int i=0;i!=tiling->get_bad_wire_w().size();i++){
      MergeGeomWire *mwire1 = (MergeGeomWire*)tiling->get_bad_wire_w().at(i);
      bad_wire_w.push_back(mwire1);
    }
    // copy cells;
    for (int i=0;i!=tiling->get_two_bad_wire_cells().size();i++){
      SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*)tiling->get_two_bad_wire_cells().at(i);
      two_bad_wire_cells.push_back(mcell1);
    }
  }
  
  // std::cout << bad_wire_u.size() << " " << bad_wire_v.size() << " " << bad_wire_w.size() << std::endl;

  
}



void WireCell2dToy::LowmemTiling::init_bad_cells(WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map){
  // form bad wires group
  form_bad_merge_wires(uplane_map, vplane_map, wplane_map);
  // create two bad wire cells // these are special ones ... 
  form_two_bad_cells();

  //  std::cout << two_bad_wire_cells.size() << std::endl; 
	    
}



void WireCell2dToy::LowmemTiling::form_two_bad_cells(){
  // form two bad wire cells ... taken from BadTiling ...  
  
  // U-V and insert Y ... 
  for (int i =0; i!=bad_wire_u.size();i++){
    const GeomWire *uwire_1 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().front();
    const GeomWire *uwire_2 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().back();
    float dis_u[3];
    float u_pitch = gds.pitch(kUwire);
    dis_u[0] = gds.wire_dist(*uwire_1) - u_pitch/2.;
    dis_u[1] = gds.wire_dist(*uwire_2) + u_pitch/2.;
    for (int j=0; j!=bad_wire_v.size();j++){
      const GeomWire *vwire_1 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().front();
      const GeomWire *vwire_2 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().back();
      float dis_v[3];
      float v_pitch = gds.pitch(kVwire);
      dis_v[0] = gds.wire_dist(*vwire_1) - v_pitch/2.;
      dis_v[1] = gds.wire_dist(*vwire_2) + v_pitch/2.;
      
      PointVector pcell;
      std::vector<Vector> pcross(4);

      bool flag1 = gds.crossing_point(dis_u[0],dis_v[0],kUwire,kVwire, pcross[0]); // check the inner point
      
      if (flag1){
	// fill the outer point
  	pcell.push_back(pcross[0]);
      }else{
  	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) - u_pitch/2.;

	  if (gds.crossing_point(dis_u[2],dis_v[0],kUwire,kVwire, pcross[0])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(k);
	  dis_v[2] = gds.wire_dist(*vwire_3) - v_pitch/2.;
	  if (gds.crossing_point(dis_u[0],dis_v[2],kUwire,kVwire, pcross[0])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
      }

      bool flag2 = gds.crossing_point(dis_u[0],dis_v[1],kUwire,kVwire, pcross[1]); 
      if (flag2){
  	pcell.push_back(pcross[1]);
      }else{
  	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) - u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_v[1],kUwire,kVwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size()-1-k);
	  dis_v[2] = gds.wire_dist(*vwire_3) + v_pitch/2.;
	  if (gds.crossing_point(dis_u[0],dis_v[2],kUwire,kVwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
      }
      
      bool flag3 = gds.crossing_point(dis_u[1],dis_v[0],kUwire,kVwire, pcross[2]); 
      if (flag3){
  	pcell.push_back(pcross[2]);
      }else{
  	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size()-1-k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) + u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_v[0],kUwire,kVwire, pcross[2])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(k);
	  dis_v[2] = gds.wire_dist(*vwire_3) - v_pitch/2.;
	  if (gds.crossing_point(dis_u[1],dis_v[2],kUwire,kVwire, pcross[2])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
      }
      
      bool flag4 = gds.crossing_point(dis_u[1],dis_v[1],kUwire,kVwire, pcross[3]);
      if (flag4){
  	pcell.push_back(pcross[3]);
      }else{
	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size()-1-k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) + u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_v[1],kUwire,kVwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[3]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size()-1-k);
	  dis_v[2] = gds.wire_dist(*vwire_3) + v_pitch/2.;
	  if (gds.crossing_point(dis_u[1],dis_v[2],kUwire,kVwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	    pcell.push_back(pcross[3]);
	    break;
	  }
	}
      }

      if (pcell.size() >=3){
	//Creat a cell and then get all the wires in ... 
	double u_max = -1e9, u_min = 1e9;
	double v_max = -1e9, v_min = 1e9;
	double w_max = -1e9, w_min = 1e9;
	for (int k=0;k!=pcell.size();k++){
	  double udist = gds.wire_dist(pcell.at(k),WirePlaneType_t(0));
	  if (udist > u_max) u_max = udist;
	  if (udist < u_min) u_min = udist;
	  double vdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(1));
	  if (vdist > v_max) v_max = vdist;
	  if (vdist < v_min) v_min = vdist;
	  double wdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(2));
	  if (wdist > w_max) w_max = wdist;
	  if (wdist < w_min) w_min = wdist;
	}
	//std::cout << u_max << " " << u_min << " " << v_max << " " << v_min << " " << w_max << " " << w_min << std::endl;
	//Create a cell
	int ident = holder.get_ncell();
	SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident); 
	mcell->AddBoundary(pcell);
	two_bad_wire_cells.push_back(mcell);
	holder.AddCell(mcell);
	
	//Insert U
	// for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	//   const GeomWire *uwire = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(k);
	//   double udist = gds.wire_dist(*uwire);
	//   if (udist>u_min && udist < u_max)
	//     mcell->AddWire(uwire,WirePlaneType_t(0));
	// }

	// for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	//   const GeomWire *vwire = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(k);
	//   double vdist = gds.wire_dist(*vwire);
	//   if (vdist>v_min && vdist < v_max)
	//     mcell->AddWire(vwire,WirePlaneType_t(1));
	// }
	// //Insert V
	// for (int k=0;k!=gds.wires_in_plane(WirePlaneType_t(2)).size();k++){
	//   const GeomWire *wwire = gds.wires_in_plane(WirePlaneType_t(2)).at(k);
	//   double wdist = gds.wire_dist(*wwire);
	//   if (wdist>w_min && wdist < w_max)
	//     mcell->AddWire(wwire,WirePlaneType_t(2));
	// }
	
	const GeomWire* uwire_min = gds.bounds(u_min,WirePlaneType_t(0)).second;
	const GeomWire* uwire_max = gds.bounds(u_max,WirePlaneType_t(0)).first;
	for (int k=uwire_min->index();k!=uwire_max->index()+1;k++){
	  const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),k);
	  mcell->AddWire(uwire,WirePlaneType_t(0));
	}
	
	//Insert V
	const GeomWire* vwire_min = gds.bounds(v_min,WirePlaneType_t(1)).second;
	const GeomWire* vwire_max = gds.bounds(v_max,WirePlaneType_t(1)).first;
	// std::cout << vwire_min->index() << " " << vwire_max->index() << std::endl;
	for (int k=vwire_min->index();k!=vwire_max->index()+1;k++){
	  const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),k);
	  mcell->AddWire(vwire,WirePlaneType_t(1));
	}

	//Insert W
	const GeomWire* wwire_min = gds.closest(w_min,WirePlaneType_t(2));//.bounds(w_min,WirePlaneType_t(2)).second;
	const GeomWire* wwire_max = gds.closest(w_max,WirePlaneType_t(2));//.bounds(w_max,WirePlaneType_t(2)).first;
	for (int k=wwire_min->index();k!=wwire_max->index()+1;k++){
	  const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),k);
	  mcell->AddWire(wwire,WirePlaneType_t(2));
	}
	
	//	std::cout << mcell->get_uwires().size() << " " << mcell->get_vwires().size() << " " << mcell->get_wwires().size() << std::endl;
      }
    }
  }


  // U-W and insert V
  for ( int i = 0; i != bad_wire_u.size() ; i++){
    const GeomWire *uwire_1 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().front();
    const GeomWire *uwire_2 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().back();
    float dis_u[2];
    float u_pitch = gds.pitch(kUwire);
    dis_u[0] = gds.wire_dist(*uwire_1) - u_pitch/2.;
    dis_u[1] = gds.wire_dist(*uwire_2) + u_pitch/2.;
    

    for (int j = 0; j != bad_wire_w.size() ; j++){
      const GeomWire *wwire_1 = ((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().front();
      const GeomWire *wwire_2 = ((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().back();
      float dis_w[2];
      float w_pitch = gds.pitch(kYwire);
      dis_w[0] = gds.wire_dist(*wwire_1) - w_pitch/2.;
      dis_w[1] = gds.wire_dist(*wwire_2) + w_pitch/2.;

      //      std::cout << dis_u[0]/units::m << " " << dis_u[1]/units::m << " " << dis_w[0]/units::m << " " << dis_w[1]/units::m << std::endl;

      PointVector pcell;
      std::vector<Vector> pcross(4);
      
      
      bool flag1 = gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, pcross[0]); 
      if (flag1){
  	pcell.push_back(pcross[0]);
      }else{
  	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) - u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_w[0],kUwire,kYwire, pcross[0])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
	// scan w-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().at(k);
	  dis_w[2] = gds.wire_dist(*wwire_3) - w_pitch/2.;
	  if (gds.crossing_point(dis_u[0],dis_w[2],kUwire,kYwire, pcross[0])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
      }

      bool flag2 = gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, pcross[1]); 
      if (flag2){
  	pcell.push_back(pcross[1]);
      }else{
  	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) - u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_w[1],kUwire,kYwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
	// scan w-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().at(((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().size()-1-k);
	  dis_w[2] = gds.wire_dist(*wwire_3) + w_pitch/2.;
	  if (gds.crossing_point(dis_u[0],dis_w[2],kUwire,kYwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
      }
      
      bool flag3 = gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, pcross[2]); 
      if (flag3){
  	pcell.push_back(pcross[2]);
      }else{
  	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size()-1-k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) + u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_w[0],kUwire,kYwire, pcross[2])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
	// scan w-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().at(k);
	  dis_w[2] = gds.wire_dist(*wwire_3) - w_pitch/2.;
	  if (gds.crossing_point(dis_u[1],dis_w[2],kUwire,kYwire, pcross[2])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
      }
      
      bool flag4 = gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, pcross[3]);
      if (flag4){
  	pcell.push_back(pcross[3]);
      }else{
	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size()-1-k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) + u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_w[1],kUwire,kYwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[3]);
	    break;
	  }
	}
	// scan w-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().at(((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().size()-1-k);
	  dis_w[2] = gds.wire_dist(*wwire_3) + w_pitch/2.;
	  if (gds.crossing_point(dis_u[1],dis_w[2],kUwire,kYwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	    pcell.push_back(pcross[3]);
	    break;
	  }
	}
      }
     
      
    
      if(pcell.size() >=3){
	//Creat a cell and then get all the wires in ... 
	double u_max = -1e9, u_min = 1e9;
	double v_max = -1e9, v_min = 1e9;
	double w_max = -1e9, w_min = 1e9;
	for (int k=0;k!=pcell.size();k++){
	  double udist = gds.wire_dist(pcell.at(k),WirePlaneType_t(0));
	  if (udist > u_max) u_max = udist;
	  if (udist < u_min) u_min = udist;
	  double vdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(1));
	  if (vdist > v_max) v_max = vdist;
	  if (vdist < v_min) v_min = vdist;
	  double wdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(2));
	  if (wdist > w_max) w_max = wdist;
	  if (wdist < w_min) w_min = wdist;
	}
	//std::cout << u_max << " " << u_min << " " << v_max << " " << v_min << " " << w_max << " " << w_min << std::endl;
	//Create a cell
	int ident = holder.get_ncell();
	SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
	mcell->AddBoundary(pcell);
	two_bad_wire_cells.push_back(mcell);
	holder.AddCell(mcell);
	//Insert U
	const GeomWire* uwire_min = gds.bounds(u_min,WirePlaneType_t(0)).second;
	const GeomWire* uwire_max = gds.bounds(u_max,WirePlaneType_t(0)).first;
	for (int k=uwire_min->index();k!=uwire_max->index()+1;k++){
	  const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),k);
	  mcell->AddWire(uwire,WirePlaneType_t(0));
	}
	
	//Insert V
	const GeomWire* vwire_min = gds.closest(v_min,WirePlaneType_t(1));//.bounds(v_min,WirePlaneType_t(1)).second;
	const GeomWire* vwire_max = gds.closest(v_max,WirePlaneType_t(1));//.bounds(v_max,WirePlaneType_t(1)).first;
	// std::cout << vwire_min->index() << " " << vwire_max->index() << std::endl;
	for (int k=vwire_min->index();k!=vwire_max->index()+1;k++){
	  const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),k);
	  mcell->AddWire(vwire,WirePlaneType_t(1));
	}

	//Insert W
	const GeomWire* wwire_min = gds.bounds(w_min,WirePlaneType_t(2)).second;
	const GeomWire* wwire_max = gds.bounds(w_max,WirePlaneType_t(2)).first;
	for (int k=wwire_min->index();k!=wwire_max->index()+1;k++){
	  const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),k);
	  mcell->AddWire(wwire,WirePlaneType_t(2));
	}


	//	std::cout << mcell->get_uwires().size() << " " << mcell->get_vwires().size() << " " << mcell->get_wwires().size() << std::endl;
      }
    } 
  }


  // W-V and add U
   // deal with w-v
  for ( int i = 0; i != bad_wire_w.size() ; i++){
    const GeomWire *wwire_1 = ((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().front();
    const GeomWire *wwire_2 = ((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().back();
    float dis_w[2];
    float w_pitch = gds.pitch(kYwire);
    dis_w[0] = gds.wire_dist(*wwire_1) - w_pitch/2.;
    dis_w[1] = gds.wire_dist(*wwire_2) + w_pitch/2.;
    

    for (int j = 0; j != bad_wire_v.size() ; j++){
      const GeomWire *vwire_1 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().front();
      const GeomWire *vwire_2 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().back();
      float dis_v[2];
      float v_pitch = gds.pitch(kVwire);
      dis_v[0] = gds.wire_dist(*vwire_1) - v_pitch/2.;
      dis_v[1] = gds.wire_dist(*vwire_2) + v_pitch/2.;
      
      PointVector pcell;
      std::vector<Vector> pcross(4);
      
      bool flag1 = gds.crossing_point(dis_w[0],dis_v[0],kYwire,kVwire, pcross[0]); 
      if (flag1){
  	pcell.push_back(pcross[0]);
      }else{
  	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().at(k);
	  dis_w[2] =  gds.wire_dist(*wwire_3) - w_pitch/2.;
	  if (gds.crossing_point(dis_w[2],dis_v[0],kYwire,kVwire, pcross[0])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(k);
	  dis_v[2] = gds.wire_dist(*vwire_3) - v_pitch/2.;
	  if (gds.crossing_point(dis_w[0],dis_v[2],kYwire,kVwire, pcross[0])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
      }

      bool flag2 = gds.crossing_point(dis_w[0],dis_v[1],kYwire,kVwire, pcross[1]); 
      if (flag2){
  	pcell.push_back(pcross[1]);
      }else{
  	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().at(k);
	  dis_w[2] =  gds.wire_dist(*wwire_3) - w_pitch/2.;
	  if (gds.crossing_point(dis_w[2],dis_v[1],kYwire,kVwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size()-1-k);
	  dis_v[2] = gds.wire_dist(*vwire_3) + v_pitch/2.;
	  if (gds.crossing_point(dis_w[0],dis_v[2],kYwire,kVwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
      }
      
      bool flag3 = gds.crossing_point(dis_w[1],dis_v[0],kYwire,kVwire, pcross[2]); 
      if (flag3){
  	pcell.push_back(pcross[2]);
      }else{
  	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().at(((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().size()-1-k);
	  dis_w[2] =  gds.wire_dist(*wwire_3) + w_pitch/2.;
	  if (gds.crossing_point(dis_w[2],dis_v[0],kYwire,kVwire, pcross[2])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(k);
	  dis_v[2] = gds.wire_dist(*vwire_3) - v_pitch/2.;
	  if (gds.crossing_point(dis_w[1],dis_v[2],kYwire,kVwire, pcross[2])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
      }
      
      bool flag4 = gds.crossing_point(dis_w[1],dis_v[1],kYwire,kVwire, pcross[3]);
      if (flag4){
  	pcell.push_back(pcross[3]);
      }else{
	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().at(((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().size()-1-k);
	  dis_w[2] =  gds.wire_dist(*wwire_3) + w_pitch/2.;
	  if (gds.crossing_point(dis_w[2],dis_v[1],kYwire,kVwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[3]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size()-1-k);
	  dis_v[2] = gds.wire_dist(*vwire_3) + v_pitch/2.;
	  if (gds.crossing_point(dis_w[1],dis_v[2],kYwire,kVwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	    pcell.push_back(pcross[3]);
	    break;
	  }
	}
      }
     
      
    
      if(pcell.size() >=3){
	//Creat a cell and then get all the wires in ... 
	double u_max = -1e9, u_min = 1e9;
	double v_max = -1e9, v_min = 1e9;
	double w_max = -1e9, w_min = 1e9;
	for (int k=0;k!=pcell.size();k++){
	  double udist = gds.wire_dist(pcell.at(k),WirePlaneType_t(0));
	  if (udist > u_max) u_max = udist;
	  if (udist < u_min) u_min = udist;
	  double vdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(1));
	  if (vdist > v_max) v_max = vdist;
	  if (vdist < v_min) v_min = vdist;
	  double wdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(2));
	  if (wdist > w_max) w_max = wdist;
	  if (wdist < w_min) w_min = wdist;
	}
	//std::cout << u_max << " " << u_min << " " << v_max << " " << v_min << " " << w_max << " " << w_min << std::endl;
	//Create a cell
	int ident = holder.get_ncell();
	SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
	// mcell->SetTimeSlice(time_slice);
	mcell->AddBoundary(pcell);
	two_bad_wire_cells.push_back(mcell);
	holder.AddCell(mcell);

	//Insert U	
	const GeomWire* uwire_min = gds.closest(u_min,WirePlaneType_t(0));//.bounds(u_min,WirePlaneType_t(0)).second;
	const GeomWire* uwire_max = gds.closest(u_max,WirePlaneType_t(0));//.bounds(u_max,WirePlaneType_t(0)).first;
	for (int k=uwire_min->index();k!=uwire_max->index()+1;k++){
	  const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),k);
	  mcell->AddWire(uwire,WirePlaneType_t(0));
	}
	
	//Insert V
	const GeomWire* vwire_min = gds.bounds(v_min,WirePlaneType_t(1)).second;
	const GeomWire* vwire_max = gds.bounds(v_max,WirePlaneType_t(1)).first;
	// std::cout << vwire_min->index() << " " << vwire_max->index() << std::endl;
	for (int k=vwire_min->index();k!=vwire_max->index()+1;k++){
	  const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),k);
	  mcell->AddWire(vwire,WirePlaneType_t(1));
	}

	//Insert W
	const GeomWire* wwire_min = gds.bounds(w_min,WirePlaneType_t(2)).second;
	const GeomWire* wwire_max = gds.bounds(w_max,WirePlaneType_t(2)).first;
	for (int k=wwire_min->index();k!=wwire_max->index()+1;k++){
	  const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),k);
	  mcell->AddWire(wwire,WirePlaneType_t(2));
	}


	// std::cout << mcell->get_uwires().size() << " " << mcell->get_vwires().size() << " " << mcell->get_wwires().size() << std::endl;
      }
    } 
  }

}


void WireCell2dToy::LowmemTiling::form_fired_merge_wires(const WireCell::Slice& slice){
  WireCell::Channel::Group group = slice.group();
  
  
  //double sum = 0;
  for (int i=0;i!=group.size();i++){
    const WireCell::GeomWire *wire = gds.by_channel_segment(group.at(i).first,0);
    float charge = group.at(i).second;
    if (charge < 11) charge = 11;
    if (wirechargemap.find(wire) == wirechargemap.end()){
      //not found
      wirechargemap[wire] = charge;
    }else{
      //wirechargemap[wire] += charge;
    }
  }
  
  // do U
  MergeGeomWire *mwire = 0;
  int ident = 0;
  int last_wire = -1;
  for (int i = 0; i!= nwire_u; i++){
    const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)0,i);
    if (wirechargemap.find(wire)!=wirechargemap.end()){
      if (mwire == 0){
	mwire = new MergeGeomWire(ident,*wire);
	//mwire->SetTimeSlice(time_slice);
	fired_wire_u.push_back(mwire);
	ident ++;
      }else{
	if (i==last_wire+1){
	  mwire->AddWire(*wire);
	}else{
	  // create a new wire
	  mwire = new MergeGeomWire(ident,*wire);
	  // mwire->SetTimeSlice(time_slice);
	  fired_wire_u.push_back(mwire);
	  ident ++;
	}
      }
      last_wire = i;
    }
  }
  
  // V
  mwire = 0;
  last_wire = -1;
  for (int i = 0; i!= nwire_v; i++){
    const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)1,i);
    if (wirechargemap.find(wire)!=wirechargemap.end()){
      if (mwire == 0){
	mwire = new MergeGeomWire(ident,*wire);
	//mwire->SetTimeSlice(time_slice);
	fired_wire_v.push_back(mwire);
	ident ++;
      }else{
	if (i==last_wire+1){
	  mwire->AddWire(*wire);
	}else{
	  // create a new wire
	  mwire = new MergeGeomWire(ident,*wire);
	  //mwire->SetTimeSlice(time_slice);
	  fired_wire_v.push_back(mwire);
	  ident ++;
	}
      }
      last_wire = i;
    }
  }
  
  // W
  mwire = 0;
  last_wire = -1;
  for (int i = 0; i!= nwire_w; i++){
    const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)2,i);
    if (wirechargemap.find(wire)!=wirechargemap.end()){
      if (mwire == 0){
	mwire = new MergeGeomWire(ident,*wire);
	//mwire->SetTimeSlice(time_slice);
	fired_wire_w.push_back(mwire);
	ident ++;
      }else{
	if (i==last_wire+1){
	  mwire->AddWire(*wire);
	}else{
	  // create a new wire
	  mwire = new MergeGeomWire(ident,*wire);
	  //mwire->SetTimeSlice(time_slice);
	  fired_wire_w.push_back(mwire);
	  ident ++;
	}
      }
      last_wire = i;
    }
  }

  
  // int sum = 0;
  // for (int i=0;i!=fired_wire_u.size();i++){
  //   MergeGeomWire *mwire = (MergeGeomWire*)fired_wire_u.at(i);
  //   sum += mwire->get_allwire().size();
  //   for (Int_t j=0;j!=mwire->get_allwire().size();j++){
  //     std::cout << i << " " << mwire->get_allwire().at(j)->index() << std::endl;
  //   }
  // }
  // for (int i=0;i!=fired_wire_v.size();i++){
  //   MergeGeomWire *mwire = (MergeGeomWire*)fired_wire_v.at(i);
  //   sum += mwire->get_allwire().size();
  // }
  // for (int i=0;i!=fired_wire_w.size();i++){
  //   MergeGeomWire *mwire = (MergeGeomWire*)fired_wire_w.at(i);
  //   sum += mwire->get_allwire().size();
  // }

  
  //std::cout << fired_wire_u.size() << " " << fired_wire_v.size() << " " << fired_wire_w.size() << " " << group.size() << " " << sum << std::endl;
  
}

void WireCell2dToy::LowmemTiling::form_bad_merge_wires(WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map){
  
  // do U
  MergeGeomWire *mwire = 0;
  int last_wire = -1;

  for (int i = 0; i!= nwire_u; i++){
    if (uplane_map.find(i)!=uplane_map.end()){
      if (time_slice >= uplane_map[i].first/nrebin && 
	  time_slice <= uplane_map[i].second/nrebin){
	const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)0,i);
	// first one 
	if (mwire == 0){
	  int ident = holder.get_nwire();
	  mwire = new MergeGeomWire(ident,*wire);
	  bad_wire_u.push_back(mwire);
	  holder.AddWire(mwire);
	}else{
	  if (i==last_wire+1){
	    mwire->AddWire(*wire);
	  }else{
	    // create a new wire
	    int ident = holder.get_nwire();
	    mwire = new MergeGeomWire(ident,*wire);
	    bad_wire_u.push_back(mwire);
	    holder.AddWire(mwire);
	  }
	}
	last_wire = i;
      }
    }
  }

  // now do V plane
  mwire = 0;
  last_wire = -1;

  for (int i = 0; i!= nwire_v; i++){
    if (vplane_map.find(i)!=vplane_map.end()){
      if (time_slice >= vplane_map[i].first/nrebin && 
	  time_slice <= vplane_map[i].second/nrebin){
	const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)1,i);
	// first one 
	if (mwire == 0){
	  int ident = holder.get_nwire();
	  mwire = new MergeGeomWire(ident,*wire);
	  bad_wire_v.push_back(mwire);
	  holder.AddWire(mwire);
	}else{
	  if (i==last_wire+1){
	    mwire->AddWire(*wire);
	  }else{
	    // create a new wire
	    int ident = holder.get_nwire();
	    mwire = new MergeGeomWire(ident,*wire);
	    bad_wire_v.push_back(mwire);
	    holder.AddWire(mwire);
	  }
	}
	last_wire = i;
      }
    }
  }
  
  // W plane
  mwire = 0;
  last_wire = -1;
  
  for (int i = 0; i!= nwire_w; i++){
    if (wplane_map.find(i)!=wplane_map.end()){
      if (time_slice >= wplane_map[i].first/nrebin && 
	  time_slice <= wplane_map[i].second/nrebin){
	const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)2,i);
	// first one 
	if (mwire == 0){
	  int ident = holder.get_nwire();
	  mwire = new MergeGeomWire(ident,*wire);
	  bad_wire_w.push_back(mwire);
	  holder.AddWire(mwire);
	}else{
	  if (i==last_wire+1){
	    mwire->AddWire(*wire);
	  }else{
	    // create a new wire
	    int ident = holder.get_nwire();
	    mwire = new MergeGeomWire(ident,*wire);
	    bad_wire_w.push_back(mwire);
	    holder.AddWire(mwire);
	  }
	}
	last_wire = i;
      }
    }
  }
  


  // int sum = 0;
  // for (int i=0;i!=bad_wire_w.size();i++){
  //   MergeGeomWire *mwire = (MergeGeomWire*)bad_wire_w.at(i);
  //   sum += mwire->get_allwire().size();
  //   for (int j=0;j!=mwire->get_allwire().size();j++){
  //     GeomWire *wire = (GeomWire*)mwire->get_allwire().at(j);
  //     std::cout << i << " " << wire->index() << std::endl;
  //   }
  // }
  // std::cout << uplane_map.size() << " " << bad_wire_u.size() << " " << sum << std::endl;
  //std::cout << vplane_map.size() << " " << bad_wire_v.size() << " " << sum << std::endl;
  //std::cout << wplane_map.size() << " " << bad_wire_w.size() << " " << sum << std::endl;

}

WireCell2dToy::LowmemTiling::~LowmemTiling(){
  // for (int i=0;i!=bad_wire_u.size();i++){
  //   delete bad_wire_u.at(i);
  // }
  // bad_wire_u.clear();
  
  // for (int i=0;i!=bad_wire_v.size();i++){
  //   delete bad_wire_v.at(i);
  // }
  // bad_wire_v.clear();

  // for (int i=0;i!=bad_wire_w.size();i++){
  //   delete bad_wire_w.at(i);
  // }
  // bad_wire_w.clear();
  
  // for (int i=0;i!=fired_wire_u.size();i++){
  //   delete fired_wire_u.at(i);
  // }
  // fired_wire_u.clear();
  
  // for (int i=0;i!=fired_wire_v.size();i++){
  //   delete fired_wire_v.at(i);
  // }
  // fired_wire_v.clear();

  // for (int i=0;i!=fired_wire_w.size();i++){
  //   delete fired_wire_w.at(i);
  // }
  // fired_wire_w.clear();
}

WireCell::GeomWireSelection WireCell2dToy::LowmemTiling::wires(const WireCell::GeomCell& cell) const{
  GeomWireSelection wires;
  return wires;
}
WireCell::GeomCellSelection WireCell2dToy::LowmemTiling::cells(const WireCell::GeomWire& wire) const{
  GeomCellSelection cells;
  return cells;
}

const WireCell::GeomCell* WireCell2dToy::LowmemTiling::cell(const WireCell::GeomWireSelection& wires) const{
  return 0;
}
