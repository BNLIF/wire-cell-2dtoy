#include "WireCell2dToy/CalcPoints.h"

using namespace WireCell;

void WireCell2dToy::calc_boundary_points_dead(WireCell::GeomDataSource& gds, WireCell::PR3DCluster* cluster){
  SMGCSelection mcells = cluster->get_mcells();
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    calc_boundary_points_dead(gds,*it);
  }
}

void WireCell2dToy::calc_sampling_points(WireCell::GeomDataSource& gds, WireCell::PR3DCluster* cluster){
  SMGCSelection mcells = cluster->get_mcells();
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    calc_sampling_points(gds,*it);
  }
}

void WireCell2dToy::calc_boundary_points_dead(WireCell::GeomDataSource& gds, WireCell::SlimMergeGeomCell* mcell){
  GeomWireSelection bad_wire_u = mcell->get_uwires();
  GeomWireSelection bad_wire_v = mcell->get_vwires();
  GeomWireSelection bad_wire_w = mcell->get_wwires();

  bool flag_u = false, flag_v = false, flag_w = false;
  std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
  for (size_t i=0;i!=bad_planes.size();i++){
    if (bad_planes.at(i)==WirePlaneType_t(0)){
      flag_u = true;
    }else if (bad_planes.at(i)==WirePlaneType_t(1)){
      flag_v = true;
    }else if (bad_planes.at(i)==WirePlaneType_t(2)){
      flag_w = true;
    }
  }


  if (flag_u && flag_v){
    const GeomWire *uwire_1 = bad_wire_u.front();
    const GeomWire *uwire_2 = bad_wire_u.back();
    float dis_u[3];
    float u_pitch = gds.pitch(kUwire);
    dis_u[0] = gds.wire_dist(*uwire_1) - u_pitch/2.;
    dis_u[1] = gds.wire_dist(*uwire_2) + u_pitch/2.;
      
    const GeomWire *vwire_1 = bad_wire_v.front();
    const GeomWire *vwire_2 = bad_wire_v.back();
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
      for (int k=0;k!=bad_wire_u.size();k++){
	const GeomWire *uwire_3 = bad_wire_u.at(k);
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
      for (int k=0;k!=bad_wire_v.size();k++){
	const GeomWire *vwire_3 = bad_wire_v.at(k);
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
      for (int k=0;k!=bad_wire_u.size();k++){
	const GeomWire *uwire_3 = bad_wire_u.at(k);
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
      for (int k=0;k!=bad_wire_v.size();k++){
	const GeomWire *vwire_3 = bad_wire_v.at(bad_wire_v.size()-1-k);
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
      for (int k=0;k!=bad_wire_u.size();k++){
	const GeomWire *uwire_3 = bad_wire_u.at(bad_wire_u.size()-1-k);
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
      for (int k=0;k!=bad_wire_v.size();k++){
	const GeomWire *vwire_3 = bad_wire_v.at(k);
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
      for (int k=0;k!=bad_wire_u.size();k++){
	const GeomWire *uwire_3 = bad_wire_u.at(bad_wire_u.size()-1-k);
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
      for (int k=0;k!=bad_wire_v.size();k++){
	const GeomWire *vwire_3 = bad_wire_v.at(bad_wire_v.size()-1-k);
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
      mcell->AddBoundary(pcell);
    }
  }else if (flag_u && flag_w){
    const GeomWire *uwire_1 = bad_wire_u.front();
    const GeomWire *uwire_2 = bad_wire_u.back();
    float dis_u[2];
    float u_pitch = gds.pitch(kUwire);
    dis_u[0] = gds.wire_dist(*uwire_1) - u_pitch/2.;
    dis_u[1] = gds.wire_dist(*uwire_2) + u_pitch/2.;
    

    const GeomWire *wwire_1 = bad_wire_w.front();
    const GeomWire *wwire_2 = bad_wire_w.back();
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
      for (int k=0;k!=bad_wire_u.size();k++){
	const GeomWire *uwire_3 = bad_wire_u.at(k);
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
      for (int k=0;k!=bad_wire_w.size();k++){
	const GeomWire *wwire_3 = bad_wire_w.at(k);
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
      for (int k=0;k!=bad_wire_u.size();k++){
	const GeomWire *uwire_3 = bad_wire_u.at(k);
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
      for (int k=0;k!=bad_wire_w.size();k++){
	const GeomWire *wwire_3 = bad_wire_w.at(bad_wire_w.size()-1-k);
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
      for (int k=0;k!=bad_wire_u.size();k++){
	const GeomWire *uwire_3 = bad_wire_u.at(bad_wire_u.size()-1-k);
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
      for (int k=0;k!=bad_wire_w.size();k++){
	const GeomWire *wwire_3 = bad_wire_w.at(k);
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
      for (int k=0;k!=bad_wire_u.size();k++){
	const GeomWire *uwire_3 = bad_wire_u.at(bad_wire_u.size()-1-k);
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
      for (int k=0;k!=bad_wire_w.size();k++){
	const GeomWire *wwire_3 = bad_wire_w.at(bad_wire_w.size()-1-k);
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
      mcell->AddBoundary(pcell);
    }
  }else if (flag_v && flag_w){
    const GeomWire *wwire_1 = bad_wire_w.front();
    const GeomWire *wwire_2 = bad_wire_w.back();
    float dis_w[2];
    float w_pitch = gds.pitch(kYwire);
    dis_w[0] = gds.wire_dist(*wwire_1) - w_pitch/2.;
    dis_w[1] = gds.wire_dist(*wwire_2) + w_pitch/2.;
    
    const GeomWire *vwire_1 = bad_wire_v.front();
    const GeomWire *vwire_2 = bad_wire_v.back();
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
      for (int k=0;k!=bad_wire_w.size();k++){
	const GeomWire *wwire_3 = bad_wire_w.at(k);
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
      for (int k=0;k!=bad_wire_v.size();k++){
	const GeomWire *vwire_3 = bad_wire_v.at(k);
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
      for (int k=0;k!=bad_wire_w.size();k++){
	const GeomWire *wwire_3 = bad_wire_w.at(k);
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
      for (int k=0;k!=bad_wire_v.size();k++){
	const GeomWire *vwire_3 = bad_wire_v.at(bad_wire_v.size()-1-k);
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
      for (int k=0;k!=bad_wire_w.size();k++){
	const GeomWire *wwire_3 = bad_wire_w.at(bad_wire_w.size()-1-k);
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
      for (int k=0;k!=bad_wire_v.size();k++){
	const GeomWire *vwire_3 = bad_wire_v.at(k);
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
      for (int k=0;k!=bad_wire_w.size();k++){
	const GeomWire *wwire_3 = bad_wire_w.at(bad_wire_w.size()-1-k);
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
      for (int k=0;k!=bad_wire_v.size();k++){
	const GeomWire *vwire_3 = bad_wire_v.at(bad_wire_v.size()-1-k);
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
      mcell->AddBoundary(pcell);
    }
  }
  
  
  
  
  
}

void WireCell2dToy::calc_sampling_points(WireCell::GeomDataSource& gds, WireCell::SlimMergeGeomCell* mcell){
  GeomWireSelection wires_u = mcell->get_uwires();
  GeomWireSelection wires_v = mcell->get_vwires();
  GeomWireSelection wires_w = mcell->get_wwires();

  GeomWireSelection max_wires = wires_u;
  WirePlaneType_t max_wire_plane_type = WirePlaneType_t(0);
  GeomWireSelection min_wires = wires_v;
  WirePlaneType_t min_wire_plane_type = WirePlaneType_t(1);
  GeomWireSelection other_wires;
  WirePlaneType_t other_wire_plane_type;
  
  if (wires_v.size() > max_wires.size()){
    max_wires = wires_v;
    max_wire_plane_type = WirePlaneType_t(1);
  }
  if (wires_w.size() > max_wires.size()){
    max_wires = wires_w;
    max_wire_plane_type = WirePlaneType_t(2);
  }
  
  if (wires_u.size() < min_wires.size()){
    min_wires = wires_u;
    min_wire_plane_type = WirePlaneType_t(0);
  }
  if (wires_w.size() < min_wires.size()){
    min_wires = wires_w;
    min_wire_plane_type = WirePlaneType_t(2);
  }

  if (max_wire_plane_type==WirePlaneType_t(0) &&
      min_wire_plane_type==WirePlaneType_t(1) ||
      max_wire_plane_type==WirePlaneType_t(1) &&
      min_wire_plane_type==WirePlaneType_t(0)){
    other_wire_plane_type = WirePlaneType_t(2);
    other_wires = wires_w;
  }else if (max_wire_plane_type==WirePlaneType_t(0) &&
	    min_wire_plane_type==WirePlaneType_t(2) ||
	    max_wire_plane_type==WirePlaneType_t(2) &&
	    min_wire_plane_type==WirePlaneType_t(0)){
    other_wire_plane_type = WirePlaneType_t(1);
    other_wires = wires_v;
  }else if (max_wire_plane_type==WirePlaneType_t(2) &&
	    min_wire_plane_type==WirePlaneType_t(1) ||
	    max_wire_plane_type==WirePlaneType_t(1) &&
	    min_wire_plane_type==WirePlaneType_t(2)){
    other_wire_plane_type = WirePlaneType_t(0);
    other_wires = wires_u;
  }

  PointVector sampling_points;
  
  float other_pitch = gds.pitch(other_wire_plane_type);
  float dis_limit[2];
  float tolerance = 0.1 * units::mm;
  const GeomWire *other_wire_1 = other_wires.front();
  const GeomWire *other_wire_2 = other_wires.back();
  dis_limit[0] = gds.wire_dist(*other_wire_1)-other_pitch/2.-tolerance;
  dis_limit[1] = gds.wire_dist(*other_wire_2)+other_pitch/2.+tolerance;

  //std::cout << dis_limit[0] << " " << dis_limit[1] << std::endl;

  
  
  int max_step = std::max(4.0,max_wires.size()/12.);
  int min_step = std::max(4.0,min_wires.size()/12.);
  //std::cout << min_step << " " << max_step << std::endl;
  GeomWireSetp max_wires_set;
  GeomWireSetp min_wires_set;
  max_wires_set.insert(max_wires.front());
  for (size_t i=0;i<max_wires.size();i+=max_step){
    max_wires_set.insert(max_wires.at(i));
  }
  max_wires_set.insert(max_wires.back());

  min_wires_set.insert(min_wires.front());
  for (size_t i=0;i<min_wires.size();i+=min_step){
    min_wires_set.insert(min_wires.at(i));
  }
  min_wires_set.insert(min_wires.back());

  for (auto it = max_wires_set.begin(); it!=max_wires_set.end(); it++){
    for (auto it1 = min_wires_set.begin(); it1!=min_wires_set.end();it1++){
      Vector point;
      gds.crossing_point(*(*it),*(*it1),point);
      float dis = gds.wire_dist(point,other_wire_plane_type);
      if (dis>=dis_limit[0]&&dis<=dis_limit[1])
	sampling_points.push_back(point);
      //std::cout << "A: " <<dis_limit[0] << " " << dis << " " << dis_limit[1] << std::endl; 
    }
  }
  
  //std::cout << sampling_points.size() << " " << wires_u.size() <<  " " << wires_v.size() << " " << wires_w.size() << " " << max_wires.size() << " " << min_wires.size() << std::endl;
  
  if (sampling_points.size()>0)
    mcell->AddSamplingPoints(sampling_points);
  
  // std::cout << min_wires.size() << " " << max_wires.size() << " " << wires_u.size() << " " << wires_v.size() << " " << wires_w.size() << std::endl;

  
  
}
