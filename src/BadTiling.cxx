#include "WireCell2dToy/BadTiling.h"

using namespace WireCell;

WireCell2dToy::BadTiling::BadTiling(int time, int scale, WireCell::ChirpMap& uplane_map, 
				    WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map, WireCell::GeomDataSource& gds){
  MergeGeomWire *mwire = 0;
  int prev_wire = -1;
  int num = 0;
  // do u wire
  for (int i=0;i!=gds.wires_in_plane(WirePlaneType_t(0)).size(); i++){
    if (uplane_map.find(i) != uplane_map.end() && i >= uplane_map[i].first/scale && i <= uplane_map[i].second/scale){
      const GeomWire* wire = gds.by_planeindex(WireCell::WirePlaneType_t(0),i);
      if (mwire == 0){
	mwire = new MergeGeomWire(num, *wire);
	wire_u.push_back(mwire);
	prev_wire = i;
	num++;
      }else{
	if ( i - prev_wire == 1){
	  mwire->AddWire(*wire);
	  prev_wire = i;
	}else{
	  mwire = new MergeGeomWire(num, *wire);
	  wire_u.push_back(mwire);
	  prev_wire = i;
	  num++;
	}
      }
    }
  }


  // do v wire
  prev_wire = -1;
  mwire = 0;
  
  for (int i=0;i!=gds.wires_in_plane(WirePlaneType_t(1)).size(); i++){
    if (vplane_map.find(i) != vplane_map.end() && i >= vplane_map[i].first/scale && i <= vplane_map[i].second/scale){
      const GeomWire* wire = gds.by_planeindex(WireCell::WirePlaneType_t(1),i);
      if (mwire == 0){
	mwire = new MergeGeomWire(num, *wire);
	wire_v.push_back(mwire);
	prev_wire = i;
	num++;
      }else{
	if ( i - prev_wire == 1){
	  mwire->AddWire(*wire);
	  prev_wire = i;
	}else{
	  mwire = new MergeGeomWire(num, *wire);
	  wire_v.push_back(mwire);
	  prev_wire = i;
	  num++;
	}
      }
    }
  }
  

  // do w wire
  prev_wire = -1;
  mwire = 0;
  
  for (int i=0;i!=gds.wires_in_plane(WirePlaneType_t(2)).size(); i++){
    if (wplane_map.find(i) != wplane_map.end() && i >= wplane_map[i].first/scale && i <= wplane_map[i].second/scale){
      const GeomWire* wire = gds.by_planeindex(WireCell::WirePlaneType_t(2),i);
      if (mwire == 0){
	mwire = new MergeGeomWire(num, *wire);
	wire_w.push_back(mwire);
	prev_wire = i;
	num++;
      }else{
	if ( i - prev_wire == 1){
	  mwire->AddWire(*wire);
	  prev_wire = i;
	}else{
	  mwire = new MergeGeomWire(num, *wire);
	  wire_w.push_back(mwire);
	  prev_wire = i;
	  num++;
	}
      }
    }
  }
 
  //Now form cells ... 

 
  

  num = 0;
  // deal with u-v
  for ( int i = 0; i != wire_u.size() ; i++){
    const GeomWire *uwire_1 = ((MergeGeomWire*)wire_u.at(i))->get_allwire().front();
    const GeomWire *uwire_2 = ((MergeGeomWire*)wire_u.at(i))->get_allwire().back();
    float dis_u[3];
    float u_pitch = gds.pitch(kUwire);
    dis_u[0] = gds.wire_dist(*uwire_1) - u_pitch/2.;
    dis_u[1] = gds.wire_dist(*uwire_2) + u_pitch/2.;
    

    for (int j = 0; j != wire_v.size() ; j++){
      const GeomWire *vwire_1 = ((MergeGeomWire*)wire_v.at(j))->get_allwire().front();
      const GeomWire *vwire_2 = ((MergeGeomWire*)wire_v.at(j))->get_allwire().back();
      float dis_v[3];
      float v_pitch = gds.pitch(kVwire);
      dis_v[0] = gds.wire_dist(*vwire_1) - v_pitch/2.;
      dis_v[1] = gds.wire_dist(*vwire_2) + v_pitch/2.;

      PointVector pcell;
      std::vector<Vector> pcross(4);
      
      bool flag1 = gds.crossing_point(dis_u[0],dis_v[0],kUwire,kVwire, pcross[0]); 
      if (flag1){
  	pcell.push_back(pcross[0]);
      }else{
  	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)wire_u.at(i))->get_allwire().at(k);
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
	for (int k=0;k!=((MergeGeomWire*)wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)wire_v.at(j))->get_allwire().at(k);
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
	for (int k=0;k!=((MergeGeomWire*)wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)wire_u.at(i))->get_allwire().at(k);
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
	for (int k=0;k!=((MergeGeomWire*)wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)wire_v.at(j))->get_allwire().at(((MergeGeomWire*)wire_v.at(j))->get_allwire().size()-1-k);
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
  	for (int k=0;k!=((MergeGeomWire*)wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)wire_u.at(i))->get_allwire().at(((MergeGeomWire*)wire_u.at(i))->get_allwire().size()-1-k);
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
	for (int k=0;k!=((MergeGeomWire*)wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)wire_v.at(j))->get_allwire().at(k);
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
	for (int k=0;k!=((MergeGeomWire*)wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)wire_u.at(i))->get_allwire().at(((MergeGeomWire*)wire_u.at(i))->get_allwire().size()-1-k);
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
	for (int k=0;k!=((MergeGeomWire*)wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)wire_v.at(j))->get_allwire().at(((MergeGeomWire*)wire_v.at(j))->get_allwire().size()-1-k);
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
     
      
    
      if(pcell.size() >=3){
	GeomCell *cell = new GeomCell(num,pcell);
	num ++;
	cell_all.push_back(cell);
      }
      
    } 
  }

  // deal with u-w
  for ( int i = 0; i != wire_u.size() ; i++){
    const GeomWire *uwire_1 = ((MergeGeomWire*)wire_u.at(i))->get_allwire().front();
    const GeomWire *uwire_2 = ((MergeGeomWire*)wire_u.at(i))->get_allwire().back();
    float dis_u[2];
    float u_pitch = gds.pitch(kUwire);
    dis_u[0] = gds.wire_dist(*uwire_1) - u_pitch/2.;
    dis_u[1] = gds.wire_dist(*uwire_2) + u_pitch/2.;
    

    for (int j = 0; j != wire_w.size() ; j++){
      const GeomWire *wwire_1 = ((MergeGeomWire*)wire_w.at(j))->get_allwire().front();
      const GeomWire *wwire_2 = ((MergeGeomWire*)wire_w.at(j))->get_allwire().back();
      float dis_w[2];
      float w_pitch = gds.pitch(kYwire);
      dis_w[0] = gds.wire_dist(*wwire_1) - w_pitch/2.;
      dis_w[1] = gds.wire_dist(*wwire_2) + w_pitch/2.;

      PointVector pcell;
      std::vector<Vector> pcross(4);
      
      
      bool flag1 = gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, pcross[0]); 
      if (flag1){
  	pcell.push_back(pcross[0]);
      }else{
  	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)wire_u.at(i))->get_allwire().at(k);
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
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)wire_w.at(j))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)wire_w.at(j))->get_allwire().at(k);
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
	for (int k=0;k!=((MergeGeomWire*)wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)wire_u.at(i))->get_allwire().at(k);
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
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)wire_w.at(j))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)wire_w.at(j))->get_allwire().at(((MergeGeomWire*)wire_w.at(j))->get_allwire().size()-1-k);
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
  	for (int k=0;k!=((MergeGeomWire*)wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)wire_u.at(i))->get_allwire().at(((MergeGeomWire*)wire_u.at(i))->get_allwire().size()-1-k);
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
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)wire_w.at(j))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)wire_w.at(j))->get_allwire().at(k);
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
	for (int k=0;k!=((MergeGeomWire*)wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)wire_u.at(i))->get_allwire().at(((MergeGeomWire*)wire_u.at(i))->get_allwire().size()-1-k);
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
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)wire_w.at(j))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)wire_w.at(j))->get_allwire().at(((MergeGeomWire*)wire_w.at(j))->get_allwire().size()-1-k);
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
	GeomCell *cell = new GeomCell(num,pcell);
	num ++;
	cell_all.push_back(cell);
      }
    } 
  }

   // deal with w-v
  for ( int i = 0; i != wire_w.size() ; i++){
    const GeomWire *wwire_1 = ((MergeGeomWire*)wire_w.at(i))->get_allwire().front();
    const GeomWire *wwire_2 = ((MergeGeomWire*)wire_w.at(i))->get_allwire().back();
    float dis_w[2];
    float w_pitch = gds.pitch(kYwire);
    dis_w[0] = gds.wire_dist(*wwire_1) - w_pitch/2.;
    dis_w[1] = gds.wire_dist(*wwire_2) + w_pitch/2.;
    

    for (int j = 0; j != wire_v.size() ; j++){
      const GeomWire *vwire_1 = ((MergeGeomWire*)wire_v.at(j))->get_allwire().front();
      const GeomWire *vwire_2 = ((MergeGeomWire*)wire_v.at(j))->get_allwire().back();
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
	for (int k=0;k!=((MergeGeomWire*)wire_w.at(i))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)wire_w.at(i))->get_allwire().at(k);
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
	for (int k=0;k!=((MergeGeomWire*)wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)wire_v.at(j))->get_allwire().at(k);
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
	for (int k=0;k!=((MergeGeomWire*)wire_w.at(i))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)wire_w.at(i))->get_allwire().at(k);
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
	for (int k=0;k!=((MergeGeomWire*)wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)wire_v.at(j))->get_allwire().at(((MergeGeomWire*)wire_v.at(j))->get_allwire().size()-1-k);
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
  	for (int k=0;k!=((MergeGeomWire*)wire_w.at(i))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)wire_w.at(i))->get_allwire().at(((MergeGeomWire*)wire_w.at(i))->get_allwire().size()-1-k);
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
	for (int k=0;k!=((MergeGeomWire*)wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)wire_v.at(j))->get_allwire().at(k);
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
	for (int k=0;k!=((MergeGeomWire*)wire_w.at(i))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)wire_w.at(i))->get_allwire().at(((MergeGeomWire*)wire_w.at(i))->get_allwire().size()-1-k);
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
	for (int k=0;k!=((MergeGeomWire*)wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)wire_v.at(j))->get_allwire().at(((MergeGeomWire*)wire_v.at(j))->get_allwire().size()-1-k);
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
	GeomCell *cell = new GeomCell(num,pcell);
	num ++;
	cell_all.push_back(cell);
      }
    } 
  }

  std::cout << "Bad: " << wire_u.size() << " " << wire_v.size() << " " << wire_w.size() << " " << cell_all.size() << std::endl;
}



WireCell2dToy::BadTiling::~BadTiling(){
}
