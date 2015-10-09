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
    float dis_u[2];
    float u_pitch = gds.pitch(kUwire);
    dis_u[0] = gds.wire_dist(*uwire_1) - u_pitch/2.;
    dis_u[1] = gds.wire_dist(*uwire_2) + u_pitch/2.;
    

    for (int j = 0; j != wire_v.size() ; j++){
      const GeomWire *vwire_1 = ((MergeGeomWire*)wire_v.at(j))->get_allwire().front();
      const GeomWire *vwire_2 = ((MergeGeomWire*)wire_v.at(j))->get_allwire().back();
      float dis_v[2];
      float v_pitch = gds.pitch(kVwire);
      dis_v[0] = gds.wire_dist(*vwire_1) - v_pitch/2.;
      dis_v[1] = gds.wire_dist(*vwire_2) + v_pitch/2.;

       PointVector pcell;
      std::vector<Vector> pcross(4);
      
      
      bool flag1 = gds.crossing_point(dis_u[0],dis_v[0],kUwire,kVwire, pcross[0]); 
      if (flag1){
	pcell.push_back(pcross[0]);
      }else{
	Point p1 = uwire_1->point1();
	Point p2 = uwire_1->point2();
	Point p3 = vwire_1->point1();
	Point p4 = vwire_1->point2();

	if (pow(p1.y-pcross[0].y,2) + pow(p1.z-pcross[0].z,2) > 
	    pow(p2.y-pcross[0].y,2) + pow(p2.z-pcross[0].z,2)){
	  pcell.push_back(p2);
	}else{
	  pcell.push_back(p1);
	}
	if (pow(p3.y-pcross[0].y,2) + pow(p3.z-pcross[0].z,2) > 
	    pow(p4.y-pcross[0].y,2) + pow(p4.z-pcross[0].z,2)){
	  pcell.push_back(p4);
	}else{
	  pcell.push_back(p3);
	}
      }
      bool flag2 = gds.crossing_point(dis_u[0],dis_v[1],kUwire,kVwire, pcross[1]); 
      if (flag2){
	pcell.push_back(pcross[1]);
      }else{
	Point p1 = uwire_1->point1();
	Point p2 = uwire_1->point2();
	Point p3 = vwire_2->point1();
	Point p4 = vwire_2->point2();

	if (pow(p1.y-pcross[1].y,2) + pow(p1.z-pcross[1].z,2) > 
	    pow(p2.y-pcross[1].y,2) + pow(p2.z-pcross[1].z,2)){
	  pcell.push_back(p2);
	}else{
	  pcell.push_back(p1);
	}
	if (pow(p3.y-pcross[1].y,2) + pow(p3.z-pcross[1].z,2) > 
	    pow(p4.y-pcross[1].y,2) + pow(p4.z-pcross[1].z,2)){
	  pcell.push_back(p4);
	}else{
	  pcell.push_back(p3);
	}
      }
      bool flag3 = gds.crossing_point(dis_u[1],dis_v[0],kUwire,kVwire, pcross[2]); 
      if (flag3){
	pcell.push_back(pcross[2]);
      }else{
	Point p1 = uwire_2->point1();
	Point p2 = uwire_2->point2();
	Point p3 = vwire_1->point1();
	Point p4 = vwire_1->point2();

	if (pow(p1.y-pcross[2].y,2) + pow(p1.z-pcross[2].z,2) > 
	    pow(p2.y-pcross[2].y,2) + pow(p2.z-pcross[2].z,2)){
	  pcell.push_back(p2);
	}else{
	  pcell.push_back(p1);
	}
	if (pow(p3.y-pcross[2].y,2) + pow(p3.z-pcross[2].z,2) > 
	    pow(p4.y-pcross[2].y,2) + pow(p4.z-pcross[2].z,2)){
	  pcell.push_back(p4);
	}else{
	  pcell.push_back(p3);
	}
      }
      bool flag4 = gds.crossing_point(dis_u[1],dis_v[1],kUwire,kVwire, pcross[3]); 
      if (flag4){
	pcell.push_back(pcross[3]);
      }else{
	Point p1 = uwire_2->point1();
	Point p2 = uwire_2->point2();
	Point p3 = vwire_2->point1();
	Point p4 = vwire_2->point2();

	if (pow(p1.y-pcross[3].y,2) + pow(p1.z-pcross[3].z,2) > 
	    pow(p2.y-pcross[3].y,2) + pow(p2.z-pcross[3].z,2)){
	  pcell.push_back(p2);
	}else{
	  pcell.push_back(p1);
	}
	if (pow(p3.y-pcross[3].y,2) + pow(p3.z-pcross[3].z,2) > 
	    pow(p4.y-pcross[3].y,2) + pow(p4.z-pcross[3].z,2)){
	  pcell.push_back(p4);
	}else{
	  pcell.push_back(p3);
	}
      }
      //std::cout << dis_u[0] << " " << dis_u[1] << " " << dis_w[0] << " " << dis_w[1] << std::endl;
      // for (int k=0;k!=4;k++){
      // 	std::cout << pcross[k].z/units::m << " " << pcross[k].y/units::m << " " << 
      // 	  pcell.at(k).z/units::m << " " << pcell.at(k).y/units::m << std::endl;
      // }
      
      if (flag1 || flag2 || flag3 || flag4){
	GeomCell *cell = new GeomCell(num,pcell);
	// for (int k=0;k!=4;k++){
	// 	std::cout << pcell.at(k).z/units::m << " " << pcell.at(k).y/units::m << " " 
	// 		  << cell->boundary().at(k).z/units::m << " " << cell->boundary().at(k).y/units::m << std::endl;
	// }
	num ++;
	cell_all.push_back(cell);
	// std::cout << std::endl;
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
	Point p1 = uwire_1->point1();
	Point p2 = uwire_1->point2();
	Point p3 = wwire_1->point1();
	Point p4 = wwire_1->point2();

	if (pow(p1.y-pcross[0].y,2) + pow(p1.z-pcross[0].z,2) > 
	    pow(p2.y-pcross[0].y,2) + pow(p2.z-pcross[0].z,2)){
	  pcell.push_back(p2);
	}else{
	  pcell.push_back(p1);
	}
	if (pow(p3.y-pcross[0].y,2) + pow(p3.z-pcross[0].z,2) > 
	    pow(p4.y-pcross[0].y,2) + pow(p4.z-pcross[0].z,2)){
	  pcell.push_back(p4);
	}else{
	  pcell.push_back(p3);
	}
      }
      bool flag2 = gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, pcross[1]); 
      if (flag2){
	pcell.push_back(pcross[1]);
      }else{
	Point p1 = uwire_1->point1();
	Point p2 = uwire_1->point2();
	Point p3 = wwire_2->point1();
	Point p4 = wwire_2->point2();

	if (pow(p1.y-pcross[1].y,2) + pow(p1.z-pcross[1].z,2) > 
	    pow(p2.y-pcross[1].y,2) + pow(p2.z-pcross[1].z,2)){
	  pcell.push_back(p2);
	}else{
	  pcell.push_back(p1);
	}
	if (pow(p3.y-pcross[1].y,2) + pow(p3.z-pcross[1].z,2) > 
	    pow(p4.y-pcross[1].y,2) + pow(p4.z-pcross[1].z,2)){
	  pcell.push_back(p4);
	}else{
	  pcell.push_back(p3);
	}
      }
      bool flag3 = gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, pcross[2]); 
      if (flag3){
	pcell.push_back(pcross[2]);
      }else{
	Point p1 = uwire_2->point1();
	Point p2 = uwire_2->point2();
	Point p3 = wwire_1->point1();
	Point p4 = wwire_1->point2();

	if (pow(p1.y-pcross[2].y,2) + pow(p1.z-pcross[2].z,2) > 
	    pow(p2.y-pcross[2].y,2) + pow(p2.z-pcross[2].z,2)){
	  pcell.push_back(p2);
	}else{
	  pcell.push_back(p1);
	}
	if (pow(p3.y-pcross[2].y,2) + pow(p3.z-pcross[2].z,2) > 
	    pow(p4.y-pcross[2].y,2) + pow(p4.z-pcross[2].z,2)){
	  pcell.push_back(p4);
	}else{
	  pcell.push_back(p3);
	}
      }
      bool flag4 = gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, pcross[3]); 
      if (flag4){
	pcell.push_back(pcross[3]);
      }else{
	Point p1 = uwire_2->point1();
	Point p2 = uwire_2->point2();
	Point p3 = wwire_2->point1();
	Point p4 = wwire_2->point2();

	if (pow(p1.y-pcross[3].y,2) + pow(p1.z-pcross[3].z,2) > 
	    pow(p2.y-pcross[3].y,2) + pow(p2.z-pcross[3].z,2)){
	  pcell.push_back(p2);
	}else{
	  pcell.push_back(p1);
	}
	if (pow(p3.y-pcross[3].y,2) + pow(p3.z-pcross[3].z,2) > 
	    pow(p4.y-pcross[3].y,2) + pow(p4.z-pcross[3].z,2)){
	  pcell.push_back(p4);
	}else{
	  pcell.push_back(p3);
	}
      }
      //std::cout << dis_u[0] << " " << dis_u[1] << " " << dis_w[0] << " " << dis_w[1] << std::endl;
      // for (int k=0;k!=4;k++){
      // 	std::cout << pcross[k].z/units::m << " " << pcross[k].y/units::m << " " << 
      // 	  pcell.at(k).z/units::m << " " << pcell.at(k).y/units::m << std::endl;
      // }
      
      if (flag1 || flag2 || flag3 || flag4){
	GeomCell *cell = new GeomCell(num,pcell);
	// for (int k=0;k!=4;k++){
	// 	std::cout << pcell.at(k).z/units::m << " " << pcell.at(k).y/units::m << " " 
	// 		  << cell->boundary().at(k).z/units::m << " " << cell->boundary().at(k).y/units::m << std::endl;
	// }
	num ++;
	cell_all.push_back(cell);
	// std::cout << std::endl;
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
	Point p1 = wwire_1->point1();
	Point p2 = wwire_1->point2();
	Point p3 = vwire_1->point1();
	Point p4 = vwire_1->point2();

	if (pow(p1.y-pcross[0].y,2) + pow(p1.z-pcross[0].z,2) > 
	    pow(p2.y-pcross[0].y,2) + pow(p2.z-pcross[0].z,2)){
	  pcell.push_back(p2);
	}else{
	  pcell.push_back(p1);
	}
	if (pow(p3.y-pcross[0].y,2) + pow(p3.z-pcross[0].z,2) > 
	    pow(p4.y-pcross[0].y,2) + pow(p4.z-pcross[0].z,2)){
	  pcell.push_back(p4);
	}else{
	  pcell.push_back(p3);
	}
      }
      bool flag2 = gds.crossing_point(dis_w[0],dis_v[1],kYwire,kVwire, pcross[1]); 
      if (flag2){
	pcell.push_back(pcross[1]);
      }else{
	Point p1 = wwire_1->point1();
	Point p2 = wwire_1->point2();
	Point p3 = vwire_2->point1();
	Point p4 = vwire_2->point2();

	if (pow(p1.y-pcross[1].y,2) + pow(p1.z-pcross[1].z,2) > 
	    pow(p2.y-pcross[1].y,2) + pow(p2.z-pcross[1].z,2)){
	  pcell.push_back(p2);
	}else{
	  pcell.push_back(p1);
	}
	if (pow(p3.y-pcross[1].y,2) + pow(p3.z-pcross[1].z,2) > 
	    pow(p4.y-pcross[1].y,2) + pow(p4.z-pcross[1].z,2)){
	  pcell.push_back(p4);
	}else{
	  pcell.push_back(p3);
	}
      }
      bool flag3 = gds.crossing_point(dis_w[1],dis_v[0],kYwire,kVwire, pcross[2]); 
      if (flag3){
	pcell.push_back(pcross[2]);
      }else{
	Point p1 = wwire_2->point1();
	Point p2 = wwire_2->point2();
	Point p3 = vwire_1->point1();
	Point p4 = vwire_1->point2();

	if (pow(p1.y-pcross[2].y,2) + pow(p1.z-pcross[2].z,2) > 
	    pow(p2.y-pcross[2].y,2) + pow(p2.z-pcross[2].z,2)){
	  pcell.push_back(p2);
	}else{
	  pcell.push_back(p1);
	}
	if (pow(p3.y-pcross[2].y,2) + pow(p3.z-pcross[2].z,2) > 
	    pow(p4.y-pcross[2].y,2) + pow(p4.z-pcross[2].z,2)){
	  pcell.push_back(p4);
	}else{
	  pcell.push_back(p3);
	}
      }
      bool flag4 = gds.crossing_point(dis_w[1],dis_v[1],kYwire,kVwire, pcross[3]); 
      if (flag4){
	pcell.push_back(pcross[3]);
      }else{
	Point p1 = wwire_2->point1();
	Point p2 = wwire_2->point2();
	Point p3 = vwire_2->point1();
	Point p4 = vwire_2->point2();

	if (pow(p1.y-pcross[3].y,2) + pow(p1.z-pcross[3].z,2) > 
	    pow(p2.y-pcross[3].y,2) + pow(p2.z-pcross[3].z,2)){
	  pcell.push_back(p2);
	}else{
	  pcell.push_back(p1);
	}
	if (pow(p3.y-pcross[3].y,2) + pow(p3.z-pcross[3].z,2) > 
	    pow(p4.y-pcross[3].y,2) + pow(p4.z-pcross[3].z,2)){
	  pcell.push_back(p4);
	}else{
	  pcell.push_back(p3);
	}
      }
      //std::cout << dis_u[0] << " " << dis_u[1] << " " << dis_w[0] << " " << dis_w[1] << std::endl;
      // for (int k=0;k!=4;k++){
      // 	std::cout << pcross[k].z/units::m << " " << pcross[k].y/units::m << " " << 
      // 	  pcell.at(k).z/units::m << " " << pcell.at(k).y/units::m << std::endl;
      // }
      
      if (flag1 || flag2 || flag3 || flag4){
	GeomCell *cell = new GeomCell(num,pcell);
	// for (int k=0;k!=4;k++){
	// 	std::cout << pcell.at(k).z/units::m << " " << pcell.at(k).y/units::m << " " 
	// 		  << cell->boundary().at(k).z/units::m << " " << cell->boundary().at(k).y/units::m << std::endl;
	// }
	num ++;
	cell_all.push_back(cell);
	// std::cout << std::endl;
      }
    } 
  }

  std::cout << "Bad: " << wire_u.size() << " " << wire_v.size() << " " << wire_w.size() << " " << cell_all.size() << std::endl;
}



WireCell2dToy::BadTiling::~BadTiling(){
}
