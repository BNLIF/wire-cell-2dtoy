#include "WireCell2dToy/LowmemTiling.h"

using namespace WireCell;

WireCell2dToy::LowmemTiling::LowmemTiling(int time_slice, int nrebin, const WireCell::Slice& slice,WireCell::GeomDataSource& gds,WireCell2dToy::WireCellHolder& holder,
					  std::vector<float>& uplane_rms, std::vector<float>& vplane_rms, std::vector<float>& wplane_rms)
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

 

  // form good wires group
  form_fired_merge_wires(slice);

  // create three good wire cells & two good wire + one bad wire cells
   
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

      bool flag1 = gds.crossing_point(dis_u[0],dis_v[0],kUwire,kVwire, pcross[0]); 
      if (flag1){
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
