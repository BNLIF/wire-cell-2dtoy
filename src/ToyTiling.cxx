#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/TotalTiling.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/Singleton.h"
#include <cmath>

using namespace WireCell;

WireCell2dToy::ToyTiling::ToyTiling()
{
}


void WireCell2dToy::ToyTiling::AddCell(WireCell::DetectorGDS& gds, int cryo, int apa, GeomCell *cell, int u_index, int v_index, int w_index, 
				       float u_charge, float v_charge, float w_charge,
				       float u_charge_err, float v_charge_err, float w_charge_err){
  
  cell_all.push_back(cell);
  const GeomWire *uwire = gds.get_apaGDS(cryo,apa)->by_planeindex(WirePlaneType_t(0),u_index); 
  const GeomWire *vwire = gds.get_apaGDS(cryo,apa)->by_planeindex(WirePlaneType_t(1),v_index); 
  const GeomWire *wwire = gds.get_apaGDS(cryo,apa)->by_planeindex(WirePlaneType_t(2),w_index); 
  
  cell->set_uwire(uwire);
  cell->set_vwire(vwire);
  cell->set_wwire(wwire);

  auto it1 = find(wire_all.begin(),wire_all.end(),uwire);
  if (it1 == wire_all.end()){
    wire_all.push_back(uwire);
    wire_u.push_back(uwire);
    wirechargemap[uwire] = u_charge;
    wirecharge_errmap[uwire] = u_charge_err;
  }
  auto it2 = find(wire_all.begin(),wire_all.end(),vwire);
  if (it2 == wire_all.end()){
    wire_all.push_back(vwire);
    wire_v.push_back(vwire);
    wirechargemap[vwire] = v_charge;
    wirecharge_errmap[vwire] = v_charge_err;
  }
  auto it3 = find(wire_all.begin(),wire_all.end(),wwire);
  if (it3 == wire_all.end()){
    wire_all.push_back(wwire);
    wire_w.push_back(wwire);
    wirechargemap[wwire] = w_charge;
    wirecharge_errmap[wwire] = w_charge_err;
  }

  GeomWireSelection wires;
  wires.push_back(uwire);
  wires.push_back(vwire);
  wires.push_back(wwire);
  cellmap[cell] = wires;
      
  if (wiremap.find(uwire) == wiremap.end()){
    //not found
    GeomCellSelection cellsel;
    cellsel.push_back(cell);
    wiremap[uwire]=cellsel;
  }else{
    //found
    wiremap[uwire].push_back(cell);
  }


  if (wiremap.find(vwire) == wiremap.end()){
    //not found
    GeomCellSelection cellsel;
    cellsel.push_back(cell);
    wiremap[vwire]=cellsel;
  }else{
    //found
    wiremap[vwire].push_back(cell);
  }
  
  if (wiremap.find(wwire) == wiremap.end()){
    //not found
    GeomCellSelection cellsel;
    cellsel.push_back(cell);
    wiremap[wwire]=cellsel;
  }else{
    //found
    wiremap[wwire].push_back(cell);
  }
 
  ave_charge = 0;
  for (int i=0;i!=wire_all.size();i++){
    ave_charge += wirechargemap[wire_all.at(i)];
  }
  if (cell_all.size()!=0){
    ave_charge = ave_charge / cell_all.size()/3.;
  }
 

}



void WireCell2dToy::ToyTiling::AddCell(WireCell::GeomDataSource& gds, GeomCell *cell, int u_index, int v_index, int w_index, 
				       float u_charge, float v_charge, float w_charge,
				       float u_charge_err, float v_charge_err, float w_charge_err){
  
  cell_all.push_back(cell);
  const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),u_index); 
  const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),v_index); 
  const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),w_index); 
  
  cell->set_uwire(uwire);
  cell->set_vwire(vwire);
  cell->set_wwire(wwire);

  auto it1 = find(wire_all.begin(),wire_all.end(),uwire);
  if (it1 == wire_all.end()){
    wire_all.push_back(uwire);
    wire_u.push_back(uwire);
    wirechargemap[uwire] = u_charge;
    wirecharge_errmap[uwire] = u_charge_err;
  }
  auto it2 = find(wire_all.begin(),wire_all.end(),vwire);
  if (it2 == wire_all.end()){
    wire_all.push_back(vwire);
    wire_v.push_back(vwire);
    wirechargemap[vwire] = v_charge;
    wirecharge_errmap[vwire] = v_charge_err;
  }
  auto it3 = find(wire_all.begin(),wire_all.end(),wwire);
  if (it3 == wire_all.end()){
    wire_all.push_back(wwire);
    wire_w.push_back(wwire);
    wirechargemap[wwire] = w_charge;
    wirecharge_errmap[wwire] = w_charge_err;
  }

  GeomWireSelection wires;
  wires.push_back(uwire);
  wires.push_back(vwire);
  wires.push_back(wwire);
  cellmap[cell] = wires;
      
  if (wiremap.find(uwire) == wiremap.end()){
    //not found
    GeomCellSelection cellsel;
    cellsel.push_back(cell);
    wiremap[uwire]=cellsel;
  }else{
    //found
    wiremap[uwire].push_back(cell);
  }


  if (wiremap.find(vwire) == wiremap.end()){
    //not found
    GeomCellSelection cellsel;
    cellsel.push_back(cell);
    wiremap[vwire]=cellsel;
  }else{
    //found
    wiremap[vwire].push_back(cell);
  }
  
  if (wiremap.find(wwire) == wiremap.end()){
    //not found
    GeomCellSelection cellsel;
    cellsel.push_back(cell);
    wiremap[wwire]=cellsel;
  }else{
    //found
    wiremap[wwire].push_back(cell);
  }
 
  ave_charge = 0;
  for (int i=0;i!=wire_all.size();i++){
    ave_charge += wirechargemap[wire_all.at(i)];
  }
  if (cell_all.size()!=0){
    ave_charge = ave_charge / cell_all.size()/3.;
  }
 

}

WireCell2dToy::ToyTiling::ToyTiling(const WireCell::Slice& slice,WireCell::DetectorGDS& gds, float rel_u , float rel_v, float rel_w, float noise_u, float noise_v, float noise_w, std::vector<float>* uplane_rms, std::vector<float>* vplane_rms, std::vector<float>* wplane_rms){
  WireCell::Channel::Group group = slice.group();
  float tolerance = 0.1 * units::mm;
  //save all the wires

  ncell = 1;

  for (int i=0;i!=group.size();i++){
    
    const GeomWireSelection& wires = gds.by_channel(group.at(i).first);
    //std::cout << "abc: " << wires.at(0)->plane() << " " << group.at(i).first << " " << group.at(i).second << std::endl;
    
    for (int i1 = 0;i1 != wires.size();i1++){
      const GeomWire *wire = wires.at(i1);
      int channel = group.at(i).first;
      float charge = group.at(i).second;
      // sum += charge;
      if (wirechargemap.find(wire) == wirechargemap.end()){
	//not found
	wirechargemap[wire] = charge;
      }else{
	wirechargemap[wire] += charge;
      }
      //fill in the error ...
      WirePlaneType_t plane = wire->plane();
      // std::cout << i << " " << plane << " " << charge << std::endl;
      
      double charge_noise;
      double rel_charge_err;
      if (plane ==0){
	charge_noise = noise_u;
	rel_charge_err = rel_u;
      }else if (plane == 1){
	charge_noise = noise_v;
	rel_charge_err = rel_v;
      }else if (plane == 2){
	charge_noise = noise_w;
	rel_charge_err = rel_w;
      }
      
      if (uplane_rms !=0 && vplane_rms !=0 && wplane_rms !=0){
	int index = wire->index();
	if (plane ==0){
	  charge_noise = uplane_rms->at(index);
	}else if (plane == 1){
	  charge_noise = vplane_rms->at(index);
	}else if (plane == 2){
	  charge_noise = wplane_rms->at(index);
	}
      }
      wirecharge_errmap[wire] = sqrt(charge_noise*charge_noise + pow(rel_charge_err,2) * charge*charge);
      
      if (i1 == 0){
	if (channelchargemap.find(channel) == channelchargemap.end()){
	  //not found
	  channelchargemap[channel] = charge;
	}else{
	  channelchargemap[channel] += charge;
	}
        channelcharge_errmap[channel] = sqrt(charge_noise*charge_noise + pow(rel_charge_err,2) * charge*charge);
	
      }
      
    }
      
    
    for (int j=0;j!=wires.size();j++){
      const GeomWire *wire = wires.at(j);
      wire_all.push_back(wire);
      if (wire->plane() == WirePlaneType_t(0)){
    	wire_u.push_back(wire);
      }else if (wire->plane() == WirePlaneType_t(1)){
    	wire_v.push_back(wire);
      }else if (wire->plane() == WirePlaneType_t(2)){
    	wire_w.push_back(wire);
      }
    }
  }

  //sort wires as APA and faces, do tiling one APA by one APA
  // for cells save n_cryos*10000 + n_apa*10 + face
  int n_cryos = gds.ncryos();
  for (int i=0;i!=n_cryos;i++){
    int n_apa = gds.napa(i);
    for (int j=0;j!=n_apa;j++){
      // loop over face
      for (int this_face = 0; this_face!=2;this_face++){
	int n_tpc = i*10000 + j*10 + this_face;
	//	std::cout << i << " " << j << " " << n_tpc << std::endl;
	GeomWireSelection temp_uwire, temp_vwire, temp_wwire;
	
	for (int k = 0; k!= wire_u.size();k++){
	  const GeomWire *wire = wire_u.at(k);
	  if (wire->cryo() == i && wire->apa() == j && wire->face() == this_face)
	    temp_uwire.push_back(wire);
	}
	for (int k = 0; k!= wire_v.size();k++){
	  const GeomWire *wire = wire_v.at(k);
	  if (wire->cryo() == i && wire->apa() == j && wire->face() == this_face)
	    temp_vwire.push_back(wire);
	}
	for (int k = 0; k!= wire_w.size();k++){
	  const GeomWire *wire = wire_w.at(k);
	  if (wire->cryo() == i && wire->apa() == j && wire->face() == this_face)
	    temp_wwire.push_back(wire);
	}
	
	//std::cout << i << " " << j << " " << this_face << " " << temp_uwire.size() << " " << temp_vwire.size() << " " << temp_wwire.size() << std::endl;
	//	det_gds.get_apaGDS(det_gds.in_which_cryo(pt), det_gds.in_which_apa(pt));
	
	const GeomDataSource *gds_apa = gds.get_apaGDS(i,j);
	
	CreateCell(tolerance, *gds_apa,this_face, n_tpc,temp_uwire, temp_vwire, temp_wwire);
	
	

      }
    }
  }
 

  // std::cout << wire_all.size() << " " << wire_u.size() << " " << wire_v.size() << " " << wire_w.size() << " " << cell_all.size() << std::endl;
  
}



WireCell2dToy::ToyTiling::ToyTiling(const WireCell::Slice& slice,WireCell::GeomDataSource& gds, float rel_u , float rel_v, float rel_w, float noise_u, float noise_v, float noise_w, std::vector<float>* uplane_rms, std::vector<float>* vplane_rms, std::vector<float>* wplane_rms){
  WireCell::Channel::Group group = slice.group();

  float tolerance = 0.1 * units::mm;

  //double sum = 0;
  for (int i=0;i!=group.size();i++){
    const WireCell::GeomWire *wire = gds.by_channel_segment(group.at(i).first,0);
    
    float charge = group.at(i).second;
    // sum += charge;

    if (wirechargemap.find(wire) == wirechargemap.end()){
      //not found
      wirechargemap[wire] = charge;
    }else{
      wirechargemap[wire] += charge;
    }

   

    //fill in the error ...
    WirePlaneType_t plane = wire->plane();

    // std::cout << i << " " << plane << " " << charge << std::endl;

    double charge_noise;
    double rel_charge_err;
    if (plane ==0){
      charge_noise = noise_u;
      rel_charge_err = rel_u;
    }else if (plane == 1){
      charge_noise = noise_v;
      rel_charge_err = rel_v;
    }else if (plane == 2){
      charge_noise = noise_w;
      rel_charge_err = rel_w;
    }

    if (uplane_rms !=0 && vplane_rms !=0 && wplane_rms !=0){
      int index = wire->index();
      if (plane ==0){
	charge_noise = uplane_rms->at(index);
      }else if (plane == 1){
	charge_noise = vplane_rms->at(index);
      }else if (plane == 2){
	charge_noise = wplane_rms->at(index);
      }
    }


    wirecharge_errmap[wire] = sqrt(charge_noise*charge_noise + pow(rel_charge_err,2) * charge*charge);

    // // hack for now
    // if (charge < 100){
    //   //std::cout << "Bad Channel, Fill in a large error " << std::endl;
    //   wirecharge_errmap[wire] = 1e4;
    // }
    // // hack for data

    
    wire_all.push_back(wire);
    if (wire->plane() == kUwire){
      wire_u.push_back(wire);
    }else if (wire->plane() == kVwire){
      wire_v.push_back(wire);
    }else if (wire->plane() == kYwire){
      wire_w.push_back(wire);
    }
    

  }


  // GeomWireSelection wire_all1;
  // GeomWireSelection wire_u1;
  // GeomWireSelection wire_v1;
  // GeomWireSelection wire_w1;



  float dis_u[3],dis_v[3],dis_w[3],
    dis_puv[5],dis_puw[5],dis_pwv[5];
  ncell = 1;

  //std::cout << "Wire Counts: " << wire_u.size() << " " << wire_v.size() << " " << wire_w.size() << std::endl;

 

  // calculate all the costant
  std::vector<float> udis,vdis,wdis;
  for (int i=0;i!=wire_u.size();i++){
    udis.push_back(gds.wire_dist(*wire_u[i]));
  }
  for (int j=0;j!=wire_v.size();j++){
    vdis.push_back(gds.wire_dist(*wire_v[j]));
  }
  for (int k=0;k!=wire_w.size();k++){
    wdis.push_back(gds.wire_dist(*wire_w[k]));
  }

  //precalculate all the offsets
  std::vector<Vector> puv_save(5), puw_save(5), pwv_save(5);
  float dis_puv_save[5], dis_puw_save[5], dis_pwv_save[5];
  
  double u_pitch, v_pitch, w_pitch;
  u_pitch = gds.pitch(kUwire);
  v_pitch = gds.pitch(kVwire);
  w_pitch = gds.pitch(kYwire);


  if (wire_u.size()>=1&&wire_v.size()>=1){
    dis_u[0] = udis.at(0) - gds.pitch(kUwire)/2.;
    dis_u[1] = dis_u[0] + gds.pitch(kUwire);
    dis_u[2] = udis.at(0);

    dis_v[0] = vdis.at(0) - gds.pitch(kVwire)/2.;
    dis_v[1] = dis_v[0] + gds.pitch(kVwire);
    dis_v[2] = vdis.at(0);

    gds.crossing_point(dis_u[0],dis_v[0],kUwire,kVwire, puv_save[0]);
    gds.crossing_point(dis_u[0],dis_v[1],kUwire,kVwire, puv_save[1]);
    gds.crossing_point(dis_u[1],dis_v[1],kUwire,kVwire, puv_save[2]);
    gds.crossing_point(dis_u[1],dis_v[0],kUwire,kVwire, puv_save[3]);
    gds.crossing_point(dis_u[2],dis_v[2],kUwire,kVwire, puv_save[4]);

    for (int k=0;k!=5;k++){
      dis_puv_save[k] = gds.wire_dist(puv_save[k],kYwire);
    }
    
  }

  if (wire_u.size()>=1&&wire_w.size()>=1){
    dis_u[0] = udis.at(0) - gds.pitch(kUwire)/2.;
    dis_u[1] = dis_u[0] + gds.pitch(kUwire);
    dis_u[2] = udis.at(0);

    dis_w[0] = wdis.at(0) - w_pitch/2.;
    dis_w[1] = dis_w[0] + w_pitch;
    dis_w[2] = wdis.at(0);

    gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, puw_save[0]);
    gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, puw_save[1]);
    gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, puw_save[2]);
    gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, puw_save[3]);
    gds.crossing_point(dis_u[2],dis_w[2],kUwire,kYwire, puw_save[4]);

    for (int k=0;k!=5;k++){
      dis_puw_save[k] = gds.wire_dist(puw_save[k],kVwire);
    }

  }

  if (wire_v.size()>=1&&wire_w.size()>=1){

    dis_w[0] = wdis.at(0) - w_pitch/2.;
    dis_w[1] = dis_w[0] + w_pitch;
    dis_w[2] = wdis.at(0);

    dis_v[0] = vdis.at(0) - gds.pitch(kVwire)/2.;
    dis_v[1] = dis_v[0] + gds.pitch(kVwire);
    dis_v[2] = vdis.at(0);
    
    gds.crossing_point(dis_v[0],dis_w[0],kVwire,kYwire, pwv_save[0]);
    gds.crossing_point(dis_v[0],dis_w[1],kVwire,kYwire, pwv_save[1]);
    gds.crossing_point(dis_v[1],dis_w[1],kVwire,kYwire, pwv_save[2]);
    gds.crossing_point(dis_v[1],dis_w[0],kVwire,kYwire, pwv_save[3]);
    gds.crossing_point(dis_v[2],dis_w[2],kVwire,kYwire, pwv_save[4]);
    
    for (int k=0;k!=5;k++){
      dis_pwv_save[k] = gds.wire_dist(pwv_save[k],kUwire);
    }
  }

  //  int counter = 0, counter1 = 0, counter2 = 0, counter3=0;

  


  for (int i=0;i!=wire_u.size();i++){
    dis_u[0] = udis.at(i) - u_pitch/2.;
    dis_u[1] = dis_u[0] + u_pitch;
    dis_u[2] = udis.at(i);

    for (int j=0;j!=wire_v.size();j++){
      //  if (wirechargemap[wire_u[i]] <100 && wirechargemap[wire_v[j]] <100 ) continue;
      dis_v[0] = vdis.at(j) - v_pitch/2.;
      dis_v[1] = dis_v[0] + v_pitch;
      dis_v[2] = vdis.at(j);
      
      //four vertices around
      //      std::vector<Vector> puv(4);
      std::vector<Vector> puv(5);
      
      //counter ++;

      if(!gds.crossing_point(dis_u[2],dis_v[2],kUwire,kVwire, puv[4])) continue;
      
      // Point pp;
      // pp.x = 0; pp.y = puv[4].y; pp.z = puv[4].z;
      // if (!gds.contained_yz(pp)) continue;

      //counter1 ++;
      dis_puv[4] = gds.wire_dist(puv[4],kYwire);
      for (int k=0;k!=4;k++){
	puv[k] = puv[4] + (puv_save[k]-puv_save[4]);
	dis_puv[k] = dis_puv[4] +(dis_puv_save[k]-dis_puv_save[4]);
      }

      // gds.crossing_point(dis_u[0],dis_v[0],kUwire,kVwire, puv[0]);
      // gds.crossing_point(dis_u[0],dis_v[1],kUwire,kVwire, puv[1]);
      // gds.crossing_point(dis_u[1],dis_v[1],kUwire,kVwire, puv[2]);
      // gds.crossing_point(dis_u[1],dis_v[0],kUwire,kVwire, puv[3]);
      //  for (int k=0;k!=4;k++){
      // 	dis_puv[k] = gds.wire_dist(puv[k],kYwire);
      // }
      // std::cout << dis_u[0] << " " << dis_u[1] << " " << dis_v[0] << " " << dis_v[1] << " " << gds.pitch(kUwire) << " " << gds.pitch(kVwire) << " " << gds.angle(WirePlaneType_t(0)) << " " << gds.angle(WirePlaneType_t(1)) << std::endl;
      // std::cout << "Xin1 " << std::endl;
      // std::cout << puv[0].y << " " << puv[0].z << std::endl;
      // std::cout << puv[1].y << " " << puv[1].z << std::endl;
      // std::cout << puv[2].y << " " << puv[2].z << std::endl;
      // std::cout << puv[3].y << " " << puv[3].z << std::endl;
      // std::cout << wire_u[i]->ident() << " " << wire_v[i]->ident()  << std::endl;
       
      for (int k=0;k!=wire_w.size();k++){
	
	// if (wirechargemap[wire_u[i]] <100 && wirechargemap[wire_w[k]] <100 ) continue;
	// if (wirechargemap[wire_w[k]] <100 && wirechargemap[wire_v[j]] <100 ) continue;
	//counter2++;

	int flag = 0;
	PointVector pcell;
	dis_w[0] = wdis.at(k) - w_pitch/2.;
  	dis_w[1] = dis_w[0] + w_pitch;//gds.wire_dist(*wire_w[k]) + gds.pitch(kYwire)/2.;	
	dis_w[2] = wdis.at(k);

	if (fabs(dis_w[0] - dis_puv[0])>3*units::cm) continue;
	
	//counter3++;
	
	for (int m = 0;m!=4;m++){
	  if (dis_puv[m] > dis_w[0]-tolerance/2. && dis_puv[m] < dis_w[1]+tolerance/2.){
	    flag = 1;
	    pcell.push_back(puv[m]);
	  }
	}

	// std::cout << "Xin " << pcell.size() << std::endl; 
	// for (int k=0;k!=pcell.size();k++){
	//   std::cout << pcell.at(k).y << " " << pcell.at(k).z  << std::endl;
	// }

	if (flag==1 ) {
	  //counter ++;
	  std::vector<Vector> puw(5);
	  // gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, puw[0]);
	  // gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, puw[1]);
	  // gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, puw[2]);
	  // gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, puw[3]);
	  
	  gds.crossing_point(dis_u[2],dis_w[2],kUwire,kYwire, puw[4]);
	  dis_puw[4] = gds.wire_dist(puw[4],kVwire);

	  for (int k1=0;k1!=4;k1++){
	    puw[k1] = puw[4] + (puw_save[k1]-puw_save[4]);
	    dis_puw[k1] = dis_puw[4] +(dis_puw_save[k1]-dis_puw_save[4]);
	    if (dis_puw[k1] > dis_v[0]-tolerance/2. && dis_puw[k1] < dis_v[1]+tolerance/2.){
	      int flag_abc = 0;
	      for (int kk = 0; kk!=pcell.size();kk++){
	      	float dis = sqrt(pow(puw[k1].y-pcell.at(kk).y,2) + pow(puw[k1].z-pcell.at(kk).z,2));
	      	if (dis < tolerance) {
	      	  flag_abc = 1;
	      	  break;
	      	}
	      }
	      if (flag_abc == 0)
	    	pcell.push_back(puw[k1]);
	    }
	  }
	  std::vector<Vector> pwv(5);
	  // gds.crossing_point(dis_v[0],dis_w[0],kVwire,kYwire, pwv[0]);
	  // gds.crossing_point(dis_v[0],dis_w[1],kVwire,kYwire, pwv[1]);
	  // gds.crossing_point(dis_v[1],dis_w[1],kVwire,kYwire, pwv[2]);
	  // gds.crossing_point(dis_v[1],dis_w[0],kVwire,kYwire, pwv[3]);
	  gds.crossing_point(dis_v[2],dis_w[2],kVwire,kYwire, pwv[4]);
	  dis_pwv[4] = gds.wire_dist(pwv[4],kUwire);

	  for (int k1=0;k1!=4;k1++){
	    pwv[k1] = pwv[4] + (pwv_save[k1]-pwv_save[4]);
	    dis_pwv[k1] = dis_pwv[4] +(dis_pwv_save[k1]-dis_pwv_save[4]);
	    //	    dis_pwv[k] = gds.wire_dist(pwv[k],kUwire);
	    if (dis_pwv[k1] > dis_u[0]-tolerance/2. && dis_pwv[k1] < dis_u[1]+tolerance/2.){
	      int flag_abc = 0;
	      for (int kk = 0; kk!=pcell.size();kk++){
	      	float dis = sqrt(pow(pwv[k1].y-pcell.at(kk).y,2) + pow(pwv[k1].z-pcell.at(kk).z,2));
	      	if (dis < tolerance) {
	      	  flag_abc = 1;
	      	  break;
	      	}
	      }
	      if (flag_abc == 0)
	      	pcell.push_back(pwv[k1]);
	    }
	  }

	  
	  if (pcell.size()>=3){
	    //order all the points by phi angle
	    
	    const GeomCell *cell = 0;
	    
	    // // old method
	    // GeomCell *cell_t = new GeomCell(ncell,pcell);
	    // cell_t->set_uwire(wire_u.at(i));
	    // cell_t->set_vwire(wire_v.at(j));
	    // cell_t->set_wwire(wire_w.at(k));
	    // cell = cell_t;

	    // new method
	    WireCell2dToy::TotalTiling& totaltiling = Singleton<WireCell2dToy::TotalTiling>::Instance();
	    GeomWireSelection temp_wires;
	    temp_wires.push_back(wire_u.at(i));
	    temp_wires.push_back(wire_v.at(j));
	    temp_wires.push_back(wire_w.at(k));
	    cell = totaltiling.cell(temp_wires);
	    if (cell == 0){
	      GeomCell *cell_t = new GeomCell(ncell,pcell);
	      cell_t->set_uwire(wire_u.at(i));
	      cell_t->set_vwire(wire_v.at(j));
	      cell_t->set_wwire(wire_w.at(k));
	      totaltiling.AddCellWire(cell_t, wire_u.at(i), wire_v.at(j), wire_w.at(k));
	      cell = cell_t;
	    }




	    Point cell_center = cell->center();
	    if (gds.contained_yz(cell_center)){
	      
	      // auto it1 = find(wire_all1.begin(),wire_all1.end(),wire_u[i]);
	      // auto it2 = find(wire_all1.begin(),wire_all1.end(),wire_v[j]);
	      // auto it3 = find(wire_all1.begin(),wire_all1.end(),wire_w[k]);
	      // if (it1 == wire_all1.end()){
	      // 	wire_all1.push_back(wire_u[i]);
	      // 	wire_u1.push_back(wire_u[i]);
	      // }
	      // if (it2 == wire_all1.end()){
	      // 	wire_all1.push_back(wire_v[j]);
	      // 	wire_v1.push_back(wire_v[j]);
	      // }
	      // if (it3 == wire_all1.end()){
	      // 	wire_all1.push_back(wire_w[k]);
	      // 	wire_w1.push_back(wire_w[k]);
	      // }
	      
	      //std::cout << i << " " << j << " " << k << " " << pcell.size() << " " << cell->center().z/units::m << " " << cell->center().y/units::m << std::endl;
	      
	      // fill cellmap
	      GeomWireSelection wiresel;
	      wiresel.push_back(wire_u[i]);
	      wiresel.push_back(wire_v[j]);
	      wiresel.push_back(wire_w[k]);	
	      cellmap[cell]=wiresel;
	      
	      //fill wiremap
	      if (wiremap.find(wire_u[i]) == wiremap.end()){
		//not found
		GeomCellSelection cellsel;
		cellsel.push_back(cell);
		wiremap[wire_u[i]]=cellsel;
	      }else{
		//found
		wiremap[wire_u[i]].push_back(cell);
	      }
	      
	      if (wiremap.find(wire_v[j]) == wiremap.end()){
		//not found
		GeomCellSelection cellsel;
		cellsel.push_back(cell);
		wiremap[wire_v[j]]=cellsel;
	      }else{
		//found
		wiremap[wire_v[j]].push_back(cell);
	      }
	      
	      if (wiremap.find(wire_w[k]) == wiremap.end()){
		//not found
		GeomCellSelection cellsel;
		cellsel.push_back(cell);
		wiremap[wire_w[k]]=cellsel;
	      }else{
		//found
		wiremap[wire_w[k]].push_back(cell);
	      }
	      
	      cell_all.push_back(cell);
	      ncell++;
	    }else{
	      delete cell;
	    }
	  }
	}
	
      } // W-loop
      
      // initialize uw and vw points
      // if (flag==1){
	// check order 
	//	pcell = cell->boundary();
	// std::cout << "Cell Count: " << pcell.size() << " " << cell->cross_section() << std::endl;
	// for (int k=0;k!=pcell.size();k++){
	//   std::cout << pcell[k].y << " " << pcell[k].z << " " << std::atan2(pcell[k].z - cell->center().z, pcell[k].y-cell->center().y) << std::endl;
	// }
      // }
    } // V-loop
  } // U-loop
  

  // wire_all.clear();
  // wire_u.clear();
  // wire_v.clear();
  // wire_w.clear();

  // wire_all = wire_all1;
  // wire_u = wire_u1;
  // wire_v = wire_v1;
  // wire_w = wire_w1;


  //std::cout << counter << " " << counter1 << " " << counter2 << " "<< counter3 << " " << wire_u.size() << " " << wire_v.size() << " " << wire_w.size() << std::endl;

  // std::cout << sum << std::endl;
  //std::cout << wire_u.size() << " " << wire_v.size() << " " << wire_w.size() << std::endl;
  
  //print the first cell and its wire
  // if (cell_all.size() > 0){
  //   GeomWireSelection qwires = cellmap[cell_all.at(0)];
  //   std::cout << qwires.at(0)->index() << " " << qwires.at(0)->plane() << " " 
  // 	      << qwires.at(1)->index() << " " << qwires.at(1)->plane() << " " 
  // 	      << qwires.at(2)->index() << " " << qwires.at(2)->plane() << " "  << std::endl;
  // }

  ave_charge = 0;
  for (int i=0;i!=wire_all.size();i++){
    ave_charge += wirechargemap[wire_all.at(i)];
  }
  if (cell_all.size()!=0){
    ave_charge = ave_charge / cell_all.size()/3.;
  }

}






void WireCell2dToy::ToyTiling::twoplane_tiling(int time, int nrebin, WireCell::GeomDataSource& gds, std::vector<float>& uplane_rms, std::vector<float>& vplane_rms, std::vector<float>& wplane_rms, WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map){
  float tolerance = 0.1 * units::mm;
  int ncell = 10000;
  


  GeomWireSelection nu_wires;
  GeomWireSelection nv_wires;
  GeomWireSelection nw_wires;
  int cut_sigma = 4;

  GeomWireSelection threeplane_wires;
  for (int i=0;i!=cell_all.size();i++){
    for (int j=0;j!=cellmap[cell_all.at(i)].size();j++){
      const GeomWire *wire = cellmap[cell_all.at(i)].at(j);
      auto it = find(threeplane_wires.begin(),threeplane_wires.end(),wire);
      if (it == threeplane_wires.end()){
	threeplane_wires.push_back(wire);
      }
    }
  }
  

  

  // calculate all the costant
  std::vector<float> udis,vdis,wdis;
  for (int i=0;i!=wire_u.size();i++){
    udis.push_back(gds.wire_dist(*wire_u[i]));
  }
  for (int j=0;j!=wire_v.size();j++){
    vdis.push_back(gds.wire_dist(*wire_v[j]));
  }
  for (int k=0;k!=wire_w.size();k++){
    wdis.push_back(gds.wire_dist(*wire_w[k]));
  }

  //  std::cout << wire_u.size() << " " << wire_v.size() << " " << wire_w.size() << std::endl;

   float dis_u[3],dis_v[3],dis_w[3],
    dis_puv[5],dis_puw[5],dis_pwv[5];

  double u_pitch, v_pitch, w_pitch;
  u_pitch = gds.pitch(kUwire);
  v_pitch = gds.pitch(kVwire);
  w_pitch = gds.pitch(kYwire);

  

  // U-V plane first
  for (int i=0;i!=wire_u.size();i++){
    const GeomWire *wire1 = wire_u.at(i);
    int channel1 = wire1->index();
    if (wirechargemap[wire1] < cut_sigma * uplane_rms.at(channel1)) continue;
    dis_u[0] = udis.at(i) - u_pitch/2.;
    dis_u[1] = dis_u[0] + u_pitch;
    dis_u[2] = udis.at(i);

    auto qt1 = find(threeplane_wires.begin(),threeplane_wires.end(),wire1);

    for (int j=0;j!=wire_v.size();j++){
      const GeomWire *wire2 = wire_v.at(j);
      int channel2 = wire2->index();
      if (wirechargemap[wire2] < cut_sigma * vplane_rms.at(channel2)) continue;
      dis_v[0] = vdis.at(j) - v_pitch/2.;
      dis_v[1] = dis_v[0] + v_pitch;
      dis_v[2] = vdis.at(j);

      auto qt2 = find(threeplane_wires.begin(),threeplane_wires.end(),wire2);

      //std::cout << i << " " << j << std::endl;

      int flag = 1;
      //check all the cells if contain both wires, go on
      for (int k=0;k!=wiremap[wire1].size();k++){
	const GeomCell *cell = wiremap[wire1].at(k);
	GeomWireSelection temp_wires = cellmap[cell];
	auto it = find(temp_wires.begin(),temp_wires.end(),wire2);
	if (it != temp_wires.end()){
	  flag = 0;
	  break;
	}
      }

      //      std::cout << i << " " << j << " " << flag << std::endl;

      if (flag == 1){
	//if (qt1!=threeplane_wires.end() && qt2!=threeplane_wires.end()) continue;

	// start to fill in
	std::vector<Vector> puv(5);
	if(!gds.crossing_point(dis_u[2],dis_v[2],kUwire,kVwire, puv[4])) continue;
	gds.crossing_point(dis_u[0],dis_v[0],kUwire,kVwire, puv[0]);
	gds.crossing_point(dis_u[0],dis_v[1],kUwire,kVwire, puv[1]);
	gds.crossing_point(dis_u[1],dis_v[1],kUwire,kVwire, puv[2]);
	gds.crossing_point(dis_u[1],dis_v[0],kUwire,kVwire, puv[3]);
	
	for (int a1 = 0;a1!=5;a1++){
	  const GeomWire *n_wire = gds.closest(puv[a1],kYwire);
	  if (n_wire == 0) continue;
	  auto it1 = find(nw_wires.begin(),nw_wires.end(),n_wire);
	  auto it2 = wplane_map.find(n_wire->index());

	  	
	  //if (it1 == nw_wires.end() && it2 == wplane_map.end()){
	  if (it1 == nw_wires.end() && 
	      (it2!=wplane_map.end() && (time <= wplane_map[n_wire->index()].second /nrebin && time >= wplane_map[n_wire->index()].first /nrebin ) )){
	    
	    nw_wires.push_back(n_wire);
	    dis_w[0] = gds.wire_dist(*n_wire) - w_pitch/2.;
	    dis_w[1] = dis_w[0] + w_pitch;
	    dis_w[2] = dis_w[0] + w_pitch/2.;
	    
	    PointVector pcell;
	    
	    for (int m = 0;m!=4;m++){
	      if (dis_puv[m] > dis_w[0]-tolerance/2. && dis_puv[m] < dis_w[1]+tolerance/2.){
		flag = 1;
		pcell.push_back(puv[m]);
	      }
	    }

	    //	    std::cout << a1 << " " << " UW" << std::endl;

	    //form a new cell and establish all the map ... 	    
	    std::vector<Vector> puw(5);
	    gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, puw[0]);
	    gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, puw[1]);
	    gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, puw[2]);
	    gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, puw[3]);
	    for (int k1=0;k1!=4;k1++){
	      dis_puw[k1] = gds.wire_dist(puw[k1],kVwire);
	      if (dis_puw[k1] > dis_v[0]-tolerance/2. && dis_puw[k1] < dis_v[1]+tolerance/2.){
		int flag_abc = 0;
		for (int kk = 0; kk!=pcell.size();kk++){
		  float dis = sqrt(pow(puw[k1].y-pcell.at(kk).y,2) + pow(puw[k1].z-pcell.at(kk).z,2));
		  if (dis < tolerance) {
		    flag_abc = 1;
		    break;
		  }
		}
		if (flag_abc == 0)
		  pcell.push_back(puw[k1]);
	      }
	    }
	    

	    //	    std::cout << a1 << " " << " VW" << std::endl;

	    std::vector<Vector> pwv(5);
	    gds.crossing_point(dis_v[0],dis_w[0],kVwire,kYwire, pwv[0]);
	    gds.crossing_point(dis_v[0],dis_w[1],kVwire,kYwire, pwv[1]);
	    gds.crossing_point(dis_v[1],dis_w[1],kVwire,kYwire, pwv[2]);
	    gds.crossing_point(dis_v[1],dis_w[0],kVwire,kYwire, pwv[3]);

	    for (int k1=0;k1!=4;k1++){
	      dis_pwv[k1] = gds.wire_dist(pwv[k1],kUwire);
	      if (dis_pwv[k1] > dis_u[0]-tolerance/2. && dis_pwv[k1] < dis_u[1]+tolerance/2.){
		int flag_abc = 0;
		for (int kk = 0; kk!=pcell.size();kk++){
		  float dis = sqrt(pow(pwv[k1].y-pcell.at(kk).y,2) + pow(pwv[k1].z-pcell.at(kk).z,2));
		  if (dis < tolerance) {
		    flag_abc = 1;
		    break;
		  }
		}
		if (flag_abc == 0)
		  pcell.push_back(pwv[k1]);
	      }
	    }

	    if (pcell.size()>=3){
	      
	      //	      std::cout << " Form Cell" << std::endl;
	      const GeomCell *cell = 0;

	      // //old 
	      // GeomCell *cell_t = new GeomCell(ncell,pcell);
	      // cell_t->set_uwire(wire_u[i]);
	      // cell_t->set_vwire(wire_v[j]);
	      // cell_t->set_wwire(n_wire);
	      // cell = cell_t;
	      
	      // new method
	      WireCell2dToy::TotalTiling& totaltiling = Singleton<WireCell2dToy::TotalTiling>::Instance();
	      GeomWireSelection temp_wires;
	      temp_wires.push_back(wire_u.at(i));
	      temp_wires.push_back(wire_v.at(j));
	      temp_wires.push_back(n_wire);
	      cell = totaltiling.cell(temp_wires);
	      if (cell == 0){
		GeomCell *cell_t = new GeomCell(ncell,pcell);
		cell_t->set_uwire(wire_u.at(i));
		cell_t->set_vwire(wire_v.at(j));
		cell_t->set_wwire(n_wire);
		totaltiling.AddCellWire(cell_t, wire_u.at(i), wire_v.at(j), n_wire);
		cell = cell_t;
	      }




	      Point cell_center = cell->center();


	      //	      std::cout << " Form Wire/Cell Map" << std::endl;

	      GeomWireSelection wiresel;
	      wiresel.push_back(wire_u[i]);
	      wiresel.push_back(wire_v[j]);
	      wiresel.push_back(n_wire);	
	      cellmap[cell]=wiresel;

	      
	      
	       //fill wiremap
	      if (wiremap.find(wire_u[i]) == wiremap.end()){
		//not found
		GeomCellSelection cellsel;
		cellsel.push_back(cell);
		wiremap[wire_u[i]]=cellsel;
	      }else{
		//found
		wiremap[wire_u[i]].push_back(cell);
	      }
	      
	      if (wiremap.find(wire_v[j]) == wiremap.end()){
		//not found
		GeomCellSelection cellsel;
		cellsel.push_back(cell);
		wiremap[wire_v[j]]=cellsel;
	      }else{
		//found
		wiremap[wire_v[j]].push_back(cell);
	      }

	       if (wiremap.find(n_wire) == wiremap.end()){
		//not found
		GeomCellSelection cellsel;
		cellsel.push_back(cell);
		wiremap[n_wire]=cellsel;
	      }else{
		//found
		wiremap[n_wire].push_back(cell);
	      }


	      
	      cell_all.push_back(cell);
	      ncell++;

	    }


	  }
	}
	
      }

    }
  }


  // U-W plane second
  for (int i=0;i!=wire_u.size();i++){
    const GeomWire *wire1 = wire_u.at(i);
    int channel1 = wire1->index();
    if (wirechargemap[wire1] < cut_sigma * uplane_rms.at(channel1)) continue;
    dis_u[0] = udis.at(i) - u_pitch/2.;
    dis_u[1] = dis_u[0] + u_pitch;
    dis_u[2] = udis.at(i);

    auto qt1 = find(threeplane_wires.begin(),threeplane_wires.end(),wire1);

    for (int j=0;j!=wire_w.size();j++){
      const GeomWire *wire2 = wire_w.at(j);
      int channel2 = wire2->index();
      if (wirechargemap[wire2] < cut_sigma * wplane_rms.at(channel2)) continue;
      dis_w[0] = wdis.at(j) - w_pitch/2.;
      dis_w[1] = dis_w[0] + w_pitch;
      dis_w[2] = wdis.at(j);

      auto qt2 = find(threeplane_wires.begin(),threeplane_wires.end(),wire2);

      int flag = 1;
      //check all the cells if contain both wires, go on
      for (int k=0;k!=wiremap[wire1].size();k++){
  	const GeomCell *cell = wiremap[wire1].at(k);
  	GeomWireSelection temp_wires = cellmap[cell];
  	auto it = find(temp_wires.begin(),temp_wires.end(),wire2);
  	if (it != temp_wires.end()){
  	  flag = 0;
  	  break;
  	}
      }
      
      
      if (flag == 1){
	//if (qt1!=threeplane_wires.end() && qt2!=threeplane_wires.end()) continue;


  	// start to fill in
  	std::vector<Vector> puw(5);
  	if(!gds.crossing_point(dis_u[2],dis_w[2],kUwire,kYwire, puw[4])) continue;
  	gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, puw[0]);
  	gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, puw[1]);
  	gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, puw[2]);
  	gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, puw[3]);
	
  	for (int a1 = 0;a1!=5;a1++){
  	  const GeomWire *n_wire = gds.closest(puw[a1],kVwire);
	  if (n_wire == 0 ) continue;
  	  auto it1 = find(nv_wires.begin(),nv_wires.end(),n_wire);
	  auto it2 = vplane_map.find(n_wire->index());

	   if (it1 == nv_wires.end() && 
	       (it2!=vplane_map.end() && (time <= vplane_map[n_wire->index()].second /nrebin && time >= vplane_map[n_wire->index()].first /nrebin ) )){

	     //  	  if (it1 == nv_wires.end() && it2 == vplane_map.end()){
  	    nv_wires.push_back(n_wire);
  	    dis_v[0] = gds.wire_dist(*n_wire) - v_pitch/2.;
  	    dis_v[1] = dis_v[0] + v_pitch;
  	    dis_v[2] = dis_v[0] + v_pitch/2.;
	    	    

  	    PointVector pcell;
  	    for (int m = 0;m!=4;m++){
  	      if (dis_puw[m] > dis_v[0]-tolerance/2. && dis_puw[m] < dis_v[1]+tolerance/2.){
  		flag = 1;
  		pcell.push_back(puw[m]);
  	      }
  	    }

  	    //form a new cell and establish all the map ... 	    
  	    std::vector<Vector> puv(5);
  	    gds.crossing_point(dis_u[0],dis_v[0],kUwire,kVwire, puv[0]);
  	    gds.crossing_point(dis_u[0],dis_v[1],kUwire,kVwire, puv[1]);
  	    gds.crossing_point(dis_u[1],dis_v[1],kUwire,kVwire, puv[2]);
  	    gds.crossing_point(dis_u[1],dis_v[0],kUwire,kVwire, puv[3]);

  	    for (int k1=0;k1!=4;k1++){
  	      dis_puv[k1] = gds.wire_dist(puv[k1],kYwire);
  	      if (dis_puv[k1] > dis_w[0]-tolerance/2. && dis_puv[k1] < dis_w[1]+tolerance/2.){
  		int flag_abc = 0;
  		for (int kk = 0; kk!=pcell.size();kk++){
  		  float dis = sqrt(pow(puv[k1].y-pcell.at(kk).y,2) + pow(puv[k1].z-pcell.at(kk).z,2));
  		  if (dis < tolerance) {
  		    flag_abc = 1;
  		    break;
  		  }
  		}
  		if (flag_abc == 0)
  		  pcell.push_back(puv[k1]);
  	      }
  	    }
	    

  	    std::vector<Vector> pwv(5);
  	    gds.crossing_point(dis_v[0],dis_w[0],kVwire,kYwire, pwv[0]);
  	    gds.crossing_point(dis_v[0],dis_w[1],kVwire,kYwire, pwv[1]);
  	    gds.crossing_point(dis_v[1],dis_w[1],kVwire,kYwire, pwv[2]);
  	    gds.crossing_point(dis_v[1],dis_w[0],kVwire,kYwire, pwv[3]);

  	    for (int k1=0;k1!=4;k1++){
  	      dis_pwv[k1] = gds.wire_dist(pwv[k1],kUwire);
  	      if (dis_pwv[k1] > dis_u[0]-tolerance/2. && dis_pwv[k1] < dis_u[1]+tolerance/2.){
  		int flag_abc = 0;
  		for (int kk = 0; kk!=pcell.size();kk++){
  		  float dis = sqrt(pow(pwv[k1].y-pcell.at(kk).y,2) + pow(pwv[k1].z-pcell.at(kk).z,2));
  		  if (dis < tolerance) {
  		    flag_abc = 1;
  		    break;
  		  }
  		}
  		if (flag_abc == 0)
  		  pcell.push_back(pwv[k1]);
  	      }
  	    }

  	    if (pcell.size()>=3){

	      const GeomCell *cell = 0;

	      // // old 
  	      // GeomCell *cell_t = new GeomCell(ncell,pcell);
  	      // cell_t->set_uwire(wire_u[i]);
	      // cell_t->set_vwire(n_wire);
	      // cell_t->set_wwire(wire_w[j]);
	      // cell = cell_t;
	      
	      // new method
	      WireCell2dToy::TotalTiling& totaltiling = Singleton<WireCell2dToy::TotalTiling>::Instance();
	      GeomWireSelection temp_wires;
	      temp_wires.push_back(wire_u.at(i));
	      temp_wires.push_back(n_wire);
	      temp_wires.push_back(wire_w.at(j));
	      cell = totaltiling.cell(temp_wires);
	      if (cell == 0){
		GeomCell *cell_t = new GeomCell(ncell,pcell);
		cell_t->set_uwire(wire_u.at(i));
		cell_t->set_vwire(n_wire);
		cell_t->set_wwire(wire_w.at(j));
		totaltiling.AddCellWire(cell_t, wire_u.at(i), n_wire, wire_w.at(j));
		cell = cell_t;
	      }



	      Point cell_center = cell->center();

  	      GeomWireSelection wiresel;
  	      wiresel.push_back(wire_u[i]);
  	      wiresel.push_back(wire_w[j]);
  	      wiresel.push_back(n_wire);	
  	      cellmap[cell]=wiresel;

  	       //fill wiremap
  	      if (wiremap.find(wire_u[i]) == wiremap.end()){
  		//not found
  		GeomCellSelection cellsel;
  		cellsel.push_back(cell);
  		wiremap[wire_u[i]]=cellsel;
  	      }else{
  		//found
  		wiremap[wire_u[i]].push_back(cell);
  	      }
	      
  	      if (wiremap.find(wire_w[j]) == wiremap.end()){
  		//not found
  		GeomCellSelection cellsel;
  		cellsel.push_back(cell);
  		wiremap[wire_w[j]]=cellsel;
  	      }else{
  		//found
  		wiremap[wire_w[j]].push_back(cell);
  	      }

  	       if (wiremap.find(n_wire) == wiremap.end()){
  		//not found
  		GeomCellSelection cellsel;
  		cellsel.push_back(cell);
  		wiremap[n_wire]=cellsel;
  	      }else{
  		//found
  		wiremap[n_wire].push_back(cell);
  	      }


	      
  	      cell_all.push_back(cell);
  	      ncell++;

  	    }


  	  }
  	}
	
      }

    }
  }


   // W-V plane first
  for (int i=0;i!=wire_w.size();i++){
    const GeomWire *wire1 = wire_w.at(i);
    int channel1 = wire1->index();
    if (wirechargemap[wire1] < cut_sigma * wplane_rms.at(channel1)) continue;
    dis_w[0] = wdis.at(i) - w_pitch/2.;
    dis_w[1] = dis_w[0] + w_pitch;
    dis_w[2] = wdis.at(i);

    auto qt1 = find(threeplane_wires.begin(),threeplane_wires.end(),wire1);

    for (int j=0;j!=wire_v.size();j++){
      const GeomWire *wire2 = wire_v.at(j);
      int channel2 = wire2->index();
      if (wirechargemap[wire2] < cut_sigma * vplane_rms.at(channel2)) continue;
      dis_v[0] = vdis.at(j) - v_pitch/2.;
      dis_v[1] = dis_v[0] + v_pitch;
      dis_v[2] = vdis.at(j);

      auto qt2 = find(threeplane_wires.begin(),threeplane_wires.end(),wire2);

      int flag = 1;
      //check all the cells if contain both wires, go on
      for (int k=0;k!=wiremap[wire1].size();k++){
  	const GeomCell *cell = wiremap[wire1].at(k);
  	GeomWireSelection temp_wires = cellmap[cell];
  	auto it = find(temp_wires.begin(),temp_wires.end(),wire2);
  	if (it != temp_wires.end()){
  	  flag = 0;
  	  break;
  	}
      }
      

      if (flag == 1){
	//	if (qt1!=threeplane_wires.end() && qt2!=threeplane_wires.end()) continue;


  	// start to fill in
  	std::vector<Vector> pwv(5);
  	if(!gds.crossing_point(dis_w[2],dis_v[2],kYwire,kVwire, pwv[4])) continue;
  	gds.crossing_point(dis_w[0],dis_v[0],kYwire,kVwire, pwv[0]);
  	gds.crossing_point(dis_w[0],dis_v[1],kYwire,kVwire, pwv[1]);
  	gds.crossing_point(dis_w[1],dis_v[1],kYwire,kVwire, pwv[2]);
  	gds.crossing_point(dis_w[1],dis_v[0],kYwire,kVwire, pwv[3]);
	
  	for (int a1 = 0;a1!=5;a1++){
  	  const GeomWire *n_wire = gds.closest(pwv[a1],kUwire);
	  if (n_wire ==0 ) continue;
  	  auto it1 = find(nu_wires.begin(),nu_wires.end(),n_wire);
	  auto it2 = uplane_map.find(n_wire->index());
  	  if (it1 == nu_wires.end() && 
	      (it2!=uplane_map.end() && (time >= uplane_map[n_wire->index()].first /nrebin && time <= uplane_map[n_wire->index()].second /nrebin ) )){
  	    nu_wires.push_back(n_wire);
  	    dis_u[0] = gds.wire_dist(*n_wire) - u_pitch/2.;
  	    dis_u[1] = dis_u[0] + u_pitch;
  	    dis_u[2] = dis_u[0] + u_pitch/2.;
	    
	    
  	    PointVector pcell;
  	    for (int m = 0;m!=4;m++){
  	      if (dis_pwv[m] > dis_u[0]-tolerance/2. && dis_pwv[m] < dis_u[1]+tolerance/2.){
  		flag = 1;
  		pcell.push_back(pwv[m]);
  	      }
  	    }

  	    //form a new cell and establish all the map ... 	    
  	    std::vector<Vector> puw(5);
  	    gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, puw[0]);
  	    gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, puw[1]);
  	    gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, puw[2]);
  	    gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, puw[3]);
  	    for (int k1=0;k1!=4;k1++){
  	      dis_puw[k1] = gds.wire_dist(puw[k1],kVwire);
  	      if (dis_puw[k1] > dis_v[0]-tolerance/2. && dis_puw[k1] < dis_v[1]+tolerance/2.){
  		int flag_abc = 0;
  		for (int kk = 0; kk!=pcell.size();kk++){
  		  float dis = sqrt(pow(puw[k1].y-pcell.at(kk).y,2) + pow(puw[k1].z-pcell.at(kk).z,2));
  		  if (dis < tolerance) {
  		    flag_abc = 1;
  		    break;
  		  }
  		}
  		if (flag_abc == 0)
  		  pcell.push_back(puw[k1]);
  	      }
  	    }
	    

  	    std::vector<Vector> puv(5);
  	    gds.crossing_point(dis_v[0],dis_u[0],kVwire,kUwire, puv[0]);
  	    gds.crossing_point(dis_v[0],dis_u[1],kVwire,kUwire, puv[1]);
  	    gds.crossing_point(dis_v[1],dis_u[1],kVwire,kUwire, puv[2]);
  	    gds.crossing_point(dis_v[1],dis_u[0],kVwire,kUwire, puv[3]);

  	    for (int k1=0;k1!=4;k1++){
  	      dis_puv[k1] = gds.wire_dist(puv[k1],kYwire);
  	      if (dis_puv[k1] > dis_w[0]-tolerance/2. && dis_puv[k1] < dis_w[1]+tolerance/2.){
  		int flag_abc = 0;
  		for (int kk = 0; kk!=pcell.size();kk++){
  		  float dis = sqrt(pow(puv[k1].y-pcell.at(kk).y,2) + pow(puv[k1].z-pcell.at(kk).z,2));
  		  if (dis < tolerance) {
  		    flag_abc = 1;
  		    break;
  		  }
  		}
  		if (flag_abc == 0)
  		  pcell.push_back(puv[k1]);
  	      }
  	    }

  	    if (pcell.size()>=3){
	      
	      const GeomCell *cell = 0;
	      // //old
  	      // GeomCell *cell_t = new GeomCell(ncell,pcell);
  	      // cell_t->set_uwire(n_wire);
	      // cell_t->set_vwire(wire_v[j]);
	      // cell_t->set_wwire(wire_w[i]);
	      // cell = cell_t;

	      // new method
	      WireCell2dToy::TotalTiling& totaltiling = Singleton<WireCell2dToy::TotalTiling>::Instance();
	      GeomWireSelection temp_wires;
	      temp_wires.push_back(n_wire);
	      temp_wires.push_back(wire_v.at(j));
	      temp_wires.push_back(wire_w.at(i));
	      cell = totaltiling.cell(temp_wires);
	      if (cell == 0){
		GeomCell *cell_t = new GeomCell(ncell,pcell);
		cell_t->set_uwire(n_wire);
		cell_t->set_vwire(wire_v.at(j));
		cell_t->set_wwire(wire_w.at(i));
		totaltiling.AddCellWire(cell_t, n_wire, wire_v.at(j), wire_w.at(i));
		cell = cell_t;
	      }
	      
	      Point cell_center = cell->center();

  	      GeomWireSelection wiresel;
  	      wiresel.push_back(wire_w[i]);
  	      wiresel.push_back(wire_v[j]);
  	      wiresel.push_back(n_wire);	
  	      cellmap[cell]=wiresel;

  	       //fill wiremap
  	      if (wiremap.find(wire_w[i]) == wiremap.end()){
  		//not found
  		GeomCellSelection cellsel;
  		cellsel.push_back(cell);
  		wiremap[wire_w[i]]=cellsel;
  	      }else{
  		//found
  		wiremap[wire_w[i]].push_back(cell);
  	      }
	      
  	      if (wiremap.find(wire_v[j]) == wiremap.end()){
  		//not found
  		GeomCellSelection cellsel;
  		cellsel.push_back(cell);
  		wiremap[wire_v[j]]=cellsel;
  	      }else{
  		//found
  		wiremap[wire_v[j]].push_back(cell);
  	      }

  	       if (wiremap.find(n_wire) == wiremap.end()){
  		//not found
  		GeomCellSelection cellsel;
  		cellsel.push_back(cell);
  		wiremap[n_wire]=cellsel;
  	      }else{
  		//found
  		wiremap[n_wire].push_back(cell);
  	      }


	      
  	      cell_all.push_back(cell);
  	      ncell++;

  	    }


  	  }
  	}
	
      }

    }
  }


  



  for (int i=0;i!=nu_wires.size();i++){
    wire_u.push_back(nu_wires.at(i));
    wire_all.push_back(nu_wires.at(i));

    wirechargemap[nu_wires.at(i)] = 0;
    wirecharge_errmap[nu_wires.at(i)] = 1e4;
  }
  for (int i=0;i!=nv_wires.size();i++){
    wire_v.push_back(nv_wires.at(i));
    wire_all.push_back(nv_wires.at(i));

    wirechargemap[nv_wires.at(i)] = 0;
    wirecharge_errmap[nv_wires.at(i)] = 1e4;
  }
  for (int i=0;i!=nw_wires.size();i++){
    wire_w.push_back(nw_wires.at(i));
    wire_all.push_back(nw_wires.at(i));

    wirechargemap[nw_wires.at(i)] = 0;
    wirecharge_errmap[nw_wires.at(i)] = 1e4;
  }
				    
  
  ave_charge = 0;
  for (int i=0;i!=wire_all.size();i++){
    ave_charge += wirechargemap[wire_all.at(i)];
  }
  if (cell_all.size()!=0){
    ave_charge = ave_charge / cell_all.size()/3.;
  }
}




WireCell2dToy::ToyTiling::~ToyTiling()
{
  //delete all the cells
  for (int i=0;i!=cell_all.size();i++){
    delete cell_all[i];
  }

  wire_u.clear();
  wire_v.clear();
  wire_w.clear();
  wire_all.clear();
  cell_all.clear();
  cellmap.clear();
  wiremap.clear();
  wirechargemap.clear();
  
}

GeomWireSelection WireCell2dToy::ToyTiling::wires(const GeomCell& cell) const
{
  if (cellmap.find(&cell) == cellmap.end()){
    //not found 
    return GeomWireSelection();
  }else{
    //found
    return cellmap.find(&cell)->second;
  }
    
}
	
GeomCellSelection WireCell2dToy::ToyTiling::cells(const GeomWire& wire) const
{
  if (wiremap.find(&wire) == wiremap.end()){
    return GeomCellSelection();
  }else{
    return wiremap.find(&wire)->second;
  }
}


const GeomCell* WireCell2dToy::ToyTiling::cell(const GeomWireSelection& wires) const
{
  if (wires.size()!=3) return 0;
  const GeomWire *wire1 = wires[0];
  const GeomWire *wire2 = wires[1];
  const GeomWire *wire3 = wires[2];

  if (wire1->plane() == wire2->plane() ||
      wire1->plane() == wire3->plane() || 
      wire2->plane() == wire3->plane()) return 0;

  // GeomCellSelection cells1 = cells(*wire1);
  // // GeomCellSelection cells2 = cells(*wire2);
  // // GeomCellSelection cells3 = cells(*wire3);
  
  // if (cells1.size() < cells(*wire2).size())
  //   cells1 = cells(*wire2);
  // if (cells1.size() < cells(*wire3).size())
  //   cells1 = cells(*wire3);
  
  // for (int i = 0; i!=cells1.size(); i++){
  //   const GeomCell *cell1 = cells1[i];
  //   GeomWireSelection twires = wires(*cell1);
  //   // for (int j =0; j!=cells2.size(); j++){
  //   //   const GeomCell *cell2 = cells2[j];
  //   //   if (*cell1==*cell2){
  //   auto it1 = find(twires.begin(),twires.end(),wire1);
  //   auto it2 = find(twires.begin(),twires.end(),wire2);
  //   auto it3 = find(twires.begin(),twires.end(),wire3);
  //   if (it1 != twires.end() && it2 != twires.end() && it3 != twires.end()){
  //     return cell1;
  //   }
  //   // for (int k=0;k!=cells3.size();k++){
  //   //   const GeomCell *cell3 = cells3[k];
  //   //   if (*cell1 == *cell3){
  //   //     //there is a problem here, not sure what to do 
  //   //     return cell1;
  //   //     //return 0;
  //   //   }
  //   // }
  //   // }
  //   // }
  // }

  GeomCellSelection cells1 = cells(*wire1);
  GeomCellSelection cells2 = cells(*wire2);
  GeomCellSelection cells3 = cells(*wire3);
  
  for (int i = 0; i!=cells1.size(); i++){
    const GeomCell *cell1 = cells1[i];
    for (int j =0; j!=cells2.size(); j++){
      const GeomCell *cell2 = cells2[j];
      if (*cell1==*cell2){
  	for (int k=0;k!=cells3.size();k++){
  	  const GeomCell *cell3 = cells3[k];
  	  if (*cell1 == *cell3){
  	    //there is a problem here, not sure what to do 
  	    return cell1;
  	    //return 0;
  	  }
  	}
      }
    }
  }

  return 0;

}



void WireCell2dToy::ToyTiling::CreateCell(float tolerance, const GeomDataSource& gds, int face, int n_tpc, GeomWireSelection& temp_wire_u, 
					  GeomWireSelection& temp_wire_v, GeomWireSelection& temp_wire_w){
   float dis_u[3],dis_v[3],dis_w[3],
    dis_puv[5],dis_puw[5],dis_pwv[5];
  
   
  // calculate all the costant
  std::vector<float> udis,vdis,wdis;
  for (int i=0;i!=temp_wire_u.size();i++){
    udis.push_back(gds.wire_dist(*temp_wire_u[i]));
  }
  for (int j=0;j!=temp_wire_v.size();j++){
    vdis.push_back(gds.wire_dist(*temp_wire_v[j]));
  }
  for (int k=0;k!=temp_wire_w.size();k++){
    wdis.push_back(gds.wire_dist(*temp_wire_w[k]));
  }

  //precalculate all the offsets
  std::vector<Vector> puv_save(5), puw_save(5), pwv_save(5);
  float dis_puv_save[5], dis_puw_save[5], dis_pwv_save[5];
  
  double u_pitch, v_pitch, w_pitch;
  u_pitch = gds.pitch(face,kUwire);
  v_pitch = gds.pitch(face,kVwire);
  w_pitch = gds.pitch(face,kYwire);

  //std::cout << u_pitch << " " << v_pitch << " " << w_pitch << std::endl;

  if (temp_wire_u.size()>=1&&temp_wire_v.size()>=1){
    dis_u[0] = udis.at(0) - gds.pitch(face,kUwire)/2.;
    dis_u[1] = dis_u[0] + gds.pitch(face,kUwire);
    dis_u[2] = udis.at(0);

    dis_v[0] = vdis.at(0) - gds.pitch(face,kVwire)/2.;
    dis_v[1] = dis_v[0] + gds.pitch(face,kVwire);
    dis_v[2] = vdis.at(0);

    gds.crossing_point(dis_u[0],dis_v[0],kUwire,kVwire, puv_save[0],face);
    gds.crossing_point(dis_u[0],dis_v[1],kUwire,kVwire, puv_save[1],face);
    gds.crossing_point(dis_u[1],dis_v[1],kUwire,kVwire, puv_save[2],face);
    gds.crossing_point(dis_u[1],dis_v[0],kUwire,kVwire, puv_save[3],face);
    gds.crossing_point(dis_u[2],dis_v[2],kUwire,kVwire, puv_save[4],face);

    for (int k=0;k!=5;k++){
      dis_puv_save[k] = gds.wire_dist(puv_save[k],kYwire,face);
    }
    
  }

  if (temp_wire_u.size()>=1&&temp_wire_w.size()>=1){
    dis_u[0] = udis.at(0) - gds.pitch(face,kUwire)/2.;
    dis_u[1] = dis_u[0] + gds.pitch(face,kUwire);
    dis_u[2] = udis.at(0);

    dis_w[0] = wdis.at(0) - w_pitch/2.;
    dis_w[1] = dis_w[0] + w_pitch;
    dis_w[2] = wdis.at(0);

    gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, puw_save[0],face);
    gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, puw_save[1],face);
    gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, puw_save[2],face);
    gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, puw_save[3],face);
    gds.crossing_point(dis_u[2],dis_w[2],kUwire,kYwire, puw_save[4],face);

    for (int k=0;k!=5;k++){
      dis_puw_save[k] = gds.wire_dist(puw_save[k],kVwire,face);
    }
  }

  if (temp_wire_v.size()>=1&&temp_wire_w.size()>=1){

    dis_w[0] = wdis.at(0) - w_pitch/2.;
    dis_w[1] = dis_w[0] + w_pitch;
    dis_w[2] = wdis.at(0);

    dis_v[0] = vdis.at(0) - gds.pitch(face,kVwire)/2.;
    dis_v[1] = dis_v[0] + gds.pitch(face,kVwire);
    dis_v[2] = vdis.at(0);
    
    gds.crossing_point(dis_v[0],dis_w[0],kVwire,kYwire, pwv_save[0],face);
    gds.crossing_point(dis_v[0],dis_w[1],kVwire,kYwire, pwv_save[1],face);
    gds.crossing_point(dis_v[1],dis_w[1],kVwire,kYwire, pwv_save[2],face);
    gds.crossing_point(dis_v[1],dis_w[0],kVwire,kYwire, pwv_save[3],face);
    gds.crossing_point(dis_v[2],dis_w[2],kVwire,kYwire, pwv_save[4],face);
    
    for (int k=0;k!=5;k++){
      dis_pwv_save[k] = gds.wire_dist(pwv_save[k],kUwire,face);
    }
  }

  for (int i=0;i!=temp_wire_u.size();i++){
    dis_u[0] = udis.at(i) - u_pitch/2.;
    dis_u[1] = dis_u[0] + u_pitch;
    dis_u[2] = udis.at(i);

    for (int j=0;j!=temp_wire_v.size();j++){
      dis_v[0] = vdis.at(j) - v_pitch/2.;
      dis_v[1] = dis_v[0] + v_pitch;
      dis_v[2] = vdis.at(j);
      
      //four vertices around
      std::vector<Vector> puv(5);
      if(!gds.crossing_point(dis_u[2],dis_v[2],kUwire,kVwire, puv[4],face)) continue;

      dis_puv[4] = gds.wire_dist(puv[4],kYwire,face);
      for (int k=0;k!=4;k++){
 	puv[k] = puv[4] + (puv_save[k]-puv_save[4]);
 	dis_puv[k] = dis_puv[4] +(dis_puv_save[k]-dis_puv_save[4]);
      }
       
      for (int k=0;k!=temp_wire_w.size();k++){
 	int flag = 0;
 	PointVector pcell;
 	dis_w[0] = wdis.at(k) - w_pitch/2.;
   	dis_w[1] = dis_w[0] + w_pitch;//gds.wire_dist(*temp_wire_w[k]) + gds.pitch(kYwire)/2.;	
 	dis_w[2] = wdis.at(k);

 	if (fabs(dis_w[0] - dis_puv[0])>3*units::cm) continue;
	
 	for (int m = 0;m!=4;m++){
 	  if (dis_puv[m] > dis_w[0]-tolerance/2. && dis_puv[m] < dis_w[1]+tolerance/2.){
 	    flag = 1;
 	    pcell.push_back(puv[m]);
 	  }
 	}

 	if (flag==1 ) {
	  std::vector<Vector> puw(5);
	  gds.crossing_point(dis_u[2],dis_w[2],kUwire,kYwire, puw[4],face);
 	  dis_puw[4] = gds.wire_dist(puw[4],kVwire,face);
	  for (int k1=0;k1!=4;k1++){
 	    puw[k1] = puw[4] + (puw_save[k1]-puw_save[4]);
 	    dis_puw[k1] = dis_puw[4] +(dis_puw_save[k1]-dis_puw_save[4]);
 	    if (dis_puw[k1] > dis_v[0]-tolerance/2. && dis_puw[k1] < dis_v[1]+tolerance/2.){
 	      int flag_abc = 0;
	      for (int kk = 0; kk!=pcell.size();kk++){
 	      	float dis = sqrt(pow(puw[k1].y-pcell.at(kk).y,2) + pow(puw[k1].z-pcell.at(kk).z,2));
 	      	if (dis < tolerance) {
 	      	  flag_abc = 1;
 	      	  break;
 	      	}
 	      }
 	      if (flag_abc == 0)
 	    	pcell.push_back(puw[k1]);
 	    }
 	  }
 	  std::vector<Vector> pwv(5);
 	  gds.crossing_point(dis_v[2],dis_w[2],kVwire,kYwire, pwv[4],face);
 	  dis_pwv[4] = gds.wire_dist(pwv[4],kUwire,face);

 	  for (int k1=0;k1!=4;k1++){
 	    pwv[k1] = pwv[4] + (pwv_save[k1]-pwv_save[4]);
 	    dis_pwv[k1] = dis_pwv[4] +(dis_pwv_save[k1]-dis_pwv_save[4]);
	    if (dis_pwv[k1] > dis_u[0]-tolerance/2. && dis_pwv[k1] < dis_u[1]+tolerance/2.){
 	      int flag_abc = 0;
 	      for (int kk = 0; kk!=pcell.size();kk++){
 	      	float dis = sqrt(pow(pwv[k1].y-pcell.at(kk).y,2) + pow(pwv[k1].z-pcell.at(kk).z,2));
 	      	if (dis < tolerance) {
 	      	  flag_abc = 1;
 	      	  break;
 	      	}
 	      }
 	      if (flag_abc == 0)
 	      	pcell.push_back(pwv[k1]);
 	    }
 	  }

	  
 	  if (pcell.size()>=3){
	    
	    const GeomCell *cell = 0;
	    
	    //old
 	    GeomCell *cell_t = new GeomCell(ncell,pcell);
 	    cell_t->set_uwire(temp_wire_u.at(i));
 	    cell_t->set_vwire(temp_wire_v.at(j));
 	    cell_t->set_wwire(temp_wire_w.at(k));
	    cell_t->set_tpc_no(n_tpc);
	    cell = cell_t;

	    //  // new method
	    // WireCell2dToy::TotalTiling& totaltiling = Singleton<WireCell2dToy::TotalTiling>::Instance();
	    // GeomWireSelection temp_wires;
	    // temp_wires.push_back(temp_wire_u.at(i));
	    // temp_wires.push_back(temp_wire_v.at(j));
	    // temp_wires.push_back(temp_wire_w.at(k));
	    // cell = totaltiling.cell(temp_wires);
	    // if (cell == 0){
	    //   GeomCell *cell_t = new GeomCell(ncell,pcell);
	    //   cell_t->set_uwire(temp_wire_u.at(i));
	    //   cell_t->set_vwire(temp_wire_v.at(j));
	    //   cell_t->set_wwire(temp_wire_w.at(k));
	    //   cell_t->set_tpc_no(n_tpc);
	    //   totaltiling.AddCellWire(cell_t, temp_wire_u.at(i), temp_wire_v.at(j), temp_wire_w.at(k));
	    //   cell = cell_t;
	    // }



 	    Point cell_center = cell->center();
 	    if (gds.contained_yz(cell_center)){
	      
 	      // fill cellmap
 	      GeomWireSelection wiresel;
 	      wiresel.push_back(temp_wire_u[i]);
 	      wiresel.push_back(temp_wire_v[j]);
 	      wiresel.push_back(temp_wire_w[k]);	
 	      cellmap[cell]=wiresel;
	      
 	      //fill wiremap
 	      if (wiremap.find(temp_wire_u[i]) == wiremap.end()){
 		//not found
		GeomCellSelection cellsel;
		cellsel.push_back(cell);
		wiremap[temp_wire_u[i]]=cellsel;
	      }else{
		//found
		wiremap[temp_wire_u[i]].push_back(cell);
	      }
	      
	      if (wiremap.find(temp_wire_v[j]) == wiremap.end()){
		//not found
		GeomCellSelection cellsel;
		cellsel.push_back(cell);
		wiremap[temp_wire_v[j]]=cellsel;
	      }else{
		//found
		wiremap[temp_wire_v[j]].push_back(cell);
	      }
	      
	      if (wiremap.find(temp_wire_w[k]) == wiremap.end()){
		//not found
		GeomCellSelection cellsel;
		cellsel.push_back(cell);
		wiremap[temp_wire_w[k]]=cellsel;
	      }else{
		//found
		wiremap[temp_wire_w[k]].push_back(cell);
	      }
	      
 	      cell_all.push_back(cell);
 	      ncell++;
	    }else{
	      delete cell;
	    }
	  }
	}
	
      } // W-loop
    }
  }
}



ClassImp(WireCell2dToy::ToyTiling);
