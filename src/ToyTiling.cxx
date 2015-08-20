#include "WireCell2dToy/ToyTiling.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include <cmath>

using namespace WireCell;

WireCell2dToy::ToyTiling::ToyTiling()
{
}


WireCell2dToy::ToyTiling::ToyTiling(const WireCell::Slice& slice,WireCell::GeomDataSource& gds, float rel_u , float rel_v, float rel_w, float noise_u, float noise_v, float noise_w){
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
    wirecharge_errmap[wire] = sqrt(charge_noise*charge_noise + pow(rel_charge_err,2) * charge*charge);

    // hack for now
    if (charge < 100){
      //std::cout << "Bad Channel, Fill in a large error " << std::endl;
      wirecharge_errmap[wire] = 1e4;
    }
    // hack for data

    
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

    dis_w[0] = wdis.at(0) - w_pitch/2.;
    dis_w[1] = dis_w[0] + w_pitch;
    dis_w[2] = wdis.at(0);
    

  
    gds.crossing_point(dis_u[0],dis_v[0],kUwire,kVwire, puv_save[0]);
    gds.crossing_point(dis_u[0],dis_v[1],kUwire,kVwire, puv_save[1]);
    gds.crossing_point(dis_u[1],dis_v[1],kUwire,kVwire, puv_save[2]);
    gds.crossing_point(dis_u[1],dis_v[0],kUwire,kVwire, puv_save[3]);
    gds.crossing_point(dis_u[2],dis_v[2],kUwire,kVwire, puv_save[4]);

    gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, puw_save[0]);
    gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, puw_save[1]);
    gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, puw_save[2]);
    gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, puw_save[3]);
    gds.crossing_point(dis_u[2],dis_w[2],kUwire,kYwire, puw_save[4]);

    gds.crossing_point(dis_v[0],dis_w[0],kVwire,kYwire, pwv_save[0]);
    gds.crossing_point(dis_v[0],dis_w[1],kVwire,kYwire, pwv_save[1]);
    gds.crossing_point(dis_v[1],dis_w[1],kVwire,kYwire, pwv_save[2]);
    gds.crossing_point(dis_v[1],dis_w[0],kVwire,kYwire, pwv_save[3]);
    gds.crossing_point(dis_v[2],dis_w[2],kVwire,kYwire, pwv_save[4]);
    
    for (int k=0;k!=5;k++){
      dis_puv_save[k] = gds.wire_dist(puv_save[k],kYwire);
      dis_puw_save[k] = gds.wire_dist(puw_save[k],kVwire);
      dis_pwv_save[k] = gds.wire_dist(pwv_save[k],kUwire);
    }
  }

  //  int counter = 0, counter1 = 0, counter2 = 0, counter3=0;

  


  for (int i=0;i!=wire_u.size();i++){
    dis_u[0] = udis.at(i) - u_pitch/2.;
    dis_u[1] = dis_u[0] + u_pitch;
    dis_u[2] = udis.at(i);

    for (int j=0;j!=wire_v.size();j++){
      if (wirechargemap[wire_u[i]] <100 && wirechargemap[wire_v[j]] <100 ) continue;
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
	
	if (wirechargemap[wire_u[i]] <100 && wirechargemap[wire_w[k]] <100 ) continue;
	if (wirechargemap[wire_w[k]] <100 && wirechargemap[wire_v[j]] <100 ) continue;
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
	    
	    
	    GeomCell *cell = new GeomCell(ncell,pcell);
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

ClassImp(WireCell2dToy::ToyTiling);
