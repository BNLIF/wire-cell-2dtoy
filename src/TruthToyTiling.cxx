#include "WireCell2dToy/TruthToyTiling.h"

using namespace WireCell;

WireCell2dToy::TruthToyTiling::TruthToyTiling(WireCell2dToy::ToyTiling& tiling, const WireCell::PointValueVector &pvv, int tbin, const GeomDataSource& gds){
  
  float sum = 0;
  for (int itruth = 0; itruth < pvv.size(); ++itruth){
    if (pvv[itruth].first.x == tbin){
      const Point& p = pvv[itruth].first; // get the point
      float charge = pvv[itruth].second;
      
      

      if (gds.contained_yz(p)){

	Point p1 = p;// hack for now, gap in the tiling ... 
	//gds.avoid_gap(p1);
	
	sum += charge;
	const GeomWire* wire_u = gds.closest(p1, static_cast<WirePlaneType_t>(0));
	const GeomWire* wire_v = gds.closest(p1, static_cast<WirePlaneType_t>(1));
	const GeomWire* wire_w = gds.closest(p1, static_cast<WirePlaneType_t>(2));
	
	GeomWireSelection wires;
	wires.push_back(wire_u);
	wires.push_back(wire_v);
	wires.push_back(wire_w);
	
	const GeomCell* cell = tiling.cell(wires);
	
	if (cell!=0){
	  
	  //Point pp = cell->center();
	  //	std::cout << pp.x << " " << pp.y << " " << pp.z << std::endl;
	  

	  if (cellchargemap.find(cell) == cellchargemap.end()){
	    //not found
	    cellchargemap[cell] = charge;
	  }else{
	    cellchargemap[cell] += charge;
	  }
	}else{
	  // Point pp = cell->center();
	  // gds.avoid_gap(p1);
	  // std::cout << p.x << " " << p.y << " " << p.z << " " << std::endl;
	  //std::cout << itruth << " " << p1.x << " " << p1.y << " " << p1.z << " " << charge << std::endl;
	  // std::cout << charge << std::endl;
	}
      }
    }
  }
  
  //  std::cout << "Xin " << sum << std::endl;

}