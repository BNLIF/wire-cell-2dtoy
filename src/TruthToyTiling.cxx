#include "WireCell2dToy/TruthToyTiling.h"

using namespace WireCell;

WireCell2dToy::TruthToyTiling::TruthToyTiling(WireCell2dToy::ToyTiling tiling, const WireCell::PointValueVector &pvv, int tbin, const GeomDataSource& gds){
  for (int itruth = 0; itruth < pvv.size(); ++itruth){
    if (pvv[itruth].first.x == tbin){
      const Point& p = pvv[itruth].first; // get the point
      float charge = pvv[itruth].second;
      
      const GeomWire* wire_u = gds.closest(p, static_cast<WirePlaneType_t>(0));
      const GeomWire* wire_v = gds.closest(p, static_cast<WirePlaneType_t>(1));
      const GeomWire* wire_w = gds.closest(p, static_cast<WirePlaneType_t>(2));

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
      }
    }
  }
}
