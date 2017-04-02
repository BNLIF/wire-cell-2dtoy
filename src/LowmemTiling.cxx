#include "WireCell2dToy/LowmemTiling.h"

using namespace WireCell;

WireCell2dToy::LowmemTiling::LowmemTiling(const WireCell::Slice& slice,WireCell::GeomDataSource& gds,
					  WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map,
					  std::vector<float>& uplane_rms, std::vector<float>& vplane_rms, std::vector<float>& wplane_rms)
  : gds(gds)
{
  // form bad wires group
  std::cout << uplane_map.size() << " " << vplane_map.size() << " " << wplane_map.size() << std::endl;

  // form good wires group
  
}


WireCell2dToy::LowmemTiling::~LowmemTiling(){

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
