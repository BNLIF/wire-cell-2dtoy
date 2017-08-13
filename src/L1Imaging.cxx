#include "WireCell2dToy/L1Imaging.h"


using namespace WireCell;

WireCell2dToy::L1Imaging::L1Imaging(const WireCell::GeomDataSource& gds, LowmemTiling& tiling)
  : gds(gds)
  , tiling(tiling)
{
}

WireCell2dToy::L1Imaging::~L1Imaging(){
}
