#include "WireCell2dToy/ChargeSolving.h"


using namespace WireCell;

WireCell2dToy::ChargeSolving::ChargeSolving(const WireCell::GeomDataSource& gds, LowmemTiling& tiling)
  : gds(gds)
  , tiling(tiling)
{
}

WireCell2dToy::ChargeSolving::~ChargeSolving(){
}
