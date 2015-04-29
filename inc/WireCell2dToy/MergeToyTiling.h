#ifndef WIRECELL_MERGETOYTILING_H
#define WIRECELL_MERGETOYTILING_H

#include "WireCell2dToy/ToyTiling.h"

namespace WireCell2dToy{
  class MergeToyTiling : public WireCell2dToy::ToyTiling {
  public:
    MergeToyTiling(WireCell2dToy::ToyTiling tiling);
    const WireCell::GeomCell* cell(const WireCell::GeomWireSelection& wires)const ;
    
    int further_merge(WireCell::GeomCellSelection &allcell, int ncell);
  };
}
#endif
