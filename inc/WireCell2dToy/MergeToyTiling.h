#ifndef WIRECELL_MERGETOYTILING_H
#define WIRECELL_MERGETOYTILING_H

#include "WireCell2dToy/ToyTiling.h"

namespace WireCell2dToy{
  class MergeToyTiling : public WireCell2dToy::ToyTiling {
  public:
    MergeToyTiling(){};
    MergeToyTiling(WireCell2dToy::ToyTiling& tiling, int time_slice=-1);
    ~MergeToyTiling();

    const WireCell::GeomCell* cell(const WireCell::GeomWireSelection& wires)const ;
    
    int further_merge(WireCell::GeomCellSelection &allcell, int ncell, int time_slice);
    int further_mergewire(WireCell::GeomWireSelection &allwire, int nwire, int time_slice);

   
  protected:
    

    WireCell::GeomCellMap cellmap1;
    WireCell::GeomWireMap wiremap1;

    WireCell::GeomWireWireMap wwmap; // wire to merged wire
    WireCell::GeomWireWiresMap wwsmap; // merged wire to wires
    
  };
}
#endif
