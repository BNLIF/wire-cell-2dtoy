#ifndef WIRECELL_MERGETOYTILING_H
#define WIRECELL_MERGETOYTILING_H

#include "WireCell2dToy/ToyTiling.h"

namespace WireCell2dToy{
  class MergeToyTiling : public WireCell2dToy::ToyTiling {
  public:
    MergeToyTiling(){};
    MergeToyTiling(WireCell2dToy::ToyTiling& tiling, int time_slice=-1, int merge_strategy = 3, int flag_remerge=0);
    ~MergeToyTiling();

    const WireCell::GeomCell* cell(const WireCell::GeomWireSelection& wires)const ;
    
    int further_merge(WireCell::GeomCellSelection &allcell, int ncell, int time_slice, double dis = 0.2 * units::mm);
    int further_mergewire(WireCell::GeomWireSelection &allwire, int nwire, int time_slice);
    bool GetRemerged(){return IsRemerged;};
    

  protected:
    bool IsRemerged;

    void form_wiremap(WireCell2dToy::ToyTiling& tiling, int time_slice);

    WireCell::GeomCellMap cellmap1;
    WireCell::GeomWireMap wiremap1;

    WireCell::GeomWireWireMap wwmap; // wire to merged wire
    WireCell::GeomWireWiresMap wwsmap; // merged wire to wires
    

   ClassDef(MergeToyTiling,1);
  };
}
#endif
