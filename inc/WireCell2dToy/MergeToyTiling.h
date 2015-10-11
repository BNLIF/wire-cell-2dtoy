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
    WireCell::GeomCellSelection& get_single_wire_cells(){return single_wire_cells;};
    WireCell::GeomCellSelection& get_three_wires_cells(){return three_wires_cells;};
    WireCell::GeomCellSelection& get_two_wires_cells(){return two_wires_cells;};

    int further_merge(WireCell::GeomCellSelection &allcell, int ncell, int time_slice, double dis = 0.2 * units::mm);
    int further_mergewire(WireCell::GeomWireSelection &allwire, int nwire, int time_slice);
    bool GetRemerged(){return IsRemerged;};
    
    static double Time2Dis(int time){return time*time_convert + dis_offset;};
    static int Dis2Time(double dis){ return round((dis/units::cm-dis_offset)/time_convert);}; //in cm

    void deghost();


  protected:
    bool IsRemerged;
    
    static double time_convert;
    static double dis_offset;


    void form_wiremap(WireCell2dToy::ToyTiling& tiling, int time_slice);

    WireCell::GeomCellMap cellmap1;
    WireCell::GeomWireMap wiremap1;

    WireCell::GeomWireWireMap wwmap; // wire to merged wire
    WireCell::GeomWireWiresMap wwsmap; // merged wire to wires
    
    WireCell::GeomCellSelection single_wire_cells;
    WireCell::GeomCellSelection three_wires_cells;
    WireCell::GeomCellSelection two_wires_cells;

   ClassDef(MergeToyTiling,1);
  };
}
#endif
