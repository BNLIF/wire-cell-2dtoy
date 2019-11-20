#ifndef WIRECELL_MERGETOYTILING_H
#define WIRECELL_MERGETOYTILING_H

#include "WCP2dToy/ToyTiling.h"

namespace WCP2dToy{
  class MergeToyTiling : public WCP2dToy::ToyTiling {
  public:
    MergeToyTiling(){};
    MergeToyTiling(const WCP::DetectorGDS& gds, WCP2dToy::ToyTiling& tiling, int time_slice = -1);
    MergeToyTiling(WCP2dToy::ToyTiling& tiling, int time_slice=-1, int merge_strategy = 3, int flag_remerge=0);
    ~MergeToyTiling();

    const WCP::GeomCell* cell(const WCP::GeomWireSelection& wires)const ;
    WCP::GeomCellSelection& get_single_wire_cells(){return single_wire_cells;};
    WCP::GeomCellSelection& get_three_wires_cells(){return three_wires_cells;};
    WCP::GeomCellSelection& get_two_wires_cells(){return two_wires_cells;};
    WCP::GeomCellCellsMap& get_not_compatible_cells_map(){ return mcmcsmap;};
    

    int further_merge(WCP::GeomCellSelection &allcell, int ncell, int time_slice, double dis = 0.2 * units::mm);
    int further_mergewire(WCP::GeomWireSelection &allwire, int nwire, int time_slice);
    bool GetRemerged(){return IsRemerged;};
    
    static double Time2Dis(int time){return time*time_convert + dis_offset;};
    static int Dis2Time(double dis){ return round((dis/units::cm-dis_offset)/time_convert);}; //in cm

    void deghost(WCP::GeomCellSelection& good_mcells); //used for MicroBooNE data ... 


  protected:
    bool IsRemerged;
    
    static double time_convert;
    static double dis_offset;


    void form_wiremap(WCP2dToy::ToyTiling& tiling, int time_slice);
    void form_wiremap(const WCP::DetectorGDS& gds, WCP2dToy::ToyTiling& tiling, int time_slice);

    WCP::GeomCellMap cellmap1;
    WCP::GeomWireMap wiremap1;

    WCP::GeomWireWireMap wwmap; // wire to merged wire
    WCP::GeomWireWiresMap wwsmap; // merged wire to wires
    
    WCP::GeomCellCellsMap mcmcsmap; //merged cell to merged cell map;

    WCP::GeomCellSelection single_wire_cells;
    WCP::GeomCellSelection three_wires_cells;
    WCP::GeomCellSelection two_wires_cells;

   ClassDef(MergeToyTiling,1);
  };
}
#endif
