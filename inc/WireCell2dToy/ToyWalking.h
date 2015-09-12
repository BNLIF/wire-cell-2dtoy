#ifndef WIRECELL_TOYWALKING_H
#define WIRECELL_TOYWALKING_H

#include "WireCellData/MergeSpaceCell.h"

namespace WireCell2dToy{
  class ToyWalking{
  public:
    ToyWalking(WireCell::MergeSpaceCell *start_cell, WireCell::Point start_point, WireCell::MergeSpaceCell *target_cell, WireCell::Point target_point, WireCell::MergeSpaceCellMap& mcells_map);
    ~ToyWalking();
    
    void Iterate(WireCell::MergeSpaceCell *curr_cell, WireCell::MergeSpaceCellSelection &curr_cells, double dis);
    
    WireCell::MergeSpaceCellSelection get_cells(){return cells;};
    double get_dist(){return dist;};
    
  protected:
    double dist;
    WireCell::MergeSpaceCellSelection cells;


    WireCell::MergeSpaceCell *start_cell;
    WireCell::MergeSpaceCell *target_cell;    
    WireCell::MergeSpaceCellMap& mcells_map;
    WireCell::Point start_point;
    WireCell::Point target_point;
    
    
  };
}

#endif
