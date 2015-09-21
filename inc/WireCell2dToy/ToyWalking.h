#ifndef WIRECELL_TOYWALKING_H
#define WIRECELL_TOYWALKING_H

#include "WireCellData/MergeSpaceCell.h"

namespace WireCell2dToy{
  class ToyWalking{
  public:
    ToyWalking(WireCell::MergeSpaceCell *start_cell, WireCell::Point start_point, WireCell::MergeSpaceCell *target_cell, WireCell::Point target_point, WireCell::MergeSpaceCellMap& mcells_map, int counter_limit=1000);
    ToyWalking(WireCell::MergeSpaceCell *start_cell, WireCell::Point start_point, WireCell::MergeSpaceCell *target_cell, WireCell::Point target_point, WireCell::MergeSpaceCellMap& mcells_map, WireCell::MergeSpaceCellSelection must_cells, int counter_limit=1000);
    ~ToyWalking();
    
    void Iterate(WireCell::MergeSpaceCell *curr_cell, WireCell::MergeSpaceCellSelection &curr_cells, double dis);
    
    WireCell::MergeSpaceCellSelection get_cells(){return cells;};
    double get_dist(){return dist;};
    int get_counter(){return counter;};
    int get_global_counter(){return global_counter;};

  protected:
    double dist;
    WireCell::MergeSpaceCellSelection cells;


    WireCell::MergeSpaceCell *start_cell;
    WireCell::MergeSpaceCell *target_cell;    
    WireCell::MergeSpaceCellMap& mcells_map;
    WireCell::Point start_point;
    WireCell::Point target_point;

    WireCell::MergeSpaceCellSelection must_cells;

    int global_counter;
    
    int counter;
    int counter_limit;
  };
}

#endif
