#ifndef WIRECELL_TOYWALKING_H
#define WIRECELL_TOYWALKING_H

#include "WCPData/MergeSpaceCell.h"

namespace WCP2dToy{
  class ToyWalking{
  public:
    ToyWalking(WCP::MergeSpaceCell *start_cell, WCP::Point start_point, WCP::MergeSpaceCell *target_cell, WCP::Point target_point, WCP::MergeSpaceCellMap& mcells_map, int counter_limit=1000);
    ToyWalking(WCP::MergeSpaceCell *start_cell, WCP::Point start_point, WCP::MergeSpaceCell *target_cell, WCP::Point target_point, WCP::MergeSpaceCellMap& mcells_map, WCP::MergeSpaceCellSelection must_cells, int counter_limit=1000);

   

    ~ToyWalking();
    
    void Iterate(WCP::MergeSpaceCell *curr_cell, WCP::MergeSpaceCellSelection &curr_cells, double dis);
    
    WCP::MergeSpaceCellSelection get_cells(){return cells;};
    double get_dist(){return dist;};
    int get_counter(){return counter;};
    int get_global_counter(){return global_counter;};

  protected:
    double dist;
    WCP::MergeSpaceCellSelection cells;


    WCP::MergeSpaceCell *start_cell;
    WCP::MergeSpaceCell *target_cell;    
    WCP::MergeSpaceCellMap& mcells_map;
    WCP::Point start_point;
    WCP::Point target_point;

    WCP::MergeSpaceCellSelection must_cells;

    int global_counter;
    
    int counter;
    int counter_limit;
  };

  class ToyNWalking{
  public:
    ToyNWalking(WCP::MergeSpaceCell *start_cell, WCP::MergeSpaceCellMap& mcells_map, WCP::MergeSpaceCellSelection& used_cell, WCP::MergeSpaceCellSelection& must_cell);
    ~ToyNWalking();
    WCP::MergeSpaceCellSelection get_cells(){return cells;};
    void Iterate(WCP::MergeSpaceCell *curr_cell);
    
  protected:
    WCP::MergeSpaceCellSelection cells;
    
    WCP::MergeSpaceCell *start_cell;
    WCP::MergeSpaceCellMap& mcells_map;
    
    WCP::MergeSpaceCellSelection& must_cell;
    WCP::MergeSpaceCellSelection& used_cell;
  };

}

#endif
