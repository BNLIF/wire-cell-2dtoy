#ifndef WIRECELL_CAVETOYTILING_H
#define WIRECELL_CAVETOYTILING_H


#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/ToyMatrix.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"

namespace WireCell2dToy{
  class CaveToyTiling : public WireCell2dToy::MergeToyTiling {
  public:

    CaveToyTiling(WireCell2dToy::MergeToyTiling& mergetiling, WireCell2dToy::ToyMatrix& toymatrix);
    ~CaveToyTiling();
    
    const WireCell::GeomCell* cell(const WireCell::GeomWireSelection& wires) const ;
    
  private:
   
      
    WireCell::GeomCellSelection cell_all_save; // save current state of merged cell before removing a cell
    WireCell::GeomWireSelection wire_all_save;

    WireCell::GeomCellMap cellmap_save;
    WireCell::GeomWireMap wiremap_save;
    
    WireCell::GeomWireWireMap wwmap_save;
    WireCell::GeomWireWiresMap wwsmap_save;
    

    
    WireCell::GeomCell* current_cell; // current ceull under consideration coming from cell_to_remove, then will be decided to put in cell_removed or cell_not_to_remove
    WireCell::GeomCellSelection cell_removed; 
    WireCell::GeomCellSelection cell_to_remove;
    WireCell::GeomCellSelection cell_not_to_remove;


    
   

   
    
  };
}
#endif
