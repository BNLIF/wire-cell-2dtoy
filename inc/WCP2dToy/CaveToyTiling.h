#ifndef WIRECELL_CAVETOYTILING_H
#define WIRECELL_CAVETOYTILING_H


#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/ToyMatrix.h"
#include "WCPData/GeomCell.h"
#include "WCPData/GeomWire.h"
#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"

namespace WCP2dToy{
  class CaveToyTiling : public WCP2dToy::MergeToyTiling {
  public:

    CaveToyTiling(WCP2dToy::ToyTiling* toytiling1, WCP2dToy::MergeToyTiling& mergetiling, WCP2dToy::ToyMatrix& toymatrix);
    ~CaveToyTiling();
    
    const WCP::GeomCell* cell(const WCP::GeomWireSelection& wires) const ;
    
  private:
   
    WCP2dToy::ToyTiling *toytiling;

    WCP::GeomCellSelection cell_all_save; // save current state of merged cell before removing a cell
    WCP::GeomWireSelection wire_all_save;

    WCP::GeomCellMap cellmap_save;
    WCP::GeomWireMap wiremap_save;
    
    WCP::GeomWireWireMap wwmap_save;
    WCP::GeomWireWiresMap wwsmap_save;
    
  

    
    WCP::GeomCell* current_cell; // current ceull under consideration coming from cell_to_remove, then will be decided to put in cell_removed or cell_not_to_remove
    WCP::GeomCellSelection cell_removed; 
    WCP::GeomCellSelection cell_to_remove;
    WCP::GeomCellSelection cell_not_to_remove;


    
   

   ClassDef(CaveToyTiling,1);
   
    
  };
}
#endif
