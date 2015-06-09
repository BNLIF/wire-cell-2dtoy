#ifndef WIRECELL_BlobTOYTILING_H
#define WIRECELL_BlobTOYTILING_H


#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/ToyMatrix.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"

namespace WireCell2dToy{
  class BlobToyTiling : public WireCell2dToy::MergeToyTiling {
  public:

    BlobToyTiling(WireCell2dToy::ToyTiling& toytiling, WireCell2dToy::MergeToyTiling& mergetiling, WireCell2dToy::ToyMatrix& toymatrix, int time_slice, int num_merge_wire);
    ~BlobToyTiling();
    
    const WireCell::GeomCell* cell(const WireCell::GeomWireSelection& wires) const ;
    
  private:
    WireCell::MergeGeomCell* FormMergeCell(WireCell::MergeGeomWire* mwireu,WireCell::MergeGeomWire* mwirev,WireCell::MergeGeomWire* mwirew,int ident_cell, int time_slice);
    
    ToyTiling *tiling;
    
    int num_wire;
    
    WireCell::GeomWireSelection wire_u1;
    WireCell::GeomWireSelection wire_v1;
    WireCell::GeomWireSelection wire_w1;


    WireCell::GeomCellSelection remain_cells;
   

   
    
  };
}
#endif
