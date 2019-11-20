#ifndef WIRECELL_BlobTOYTILING_H
#define WIRECELL_BlobTOYTILING_H


#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/ToyMatrix.h"
#include "WCPData/GeomCell.h"
#include "WCPData/GeomWire.h"
#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"

#include "Rtypes.h"

namespace WCP2dToy{
  class BlobToyTiling : public WCP2dToy::MergeToyTiling {
  public:

    BlobToyTiling(WCP2dToy::ToyTiling& toytiling, WCP2dToy::MergeToyTiling& mergetiling, WCP2dToy::ToyMatrix& toymatrix, int time_slice, int num_merge_wire);
    ~BlobToyTiling();
    
    const WCP::GeomCell* cell(const WCP::GeomWireSelection& wires) const ;
    
  private:
    WCP::MergeGeomCell* FormMergeCell(WCP::MergeGeomWire* mwireu,WCP::MergeGeomWire* mwirev,WCP::MergeGeomWire* mwirew,int ident_cell, int time_slice);
    
    ToyTiling *tiling;
    
    int num_wire;
    
    WCP::GeomWireSelection wire_u1;
    WCP::GeomWireSelection wire_v1;
    WCP::GeomWireSelection wire_w1;


    WCP::GeomCellSelection remain_cells;
   

   
    ClassDef(BlobToyTiling,1);
  };
}
#endif
