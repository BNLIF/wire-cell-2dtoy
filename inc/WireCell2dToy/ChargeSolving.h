#ifndef WIRECELL_CHARGESOLVING_H
#define WIRECELL_CHARGESOLVING_H

#include "TMatrixD.h"
#include "TVectorD.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"
#include "LowmemTiling.h"

namespace WireCell2dToy{
  class ChargeSolving {
  public:
    ChargeSolving(const WireCell::GeomDataSource& gds, LowmemTiling& tiling);
    virtual ~ChargeSolving();


    
  protected:
    void divide_groups();

    std::vector<WireCell::GeomCellSelection> final_cells_vec;
    std::vector<WireCell::GeomWireSelection> final_wires_vec;
    
    const WireCell::GeomDataSource& gds;
    LowmemTiling& tiling;

   
    
    
    
  };

}

#endif
