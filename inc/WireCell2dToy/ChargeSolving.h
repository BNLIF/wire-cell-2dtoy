#ifndef WIRECELL_CHARGESOLVING_H
#define WIRECELL_CHARGESOLVING_H

#include "TMatrixD.h"
#include "TVectorD.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"
#include "WireCell2dToy/LowmemTiling.h"
#include "WireCell2dToy/MatrixSolver.h"

namespace WireCell2dToy{
  class ChargeSolving {
  public:
    ChargeSolving(const WireCell::GeomDataSource& gds, LowmemTiling& tiling);
    virtual ~ChargeSolving();

    int get_ndirect_solved(){return ndirect_solved;};
    int get_nL1_solved(){return nL1_solved;};
    
    
  protected:
    int ndirect_solved;
    int nL1_solved;
    
    void divide_groups();

    std::vector<WireCell::GeomCellSelection> final_cells_vec;
    std::vector<WireCell::GeomWireSelection> final_wires_vec;
    std::vector<WireCell2dToy::MatrixSolver*> group_matrices;
    
    const WireCell::GeomDataSource& gds;
    LowmemTiling& tiling;

   
    
    
    
  };

}

#endif
