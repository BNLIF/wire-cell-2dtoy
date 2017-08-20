#ifndef WIRECELL2DTOY_MATRIXSOLVER_H
#define WIRECELL2DTOY_MATRIXSOLVER_H

#include "TMatrixD.h"
#include "TVectorD.h"

#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"

namespace WireCell2dToy{
  class MatrixSolver{
  public:
    MatrixSolver(WireCell::GeomCellSelection& cells, WireCell::GeomWireSelection& wires, WireCell::GeomCellMap& cell_wire_map, WireCell::GeomWireMap& wire_cell_map, WireCell::WireChargeMap& wire_charge_map, WireCell::WireChargeMap& wire_charge_error_map);
    virtual ~MatrixSolver();
  private:

    TMatrixD *MA, *MB, *MAT, *MBT;
    TMatrixD *Vy, *VBy, *Vx, *VBy_inv, *Vx_inv;

    TMatrixD *MC, *MC_inv;

    TVectorD *Wy, *Cx, *dCx;
    
    TVectorD *MWy_pred, *MWy;
    
    WireCell::CellIndexMap mcimap;
    WireCell::WireIndexMap mwimap, swimap;
    
    
    int mcindex;
    int mwindex;
    int swindex;
    
  };
}

#endif
