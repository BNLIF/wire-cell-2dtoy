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

    int get_solve_flag(){return solve_flag;};
  private:

    void Direct_Solve();
    void L1_Solve();
    
    int solve_flag;
    double direct_chi2;
    int direct_ndf;

    double L1_chi2_base, L1_chi2_penalty;
    int L1_ndf;
    
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
