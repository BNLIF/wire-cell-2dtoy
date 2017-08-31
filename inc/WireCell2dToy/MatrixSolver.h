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

    double get_mcell_charge(WireCell::MergeGeomCell *mcell);

    double get_direct_chi2(){return direct_chi2;};
    int get_direct_ndf(){return direct_ndf;};

    double get_L1_chi2_base(){return L1_chi2_base;};
    double get_L1_chi2_penalty(){return L1_chi2_penalty;};
    int get_L1_ndf(){return L1_ndf;};
    
    int get_solve_flag(){return solve_flag;};

    void L1_Solve(std::map<const WireCell::GeomCell*, double>& cell_weight_map);

    WireCell::GeomCellSelection get_all_cells();
    
  private:
    void Direct_Solve();
        
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
