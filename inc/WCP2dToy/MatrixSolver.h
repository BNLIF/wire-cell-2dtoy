#ifndef WIRECELL2DTOY_MATRIXSOLVER_H
#define WIRECELL2DTOY_MATRIXSOLVER_H

#include "TMatrixD.h"
#include "TVectorD.h"


#include "WCPData/GeomCell.h"
#include "WCPData/GeomWire.h"
#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"


#include <Eigen/Dense>

namespace WCP2dToy{
  class MatrixSolver{
  public:
    MatrixSolver(WCP::GeomCellSelection& cells, WCP::GeomWireSelection& wires, WCP::GeomCellMap& cell_wire_map, WCP::GeomWireMap& wire_cell_map, WCP::WireChargeMap& wire_charge_map, WCP::WireChargeMap& wire_charge_error_map);
    virtual ~MatrixSolver();

    double get_mcell_charge(WCP::MergeGeomCell *mcell);

    double get_direct_chi2(){return direct_chi2;};
    int get_direct_ndf(){return direct_ndf;};

    double get_L1_chi2_base(){return L1_chi2_base;};
    double get_L1_chi2_penalty(){return L1_chi2_penalty;};
    int get_L1_ndf(){return L1_ndf;};
    
    int get_solve_flag(){return solve_flag;};

    void L1_Solve(std::map<const WCP::GeomCell*, double, WCP::GeomCellComparep >& cell_weight_map);

    WCP::GeomCellSelection get_all_cells();

    // added to expose more information out ...
    double get_lambda(){return lambda;};
    double get_TOL(){return TOL;};
    int get_mc_index(const WCP::GeomCell* cell);
    int get_mw_index(const WCP::GeomWire* wire);
    double get_W_value(int index);
    double get_G_value(int index_mw, int index_mc);

    void set_id(int id1){id = id1;};
    int get_id(){return id;};
    
  private:
    //    void Direct_Solve();
    int id;
        
    int solve_flag;
    double direct_chi2;
    int direct_ndf;

    double L1_chi2_base, L1_chi2_penalty;
    int L1_ndf;
    
    /* TMatrixD *MA, *MB, *MAT, *MBT; */
    /* TMatrixD *Vy, *VBy, *Vx, *VBy_inv, *Vx_inv; */
    /* TMatrixD *MC, *MC_inv; */

    /* TMatrixD *UMA;  */
    /* TVectorD *UMWy; */ 
    double lambda;
    double TOL;
    float scale_factor;
    
    TVectorD *Wy, *Cx, *dCx;
    TVectorD *MWy_pred, *MWy;

    Eigen::VectorXd W;
    Eigen::MatrixXd G;
    
    WCP::CellIndexMap mcimap;
    WCP::WireIndexMap mwimap, swimap;
    
    
    int mcindex;
    int mwindex;
    int swindex;
    
  };
}

#endif
