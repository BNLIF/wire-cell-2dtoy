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

    double get_mcell_charge(WireCell::MergeGeomCell *mcell){return ccmap[mcell];};
    void Update_ndf_chi2();
    int get_ndf();
    double get_chi2();

    int get_L1_ndf(int i){ return L1_ndf.at(i);};
    int get_direct_ndf(int i){ return direct_ndf.at(i);};
    double get_L1_chi2_base(int i){ return L1_chi2_base.at(i);};
    double get_L1_chi2_penalty(int i){return L1_chi2_penalty.at(i);};
    double get_direct_chi2(int i){return direct_chi2.at(i);};
    
  protected:
    int ndf;
    double chi2;

    int nL1_solved;
    std::vector<int> L1_ndf;
    std::vector<double> L1_chi2_base;
    std::vector<double> L1_chi2_penalty;
    
    int ndirect_solved;
    std::vector<int> direct_ndf;
    std::vector<double> direct_chi2;
    
    void divide_groups();
    void init_cell_weight_map();
    
    std::vector<WireCell::GeomCellSelection> final_cells_vec;
    std::vector<WireCell::GeomWireSelection> final_wires_vec;
    std::vector<WireCell2dToy::MatrixSolver*> group_matrices;

    WireCell::CellChargeMap ccmap;
    
    const WireCell::GeomDataSource& gds;
    LowmemTiling& tiling;

    std::map<const WireCell::GeomCell*, double> cell_weight_map;
    
    
    
  };

}

#endif
