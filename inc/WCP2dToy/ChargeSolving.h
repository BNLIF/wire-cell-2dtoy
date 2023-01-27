#ifndef WIRECELL_CHARGESOLVING_H
#define WIRECELL_CHARGESOLVING_H

#include "TMatrixD.h"
#include "TVectorD.h"
#include "WCPData/GeomCell.h"
#include "WCPData/GeomWire.h"
#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"
#include "WCP2dToy/LowmemTiling.h"
#include "WCP2dToy/MatrixSolver.h"

//#include <set>

namespace WCP2dToy{
  class ChargeSolving {
  public:
    ChargeSolving(const WCP::GeomDataSource& gds, LowmemTiling& tiling);
    virtual ~ChargeSolving();

    void L1_resolve(float weight, float reduce_weight_factor=3);
    
    int get_ndirect_solved(){return ndirect_solved;};
    int get_nL1_solved(){return nL1_solved;};

    double get_mcell_charge(const WCP::GeomCell *mcell){return ccmap[mcell];};
    void Update_ndf_chi2();
    int get_ndf();
    double get_chi2();

    int get_L1_ndf(int i){ return L1_ndf.at(i);};
    int get_direct_ndf(int i){ return direct_ndf.at(i);};
    double get_L1_chi2_base(int i){ return L1_chi2_base.at(i);};
    double get_L1_chi2_penalty(int i){return L1_chi2_penalty.at(i);};
    double get_direct_chi2(int i){return direct_chi2.at(i);};

    void clear_connectivity(){front_connectivity.clear(); back_connectivity.clear();};
    void add_front_connectivity(const WCP::GeomCell *cell){ front_connectivity.insert(cell);};
    void add_back_connectivity(const WCP::GeomCell *cell){back_connectivity.insert(cell); };

    
    // information
    void update_information();
    float get_weight(WCP::GeomCell *cell);
    int get_matrix_id(WCP::GeomCell *cell);
    int get_mc_id(WCP::GeomCell *cell);
    double get_mc_lambda(WCP::GeomCell *cell);
    double get_mc_TOL(WCP::GeomCell *cell);

    int get_mw_id(WCP::GeomWire *mwire);
    double get_W(WCP::GeomWire *mwire);
    double get_G(WCP::GeomCell *cell, WCP::GeomWire *mwire);
    
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
    void init_cell_weight_map(float weight=1);
    void update_cell_weight_map(float weight, float reduce_weight_factor);

    
    
    std::vector<WCP::GeomCellSelection> final_cells_vec;
    std::vector<WCP::GeomWireSelection> final_wires_vec;
    std::vector<WCP2dToy::MatrixSolver*> group_matrices;

    WCP::CellChargeMap ccmap;
    
    const WCP::GeomDataSource& gds;
    LowmemTiling& tiling;

    std::map<const WCP::GeomCell*, double, WCP::GeomCellComparep > cell_weight_map;
    
    std::set<const WCP::GeomCell*, WCP::GeomCellComparep > front_connectivity;
    std::set<const WCP::GeomCell*, WCP::GeomCellComparep > back_connectivity;

    //information ...
    std::map<const WCP::GeomCell*, int> map_mc_matrix_id;
    std::map<const WCP::GeomCell*, int> map_mc_index;
    std::map<const WCP::GeomCell*, double> map_mc_TOL;
    std::map<const WCP::GeomCell*, double> map_mc_lambda;
    
    std::map<const WCP::GeomWire*, int> map_mw_index;
    std::map<const WCP::GeomWire*, double> map_W;
    std::map<std::pair<const WCP::GeomCell*, const WCP::GeomWire*>, double> map_G;
    
  };

}

#endif
