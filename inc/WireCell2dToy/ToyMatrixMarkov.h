#ifndef WIRECELL_TOYMATRIXMARKOV_H
#define WIRECELL_TOYMATRIXMARKOV_H

#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/ToyMatrixKalman.h"
#include "WireCellData/MergeGeomCell.h"

typedef std::pair<int, double> CellRankPair;
struct CellRankCompare {
  bool operator() (const CellRankPair& a, const CellRankPair& b) const {
    if (a.second < b.second) {
      return false;
    }else{
      return true;
    }
  }
};
typedef std::set<CellRankPair, CellRankCompare> CellRankSet;

namespace WireCell2dToy{
  class ToyMatrixMarkov {
  public:
    ToyMatrixMarkov(WireCell2dToy::ToyMatrix *toymatrix1,WireCell::GeomCellSelection *allmcell1);
    virtual ~ToyMatrixMarkov();

    std::vector<int>& Get_cur_cell_status(){ return cur_cell_status;};
    std::vector<int>& Get_cur_cell_pol(){ return cur_cell_pol;};
    std::vector<double>& Get_cell_res(){ return cell_res;};
    std::vector<double>& Get_cell_charge(){ return cell_charge;};
    std::vector<double>& Get_cell_prob(){ return cell_prob;};
    

  protected:
    void make_guess();
    void Iterate(WireCell2dToy::ToyMatrixKalman &toykalman);
    
    WireCell2dToy::ToyMatrix *toymatrix;
    WireCell2dToy::ToyMatrixKalman *toymatrixkalman;
    WireCell::GeomCellSelection *allmcell;
    int ncount; // count the number of iterations

    int first_flag;

    //store the MC information
    int mcindex;
    
    std::vector<int>    cur_cell_status; // which one is off from beginning
    std::vector<int>    cur_cell_index;
    std::vector<int>    cur_cell_pol; // on or off
    std::vector<double> cell_res; // residual
    std::vector<double> cell_charge;
    std::vector<double> cell_prob; // probability to go off
    CellRankSet cell_set;

    std::vector<int>    next_cell_status; // to be predicted which one is off from beginning
   
    double_t cur_chi2; //current chi2
    double_t cur_dof;
    
    double_t next_chi2; // chi2 for the next cell
    double_t next_dof;
  };
}

#endif
