#ifndef WIRECELL_TOYMATRIXMARKOV_H
#define WIRECELL_TOYMATRIXMARKOV_H

#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/ToyMatrixKalman.h"
#include "WCPData/MergeGeomCell.h"

typedef std::pair<int, double> CellRankPair;
struct CellRankCompare {
  bool operator() (const CellRankPair& a, const CellRankPair& b) const {
    if (a.second == b.second){
      return a.first < b.first;
    }
    return a.second > b.second;
  }
};
typedef std::set<CellRankPair, CellRankCompare> CellRankSet;

namespace WCP2dToy{
  class ToyMatrixMarkov {
  public:
    ToyMatrixMarkov(WCP2dToy::ToyMatrix *toymatrix1, WCP::GeomCellSelection *allmcell1, int recon_t1=1500, int recon_t2=2000);
    ToyMatrixMarkov(WCP2dToy::ToyMatrix &toymatrix1, WCP2dToy::MergeToyTiling &mergecur, WCP::GeomCellSelection *allmcell1, WCP::GeomCellSelection &cells, int recon_t1=1500, int recon_t2=2000);
    
    ToyMatrixMarkov(WCP2dToy::ToyMatrix *toybefore, WCP2dToy::ToyMatrix *toycur, WCP2dToy::ToyMatrix *toyafter, WCP2dToy::MergeToyTiling *mergebefore, WCP2dToy::MergeToyTiling *mergecur, WCP2dToy::MergeToyTiling *mergeafter, WCP::GeomCellSelection *allmcell1, int recon_t1=1500, int recon_t2=2000, double penalty = 5, double penalty_ncpt = 0);

    ToyMatrixMarkov(WCP2dToy::ToyMatrix &toybefore, WCP2dToy::ToyMatrix &toycur, WCP2dToy::ToyMatrix &toyafter, WCP2dToy::MergeToyTiling &mergebefore, WCP2dToy::MergeToyTiling &mergecur, WCP2dToy::MergeToyTiling &mergeafter, WCP::GeomCellSelection *allmcell1, int recon_t1=1500, int recon_t2=2000, double penalty = 5, double penalty_ncpt = 0);

    virtual ~ToyMatrixMarkov();

    std::vector<int>& Get_cur_cell_status(){ return cur_cell_status;};
    std::vector<int>& Get_cur_cell_pol(){ return cur_cell_pol;};
    std::vector<double>& Get_cell_res(){ return cell_res;};
    std::vector<double>& Get_cell_charge(){ return cell_charge;};
    std::vector<double>& Get_cell_prob(){ return cell_prob;};
    

  protected:
    void find_subset(WCP2dToy::ToyMatrixKalman &toymatrix,WCP2dToy::ToyMatrix &toymatrix1, std::vector<int>& vec);
    void make_guess();
    void Iterate(WCP2dToy::ToyMatrixKalman &toykalman);
    
    WCP2dToy::ToyMatrix *toymatrix;
    WCP2dToy::ToyMatrixKalman *toymatrixkalman;
    WCP::GeomCellSelection *allmcell;
    int ncount; // count the number of iterations

    double penalty_ncpt;

    int first_flag;
    float  use_time_threshold;

    //store the MC information
    int mcindex;

    std::vector<int> use_time; //to be added

    
    std::vector<int>    cur_cell_status; // which one is off from beginning
    std::vector<int>    cur_cell_index;
    std::vector<int>    cur_cell_pol; // on or off
    std::vector<double> cell_res; // residual
    std::vector<double> cell_charge;
    std::vector<double> cell_prob; // probability to go off
    CellRankSet cell_set;

    std::vector<int>    next_cell_status; // to be predicted which one is off from beginning
   
    double cur_chi2; //current chi2
    double cur_chi2_save;
    double cur_dof;
    
    double next_chi2; // chi2 for the next cell
    double next_dof;

    std::map<int,double> cell_penal;
    std::map<int,std::vector<int>> cells_ncpt;
    
    int recon_threshold1;
    int recon_threshold2;


    int nlevel;
  };
}

#endif
