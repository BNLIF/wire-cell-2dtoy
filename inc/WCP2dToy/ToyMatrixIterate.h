#ifndef WIRECELL_TOYMATRIXITERATE_H
#define WIRECELL_TOYMATRIXITERATE_H

#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/ToyMatrixKalman.h"

namespace WCP2dToy {
  class ToyMatrixIterate {
  public:
    ToyMatrixIterate(WCP2dToy::ToyMatrix &toymatrix, int recon_t = 2000, float limit = 1e6);
    ToyMatrixIterate(WCP2dToy::ToyMatrix &toymatrix, WCP2dToy::MergeToyTiling &mergecur, WCP::GeomCellSelection &cells, int recon_t = 2000);
    ToyMatrixIterate(WCP2dToy::ToyMatrix *toybefore, WCP2dToy::ToyMatrix *toycur, WCP2dToy::ToyMatrix *toyafter, WCP2dToy::MergeToyTiling *mergebefore, WCP2dToy::MergeToyTiling *mergecur, WCP2dToy::MergeToyTiling *mergeafter, 
		     int recon_t = 2000, float limit = 1e6, double penalty = 5, double penalty_ncpt = 0);    

    ToyMatrixIterate(WCP2dToy::ToyMatrix &toymatrix, std::vector<int>& already_removed, int recon_t = 2000, int limit = 1e6); // not very useful any way ... 

    void UseTime(WCP2dToy::ToyMatrix &toybefore, WCP2dToy::ToyMatrix &toycur, WCP2dToy::ToyMatrix &toyafter, WCP2dToy::MergeToyTiling &mergebefore, WCP2dToy::MergeToyTiling &mergecur, WCP2dToy::MergeToyTiling &mergeafter); //previous use time, not popular now ... 


    virtual ~ToyMatrixIterate();
    int Get_ncount(){return ncount;};
    int Get_timeflag(){return time_flag;};
    
  protected:
    void Iterate(WCP2dToy::ToyMatrixKalman &toymatrixkalman,WCP2dToy::ToyMatrix &toymatrix);
    void Iterate_simple(WCP2dToy::ToyMatrixKalman &toymatrixkalman,WCP2dToy::ToyMatrix &toymatrix);
    void Iterate_simple1(WCP2dToy::ToyMatrixKalman &toymatrixkalman,WCP2dToy::ToyMatrix &toymatrix);
    void find_subset(WCP2dToy::ToyMatrixKalman &toymatrixkalman,WCP2dToy::ToyMatrix &toymatrix,  std::vector<int>& vec);
    int ncount;
    int prev_ncount;
    int nlevel;
    double estimated_loop;

    double penalty_ncpt;

    int time_flag;
    int recon_threshold;

    std::map<int,double> cell_penal;
    std::map<int,std::vector<int>> cells_ncpt;
    
    WCP2dToy::ToyMatrixKalman *toymatrixkalman;
 };
}
#endif
