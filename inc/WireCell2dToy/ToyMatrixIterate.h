#ifndef WIRECELL_TOYMATRIXITERATE_H
#define WIRECELL_TOYMATRIXITERATE_H

#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/ToyMatrixKalman.h"

namespace WireCell2dToy {
  class ToyMatrixIterate {
  public:
    ToyMatrixIterate(WireCell2dToy::ToyMatrix &toymatrix, int recon_t = 2000, float limit = 1e6);
    ToyMatrixIterate(WireCell2dToy::ToyMatrix &toymatrix, WireCell2dToy::MergeToyTiling &mergecur, WireCell::GeomCellSelection &cells, int recon_t = 2000);
    ToyMatrixIterate(WireCell2dToy::ToyMatrix *toybefore, WireCell2dToy::ToyMatrix *toycur, WireCell2dToy::ToyMatrix *toyafter, WireCell2dToy::MergeToyTiling *mergebefore, WireCell2dToy::MergeToyTiling *mergecur, WireCell2dToy::MergeToyTiling *mergeafter, 
		     int recon_t = 2000, float limit = 1e6, double penalty = 5, double penalty_ncpt = 0);    

    ToyMatrixIterate(WireCell2dToy::ToyMatrix &toymatrix, std::vector<int>& already_removed, int recon_t = 2000, int limit = 1e6); // not very useful any way ... 

    void UseTime(WireCell2dToy::ToyMatrix &toybefore, WireCell2dToy::ToyMatrix &toycur, WireCell2dToy::ToyMatrix &toyafter, WireCell2dToy::MergeToyTiling &mergebefore, WireCell2dToy::MergeToyTiling &mergecur, WireCell2dToy::MergeToyTiling &mergeafter); //previous use time, not popular now ... 


    virtual ~ToyMatrixIterate();
    int Get_ncount(){return ncount;};
    int Get_timeflag(){return time_flag;};
    
  protected:
    void Iterate(WireCell2dToy::ToyMatrixKalman &toymatrixkalman,WireCell2dToy::ToyMatrix &toymatrix);
    void Iterate_simple(WireCell2dToy::ToyMatrixKalman &toymatrixkalman,WireCell2dToy::ToyMatrix &toymatrix);
    void Iterate_simple1(WireCell2dToy::ToyMatrixKalman &toymatrixkalman,WireCell2dToy::ToyMatrix &toymatrix);
    void find_subset(WireCell2dToy::ToyMatrixKalman &toymatrixkalman,WireCell2dToy::ToyMatrix &toymatrix,  std::vector<int>& vec);
    int ncount;
    int prev_ncount;
    int nlevel;
    double estimated_loop;

    double penalty_ncpt;

    int time_flag;
    int recon_threshold;

    std::map<int,double> cell_penal;
    std::map<int,std::vector<int>> cells_ncpt;
    
    WireCell2dToy::ToyMatrixKalman *toymatrixkalman;
 };
}
#endif
