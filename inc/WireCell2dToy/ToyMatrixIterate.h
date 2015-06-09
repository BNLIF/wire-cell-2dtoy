#ifndef WIRECELL_TOYMATRIXITERATE_H
#define WIRECELL_TOYMATRIXITERATE_H

#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/ToyMatrixKalman.h"


namespace WireCell2dToy {
  class ToyMatrixIterate {
  public:
    ToyMatrixIterate(WireCell2dToy::ToyMatrix &toymatrix);
    ToyMatrixIterate(WireCell2dToy::ToyMatrix &toymatrix, std::vector<int>& already_removed);

    void UseTime(WireCell2dToy::ToyMatrix &toybefore, WireCell2dToy::ToyMatrix &toycur, WireCell2dToy::ToyMatrix &toyafter, WireCell2dToy::MergeToyTiling &mergebefore, WireCell2dToy::MergeToyTiling &mergecur, WireCell2dToy::MergeToyTiling &mergeafter);

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

    int time_flag;
    
    WireCell2dToy::ToyMatrixKalman *toymatrixkalman;
 };
}
#endif
