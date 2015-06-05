#ifndef WIRECELL_TOYMATRIXITERATE_H
#define WIRECELL_TOYMATRIXITERATE_H

#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/ToyMatrixKalman.h"


namespace WireCell2dToy {
  class ToyMatrixIterate {
  public:
    ToyMatrixIterate(WireCell2dToy::ToyMatrix &toymatrix);
    virtual ~ToyMatrixIterate();
    int Get_ncount(){return ncount;};
  protected:
    void Iterate(WireCell2dToy::ToyMatrixKalman &toymatrixkalman,WireCell2dToy::ToyMatrix &toymatrix);
    int ncount;
    int nlevel;
    WireCell2dToy::ToyMatrixKalman *toymatrixkalman;
 };
}
#endif
