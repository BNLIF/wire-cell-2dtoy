#ifndef WIRECELL_TOYMATRIXMARKOV_H
#define WIRECELL_TOYMATRIXMARKOV_H

#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/ToyMatrixKalman.h"

namespace WireCell2dToy{
  class ToyMatrixMarkov {
  public:
    ToyMatrixMarkov(WireCell2dToy::ToyMatrix *toymatrix1);
    virtual ~ToyMatrixMarkov();

  protected:
    void make_first_guess();
    void Iterate(WireCell2dToy::ToyMatrixKalman &toykalman);

    WireCell2dToy::ToyMatrix *toymatrix;
    WireCell2dToy::ToyMatrixKalman *toymatrixkalman;
    int ncount; // count the number of iterations
    
  };
}

#endif
