#ifndef WIRECELL_TOYMATRIXEXCLUSIVE_H
#define WIRECELL_TOYMATRIXEXCLUSIVE_H

#include "WireCell2dToy/ToyMatrix.h"
#include "TMatrixDEigen.h"


namespace WireCell2dToy {
  class ToyMatrixExclusive {
  public:
    ToyMatrixExclusive(WireCell2dToy::ToyMatrix &toymatrix);
    virtual ~ToyMatrixExclusive();
    

  protected: 
    int numz; // number of zeros in the matrix's eigenvalues
    int num_size;

    TMatrixDEigen *Eigen;
    TVectorD *EigenValue;
    TMatrixD *trans;
    TMatrixD *transT;
  };
}
#endif
