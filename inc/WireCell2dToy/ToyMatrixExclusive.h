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
    int ncount;
    int Solve(std::vector<int>& flag, WireCell2dToy::ToyMatrix &toymatrix);
    void move(std::vector<int>& flag, int index, int max);

    int numz; // number of zeros in the matrix's eigenvalues
   
    int mcindex;
    int mwindex;
    int swindex;

    TMatrixDEigen *Eigen;
    TVectorD *EigenValue;
    TMatrixD *trans; 
    /* TMatrixD *transT; */

    TMatrixD *MA, *MAT;
    TMatrixD *MC, *MC_inv;
    TMatrixD *Vx, *Vx_inv;
    //temporary solution
    TVectorD *Cxt, *dCxt;
    double chi2t;
    
    //final solution
    //TVectorD *Cx, *dCx;
    double chi2;
    
  };
}
#endif
