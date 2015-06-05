#ifndef WIRECELL_TOYMATRIXKALMAN_H
#define WIRECELL_TOYMATRIXKALMAN_H

#include "WireCell2dToy/ToyMatrix.h"
#include "TMatrixDEigen.h"

namespace WireCell2dToy {
  class ToyMatrixKalman {
  public:
    ToyMatrixKalman(WireCell2dToy::ToyMatrix &toymatrix);
    ToyMatrixKalman(std::vector<int>& already_removed1, std::vector<int>& no_need_remove1 ,WireCell2dToy::ToyMatrix &toymatrix);
    virtual ~ToyMatrixKalman();

    int Get_numz(){return numz;};
    int Get_ncount(){return ncount;};
    std::vector<int>& Get_already_removed(){return already_removed;};
    std::vector<int>& Get_no_need_remove(){return no_need_remove;};
    int Get_mcindex(){return mcindex;};

  protected:
    void init(WireCell2dToy::ToyMatrix &toymatrix);
    
    std::vector<int> already_removed;
    std::vector<int> no_need_remove;
    
    int numz; //number of zeros in the matrix
    int ncount; // how many solve was done
    
    int mcindex;
    int mwindex;
    int swindex;

    TMatrixDEigen *Eigen;
    TVectorD *EigenValue;

    TMatrixD *MA, *MAT;
    TMatrixD *MC, *MC_inv;
  };
}
#endif
