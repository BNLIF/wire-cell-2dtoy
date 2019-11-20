#ifndef WIRECELL_TOYMATRIXKALMAN_H
#define WIRECELL_TOYMATRIXKALMAN_H

#include "WCP2dToy/ToyMatrix.h"
#include "TMatrixDEigen.h"

namespace WCP2dToy {
  class ToyMatrixKalman {
  public:
    ToyMatrixKalman(WCP2dToy::ToyMatrix &toymatrix, int flag_no_need=1, double chi2_penalty = 0);
    ToyMatrixKalman(std::vector<int>& already_removed1, std::vector<int>& no_need_remove1 ,WCP2dToy::ToyMatrix &toymatrix, int check, int flag_no_need=1, double chi2_penalty = 0);
    virtual ~ToyMatrixKalman();

    int Get_numz(){return numz;};
    int Get_ncount(){return ncount;};
    std::vector<int>& Get_already_removed(){return already_removed;};
    std::vector<int>& Get_no_need_remove(){return no_need_remove;};
    int Get_mcindex(){return mcindex;};
    int Cal_numz(WCP2dToy::ToyMatrix &toymatrix);
    void init(WCP2dToy::ToyMatrix &toymatrix);
    void Set_penalty(double val){chi2_penalty = val;}

  protected:
    int flag_no_need;
    double chi2_penalty;

    std::vector<int> already_removed;
    std::vector<int> no_need_remove;
    
    int numz; //number of zeros in the matrix
    int ncount; // how many solve was done
    
    int mcindex;
    int mwindex;
    int swindex;

    TMatrixDEigen *Eigen;
    TVectorD *EigenValue;

    int check_flag;

    /* TMatrixD *MA, *MAT; */
    /* TMatrixD *MC, *MC_inv; */
  };
}
#endif
