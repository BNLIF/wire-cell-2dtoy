#include "WireCell2dToy/ToyMatrixExclusive.h"
using namespace WireCell;

#include "TMath.h"

typedef std::pair<int, double> EigenPair;
struct EigenCompare {
  bool operator() (const EigenPair& a, const EigenPair& b) const {
    if (a.second < b.second) {
      return false;
    }else{
      return true;
    }
  }
};
typedef std::set<EigenPair, EigenCompare> EigenSet;


WireCell2dToy::ToyMatrixExclusive::ToyMatrixExclusive(WireCell2dToy::ToyMatrix &toymatrix){
  Eigen = new TMatrixDEigen(*(toymatrix.Get_MC()));
  EigenValue = new TVectorD(Eigen->GetEigenValuesRe());
  trans = new TMatrixD(Eigen->GetEigenVectors());
  transT= new TMatrixD(*trans);
  transT->Transpose(*trans);

  num_size = toymatrix.Get_mcindex();
  numz = 0;
  TVectorD temp(num_size);
  for (int i=0;i!=num_size;i++){
    if (fabs((*EigenValue)[i])<1e-5){
      // for (int j=0;j!=num_size;j++){
      // 	temp[j] += (*trans)(j,i);
      // }
      numz ++;
    }
  }

  


  // EigenSet eset;
  // for (int i=0;i!=num_size;i++){
  //   EigenPair abc(i,fabs(temp[i]));
  //   eset.insert(abc);
  //   //std::cout << i << " " << temp[i] << std::endl;
  // }
  // for (auto it = eset.begin(); it != eset.end(); it++){
  //   std::cout << it->first << " " << it->second << std::endl;
  // }
  //  trans->Print();
  //std::cout << (*trans)(44,45) << std::endl; // second number is column, first number is row  
  // std::cout << num_size << " " << numz << " " << TMath::Factorial(num_size)/TMath::Factorial(numz)/TMath::Factorial(num_size-numz) << std::endl;
  

}

WireCell2dToy::ToyMatrixExclusive::~ToyMatrixExclusive(){
  delete Eigen;
  delete EigenValue;
  delete trans;
  delete transT;
}



