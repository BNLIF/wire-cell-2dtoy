#include "WireCell2dToy/ToyMatrixMarkov.h"

using namespace WireCell;
#include "TMath.h"

WireCell2dToy::ToyMatrixMarkov::ToyMatrixMarkov(WireCell2dToy::ToyMatrix *toymatrix1){
  ncount = 0;
  toymatrix = toymatrix1; //save the matrix in here
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(*toymatrix);  // hold the current results 
  make_first_guess();
}

WireCell2dToy::ToyMatrixMarkov::~ToyMatrixMarkov(){
}

void WireCell2dToy::ToyMatrixMarkov::make_first_guess(){
  Iterate(*toymatrixkalman);
}

void WireCell2dToy::ToyMatrixMarkov::Iterate(WireCell2dToy::ToyMatrixKalman &toykalman){
  if (toykalman.Get_numz()!=0&&toymatrix->Get_Chi2()<0){
    for (int i=0;i!=toykalman.Get_mcindex();i++){
      auto it1 = find(toykalman.Get_already_removed().begin(),toykalman.Get_already_removed().end(),i);
      auto it2 = find(toykalman.Get_no_need_remove().begin(),toykalman.Get_no_need_remove().end(),i);
      if (it1 == toykalman.Get_already_removed().end() && it2 == toykalman.Get_no_need_remove().end()){
	std::vector<int> already_removed = toykalman.Get_already_removed();
	already_removed.push_back(i);
	WireCell2dToy::ToyMatrixKalman kalman(already_removed,toykalman.Get_no_need_remove(),*toymatrix,0);
	Iterate(kalman);
	toykalman.Get_no_need_remove().push_back(i);
	if (toymatrix->Get_Chi2()>0) break;
      }
    }
  }
}
