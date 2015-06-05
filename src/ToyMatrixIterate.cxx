#include "WireCell2dToy/ToyMatrixIterate.h"

using namespace WireCell;

#include "TMath.h"

WireCell2dToy::ToyMatrixIterate::ToyMatrixIterate(WireCell2dToy::ToyMatrix &toymatrix){
  ncount = 0;
  nlevel = 0;
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(toymatrix);
  Iterate(*toymatrixkalman,toymatrix);
}


WireCell2dToy::ToyMatrixIterate::~ToyMatrixIterate(){
  delete toymatrixkalman;
}

void WireCell2dToy::ToyMatrixIterate::Iterate(WireCell2dToy::ToyMatrixKalman &toymatrix,WireCell2dToy::ToyMatrix &toymatrix1){
  nlevel ++;
  if (toymatrix.Get_numz()!=0){
    for (int i=0;i!=toymatrix.Get_mcindex();i++){
      auto it1 = find(toymatrix.Get_already_removed().begin(),toymatrix.Get_already_removed().end(),i);
      auto it2 = find(toymatrix.Get_no_need_remove().begin(),toymatrix.Get_no_need_remove().end(),i);
      if (it1 == toymatrix.Get_already_removed().end() && it2 == toymatrix.Get_no_need_remove().end()){
	std::vector<int> already_removed = toymatrix.Get_already_removed();
	already_removed.push_back(i);
	WireCell2dToy::ToyMatrixKalman *kalman = new ToyMatrixKalman(already_removed,toymatrix.Get_no_need_remove(),toymatrix1);
	//	std::cout << nlevel << std::endl;
	//if (nlevel<5){
	Iterate(*kalman,toymatrix1);
	ncount += kalman->Get_ncount();
	std::cout << ncount << std::endl;
	nlevel --;
	//}else{
	//std::cout << nlevel << std::endl;
	//}
	
	delete kalman;
      }
    }
    
    // ncount += toymatrix.Get_ncount();
    
  }
  
  
}
