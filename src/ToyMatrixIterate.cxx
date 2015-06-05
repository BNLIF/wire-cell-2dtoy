#include "WireCell2dToy/ToyMatrixIterate.h"

using namespace WireCell;

#include "TMath.h"

WireCell2dToy::ToyMatrixIterate::ToyMatrixIterate(WireCell2dToy::ToyMatrix &toymatrix){
  prev_ncount = -1;
  ncount = 0;
  nlevel = 0;

 
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(toymatrix);
  
  std::cout << toymatrixkalman->Get_numz() << std::endl;

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
	WireCell2dToy::ToyMatrixKalman kalman(already_removed,toymatrix.Get_no_need_remove(),toymatrix1);

	

	//	std::cout << nlevel << std::endl;
	//if (nlevel<5){
	Iterate(kalman,toymatrix1);
	//move to "no need to remove"??
	toymatrix.Get_no_need_remove().push_back(i);

	ncount += kalman.Get_ncount();


	if (ncount != prev_ncount)
	  std::cout << ncount << 
	    " " << already_removed.at(0) << " " << already_removed.at(1) << " " << already_removed.at(2) << " " << already_removed.at(3) << " " << already_removed.at(4) << " " << already_removed.at(5) << 
	    // " " << already_removed.at(6) << " " << already_removed.at(7) << " " << already_removed.at(8) << " " << already_removed.at(9) << " " << already_removed.at(10) << " " << already_removed.at(11) << 
	    //" " << already_removed.at(12) << " " << already_removed.at(13) << " " << already_removed.at(14) << " " << already_removed.at(15) << " " << already_removed.at(16) << " " << already_removed.at(17) << 
	    //" " << already_removed.at(18) << " " << already_removed.at(19) << " " << already_removed.at(20) << " " << already_removed.at(21) << " " << already_removed.at(22) << " " << already_removed.at(23) << 
	    "  " <<  toymatrix.Get_no_need_remove().size() << std::endl;

	prev_ncount = ncount;
	nlevel --;

	//}else{
	//std::cout << nlevel << std::endl;
	//}
	
	//	delete kalman;
      }
    }
    
    // ncount += toymatrix.Get_ncount();
    
  }
  
  
}
