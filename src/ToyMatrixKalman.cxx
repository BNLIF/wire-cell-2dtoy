#include "WireCell2dToy/ToyMatrixKalman.h"
using namespace WireCell;

#include "TMath.h"
void WireCell2dToy::ToyMatrixKalman::init(WireCell2dToy::ToyMatrix &toymatrix){
  ncount = 0;
  numz = 0;
  
  mcindex = toymatrix.Get_mcindex();
  mwindex = toymatrix.Get_mwindex();
  swindex = toymatrix.Get_swindex();

  //  std::cout << " Xin " << already_removed.size() << " " << no_need_remove.size() << std::endl;
  
  int n_removed = already_removed.size();

  MA = new TMatrixD(mwindex,mcindex-n_removed);
  MAT = new TMatrixD(mcindex-n_removed,mwindex);

  MC = new TMatrixD(mcindex-n_removed,mcindex-n_removed);
  MC_inv = new TMatrixD(mcindex-n_removed,mcindex-n_removed);
  
  const TMatrixD *MA_big = toymatrix.Get_MA();
  const TMatrixD *VBy_inv = toymatrix.Get_VBy_inv();

  int index2 = 0;
  for (int i=0;i!=mcindex;i++){
    auto it1 = find(already_removed.begin(),already_removed.end(),i);
    if (it1==already_removed.end()){
      for (int j=0;j!=mwindex;j++){
	(*MA)(j,index2) = (*MA_big)(j,i);
      }
      index2 ++;
    }
  }
  MAT->Transpose(*MA);
  *MC = (*MAT) * (*VBy_inv) * (*MA);

  Eigen = new TMatrixDEigen(*MC);
  EigenValue = new TVectorD(Eigen->GetEigenValuesRe());
    
  for (int i=0;i!=mcindex-n_removed;i++){
    if (fabs((*EigenValue)[i])<1e-5){
      numz ++;
    }
  }
  //  std::cout << " Xin " << numz<< std::endl;


  if (numz!=0){
    for (int i=0;i!=mcindex;i++){ //loop all possibility
      auto it1 = find(already_removed.begin(),already_removed.end(),i);
      auto it2 = find(no_need_remove.begin(),no_need_remove.end(),i);
      if (it1==already_removed.end() && it2 == no_need_remove.end()){ // not in removed
	TMatrixD *MA1 = new TMatrixD(mwindex,mcindex-n_removed-1);
	TMatrixD *MA1T = new TMatrixD(mcindex-n_removed-1,mwindex);
	TMatrixD *MC1 = new TMatrixD(mcindex-n_removed-1,mcindex-n_removed-1);
	
	int index1 = 0;
	for (int k=0;k!=mcindex;k++){
	  auto it1 = find(already_removed.begin(),already_removed.end(),k);
	  if (it1==already_removed.end() && k!=i ){
	    for (int j=0;j!=mwindex;j++){
	      (*MA1)(j,index1) = (*MA_big)(j,k);
	    }
	    index1 ++;
	  }
	}
	MA1T->Transpose(*MA1);
	*MC1 = (*MA1T) * (*VBy_inv) * (*MA1);
	
	TMatrixDEigen *Eigen1 = new TMatrixDEigen(*MC1);
	TVectorD *EigenValue1 = new TVectorD(Eigen1->GetEigenValuesRe());
	
	int numz1 = 0;
       	for (int k=0;k!=mcindex-n_removed-1;k++){
       	  if (fabs((*EigenValue1)[k])<1e-5){
       	    numz1 ++;
       	  }
       	}
	
	//	std::cout << i << " " << numz << " " << numz1 << std::endl;
	if (numz == numz1) no_need_remove.push_back(i);
	
	delete MA1, MA1T, MC1;
	delete Eigen1, EigenValue1;
      }
    }
  }else{
    // *MC_inv = *MC;
    // MC_inv->Invert();
    
    // const TMatrixD *MB = toymatrix.Get_MB();
    // const TVectorD *Wy = toymatrix.Get_Wy();
    
    // TVectorD *Cxt = new TVectorD(mcindex-n_removed);
    // TVectorD *dCxt = new TVectorD(mcindex-n_removed);
    // TMatrixD *Vx = new TMatrixD(mcindex-n_removed,mcindex-n_removed);
    // TMatrixD *Vx_inv = new TMatrixD(mcindex-n_removed,mcindex-n_removed);

    // *Cxt = (*MC_inv) * (*MAT) * (*VBy_inv) * (*MB) * (*Wy);
    // *Vx_inv = (*MAT) * (*VBy_inv) * (*MA);
    // *Vx = *Vx_inv;
    // Vx->Invert();
    
    // for (int i=0;i!=mcindex-n_removed;i++){
    //   (*dCxt)[i] = sqrt( (*Vx)(i,i)) * 1000.;
    // }

    // TVectorD sol = (*MB) * (*Wy) - (*MA) * (*Cxt);
    // TVectorD sol1 =  (*VBy_inv) * sol;
    // double chi2 = sol * (sol1)/1e6;
    
    
    // if (chi2 < toymatrix.Get_Chi2() || toymatrix.Get_Chi2()==-1){
    //   //copy the Cx etc      
    //   toymatrix.Set_chi2(chi2);
    //   int index = 0;
    //   for (int i=0;i!=mcindex;i++){
    // 	auto it = find(already_removed.begin(),already_removed.end(),i);
    // 	if (it==already_removed.end()){
    // 	  toymatrix.Set_value((*Cxt)[index],i);
    // 	  toymatrix.Set_error((*dCxt)[index],i);
    // 	  index ++;
    // 	}else{
    // 	  toymatrix.Set_value(0,i);
    // 	  toymatrix.Set_error(0,i);
    // 	}
    //   }
    // }

    //std::cout << ncount << " " << chi2 << " " << toymatrix.Get_Chi2() << std::endl;
    ncount ++;
    // delete dCxt, Cxt, Vx, Vx_inv;
  }
  
}

WireCell2dToy::ToyMatrixKalman::ToyMatrixKalman(WireCell2dToy::ToyMatrix &toymatrix){
  init(toymatrix);
       
  //  std::cout << numz << std::endl;
  //  std::cout << already_removed.size() << " " << no_need_remove.size() << std::endl;
}


WireCell2dToy::ToyMatrixKalman::ToyMatrixKalman(std::vector<int>& already_removed1, std::vector<int>& no_need_remove1, WireCell2dToy::ToyMatrix &toymatrix){
  already_removed = already_removed1;
  no_need_remove = no_need_remove1;
  init(toymatrix);
       
  //  std::cout << numz << std::endl;
  //  std::cout << already_removed.size() << " " << no_need_remove.size() << std::endl;
}



WireCell2dToy::ToyMatrixKalman::~ToyMatrixKalman(){
  delete MA, MAT;
  delete MC, MC_inv;
  delete Eigen, EigenValue;
}
