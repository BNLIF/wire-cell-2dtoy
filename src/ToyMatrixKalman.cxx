#include "WireCell2dToy/ToyMatrixKalman.h"
using namespace WireCell;

#include "TMath.h"

int WireCell2dToy::ToyMatrixKalman::Cal_numz(WireCell2dToy::ToyMatrix &toymatrix){
  
  //std::cout << no_need_remove.size() << std::endl;
  if (no_need_remove.size()!=0){
  TMatrixD MA2(mwindex,no_need_remove.size());
  TMatrixD MA2T(no_need_remove.size(),mwindex);
  TMatrixD MC2(no_need_remove.size(),no_need_remove.size());
	    
  const TMatrixD *MA_big = toymatrix.Get_MA();
  const TMatrixD *VBy_inv = toymatrix.Get_VBy_inv();

  int index2 = 0;
  for (int k=0;k!=mcindex;k++){ //loop all possibility
    auto it3 = find(no_need_remove.begin(),no_need_remove.end(),k);
    if (it3 != no_need_remove.end()){ // not in removed
      for (int j=0;j!=mwindex;j++){
	MA2(j,index2) = (*MA_big)(j,k);
      }
      index2 ++;
    }
  }
	    
  MA2T.Transpose(MA2);
  MC2 = (MA2T) * (*VBy_inv) * (MA2);
  TMatrixDEigen Eigen2(MC2);
  TVectorD EigenValue2(Eigen2.GetEigenValuesRe());
	    
  int numz2 = 0;
  for (int k=0;k!=no_need_remove.size();k++){
    if (fabs(EigenValue2[k])<1e-5){
      numz2 ++;
    }
  }
  return numz2;
  }else{
    return 0;
  }
}

void WireCell2dToy::ToyMatrixKalman::init(WireCell2dToy::ToyMatrix& toymatrix){
  ncount = 0;
  numz = 0;
  
  mcindex = toymatrix.Get_mcindex();
  mwindex = toymatrix.Get_mwindex();
  swindex = toymatrix.Get_swindex();

  
  
  int n_removed = already_removed.size();

  // MA = new TMatrixD(mwindex,mcindex-n_removed);
  // MAT = new TMatrixD(mcindex-n_removed,mwindex);

  // MC = new TMatrixD(mcindex-n_removed,mcindex-n_removed);
  // MC_inv = new TMatrixD(mcindex-n_removed,mcindex-n_removed);

  //std::cout << " Xin " << already_removed.size() << " " << no_need_remove.size() << std::endl;

  TMatrixD MA(mwindex,mcindex-n_removed);
  TMatrixD MAT(mcindex-n_removed,mwindex);
  
  TMatrixD MC(mcindex-n_removed,mcindex-n_removed);
  TMatrixD MC_inv(mcindex-n_removed,mcindex-n_removed);
  
  
  
  const TMatrixD *MA_big = toymatrix.Get_MA();
  const TMatrixD *VBy_inv = toymatrix.Get_VBy_inv();
  
  int index2 = 0;
  for (int i=0;i!=mcindex;i++){
    auto it1 = find(already_removed.begin(),already_removed.end(),i);
    if (it1==already_removed.end()){
      for (int j=0;j!=mwindex;j++){
  	MA(j,index2) = (*MA_big)(j,i);
      }
      index2 ++;
    }
  }
  MAT.Transpose(MA);

  
  MC = (MAT) * (*VBy_inv) * (MA);
  


  TMatrixDEigen Eigen(MC);
  TVectorD EigenValue(Eigen.GetEigenValuesRe());
    
  for (int i=0;i!=mcindex-n_removed;i++){
    if (fabs(EigenValue[i])<1e-5){
      numz ++;
    }
  }
  

  if (numz!=0){
    for (int i=0;i!=mcindex;i++){ //loop all possibility
      //std::cout << i << " " << mcindex << std::endl;

      auto it1 = find(already_removed.begin(),already_removed.end(),i);
      auto it2 = find(no_need_remove.begin(),no_need_remove.end(),i);
      if (it1==already_removed.end() && it2 == no_need_remove.end()){ // not in removed
    	TMatrixD MA1(mwindex,mcindex-n_removed-1);
    	TMatrixD MA1T(mcindex-n_removed-1,mwindex);
    	TMatrixD MC1(mcindex-n_removed-1,mcindex-n_removed-1);
	
    	int index1 = 0;
    	for (int k=0;k!=mcindex;k++){
	  

    	  auto it1 = find(already_removed.begin(),already_removed.end(),k);
    	  if (it1==already_removed.end() && k!=i ){
    	    for (int j=0;j!=mwindex;j++){
    	      MA1(j,index1) = (*MA_big)(j,k);
    	    }
    	    index1 ++;
    	  }
    	}
     	MA1T.Transpose(MA1);
	
    	MC1 = (MA1T) * (*VBy_inv) * (MA1);
	
    	TMatrixDEigen Eigen1(MC1);
    	TVectorD EigenValue1(Eigen1.GetEigenValuesRe());
	
    	int numz1 = 0;
    	for (int k=0;k!=mcindex-n_removed-1;k++){
    	  if (fabs(EigenValue1[k])<1e-5){
    	    numz1 ++;
    	  }
    	}
	
    	if (numz == numz1) {
	  no_need_remove.push_back(i);
	  // if (check_flag!=0){
	  // //   TMatrixD MA2(mwindex,no_need_remove.size()+1);
	  // //   TMatrixD MA2T(no_need_remove.size()+1,mwindex);
	  // //   TMatrixD MC2(no_need_remove.size()+1,no_need_remove.size()+1);
	    
	  // //   int index2 = 0;
	  // //   for (int k=0;k!=mcindex;k++){ //loop all possibility
	  // //     auto it3 = find(no_need_remove.begin(),no_need_remove.end(),k);
	  // //     if (it3 != no_need_remove.end() || k==i){ // not in removed
	  // //   	for (int j=0;j!=mwindex;j++){
	  // //   	  MA2(j,index2) = (*MA_big)(j,k);
	  // //   	}
	  // //   	index2 ++;
	  // //     }
	  // //   }
	    
	  // //   MA2T.Transpose(MA2);
	  // //   MC2 = (MA2T) * (*VBy_inv) * (MA2);
	  // //   TMatrixDEigen Eigen2(MC2);
	  // //   TVectorD EigenValue2(Eigen2.GetEigenValuesRe());
	    
	  // //   int numz2 = 0;
	  // //   for (int k=0;k!=no_need_remove.size()+1;k++){
	  // //     if (fabs(EigenValue2[k])<1e-5){
	  // //   	numz2 ++;
	  // //     }
	  // //   }
	  // //   if (numz2==0){
	    
	  //   // }else{
	  //   // //   //std::cout << numz2 << std::endl;
	  //   // }
	  // }else{
	  //   no_need_remove.push_back(i);
	  // }
	}
      }
    }
  }else{
    MC_inv = MC;
    MC_inv.Invert();
    
    const TMatrixD *MB = toymatrix.Get_MB();
    const TVectorD *Wy = toymatrix.Get_Wy();
    
    TVectorD Cxt(mcindex-n_removed);
    TVectorD dCxt(mcindex-n_removed);
    TMatrixD Vx(mcindex-n_removed,mcindex-n_removed);
    TMatrixD Vx_inv(mcindex-n_removed,mcindex-n_removed);
    
    Cxt = (MC_inv) * (MAT) * (*VBy_inv) * (*MB) * (*Wy);
    Vx_inv = (MAT) * (*VBy_inv) * (MA);
    Vx = Vx_inv;
    Vx.Invert();
    
    
    for (int i=0;i!=mcindex-n_removed;i++){
      dCxt[i] = sqrt( Vx(i,i)) * 1000.;
    }

    TVectorD sol = (*MB) * (*Wy) - (MA) * (Cxt);
    TVectorD sol1 =  (*VBy_inv) * sol;
    double chi2 = sol * (sol1)/1e6;
    
    for (int i=0;i!=mcindex-already_removed.size();i++){
      if (Cxt[i] <0){
	chi2 += 10*pow(Cxt[i]/dCxt[i],2);
      }
    }

    // std::cout << numz << " " << chi2 << " ";
    // for (int i=0;i!=already_removed.size();i++){
    //   std::cout << already_removed.at(i) << " ";
    // }
    // std::cout << std::endl;

    if (chi2 < toymatrix.Get_Chi2() || toymatrix.Get_Chi2()==-1){
      //copy the Cx etc      
      toymatrix.Set_chi2(chi2);
      toymatrix.Set_ndf(mwindex - (mcindex - already_removed.size()) );
      toymatrix.Set_Solve_Flag(2);
      int index = 0;
      for (int i=0;i!=mcindex;i++){
    	auto it = find(already_removed.begin(),already_removed.end(),i);
    	if (it==already_removed.end()){
    	  toymatrix.Set_value(Cxt[index],i);
    	  toymatrix.Set_error(dCxt[index],i);
    	  index ++;
    	}else{
    	  toymatrix.Set_value(0,i);
    	  toymatrix.Set_error(0,i);
    	}
      }
   }

    ncount ++;
  }
  
}

WireCell2dToy::ToyMatrixKalman::ToyMatrixKalman(WireCell2dToy::ToyMatrix& toymatrix){
  check_flag = 0;
  init(toymatrix);
       
  //  std::cout << numz << std::endl;
  //  std::cout << already_removed.size() << " " << no_need_remove.size() << std::endl;
}


WireCell2dToy::ToyMatrixKalman::ToyMatrixKalman(std::vector<int>& already_removed1, std::vector<int>& no_need_remove1, WireCell2dToy::ToyMatrix& toymatrix, int check){
  already_removed = already_removed1;
  no_need_remove = no_need_remove1;
  check_flag = check;
  init(toymatrix);
       
  //  std::cout << numz << std::endl;
  //  std::cout << already_removed.size() << " " << no_need_remove.size() << std::endl;
}



WireCell2dToy::ToyMatrixKalman::~ToyMatrixKalman(){
  // already_removed.clear();
  // no_need_remove.clear();
  
  // delete MA, MAT;
  // delete MC, MC_inv;
  //delete Eigen, EigenValue;
}
