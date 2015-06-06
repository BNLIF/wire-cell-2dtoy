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
  // transT= new TMatrixD(*trans);
  // transT->Transpose(*trans);

  mcindex = toymatrix.Get_mcindex();
  mwindex = toymatrix.Get_mwindex();
  swindex = toymatrix.Get_swindex();
  numz = 0;
  TVectorD temp(mcindex);
  for (int i=0;i!=mcindex;i++){
    //std::cout << (*EigenValue)[i] << std::endl;
    if (fabs((*EigenValue)[i])<1e-5){
      for (int j=0;j!=mcindex;j++){
      	temp[j] += (*trans)(j,i);
      }
      numz ++;
    }
  }

  EigenSet eset;
  for (int i=0;i!=mcindex;i++){
    EigenPair abc(i,fabs(temp[i]));
    eset.insert(abc);
  }
  chi2 = 1e9;
  std::vector<int> flag;
  for (auto it = eset.begin(); it != eset.end(); it++){
    if (flag.size()!=numz){
      flag.push_back(it->first);
    }else{
      break;
    }
    //   std::cout << it->first << " " << it->second << std::endl;
  }
  
  //  trans->Print();
  //std::cout << (*trans)(44,45) << std::endl; // second number is column, first number is row  
  // std::cout << num_size << " " << numz << " " << TMath::Factorial(num_size)/TMath::Factorial(numz)/TMath::Factorial(num_size-numz) << std::endl;
  


  
  // for (int i=0;i!=numz;i++){
  //   flag.push_back(i);
  // }

  // numz = 26;
  // //flag.push_back(1);  
  // flag.push_back(2);  flag.push_back(3);  flag.push_back(4);
  // // flag.push_back(5);  flag.push_back(6);  
  // flag.push_back(10); flag.push_back(13); flag.push_back(14);
  // flag.push_back(15); flag.push_back(16); flag.push_back(17);  flag.push_back(19);
  // flag.push_back(21); flag.push_back(22); // flag.push_back(24);
  // flag.push_back(25); flag.push_back(26); flag.push_back(27); flag.push_back(28); //flag.push_back(29);
  // flag.push_back(30); flag.push_back(31); flag.push_back(32);  flag.push_back(34);
  // flag.push_back(35); flag.push_back(37); flag.push_back(38); //flag.push_back(39);
  // flag.push_back(41); flag.push_back(42);  flag.push_back(44);
  // //flag.push_back(45); 
  
 
  MA = new TMatrixD(mwindex,mcindex-numz);
  MAT = new TMatrixD(mcindex-numz,mwindex);
  Vx = new TMatrixD(mcindex-numz,mcindex-numz);
  Vx_inv = new TMatrixD(mcindex-numz,mcindex-numz);
  MC = new TMatrixD(mcindex-numz,mcindex-numz);
  MC_inv = new TMatrixD(mcindex-numz,mcindex-numz);

  Cxt = new TVectorD(mcindex-numz);
  dCxt = new TVectorD(mcindex-numz);

  // Cx = new TVectorD(mcindex);
  // dCx = new TVectorD(mcindex);
  
  std::cout << numz << " " << flag.size() << " " << mcindex-numz-1 << std::endl;
  ncount = 0;
  int flag1 = 1;
  int ncount1 = 0;
  
  double limit = TMath::Factorial(mcindex)/TMath::Factorial(numz)/TMath::Factorial(mcindex-numz);
  
  std::cout << limit << std::endl;

  if (flag.size()!=0){

    while( ncount1 < 100 && ncount <1e6 && ncount < fabs(2*limit)){
      if (ncount%100000==0 && ncount !=0) std::cout << ncount << std::endl;
      for (int i=0;i!=flag.size();i++){
	for (int j=0; j!=mcindex-numz-1; j++){
	  if (Solve(flag,toymatrix)==1 && chi2 < mwindex + 5*sqrt(mwindex)) flag1 = 0;
	  move(flag,i,mcindex);
	}
      }
      if (flag1 == 0) ncount1 ++;
    }
    

  }else{
    
    Solve(flag,toymatrix);
  }
  
}



WireCell2dToy::ToyMatrixExclusive::~ToyMatrixExclusive(){
  delete Eigen;
  delete EigenValue;
  delete MA, MAT;
  delete MC, MC_inv;
  delete Vx, Vx_inv;
  delete Cxt, dCxt;
  //delete Cx, dCx;
  delete trans;
  // delete transT;
}


int WireCell2dToy::ToyMatrixExclusive::Solve(std::vector<int>& flag, WireCell2dToy::ToyMatrix &toymatrix){
  const TMatrixD *MA_big = toymatrix.Get_MA();
  const TMatrixD *VBy = toymatrix.Get_VBy();
  const TMatrixD *VBy_inv = toymatrix.Get_VBy_inv();
  const TMatrixD *MB = toymatrix.Get_MB();
  const TVectorD *Wy = toymatrix.Get_Wy();
  
  int index = 0;
  for (int i=0;i!=mcindex;i++){
    auto it = find(flag.begin(),flag.end(),i);
    if (it==flag.end()){
      for (int j=0;j!=mwindex;j++){
	(*MA)(j,index) = (*MA_big)(j,i);
      }
      // std::cout << index << std::endl;
      index ++;
    }
  }

  //MA->Print();

  MAT->Transpose(*MA);
  *MC = (*MAT) * (*VBy_inv) * (*MA);
  double det = MC->Determinant();
  ncount ++; 
 
  if (fabs(det)<1e-5){
    return 0;
  }else{
    *MC_inv = *MC;
    MC_inv->Invert();
    *Cxt = (*MC_inv) * (*MAT) * (*VBy_inv) * (*MB) * (*Wy);
    *Vx_inv = (*MAT) * (*VBy_inv) * (*MA);
    *Vx = *Vx_inv;
    Vx->Invert();
    
    for (int i=0;i!=mcindex-numz;i++){
      (*dCxt)[i] = sqrt( (*Vx)(i,i)) * 1000.;
    }
    
    
    
    

    TVectorD sol = (*MB) * (*Wy) - (*MA) * (*Cxt);

    //sol.Print();
    
    TVectorD sol1 =  (*VBy_inv) * sol;
    chi2t = sol * (sol1)/1e6;
    
    for (int i=0;i!=mcindex-numz;i++){
      if ((*Cxt)[i] <0){
	chi2t += 10*pow((*Cxt)[i]/(*dCxt)[i],2);
      }
    }

    if (chi2t < chi2){
      //copy the Cx etc      
      chi2 = chi2t;
      toymatrix.Set_chi2(chi2);
      toymatrix.Set_ndf(mwindex - (mcindex - flag.size()) );
      toymatrix.Set_Solve_Flag(3);
      index = 0;
      for (int i=0;i!=mcindex;i++){
	auto it = find(flag.begin(),flag.end(),i);
	if (it==flag.end()){
	  toymatrix.Set_value((*Cxt)[index],i);
	  toymatrix.Set_error((*dCxt)[index],i);
	  index ++;
	}else{
	  toymatrix.Set_value(0,i);
	  toymatrix.Set_error(0,i);
	}
      }
      // for (int i=0;i!=numz;i++){
      // 	std::cout << flag[i] << std::endl;
      // }
      std::cout << "xin: " << chi2 << " " << ncount << std::endl;
    }
    //std::cout << det << std::endl;
    return 1;
  }

}

void WireCell2dToy::ToyMatrixExclusive::move(std::vector<int>& flag, int index, int max){
  int temp = flag.at(index);
  auto it = find(flag.begin(),flag.end(),temp);
  while(it!=flag.end()){
    temp ++;
    if (temp >=max) temp = 0;
    it = find(flag.begin(),flag.end(),temp);
  }
  flag.at(index) = temp;
}
