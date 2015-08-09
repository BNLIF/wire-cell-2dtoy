#include "WireCell2dToy/ToyMatrixIterate.h"

using namespace WireCell;

#include "TMath.h"

WireCell2dToy::ToyMatrixIterate::ToyMatrixIterate(WireCell2dToy::ToyMatrix &toymatrix, std::vector<int>& already_removed, int recon_t){
  prev_ncount = -1;
  ncount = 0;
  nlevel = 0;
  recon_threshold = recon_t;
  
  std::vector<int> no_need_remove;

  
  //std::cout << "Number of zeros: " << toymatrixkalman->Get_numz() << std::endl;
  //int numz = toymatrix.Get_mcindex() - toymatrix.Get_mwindex();
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, no_need_remove, toymatrix, 0,0);  
  
  //estimated_loop = TMath::Factorial(toymatrix.Get_mcindex())/TMath::Factorial(toymatrix.Get_mcindex()-numz)/TMath::Factorial(numz)/25.;
  //std::cout << estimated_loop << std::endl;
  estimated_loop = TMath::Factorial(toymatrix.Get_mcindex())/TMath::Factorial(toymatrix.Get_mcindex()-toymatrixkalman->Get_numz())/TMath::Factorial(toymatrixkalman->Get_numz())/25.;
  
  // if(toymatrixkalman->Get_numz()==toymatrix.Get_mcindex()){
  //   delete toymatrixkalman;
  //   toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, no_need_remove, toymatrix, 0,1); 
  // }
  
  std::cout << estimated_loop << " " << toymatrix.Get_mcindex() << " " << toymatrixkalman->Get_numz() << std::endl;

  if (estimated_loop < 1e6 && toymatrixkalman->Get_numz()!=toymatrix.Get_mcindex()){
    time_flag = 0;
    delete toymatrixkalman;
    toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, no_need_remove, toymatrix, 0,1); 
    
    
    Iterate(*toymatrixkalman,toymatrix);
  }


}

WireCell2dToy::ToyMatrixIterate::ToyMatrixIterate(WireCell2dToy::ToyMatrix &toymatrix, int recon_t, float limit){
  prev_ncount = -1;
  ncount = 0;
  nlevel = 0;
  recon_threshold = recon_t;

 
  
  //std::cout << "Number of zeros: " << toymatrixkalman->Get_numz() << std::endl;
  //int numz = toymatrix.Get_mcindex() - toymatrix.Get_mwindex();
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(toymatrix,0);  
  // estimated_loop = TMath::Factorial(toymatrix.Get_mcindex())/TMath::Factorial(toymatrix.Get_mcindex()-numz)/TMath::Factorial(numz)/25.;
  // std::cout << estimated_loop << std::endl;
  estimated_loop = TMath::Factorial(toymatrix.Get_mcindex())/TMath::Factorial(toymatrix.Get_mcindex()-toymatrixkalman->Get_numz())/TMath::Factorial(toymatrixkalman->Get_numz())/25.;
    

  std::cout << estimated_loop << " " << toymatrix.Get_mcindex() << " " << toymatrixkalman->Get_numz() << std::endl;
  
  // if(toymatrixkalman->Get_numz()==toymatrix.Get_mcindex()){
  //   delete toymatrixkalman;
  //   toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(toymatrix,1); 
  // }
  
  if (estimated_loop < limit && toymatrixkalman->Get_numz()!=toymatrix.Get_mcindex()){
    time_flag = 0;
    delete toymatrixkalman;
    toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(toymatrix,1);  
    //  std::cout << toymatrixkalman->Get_numz() << std::endl;
    Iterate(*toymatrixkalman,toymatrix);
  }
  //if not use time information ... 
}


WireCell2dToy::ToyMatrixIterate::ToyMatrixIterate(WireCell2dToy::ToyMatrix &toycur, WireCell2dToy::MergeToyTiling &mergecur, WireCell::GeomCellSelection &cells, int recon_t ){
  std::vector<int> already_removed; //dummy
  std::vector<int> no_need_remove; //to be added
 
  GeomCellSelection allmcell_c =  mergecur.get_allcell();
  for (int i = 0;i!=allmcell_c.size();i++){
    MergeGeomCell *mcell_c = (MergeGeomCell*)allmcell_c[i];
    int index_c = toycur.Get_mcindex(mcell_c);
    
    auto it = find(cells.begin(),cells.end(),mcell_c);
    if (it != cells.end()){
      no_need_remove.push_back(index_c);
    }
  }

  //std::cout << "Xin: " << no_need_remove.size() << " " << cells.size() << std::endl;

  
   for (int i=0;i!=allmcell_c.size();i++){
    auto it = find(no_need_remove.begin(),no_need_remove.end(),i);
    if (it==no_need_remove.end())
      already_removed.push_back(i);
  }
  no_need_remove.clear();
  

  
  
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, no_need_remove, toycur,1);
  std::cout << "With Time: " << toymatrixkalman->Get_numz() << " " << allmcell_c.size() << " " <<  already_removed.size() << std::endl;
  // Find a sub-set that is not degenerated
  // put things into no_need_remove
  find_subset(*toymatrixkalman,toycur,no_need_remove);
  already_removed.clear();
  prev_ncount = -1;
  ncount = 0;
  nlevel = 0;
  // recalculate
  delete toymatrixkalman;
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, no_need_remove, toycur,1);
  toycur.Set_Solve_Flag(0);
  toycur.Set_chi2(-1);

  //std::cout << "With Time: " << toymatrixkalman->Get_numz() << std::endl;
  
  
  estimated_loop = TMath::Factorial(toycur.Get_mcindex()-no_need_remove.size())/TMath::Factorial(toycur.Get_mcindex()-no_need_remove.size()-toymatrixkalman->Get_numz())/TMath::Factorial(toymatrixkalman->Get_numz())/5;
  std::cout << "With Cluster: " << estimated_loop << " " << toycur.Get_mcindex() << " " << no_need_remove.size() << " " << toymatrixkalman->Get_numz() <<std::endl;

  if (estimated_loop < 2e5 && toymatrixkalman->Get_numz()!=toycur.Get_mcindex()){
    time_flag = 1;
    delete toymatrixkalman;
    toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, no_need_remove, toycur, 0,1); 
    
    Iterate(*toymatrixkalman,toycur);
  }

  
  // if no good, release ... 
  if (toycur.Get_Chi2() > 5*(toycur.Get_ndf()+0.1)){
    no_need_remove.clear();
    already_removed.clear();
    
    
    toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, no_need_remove, toycur, 0,0);  
    estimated_loop = TMath::Factorial(toycur.Get_mcindex())/TMath::Factorial(toycur.Get_mcindex()-toymatrixkalman->Get_numz())/TMath::Factorial(toymatrixkalman->Get_numz())/25.;
    std::cout << "Try again: " << estimated_loop << " " << toycur.Get_mcindex() << " " << toymatrixkalman->Get_numz() << std::endl;
    
    if (estimated_loop < 5e5 && toymatrixkalman->Get_numz()!=toycur.Get_mcindex()){
      time_flag = 0;
      delete toymatrixkalman;
      toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, no_need_remove, toycur, 0,1); 
      Iterate(*toymatrixkalman,toycur);
    }else{
      toycur.Set_Solve_Flag(0);
    }
  }


}




void WireCell2dToy::ToyMatrixIterate::UseTime(WireCell2dToy::ToyMatrix &toybefore, WireCell2dToy::ToyMatrix &toycur, WireCell2dToy::ToyMatrix &toyafter, WireCell2dToy::MergeToyTiling &mergebefore, WireCell2dToy::MergeToyTiling &mergecur, WireCell2dToy::MergeToyTiling &mergeafter){
  
  GeomCellSelection allmcell_p = mergebefore.get_allcell();
  GeomCellSelection allmcell_c = mergecur.get_allcell();
  GeomCellSelection allmcell_n = mergeafter.get_allcell();
  
  //form two vectors 
  std::vector<int> already_removed; //dummy
  std::vector<int> no_need_remove; //to be added
  for (int i=0;i!=allmcell_c.size();i++){
    MergeGeomCell *mcell_c = (MergeGeomCell*)allmcell_c[i];
    int index_c = toycur.Get_mcindex(mcell_c);
    
    if (toybefore.Get_Solve_Flag()!=0 ){
      for (int j=0;j!=allmcell_p.size();j++){
	MergeGeomCell *mcell_p = (MergeGeomCell*)allmcell_p[j];
	int index_p = toybefore.Get_mcindex(mcell_p);
	double charge = toybefore.Get_Cell_Charge(mcell_p,1);
	if ( charge > recon_threshold && mcell_c->Overlap(*mcell_p)){
	  auto it = find(no_need_remove.begin(),no_need_remove.end(),index_c);
	  if (it == no_need_remove.end()){
	    no_need_remove.push_back(index_c);
	  }
	}
      }
    }
    if (toyafter.Get_Solve_Flag()!=0){
      for (int j=0;j!=allmcell_n.size();j++){
	MergeGeomCell *mcell_n = (MergeGeomCell*)allmcell_n[j];
	int index_n = toyafter.Get_mcindex(mcell_n);
	double charge = toyafter.Get_Cell_Charge(mcell_n,1);
	if ( charge > recon_threshold && mcell_c->Overlap(*mcell_n)){
	  auto it = find(no_need_remove.begin(),no_need_remove.end(),index_c);
	  if (it == no_need_remove.end()){
	    no_need_remove.push_back(index_c);
	  }
	}
      }
    }
  }
  


  for (int i=0;i!=allmcell_c.size();i++){
    auto it = find(no_need_remove.begin(),no_need_remove.end(),i);
    if (it==no_need_remove.end())
      already_removed.push_back(i);
  }
  no_need_remove.clear();
  

  
  delete toymatrixkalman;
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, no_need_remove, toycur,1);
  std::cout << "With Time: " << toymatrixkalman->Get_numz() << " " << allmcell_c.size() << " " <<  already_removed.size() << std::endl;
  // Find a sub-set that is not degenerated
  // put things into no_need_remove
  find_subset(*toymatrixkalman,toycur,no_need_remove);
  already_removed.clear();
  prev_ncount = -1;
  ncount = 0;
  nlevel = 0;
  // recalculate
  delete toymatrixkalman;
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, no_need_remove, toycur,1);
  toycur.Set_Solve_Flag(0);
  toycur.Set_chi2(-1);

  //std::cout << "With Time: " << toymatrixkalman->Get_numz() << std::endl;
  
  
  estimated_loop = TMath::Factorial(toycur.Get_mcindex()-no_need_remove.size())/TMath::Factorial(toycur.Get_mcindex()-no_need_remove.size()-toymatrixkalman->Get_numz())/TMath::Factorial(toymatrixkalman->Get_numz())/5;
  std::cout << "With Time: " << estimated_loop << " " << toycur.Get_mcindex() << " " << no_need_remove.size() << " " << toymatrixkalman->Get_numz() <<std::endl;
  
  
  if (estimated_loop < 1e6){
    time_flag = 1;
    // std::cout << toycur.Get_Solve_Flag() << " " << toycur.Get_Chi2() << " " << toycur.Get_ndf() << std::endl;
    Iterate(*toymatrixkalman,toycur);
    //std::cout << toycur.Get_Solve_Flag() << " " << toycur.Get_Chi2() << " " << toycur.Get_ndf() << std::endl;
  }else{
    time_flag = 0;
    //start with the smallest area? need a reasonable chi-square, and a limit of iterations
    Iterate_simple(*toymatrixkalman,toycur);
  }
  
 

  if (toycur.Get_Chi2() > 3*(toycur.Get_ndf()+0.1)){
    no_need_remove.clear();
    already_removed.clear();
    delete toymatrixkalman;
    toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, no_need_remove, toycur,0);
    Iterate_simple1(*toymatrixkalman,toycur);
  }
  

}


WireCell2dToy::ToyMatrixIterate::~ToyMatrixIterate(){
  if (estimated_loop < 1e6){
    delete toymatrixkalman;
  }
}



void WireCell2dToy::ToyMatrixIterate::Iterate(WireCell2dToy::ToyMatrixKalman &toymatrix,WireCell2dToy::ToyMatrix &toymatrix1){
  nlevel ++;
  if (toymatrix.Get_numz()!=0 && toymatrix.Cal_numz(toymatrix1)==0){
    for (int i=0;i!=toymatrix.Get_mcindex();i++){
      auto it1 = find(toymatrix.Get_already_removed().begin(),toymatrix.Get_already_removed().end(),i);
      auto it2 = find(toymatrix.Get_no_need_remove().begin(),toymatrix.Get_no_need_remove().end(),i);
      if (it1 == toymatrix.Get_already_removed().end() && it2 == toymatrix.Get_no_need_remove().end()){
	std::vector<int> already_removed = toymatrix.Get_already_removed();
	already_removed.push_back(i);
	WireCell2dToy::ToyMatrixKalman kalman(already_removed,toymatrix.Get_no_need_remove(),toymatrix1,0);

	

	//	std::cout << nlevel << std::endl;
	//if (nlevel<5){
	Iterate(kalman,toymatrix1);
	//move to "no need to remove"??
	toymatrix.Get_no_need_remove().push_back(i);
	ncount += kalman.Get_ncount();


	//if (ncount != prev_ncount)
	//  std::cout << ncount << " " << already_removed.size() << " " << toymatrix.Get_no_need_remove().size() << std::endl;
	//     " " << already_removed.at(0) << " " << already_removed.at(1) << " " << already_removed.at(2) << " " << already_removed.at(3) << " " << already_removed.at(4) << " " << already_removed.at(5) << 
	//     // " " << already_removed.at(6) << " " << already_removed.at(7) << " " << already_removed.at(8) << " " << already_removed.at(9) << " " << already_removed.at(10) << " " << already_removed.at(11) << 
	//     //" " << already_removed.at(12) << " " << already_removed.at(13) << " " << already_removed.at(14) << " " << already_removed.at(15) << " " << already_removed.at(16) << " " << already_removed.at(17) << 
	//     //" " << already_removed.at(18) << " " << already_removed.at(19) << " " << already_removed.at(20) << " " << already_removed.at(21) << " " << already_removed.at(22) << " " << already_removed.at(23) << 
	//     "  " <<  toymatrix.Get_no_need_remove().size() << std::endl;

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



void WireCell2dToy::ToyMatrixIterate::find_subset(WireCell2dToy::ToyMatrixKalman &toymatrix,WireCell2dToy::ToyMatrix &toymatrix1, std::vector<int>& vec){
  //std::cout << toymatrix.Get_already_removed().size() << " " << toymatrix.Get_no_need_remove().size() << " " << vec.size() << std::endl;

  if (toymatrix.Get_numz()!=0){
    for (int i=0;i!=toymatrix.Get_mcindex();i++){
      auto it1 = find(toymatrix.Get_already_removed().begin(),toymatrix.Get_already_removed().end(),i);
      auto it2 = find(toymatrix.Get_no_need_remove().begin(),toymatrix.Get_no_need_remove().end(),i);
      if (it1 == toymatrix.Get_already_removed().end() && it2 == toymatrix.Get_no_need_remove().end()){
	std::vector<int> already_removed = toymatrix.Get_already_removed();
	already_removed.push_back(i);
	WireCell2dToy::ToyMatrixKalman kalman(already_removed,toymatrix.Get_no_need_remove(),toymatrix1,0);
	
	if (kalman.Get_numz()==toymatrix.Get_numz()){
	}else{
	  find_subset(kalman,toymatrix1,vec);
	  //move to "no need to remove"??
	  toymatrix.Get_no_need_remove().push_back(i);
	  if (vec.size()!=0) break;
	}
      }
    }
  }else{
    for (int i=0;i!=toymatrix.Get_mcindex();i++){
      auto it1 = find(toymatrix.Get_already_removed().begin(),toymatrix.Get_already_removed().end(),i);
      if (it1 == toymatrix.Get_already_removed().end()){
	vec.push_back(i);
      }
    }
  }
}





void WireCell2dToy::ToyMatrixIterate::Iterate_simple(WireCell2dToy::ToyMatrixKalman &toymatrix,WireCell2dToy::ToyMatrix &toymatrix1){
  nlevel ++;
  
  int ncount_cut = 50000;
  
  //std::cout << "Simple " << nlevel << " " << ncount << std::endl;
  if (toymatrix.Get_numz()!=0&& toymatrix.Cal_numz(toymatrix1)==0&& (ncount < ncount_cut&&(toymatrix1.Get_Chi2()<0 || toymatrix1.Get_Chi2()>1.5*(toymatrix1.Get_ndf()+0.1)))){
    for (int i=0;i!=toymatrix.Get_mcindex();i++){
      auto it1 = find(toymatrix.Get_already_removed().begin(),toymatrix.Get_already_removed().end(),i);
      auto it2 = find(toymatrix.Get_no_need_remove().begin(),toymatrix.Get_no_need_remove().end(),i);
      if (it1 == toymatrix.Get_already_removed().end() && it2 == toymatrix.Get_no_need_remove().end()){
	std::vector<int> already_removed = toymatrix.Get_already_removed();
	already_removed.push_back(i);
	WireCell2dToy::ToyMatrixKalman kalman(already_removed,toymatrix.Get_no_need_remove(),toymatrix1,1);

	Iterate_simple(kalman,toymatrix1);
	toymatrix.Get_no_need_remove().push_back(i);
	ncount += kalman.Get_ncount();

	if (ncount != prev_ncount && ncount%10000==0){
	//if (ncount != prev_ncount){
	  std::cout << "Simple: " << ncount << " " << toymatrix1.Get_Chi2() << std::endl;
	  //if (ncount != prev_ncount){
	  // std::cout << "Simple: " << nlevel << " " << toymatrix.Get_numz() << " " << kalman.Get_numz() << " " << ncount << " " << toymatrix1.Get_Chi2() <<  " ";
	  // for (int kk=0;kk!=already_removed.size();kk++){
	  //   std::cout << already_removed.at(kk) << " " ;
	  // }
	  // std::cout << std::endl;
	  // for (int kk = 0; kk!=toymatrix.Get_no_need_remove().size();kk++){
	  //   std::cout << toymatrix.Get_no_need_remove().at(kk) << " ";
	  // }
	  // std::cout << std::endl;
	}

	prev_ncount = ncount;
	nlevel --;

	if (!(ncount < ncount_cut &&(toymatrix1.Get_Chi2()<0 || toymatrix1.Get_Chi2()>1.5*(toymatrix1.Get_ndf()+0.1))&& toymatrix.Cal_numz(toymatrix1)==0)) 
	  break;

      }
    }
  }
}

void WireCell2dToy::ToyMatrixIterate::Iterate_simple1(WireCell2dToy::ToyMatrixKalman &toymatrix,WireCell2dToy::ToyMatrix &toymatrix1){
  nlevel ++;
  int ncount_cut = 100000;
  //std::cout << "Simple1 " << nlevel << " " << ncount << std::endl;
  if (toymatrix.Get_numz()!=0&& toymatrix.Cal_numz(toymatrix1)==0&& (toymatrix1.Get_Chi2()<0 || toymatrix1.Get_Chi2()>3*(toymatrix1.Get_ndf()+0.1)) && ncount < ncount_cut){
    for (int i=0;i!=toymatrix.Get_mcindex();i++){
      auto it1 = find(toymatrix.Get_already_removed().begin(),toymatrix.Get_already_removed().end(),i);
      auto it2 = find(toymatrix.Get_no_need_remove().begin(),toymatrix.Get_no_need_remove().end(),i);
      if (it1 == toymatrix.Get_already_removed().end() && it2 == toymatrix.Get_no_need_remove().end()){
	std::vector<int> already_removed = toymatrix.Get_already_removed();
	already_removed.push_back(i);
	WireCell2dToy::ToyMatrixKalman kalman(already_removed,toymatrix.Get_no_need_remove(),toymatrix1,0);

	Iterate_simple1(kalman,toymatrix1);
	toymatrix.Get_no_need_remove().push_back(i);
	ncount += kalman.Get_ncount();

	if (ncount != prev_ncount && ncount%10000==0)
	//if (ncount != prev_ncount )
	std::cout << "Simple1: " << ncount << " " << nlevel << " " << toymatrix1.Get_Chi2() << std::endl;
	
	prev_ncount = ncount;
	nlevel --;

	if (!((toymatrix1.Get_Chi2()<0 || toymatrix1.Get_Chi2()>3*(toymatrix1.Get_ndf()+0.1)) && ncount < ncount_cut&& toymatrix.Cal_numz(toymatrix1)==0))
	  break;

      }
    }
  }
}
