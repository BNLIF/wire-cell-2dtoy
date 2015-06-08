#include "WireCell2dToy/ToyMatrixMarkov.h"

using namespace WireCell;
#include "TMath.h"

WireCell2dToy::ToyMatrixMarkov::ToyMatrixMarkov(WireCell2dToy::ToyMatrix *toymatrix1,WireCell::GeomCellSelection *allmcell1){
  ncount = 0;
  toymatrix = toymatrix1; //save the matrix in here
  allmcell = allmcell1;
  mcindex = toymatrix->Get_mcindex();
  
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(*toymatrix);  // hold the current results 
  make_first_guess();

  
}

WireCell2dToy::ToyMatrixMarkov::~ToyMatrixMarkov(){
}

void WireCell2dToy::ToyMatrixMarkov::make_first_guess(){
  Iterate(*toymatrixkalman);
  //need to save the results???
  
  cur_chi2 = toymatrix->Get_Chi2();
  cur_dof = toymatrix->Get_ndf();

  toymatrix->Update_pred();

  for (int i=0;i!= (*allmcell).size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*)(*allmcell).at(i);
    double charge = toymatrix->Get_Cell_Charge(mcell,1);
    double charge_err = toymatrix->Get_Cell_Charge(mcell,2);
    
    if (charge ==0 && charge_err ==0 ){
      cur_cell_status.push_back(0); // removed for matrix calculation
    }else{
      cur_cell_status.push_back(1); // not removed
    }
    
    if (charge > 1500){ // hard coded to be fixed later
      cur_cell_pol.push_back(1); // on
    }else{
      cur_cell_pol.push_back(0); //off
    }

    //calculate residual for each cell, part of toymatrix
    // with residual and polarity, one can calculate probability
    cell_res.push_back(toymatrix->Get_residual(mcell));
    cell_charge.push_back(charge);
    // std::cout << cur_cell_status.at(i) << " " << cur_cell_pol.at(i) <<
    //   " " << cell_res.at(i) << " " << std::endl;
  }
  
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
