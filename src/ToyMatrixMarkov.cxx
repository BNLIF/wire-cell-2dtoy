#include "WireCell2dToy/ToyMatrixMarkov.h"

using namespace WireCell;
#include "TMath.h"




WireCell2dToy::ToyMatrixMarkov::ToyMatrixMarkov(WireCell2dToy::ToyMatrix *toymatrix1,WireCell::GeomCellSelection *allmcell1){
  ncount = 0;
  first_flag = 0;
  toymatrix = toymatrix1; //save the matrix in here
  allmcell = allmcell1;
  mcindex = toymatrix->Get_mcindex();
  
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(*toymatrix);  // hold the current results 
  make_guess();

  

  
}

WireCell2dToy::ToyMatrixMarkov::~ToyMatrixMarkov(){
}

void WireCell2dToy::ToyMatrixMarkov::make_guess(){
  if (first_flag ==0){
    Iterate(*toymatrixkalman);
    first_flag = 1;
  }
  //need to save the results???
  
  cur_chi2 = toymatrix->Get_Chi2();
  cur_dof = toymatrix->Get_ndf();

  toymatrix->Update_pred();

  double_t max_res = 0;

  for (int i=0;i!= (*allmcell).size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*)(*allmcell).at(i);
    double charge = toymatrix->Get_Cell_Charge(mcell,1);
    double charge_err = toymatrix->Get_Cell_Charge(mcell,2);
    
    cur_cell_index.push_back(toymatrix->Get_mcindex(mcell));

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

    if (toymatrix->Get_residual(mcell) > max_res) 
      max_res = toymatrix->Get_residual(mcell);

    cell_charge.push_back(charge);

    if (cur_cell_status.at(i)==1){
      double ratio;
      
      if (cell_res.at(i)!=0){
	ratio = cell_charge.at(i)/cell_res.at(i);
      }else{
	ratio = cell_charge.at(i)/10.;
      }
      
      if (ratio > 10.){
	cell_prob.push_back(0.9); 
      }else if (ratio <=10 && ratio>=-10.){
	cell_prob.push_back(0.7);
      }else{
	cell_prob.push_back(0.05);
      }
    }else{
      if (cell_res.at(i) > 1000){
	cell_prob.push_back(0.5);
      }else{
	cell_prob.push_back(0.1);
      }
    }
    //   std::cout << cur_cell_status.at(i) << " " << cur_cell_pol.at(i) <<
    //   " " << cell_res.at(i) << " " << std::endl;
  }

  //rank stuff ... 
  // first probability, second bias
  for (int i=0;i!= (*allmcell).size();i++){
    CellRankPair a(i,cell_prob.at(i)+cell_res.at(i)/max_res/10.);
    cell_set.insert(a);
  }

  // for (auto it= cell_set.begin(); it!=cell_set.end();it++){
  //   std::cout << (*it).first << " " << (*it).second << std::endl;
  // }
  
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
