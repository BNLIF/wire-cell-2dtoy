#include "WireCell2dToy/ToyMetric.h"

using namespace WireCell;

WireCell2dToy::ToyMetric::ToyMetric(){
  rm_cell_true = 0;
  rm_cell_false = 0;
  el_cell_true = 0;
  el_cell_false = 0;

  solve_condition[0] = 0;
  solve_condition[1] = 0;
  solve_condition[2] = 0;
}

void WireCell2dToy::ToyMetric::AddSolve(int cond){
  if (cond ==0 ){
    solve_condition[0] ++;
  }else if (cond==1){
    solve_condition[1] ++;
  }else{
    solve_condition[2]++;
  }
}

WireCell2dToy::ToyMetric::~ToyMetric(){
}

void WireCell2dToy::ToyMetric::Add(GeomCellSelection &allmcell,WireCell2dToy::ToyMatrix& toymatrix, CellChargeMap& ccmap){
  for (int j=0;j!=allmcell.size();j++){
    MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
    double charge = toymatrix.Get_Cell_Charge(mcell,1);
    double charge_err = toymatrix.Get_Cell_Charge(mcell,2);
    bool pass_threshold = false;
    //std::cout << charge << " " << charge_err << std::endl;
    if (charge + charge_err > threshold) pass_threshold = true;
    bool contain_truth = mcell->CheckContainTruthCell(ccmap);
    
    if (pass_threshold == true && contain_truth == true){
      rm_cell_true ++;
      charge_rm_cell_true += mcell->GetTruthCharge();
    }else if (pass_threshold == true && contain_truth == false){
      rm_cell_false ++;
      charge_rm_cell_false += charge;
    }else if (pass_threshold == false && contain_truth == true){
      el_cell_true ++ ;
      charge_el_cell_true += mcell->GetTruthCharge();
    }else if (pass_threshold == false && contain_truth == false){
      el_cell_false ++;
      charge_el_cell_false += charge;
    }
    
  }
  
}

void WireCell2dToy::ToyMetric::Print(){
  std::cout << "Remaining Cells Containing Truth     : " << rm_cell_true << " " << charge_rm_cell_true << std::endl;
  std::cout << "Remaining Cells Not Containing Truth : " << rm_cell_false << " " << charge_rm_cell_false << std::endl;
  std::cout << "Eliminated Cells Containing Truth    : " << el_cell_true << " " << charge_el_cell_true << std::endl;
  std::cout << "Eliminated Cells Not Containing Truth: " << el_cell_false << " " << charge_el_cell_false << std::endl;

  std::cout << "Not       Solved Case: " << solve_condition[0] << std::endl;
  std::cout << "Direct    Solved Case: " << solve_condition[1] << std::endl;
  std::cout << "Iterative Solved Case: " << solve_condition[2] << std::endl;
}
