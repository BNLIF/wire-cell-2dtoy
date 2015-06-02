#include "WireCell2dToy/ToyMetric.h"

using namespace WireCell;

WireCell2dToy::ToyMetric::ToyMetric(){
  rm_cell_true = 0;
  rm_cell_false = 0;
  el_cell_true = 0;
  el_cell_false = 0;
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
    }else if (pass_threshold == true && contain_truth == false){
      rm_cell_false ++;
    }else if (pass_threshold == false && contain_truth == true){
      el_cell_true ++ ;
    }else if (pass_threshold == false && contain_truth == false){
      el_cell_false ++;
    }
    
  }
  
}

void WireCell2dToy::ToyMetric::Print(){
  std::cout << "Remaining Cells Containing Truth     : " << rm_cell_true << std::endl;
  std::cout << "Remaining Cells Not Containing Truth : " << rm_cell_false << std::endl;
  std::cout << "Eliminated Cells Containing Truth    : " << el_cell_true << std::endl;
  std::cout << "Eliminated Cells Not Containing Truth: " << el_cell_false << std::endl;
}
