#include "WCP2dToy/BlobMetric.h"

using namespace WCP;

WCP2dToy::BlobMetric::BlobMetric(){
  rm_cell_true = 0;
  rm_cell_false = 0;
  el_cell_true = 0;
  el_cell_false = 0;

  charge_rm_cell_true = 0;
  charge_rm_cell_false = 0;
  charge_el_cell_true = 0;
  charge_el_cell_false = 0;
}


void WCP2dToy::BlobMetric::Print(){
  std::cout << "Summary of Simple DeBlob " << std::endl;
  std::cout << "Remaining Cells Containing Truth     : " << rm_cell_true << "  charge:" << charge_rm_cell_true << " true charge:" << Tcharge_rm_cell_true << std::endl;
  std::cout << "Remaining Cells Not Containing Truth : " << rm_cell_false << "  charge:" << charge_rm_cell_false << std::endl;
  std::cout << "Eliminated Cells Containing Truth    : " << el_cell_true << "  charge:" << charge_el_cell_true << " true charge:" << Tcharge_el_cell_true << std::endl;
  std::cout << "Eliminated Cells Not Containing Truth: " << el_cell_false << "  charge:" << charge_el_cell_false << std::endl;

  
}

WCP2dToy::BlobMetric::~BlobMetric(){
}


void WCP2dToy::BlobMetric::Add(WCP2dToy::SimpleBlobToyTiling &blobtiling, WCP::CellChargeMap& ccmap){
  GeomCellSelection sbcells = blobtiling.Get_SB_Cells();
  GeomCellSelection cells = blobtiling.Get_Cells();

  bool contain_truth;
  bool pass_threshold;

  for (int i = 0; i!=sbcells.size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*)sbcells.at(i);
    for (int j=0;j!=mcell->get_allcell().size();j++){
      const GeomCell *cell = mcell->get_allcell().at(j);
      
      if (ccmap.find(cell)==ccmap.end()){
	contain_truth = false;
      }else{
	contain_truth = true;
      }
      
      auto it = find(cells.begin(),cells.end(),cell);
      if (it == cells.end()){
	pass_threshold = false;
      }else{
	pass_threshold = true;
      }

      if (pass_threshold == true && contain_truth == true){
	rm_cell_true ++;
	charge_rm_cell_true += blobtiling.Get_Cell_Charge(cell,1);
	Tcharge_rm_cell_true += ccmap[cell];
      }else if (pass_threshold == true && contain_truth == false){
	rm_cell_false ++;
	charge_rm_cell_false += blobtiling.Get_Cell_Charge(cell,1);
      }else if (pass_threshold == false && contain_truth == true){
	el_cell_true ++ ;
	charge_el_cell_true += 0;
	Tcharge_el_cell_true += ccmap[cell];
      }else if (pass_threshold == false && contain_truth == false){
	el_cell_false ++;
	charge_el_cell_false += 0;
      }
      
      
    }
  }
  
}
