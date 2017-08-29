#include "WireCell2dToy/ChargeSolving.h"


using namespace WireCell;

WireCell2dToy::ChargeSolving::ChargeSolving(const WireCell::GeomDataSource& gds, LowmemTiling& tiling)
  : gds(gds)
  , tiling(tiling)
{
  ndirect_solved = 0;
  nL1_solved = 0;

  // initialize the cell weight map
  init_cell_weight_map();
  
  // divide into small groups
  divide_groups();

}

WireCell2dToy::ChargeSolving::~ChargeSolving(){
}

void WireCell2dToy::ChargeSolving::init_cell_weight_map(){
  GeomCellSelection& one_wire_cells = tiling.get_one_good_wire_cells();
  GeomCellSelection& two_wire_cells = tiling.get_two_good_wire_cells();
  GeomCellSelection& three_wire_cells = tiling.get_three_good_wire_cells();

  GeomCellMap cell_wire_map = tiling.get_cell_wires_map();
  GeomWireMap wire_cell_map = tiling.get_wire_cells_map();

  for (auto it = three_wire_cells.begin(); it!=three_wire_cells.end(); it++){
    const GeomCell *mcell = (*it);
    cell_weight_map[mcell] = 1;
  }

  for (auto it = two_wire_cells.begin(); it!=two_wire_cells.end(); it++){
    const GeomCell *mcell = (*it);
    cell_weight_map[mcell] = 1;
  }

  for (auto it = one_wire_cells.begin(); it!=one_wire_cells.end(); it++){
    const GeomCell *mcell = (*it);
    cell_weight_map[mcell] = 1;
  }
  
}



void WireCell2dToy::ChargeSolving::divide_groups(){

  GeomCellSelection& one_wire_cells = tiling.get_one_good_wire_cells();
  GeomCellSelection& two_wire_cells = tiling.get_two_good_wire_cells();
  GeomCellSelection& three_wire_cells = tiling.get_three_good_wire_cells();
  
  WireChargeMap& wire_charge_map = tiling.get_wire_charge_map();
  WireChargeMap& wire_charge_error_map = tiling.get_wire_charge_error_map();
  
  GeomWireSelection good_wires = tiling.get_all_good_wires();
  GeomWireSelection bad_wires = tiling.get_all_bad_wires();

  GeomCellMap cell_wire_map = tiling.get_cell_wires_map();
  GeomWireMap wire_cell_map = tiling.get_wire_cells_map();
  
  // original cell/wire holders
  GeomCellSetp all_cells;
  GeomWireSetp all_wires;

  std::copy(good_wires.begin(), good_wires.end(), std::inserter(all_wires, all_wires.end()));
  std::copy(three_wire_cells.begin(), three_wire_cells.end(), std::inserter(all_cells, all_cells.end()));
  std::copy(two_wire_cells.begin(), two_wire_cells.end(), std::inserter(all_cells, all_cells.end()));
  std::copy(one_wire_cells.begin(), one_wire_cells.end(), std::inserter(all_cells, all_cells.end()));

  
  // final cell/wire holders
  // std::cout << all_wires.size() << " " << all_cells.size() << std::endl;

  while(all_wires.size()){
    //temporary cell/wire holders
    GeomCellSelection temp_cells;
    GeomWireSelection temp_wires;
    
    GeomCellSelection grouped_cells;
    GeomWireSelection grouped_wires;
    
    // start one good wire, move it in
    auto it5 = all_wires.begin();
    const GeomWire* mwire = *it5;
    all_wires.erase(mwire);
    //all_wires.pop_front();
    temp_wires.push_back(mwire);
    grouped_wires.push_back(mwire);
    
    while(temp_wires.size()!=0){
      // find all the cells associated with it, move in 
      temp_cells.clear();
      for (auto it = temp_wires.begin(); it!= temp_wires.end(); it++){
	const GeomWire *mwire = *it;
	for (auto it1 = wire_cell_map[mwire].begin(); it1!= wire_cell_map[mwire].end(); it1++){
	  const GeomCell *mcell = *it1;
	  auto it2 = all_cells.find(mcell);
	  if (it2 != all_cells.end()){
	    // fill in the temp cells ... 
	    temp_cells.push_back(mcell);
	    // remove from the all_cells
	    all_cells.erase(it2);
	    // fill in grouped_wires
	    grouped_cells.push_back(mcell);
	  }
	}
      }
      temp_wires.clear();
      // find all the wires associated with these cells, move in
      for (auto it = temp_cells.begin(); it!= temp_cells.end(); it++){
	const GeomCell *mcell = *it;
	for (auto it1 = cell_wire_map[mcell].begin(); it1!= cell_wire_map[mcell].end(); it1++){
	  const GeomWire *mwire = *it1;
	  auto it2 = all_wires.find(mwire);
	  if (it2!= all_wires.end()){
	    temp_wires.push_back(mwire);
	    all_wires.erase(it2);
	    grouped_wires.push_back(mwire);
	  }
	}
      }
      
    }

    // std::cout << "Xin " << grouped_wires.size() << " " << grouped_cells.size() << std::endl;
    final_cells_vec.push_back(grouped_cells);
    final_wires_vec.push_back(grouped_wires);

    MatrixSolver *matrix = new MatrixSolver(grouped_cells, grouped_wires, cell_wire_map, wire_cell_map, wire_charge_map, wire_charge_error_map);

    // do L1 solve ... 
    if (matrix->get_solve_flag()==0){
      matrix->L1_Solve(cell_weight_map);
    }
    
    if (matrix->get_solve_flag()==1){
      ndirect_solved ++;
    }else if (matrix->get_solve_flag()==2){
      nL1_solved++;
    }

    for (auto it = grouped_cells.begin(); it!=grouped_cells.end(); it++){
      MergeGeomCell *mcell = (MergeGeomCell*)(*it);
      //test
      // if (matrix->get_solve_flag()==1){
      // 	ccmap[mcell] = matrix->get_mcell_charge(mcell) ;
      // }else{
      // 	ccmap[mcell] = 10000;
      // }
      
      ccmap[mcell] = matrix->get_mcell_charge(mcell) ;
      //   std::cout << matrix->get_mcell_charge(mcell) << std::endl;
    }
    
    group_matrices.push_back(matrix);
    
  }

  
  
  // std::cout << one_wire_cells.size() << " " << two_wire_cells.size() << " " << three_wire_cells.size() << " " << good_wires.size() << " " << bad_wires.size() << " " << cell_wire_map.size() << " " << wire_cell_map.size() << std::endl;
  // divide all wires into seperate groups ...

}

double WireCell2dToy::ChargeSolving::get_chi2(){
  return chi2;
}

void WireCell2dToy::ChargeSolving::Update_ndf_chi2(){
  ndf = 0;
  L1_ndf.clear();
  direct_ndf.clear();
  for (auto it = group_matrices.begin(); it!=group_matrices.end(); it++){
    MatrixSolver* matrix = (*it);
    if (matrix->get_solve_flag()==1){
      ndf += matrix->get_direct_ndf();
      direct_ndf.push_back(matrix->get_direct_ndf());
      //std::cout << "1 " << matrix->get_direct_ndf() << std::endl;
    }else if (matrix->get_solve_flag()==2){
      ndf += matrix->get_L1_ndf();
      L1_ndf.push_back(matrix->get_L1_ndf());
      //std::cout << "2 " << matrix->get_L1_ndf() << std::endl;
    }
    // std::cout << ndf << std::endl;
  }

  chi2 = 0;
  L1_chi2_base.clear();
  L1_chi2_penalty.clear();
  direct_chi2.clear();
  for (auto it = group_matrices.begin(); it!=group_matrices.end(); it++){
    MatrixSolver* matrix = (*it);
    if (matrix->get_solve_flag()==1){
      chi2 += matrix->get_direct_chi2();
      direct_chi2.push_back(matrix->get_direct_chi2());
      //std::cout << "1 " << matrix->get_direct_chi2() << std::endl;
    }else if (matrix->get_solve_flag()==2){
      chi2 += matrix->get_L1_chi2_base();
      L1_chi2_base.push_back(matrix->get_L1_chi2_base());
      L1_chi2_penalty.push_back(matrix->get_L1_chi2_penalty());
      //std::cout << "2 " << matrix->get_L1_chi2_base() << std::endl;
    }
    //    std::cout << chi2 << std::endl;
  }
  
  
}

int WireCell2dToy::ChargeSolving::get_ndf(){
  
  
  return ndf;
}


