#include "WireCell2dToy/ToyMatrixIterate_Only.h"

using namespace WireCell;

WireCell2dToy::ToyMatrixIterate_Only::~ToyMatrixIterate_Only(){

}

WireCell2dToy::ToyMatrixIterate_Only::ToyMatrixIterate_Only(WireCell2dToy::ToyMatrix &toymatrix, WireCell2dToy::MergeToyTiling* mergetiling)
  : toymatrix(toymatrix)
  , mergetiling(mergetiling)
{
  wirechargemap = mergetiling->wcmap();
  cellmap = mergetiling->cmap();
  
  // get all the cells,
  all_cells = mergetiling->get_allcell();
 
  ncount = 0;
  nlevel = 0;

  WireCell::GeomCellSelection good_cells;
  WireCell::GeomWireSelection used_wires;
  WireCell::GeomCellSelection bad_cells;
  WireCell::GeomCellSelection tried_cells;

  std::cout << "Test: " << all_cells.size() << std::endl;

  Iterate(all_cells,good_cells,used_wires,bad_cells,tried_cells);
 
}

void WireCell2dToy::ToyMatrixIterate_Only::Iterate(GeomCellSelection remaining_cells, GeomCellSelection good_cells, GeomWireSelection used_wires, GeomCellSelection bad_cells, GeomCellSelection& tried_cells){
  nlevel ++;

  for (int i=0;i!=remaining_cells.size();i++){
    GeomCellSelection remaining_cells_1 = remaining_cells;
    GeomCellSelection good_cells_1 = good_cells;
    GeomWireSelection used_wires_1 = used_wires;
    GeomCellSelection bad_cells_1 = bad_cells;

    MergeGeomCell *mcell = (MergeGeomCell*)remaining_cells.at(i);
    
    auto it1 = find(tried_cells.begin(),tried_cells.end(),mcell);
    if (it1 == tried_cells.end()){
      // put it into tried cells
      tried_cells.push_back(mcell);

      // put it into good_cells
      auto it2 = find(good_cells_1.begin(),good_cells_1.end(),mcell);
      if (it2 == good_cells_1.end()){
	good_cells_1.push_back(mcell);
	// put its wires into used_wires
	GeomWireSelection wires = cellmap[mcell];
	for (int j=0;j!=wires.size();j++){
	  if (wirechargemap[wires.at(j)] > 10){
	    auto it3 = find(used_wires_1.begin(),used_wires_1.end(),wires.at(j));
	    if (it3 == used_wires_1.end())
	      used_wires_1.push_back(wires.at(j));
	  }
	}
      }
      
      // remove it from remaining_cells
      auto it3 = find(remaining_cells_1.begin(),remaining_cells_1.end(),mcell);
      remaining_cells_1.erase(it3);

      GeomCellSelection to_be_removed;
      // fill bad_cells and further reduce remaining_cells
      for (int j=0;j!=remaining_cells_1.size();j++){
	GeomWireSelection wires = cellmap[remaining_cells_1.at(j)];
	for (int k=0;k!=wires.size();k++){
	  if (wirechargemap[wires.at(k)]>10){
	    auto it4 = find(used_wires.begin(),used_wires.end(),wires.at(k));
	    if (it4 != used_wires.end()){
	      to_be_removed.push_back(remaining_cells_1.at(j));
	      break;
	    }
	  }
	}
      }
      for (int j=0;j!=to_be_removed.size();j++){
	auto it4 = find(remaining_cells_1.begin(),remaining_cells_1.end(),to_be_removed.at(j));
	remaining_cells_1.erase(it4);
	auto it5 = find(bad_cells_1.begin(),bad_cells_1.end(),to_be_removed.at(j));
	if (it5 == bad_cells_1.end())
	  bad_cells_1.push_back(to_be_removed.at(j));
      }
      
      if (remaining_cells_1.size() == 0 ){
	// if bad_cells + good_cells == all_cells, solve
	std::vector<int> already_removed; 
	std::vector<int> no_need_remove; 
	for (int j=0;j!=good_cells_1.size();j++){
	  int index = toymatrix.Get_mcindex(good_cells_1.at(j));
	  no_need_remove.push_back(index);
	}
	for (int j=0; j!= bad_cells_1.size();j++){
	  int index = toymatrix.Get_mcindex(bad_cells_1.at(j));
	  already_removed.push_back(index);
	}
	 WireCell2dToy::ToyMatrixKalman toykalman(already_removed, no_need_remove, toymatrix,0,0);
	 ncount ++;
	 
	 std::cout << "Test: " << all_cells.size() << " " << bad_cells_1.size() + good_cells_1.size() << " " << toykalman.Get_numz() << " " << nlevel << " " << ncount << std::endl;
      }else{
	// otherwise, go deep
	Iterate(remaining_cells_1, good_cells_1, used_wires_1, bad_cells_1, tried_cells);
      }
      

    }
  }
  
  
  nlevel -- ;

}


