#include "WireCell2dToy/ToyWalking.h"

using namespace WireCell;

WireCell2dToy::ToyWalking::ToyWalking(WireCell::MergeSpaceCell *start_cell, WireCell::MergeSpaceCell *target_cell, WireCell::MergeSpaceCellMap& mcells_map)
  : start_cell(start_cell)
  , target_cell(target_cell)
  , mcells_map(mcells_map)
{
  dist = 1e9;
  MergeSpaceCellSelection curr_cells;
  Iterate(start_cell,curr_cells,0);
  
}

void WireCell2dToy::ToyWalking::Iterate(MergeSpaceCell *curr_cell, MergeSpaceCellSelection &curr_cells, double dis){
  curr_cells.push_back(curr_cell);

  if (curr_cells.size()>=2){
    dis += pow(curr_cell->Get_Center().x - curr_cells.at(curr_cells.size()-2)->Get_Center().x,2)
       + pow(curr_cell->Get_Center().y - curr_cells.at(curr_cells.size()-2)->Get_Center().y,2)
      + pow(curr_cell->Get_Center().z - curr_cells.at(curr_cells.size()-2)->Get_Center().z,2);
  }

  if (curr_cell == target_cell){
    //end of it, calculate distance ... 
    dist = dis;
    cells.clear();
    cells = curr_cells;
  }else{
    if (dis < dist){
      // go in
      for (int i=0;i!=mcells_map[curr_cell].size();i++){
	auto it = find(curr_cells.begin(),curr_cells.end(),mcells_map[curr_cell].at(i));
	if (it == curr_cells.end())
	  Iterate(mcells_map[curr_cell].at(i), curr_cells,dis);
      }
    }
  }

  //move out the last element;
  curr_cells.pop_back();
}


WireCell2dToy::ToyWalking::~ToyWalking(){
}
