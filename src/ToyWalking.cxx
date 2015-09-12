#include "WireCell2dToy/ToyWalking.h"

using namespace WireCell;

WireCell2dToy::ToyWalking::ToyWalking(WireCell::MergeSpaceCell *start_cell, Point start_point, WireCell::MergeSpaceCell *target_cell, Point target_point, WireCell::MergeSpaceCellMap& mcells_map)
  : start_cell(start_cell)
  , target_cell(target_cell)
  , mcells_map(mcells_map)
  , start_point(start_point)
  , target_point(target_point)
{
  dist = 1e9;
  MergeSpaceCellSelection curr_cells;
  Iterate(start_cell,curr_cells,0);
  
}

void WireCell2dToy::ToyWalking::Iterate(MergeSpaceCell *curr_cell, MergeSpaceCellSelection &curr_cells, double dis){
  curr_cells.push_back(curr_cell);

  float dis1=0; 
  if (curr_cells.size()>=2){
    // Point p1,p2;
    // if (curr_cell == target_cell){
    //   p1 = target_point;
    // }else{
    //   p1 = curr_cell->Get_Center();
    // }
    // if (curr_cells.at(curr_cells.size()-2) == start_cell){
    //   p2 = start_point;
    // }else{
    //   p2 = curr_cells.at(curr_cells.size()-2)->Get_Center();
    // }
    
    // dis1 = pow(p1.x - p2.x,2)
    //    + pow(p1.y - p2.y,2)
    //   + pow(p1.z - p2.z,2);
    dis1 = 1;

    dis += dis1;
  }
  
  // std::cout << curr_cell->Get_Center().x/units::cm << " " << target_cell->Get_Center().x/units::cm << " " << mcells_map[curr_cell].size() << std::endl;

  if (curr_cell == target_cell){
    //end of it, calculate distance ... 
    if (dis < dist){
      dist = dis;
      cells.clear();
      cells = curr_cells;
    }
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

  dis -= dis1;
  
  //move out the last element;
  curr_cells.pop_back();
}


WireCell2dToy::ToyWalking::~ToyWalking(){
}
