#include "WireCell2dToy/ToyWalking.h"

using namespace WireCell;

WireCell2dToy::ToyWalking::ToyWalking(WireCell::MergeSpaceCell *start_cell, Point start_point, WireCell::MergeSpaceCell *target_cell, Point target_point, WireCell::MergeSpaceCellMap& mcells_map, int counter_limit)
  : start_cell(start_cell)
  , target_cell(target_cell)
  , mcells_map(mcells_map)
  , start_point(start_point)
  , target_point(target_point)
  , counter_limit(counter_limit)
{
  dist = 1e9;
  MergeSpaceCellSelection curr_cells;
  Iterate(start_cell,curr_cells,0);
  counter = 0;
  global_counter = 0;
}

void WireCell2dToy::ToyWalking::Iterate(MergeSpaceCell *curr_cell, MergeSpaceCellSelection &curr_cells, double dis){
  curr_cells.push_back(curr_cell);
  global_counter ++;
  
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
  
  // std::cout << curr_cell->Get_Center().x/units::cm << " " 
  // 	    << curr_cell->Get_Center().y/units::cm << " " 
  // 	    << curr_cell->Get_Center().z/units::cm << " " 
  // 	    << mcells_map[curr_cell].size() << " " 
  // 	    << counter << std::endl;

  if (curr_cell == target_cell){
    counter ++;
    //    std::cout << counter << std::endl;
    //end of it, calculate distance ... 
    if (dis < dist){
      dist = dis;
      cells.clear();
      cells = curr_cells;
    }
  }else{
    if (dis < dist && counter < counter_limit && global_counter < 50*counter_limit){
      // go in
      for (int i=0;i!=mcells_map[curr_cell].size();i++){
	auto it = find(curr_cells.begin(),curr_cells.end(),mcells_map[curr_cell].at(i));
	//	auto it1 = find(tried_cells.begin(),tried_cells.end(),mcells_map[curr_cell].at(i));
	if (it == curr_cells.end() )
	  Iterate(mcells_map[curr_cell].at(i), curr_cells,dis);
      }
    }else{
      counter ++;
    }
  }

  dis -= dis1;
  //move out the last element;
  curr_cells.pop_back();
  // tried_cells.push_back(curr_cell);
}


WireCell2dToy::ToyWalking::~ToyWalking(){
}
