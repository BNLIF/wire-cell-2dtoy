#include "WireCell2dToy/ToyWalking.h"

using namespace WireCell;

WireCell2dToy::ToyNWalking::ToyNWalking(WireCell::MergeSpaceCell *start_cell, WireCell::MergeSpaceCellMap& mcells_map, WireCell::MergeSpaceCellSelection& used_cell, WireCell::MergeSpaceCellSelection& must_cell)
  : start_cell(start_cell)
  , mcells_map(mcells_map)
  , used_cell(used_cell)
  , must_cell(must_cell)
{
  Iterate(start_cell);
}

void WireCell2dToy::ToyNWalking::Iterate(MergeSpaceCell *curr_cell){
  cells.push_back(curr_cell);
  
  for (int i=0;i!=mcells_map[curr_cell].size();i++){
    MergeSpaceCell *next_cell = mcells_map[curr_cell].at(i);
    auto it3 = find(cells.begin(), cells.end(), next_cell);
    if (it3 == cells.end()){
      auto it2 = find(used_cell.begin(), used_cell.end(), next_cell);
      if (it2 == used_cell.end()){
	auto it1 = find(must_cell.begin(), must_cell.end(), next_cell);
	if (it1 != must_cell.end())
	  Iterate(next_cell);
      }
    }
  }
}

WireCell2dToy::ToyNWalking::~ToyNWalking(){
}


WireCell2dToy::ToyWalking::ToyWalking(WireCell::MergeSpaceCell *start_cell, Point start_point, WireCell::MergeSpaceCell *target_cell, Point target_point, WireCell::MergeSpaceCellMap& mcells_map, int counter_limit)
  : start_cell(start_cell)
  , target_cell(target_cell)
  , mcells_map(mcells_map)
  , start_point(start_point)
  , target_point(target_point)
  , counter_limit(counter_limit)
{
  dist = 1e9;
  counter = 0;
  global_counter = 0;
  must_cells.clear();
  MergeSpaceCellSelection curr_cells;
  Iterate(start_cell,curr_cells,0);
  
}


WireCell2dToy::ToyWalking::ToyWalking(WireCell::MergeSpaceCell *start_cell, Point start_point, WireCell::MergeSpaceCell *target_cell, Point target_point, WireCell::MergeSpaceCellMap& mcells_map, WireCell::MergeSpaceCellSelection must_cells, int counter_limit)
  : start_cell(start_cell)
  , target_cell(target_cell)
  , mcells_map(mcells_map)
  , start_point(start_point)
  , target_point(target_point)
  , counter_limit(counter_limit)
  , must_cells(must_cells)
{
  dist = 1e9;
  counter = 0;
  global_counter = 0;
  MergeSpaceCellSelection curr_cells;
  Iterate(start_cell,curr_cells,0);
  
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
    if (dis < dist && counter < counter_limit && global_counter < 25*counter_limit){
      // go in
      //rank everything according to size ... 

      //std::cout << target_cell->Get_Center().x/units::cm << " " << curr_cell->Get_Center().x/units::cm  << " " << curr_cell->Get_Center().y/units::cm << " " << curr_cell->Get_Center().z/units::cm << " " << curr_cells.size() << std::endl;

      if (target_cell->Get_Center().x > curr_cell->Get_Center().x){
	MergeSpaceCellSet1 msc_set;
	for (int i=0;i!=mcells_map[curr_cell].size();i++){
	  msc_set.insert(mcells_map[curr_cell].at(i));
	}
	
	for (auto it2 = msc_set.begin(); it2!= msc_set.end(); it2++){
	  auto it = find(curr_cells.begin(),curr_cells.end(),*it2);
	  auto it1 = find(must_cells.begin(),must_cells.end(),*it2);
	  if (it == curr_cells.end() && (it1 != must_cells.end() || must_cells.size()==0)){ 
	    // std::cout << "abc1: " << (*it2)->Get_Center().x/units::cm << " " << (*it2)->Get_Center().y/units::cm << " " << (*it2)->Get_Center().z/units::cm << " " <<  std::endl;
	    Iterate(*it2, curr_cells,dis);
	  
	  }
	}
      }else{
	MergeSpaceCellSet msc_set;
	for (int i=0;i!=mcells_map[curr_cell].size();i++){
	  msc_set.insert(mcells_map[curr_cell].at(i));
	}

	//std::cout << mcells_map[curr_cell].size() << " " << msc_set.size() << std::endl;
	
	for (auto it2 = msc_set.begin(); it2!= msc_set.end(); it2++){
	  auto it = find(curr_cells.begin(),curr_cells.end(),*it2);
	  auto it1 = find(must_cells.begin(),must_cells.end(),*it2);
	  if (it == curr_cells.end() && (it1 != must_cells.end() || must_cells.size()==0)){
	    //std::cout << "abc2: " << (*it2)->Get_Center().x/units::cm << " " << (*it2)->Get_Center().y/units::cm << " " << (*it2)->Get_Center().z/units::cm << " " <<  std::endl;
	    Iterate(*it2, curr_cells,dis);
	  }
	}
      }
      
      //
    }else{
      //counter ++;
    }
  }

  dis -= dis1;
  //move out the last element;
  curr_cells.pop_back();
  // tried_cells.push_back(curr_cell);
}


WireCell2dToy::ToyWalking::~ToyWalking(){
}
