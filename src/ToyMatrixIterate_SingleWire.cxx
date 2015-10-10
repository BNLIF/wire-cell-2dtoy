#include "WireCell2dToy/ToyMatrixIterate_SingleWire.h"

using namespace WireCell;


WireCell2dToy::ToyMatrixIterate_SingleWire::ToyMatrixIterate_SingleWire(WireCell2dToy::ToyMatrix &toymatrix, WireCell2dToy::MergeToyTiling* mergetiling)
  : toymatrix(toymatrix)
  , mergetiling(mergetiling)
{
  // figure out single-wire cell, and no need to remove these
  WireCell::GeomCellSelection all_cells = mergetiling->get_allcell();
  WireCell::GeomCellSelection single_wire_cells = mergetiling->get_single_wire_cells();
  
  wirechargemap = mergetiling->wcmap();

  GeomCellSelection cells;
  for (int i=0;i!=all_cells.size();i++){
    auto it = find(single_wire_cells.begin(),single_wire_cells.end(),all_cells.at(i));
    if (it == single_wire_cells.end())
      cells.push_back(all_cells.at(i));
  }
  GeomCellMap cellmap = mergetiling->cmap();
  GeomWireMap wiremap = mergetiling->wmap();
  
  std::cout << "Test: " << cells.size() << std::endl;
  ncount = 0;
  GeomCellSelection tried_cell;
  Iterate(cells,single_wire_cells, tried_cell, cellmap,wiremap);
  

  
  
}

void WireCell2dToy::ToyMatrixIterate_SingleWire::Iterate(WireCell::GeomCellSelection cells, WireCell::GeomCellSelection single_cells, WireCell::GeomCellSelection tried_cell, WireCell::GeomCellMap cellmap, WireCell::GeomWireMap wiremap){
  GeomCellSelection tried_cells = tried_cell;
  
  for (int i1=0;i1!=cells.size();i1++){
     // stuff that will get passed ... 
    GeomCellSelection cells_1 = cells;
    GeomCellMap cellmap_1 = cellmap;
    GeomWireMap wiremap_1 = wiremap;
    GeomCellSelection single_wire_cells = single_cells;


    auto it = find(tried_cells.begin(),tried_cells.end(),cells.at(i1));
    if (it == tried_cells.end()){
      //remove the cell
      tried_cells.push_back(cells.at(i1));
      auto it1 = find(cells_1.begin(),cells_1.end(),cells.at(i1));
      cells_1.erase(it1);
      
      // Update the maps
      cellmap_1.erase(cells.at(i1));
      for (auto it2 = wiremap_1.begin(); it2!=wiremap_1.end();it2++){
	auto it3 = find(it2->second.begin(),it2->second.end(),cells.at(i1));
	if (it3 != it2->second.end())
	  it2->second.erase(it3);
      }
      
      
      int flag1 = 1;
      while(flag1 == 1){
	flag1 = 0;
	GeomCellSelection to_be_removed_cells;
	
	for (int i=0;i!=cells_1.size();i++){
	  MergeGeomCell *mcell = (MergeGeomCell*) cells_1.at(i);
	  GeomWireSelection wires = cellmap_1[mcell];
	  int flag = 0;
	  for (int j=0;j!=wires.size();j++){
	    MergeGeomWire *mwire = (MergeGeomWire*)wires.at(j);
	    if (wirechargemap[mwire] > 10){
	      if (wiremap_1[mwire].size()==1){
		flag = 1;
		break;
	      }
	    }
	  }
	  if (flag==1)
	    single_wire_cells.push_back(mcell);
	}
	
	for (int i=0;i!=single_wire_cells.size();i++){
	  MergeGeomCell *mcell = (MergeGeomCell*)single_wire_cells.at(i);
	  
	  GeomWireSelection wires = cellmap_1[mcell];
	  for (int j=0;j!=wires.size();j++){
	    MergeGeomWire *mwire = (MergeGeomWire*)wires.at(j);
	    GeomCellSelection tmp_cells = wiremap_1[mwire];
	    for (int k=0;k!=tmp_cells.size();k++){
	      auto it4 = find(single_wire_cells.begin(),single_wire_cells.end(),
			      tmp_cells.at(k));
	      if (it4 == single_wire_cells.end()){
		auto it5 = find(to_be_removed_cells.begin(), to_be_removed_cells.end(), tmp_cells.at(k));
		
		if (it5 == to_be_removed_cells.end())
		  to_be_removed_cells.push_back(tmp_cells.at(k));
	      }
	    }
	  }
	}
	
	if (to_be_removed_cells.size() !=0) flag1 = 1; 
	
	//Now remove these cells, and the update the two maps
	for (int i=0;i!=to_be_removed_cells.size();i++){
	  MergeGeomCell *mcell = (MergeGeomCell*)to_be_removed_cells.at(i);
	  // put into tried pool
	  auto it7 = find(tried_cells.begin(),tried_cells.end(),mcell);
	  if (it7 == tried_cells.end())
	    tried_cells.push_back(mcell);
	  //remove it
	  auto it4 = find(cells_1.begin(),cells_1.end(),mcell);
	  cells_1.erase(it4);
	  //Upate map
	  cellmap_1.erase(mcell);
	  for (auto it5 = wiremap_1.begin(); it5!=wiremap_1.end();it5++){
	    auto it6 = find(it5->second.begin(),it5->second.end(),mcell);
	    if (it6 != it5->second.end())
	      it5->second.erase(it6);
	  }
	}
      }
      // go to next level
      if (cells_1.size() >0){
	Iterate(cells_1,single_wire_cells, tried_cells,cellmap_1,wiremap_1);
      }else{
	std::vector<int> already_removed; 
	std::vector<int> no_need_remove; 
	
	for (int i=0;i!=mergetiling->get_allcell().size();i++){
	  MergeGeomCell *mcell = (MergeGeomCell*)mergetiling->get_allcell().at(i);
	  int index = toymatrix.Get_mcindex(mcell);
	  if (cellmap.find(mcell) == cellmap.end()){
	    already_removed.push_back(index);
	  }else{
	    no_need_remove.push_back(index);
	  }
	}
	
	WireCell2dToy::ToyMatrixKalman toykalman(already_removed, no_need_remove, toymatrix,0,0);

	ncount ++;
	std::cout << "Test: " << cellmap.size() << " " << wiremap.size() << " " << toykalman.Get_numz() << " " << ncount << std::endl;
      }
    }
  }
  
}


WireCell2dToy::ToyMatrixIterate_SingleWire::~ToyMatrixIterate_SingleWire(){
}
