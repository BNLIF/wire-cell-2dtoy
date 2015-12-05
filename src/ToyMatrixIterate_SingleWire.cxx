#include "WireCell2dToy/ToyMatrixIterate_SingleWire.h"

using namespace WireCell;

WireCell2dToy::ToyMatrixIterate_SingleWire::ToyMatrixIterate_SingleWire(WireCell2dToy::ToyMatrix *toybefore, WireCell2dToy::ToyMatrix *toycur, WireCell2dToy::ToyMatrix *toyafter, WireCell2dToy::MergeToyTiling *mergebefore, WireCell2dToy::MergeToyTiling *mergecur, WireCell2dToy::MergeToyTiling *mergeafter, int recon_t, float limit, double penalty, double penalty_ncpt)
  : toymatrix(*toycur)
  , mergetiling(mergecur)
  , limit(limit)
  , penalty_ncpt(penalty_ncpt)
{
  int recon_threshold = recon_t;
  
  // find the penalties for the current set ... 
  GeomCellSelection allmcell_p; 
  if (mergebefore !=0 ){
    allmcell_p = mergebefore->get_allcell();
  }
  GeomCellSelection allmcell_c = mergecur->get_allcell();
  GeomCellSelection allmcell_n; 
  if (mergeafter != 0){
    allmcell_n = mergeafter->get_allcell();
  }

  for (int i=0;i!=allmcell_c.size();i++){
    MergeGeomCell *mcell_c = (MergeGeomCell*)allmcell_c[i];
    int index_c = toycur->Get_mcindex(mcell_c);
    cell_penal[index_c] = 0;

    int flag_before =0;
    if (toybefore!=0){
      if (toybefore->Get_Solve_Flag()!=0){
	 for (int j=0;j!=allmcell_p.size();j++){
	   MergeGeomCell *mcell_p = (MergeGeomCell*)allmcell_p[j];
	   double charge = toybefore->Get_Cell_Charge(mcell_p,1);
	   if ( charge > recon_threshold && mcell_c->Overlap(*mcell_p)){
	     flag_before = 1;
	     break;
	   }
	 }
      }
    }
    if (flag_before == 1){
      cell_penal[index_c] += penalty;
    }
    
    int flag_after = 0;
    if (toyafter !=0){
      if (toyafter->Get_Solve_Flag()!=0){
	for (int j=0;j!=allmcell_n.size();j++){
	  MergeGeomCell *mcell_n = (MergeGeomCell*)allmcell_n[j];
	  double charge = toyafter->Get_Cell_Charge(mcell_n,1);
	  if ( charge > recon_threshold && mcell_c->Overlap(*mcell_n)){
	    flag_after = 1;
	    break;
	  }
	}
      }
    }
    if (flag_after == 1){
      cell_penal[index_c] += penalty;
    }
  }

   // fill in the cell_connect map;
  GeomCellCellsMap& cells_map = mergecur->get_not_compatible_cells_map();
  GeomCellSelection cells = mergecur->get_allcell();
  for (int i=0;i!=cells.size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*)cells.at(i);
    int index = toycur->Get_mcindex(mcell);
    if (cells_map[mcell].size() > 0){
      std::vector<int> cells_indices;
      for (int j=0;j!=cells_map[mcell].size();j++){
	int index_c = toycur->Get_mcindex(cells_map[mcell].at(j));
	cells_indices.push_back(index_c);
      }
      cells_ncpt[index] = cells_indices;
    }
  }
  
  //std::cout << "Check " << cells_ncpt.size() << std::endl;


   // figure out single-wire cell, and no need to remove these
  WireCell::GeomCellSelection all_cells = mergetiling->get_allcell();
  WireCell::GeomCellSelection single_wire_cells = mergetiling->get_single_wire_cells();
  
  wirechargemap = mergetiling->wcmap();

  cells.clear();
  for (int i=0;i!=all_cells.size();i++){
    auto it = find(single_wire_cells.begin(),single_wire_cells.end(),all_cells.at(i));
    if (it == single_wire_cells.end())
      cells.push_back(all_cells.at(i));
  }
  GeomCellMap cellmap = mergetiling->cmap();
  GeomWireMap wiremap = mergetiling->wmap();
  
  //  std::cout << "Test: " << cells.size() << std::endl;
  ncount = 0;
  nlevel = 0;
  if (cells.size() != 0){
    GeomCellSelection tried_cell;
    Iterate(cells,single_wire_cells, tried_cell, cellmap,wiremap);
  }else{
    std::vector<int> already_removed; 
    std::vector<int> no_need_remove; 
	
    for (int i=0;i!=mergetiling->get_allcell().size();i++){
      MergeGeomCell *mcell = (MergeGeomCell*)mergetiling->get_allcell().at(i);
      int index = toymatrix.Get_mcindex(mcell);
      no_need_remove.push_back(index);
    }

    double chi2_p = 0;
    for (int j = 0; j!=already_removed.size(); j++){
      chi2_p += cell_penal[already_removed.at(j)];
    }
    
    if (penalty_ncpt>0){
      //add the penalty due to not compatible cells (ncpt)
      //int flag_test = 0;
      for (auto it = cells_ncpt.begin(); it!= cells_ncpt.end(); it++){
	int index1 = it->first;
	if (find(already_removed.begin(),already_removed.end(),index1) != already_removed.end()) continue;
	std::vector<int> indices = it->second;
	
	for (int j=0;j!=indices.size();j++){
	  int index2 = indices.at(j);
	  if (find(already_removed.begin(),already_removed.end(),index2) == already_removed.end()){
	    chi2_p += penalty_ncpt;
	    // flag_test = 1;
	    // break;
	  }
	}
	// if (flag_test ==1) break;
      }
      // if (flag_test == 1){
      //   chi2_p += penalty_ncpt;
      // }
    }
    
    WireCell2dToy::ToyMatrixKalman toykalman(already_removed, no_need_remove, toymatrix,0,0,chi2_p);
    //std::cout << "Test: " << cellmap.size() << " " << wiremap.size() << " " << toykalman.Get_numz() << " " << ncount << std::endl;
  }



}


WireCell2dToy::ToyMatrixIterate_SingleWire::ToyMatrixIterate_SingleWire(WireCell2dToy::ToyMatrix &toymatrix, WireCell2dToy::MergeToyTiling* mergetiling)
  : toymatrix(toymatrix)
  , mergetiling(mergetiling)
  , penalty_ncpt(0)
{
  for (int i=0;i!=toymatrix.Get_mcindex();i++){
    cell_penal[i] = 0;
  }
  limit = 1e5;

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
  
  //  std::cout << "Test: " << cells.size() << std::endl;
  ncount = 0;
  nlevel = 0;
  if (cells.size() != 0){
    GeomCellSelection tried_cell;
    Iterate(cells,single_wire_cells, tried_cell, cellmap,wiremap);
  }else{
    std::vector<int> already_removed; 
    std::vector<int> no_need_remove; 
	
    for (int i=0;i!=mergetiling->get_allcell().size();i++){
      MergeGeomCell *mcell = (MergeGeomCell*)mergetiling->get_allcell().at(i);
      int index = toymatrix.Get_mcindex(mcell);
      no_need_remove.push_back(index);
    }
	
    WireCell2dToy::ToyMatrixKalman toykalman(already_removed, no_need_remove, toymatrix,0,0);
    //std::cout << "Test: " << cellmap.size() << " " << wiremap.size() << " " << toykalman.Get_numz() << " " << ncount << std::endl;
  }

  
  
}

void WireCell2dToy::ToyMatrixIterate_SingleWire::Iterate(WireCell::GeomCellSelection cells, WireCell::GeomCellSelection single_cells, WireCell::GeomCellSelection tried_cell, WireCell::GeomCellMap cellmap, WireCell::GeomWireMap wiremap){
  if (ncount > limit) return;

  

  nlevel ++;
  for (int i1=0;i1!=cells.size();i1++){
     // stuff that will get passed ... 
    GeomCellSelection cells_1 = cells;
    GeomCellMap cellmap_1 = cellmap;
    GeomWireMap wiremap_1 = wiremap;
    GeomCellSelection single_wire_cells = single_cells;
    GeomCellSelection tried_cells = tried_cell;  
    auto it = find(tried_cells.begin(),tried_cells.end(),cells.at(i1));
    if (it == tried_cells.end()){
      
      // std::cout << "Start: " << nlevel << " " << tried_cells.size() << " " << 
      // 	cells_1.size() << " " << single_wire_cells.size() << " " <<
      // 	cellmap_1.size() << " " << wiremap_1.size() << std::endl;

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
	  if (flag==1){
	    auto it4 = find(single_wire_cells.begin(),single_wire_cells.end(),mcell);
	    if (it4 == single_wire_cells.end())
	      single_wire_cells.push_back(mcell);
	    
	  }
	}
	
	for (int i=0;i!=single_wire_cells.size();i++){
	  MergeGeomCell *mcell = (MergeGeomCell*)single_wire_cells.at(i);
	  
	  auto it6 = find(cells_1.begin(),cells_1.end(),mcell);
	  if (it6 != cells_1.end())
	    cells_1.erase(it6);

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
	  // // put into tried pool
	  // auto it7 = find(tried_cells.begin(),tried_cells.end(),mcell);
	  // if (it7 == tried_cells.end())
	  //   tried_cells.push_back(mcell);
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

	// std::cout << "Middle: " << nlevel << " " << tried_cells.size() << " " << 
	//   cells_1.size() << " " << single_wire_cells.size() << " " <<
	//   cellmap_1.size() << " " << wiremap_1.size() << std::endl;

      }
      // go to next level
      if (cells_1.size() >0){

	// std::cout << "End: " << nlevel << " " << tried_cells.size() << " " << 
	//   cells_1.size() << " " << single_wire_cells.size() << " " <<
	//   cellmap_1.size() << " " << wiremap_1.size() << std::endl;

	Iterate(cells_1,single_wire_cells, tried_cells,cellmap_1,wiremap_1);
      }else{
	std::vector<int> already_removed; 
	std::vector<int> no_need_remove; 
	
	for (int i=0;i!=mergetiling->get_allcell().size();i++){
	  MergeGeomCell *mcell = (MergeGeomCell*)mergetiling->get_allcell().at(i);
	  int index = toymatrix.Get_mcindex(mcell);
	  if (cellmap_1.find(mcell) == cellmap_1.end()){
	    already_removed.push_back(index);
	  }else{
	    no_need_remove.push_back(index);
	  }
	}
	

	double chi2_p = 0;
	for (int j = 0; j!=already_removed.size(); j++){
	  chi2_p += cell_penal[already_removed.at(j)];
	}
	
	if (penalty_ncpt>0){
	  //add the penalty due to not compatible cells (ncpt)
	  //int flag_test = 0;
	  for (auto it = cells_ncpt.begin(); it!= cells_ncpt.end(); it++){
	    int index1 = it->first;
	    if (find(already_removed.begin(),already_removed.end(),index1) != already_removed.end()) continue;
	    std::vector<int> indices = it->second;
	    
	    for (int j=0;j!=indices.size();j++){
	      int index2 = indices.at(j);
	      if (find(already_removed.begin(),already_removed.end(),index2) == already_removed.end()){
		chi2_p += penalty_ncpt;
		// flag_test = 1;
		// break;
	      }
	    }
	    // if (flag_test ==1) break;
	  }
	  // if (flag_test == 1){
	  //   chi2_p += penalty_ncpt;
	  // }
	}
	
	
	WireCell2dToy::ToyMatrixKalman toykalman(already_removed, no_need_remove, toymatrix,0,0,chi2_p);
	if (ncount > limit) return;
	ncount ++;

	//std::cout << "Test: " << cellmap.size() << " " << wiremap.size() << " " << toykalman.Get_numz() << " " << ncount << std::endl;
      }
    }
  }

  nlevel --;
  
}


WireCell2dToy::ToyMatrixIterate_SingleWire::~ToyMatrixIterate_SingleWire(){
}
