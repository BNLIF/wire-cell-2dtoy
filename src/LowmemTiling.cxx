#include "WireCell2dToy/LowmemTiling.h"

using namespace WireCell;

WireCell2dToy::LowmemTiling::LowmemTiling(int time_slice, int nrebin, WireCell::GeomDataSource& gds,WireCell2dToy::WireCellHolder& holder)
  : gds(gds)
  , nrebin(nrebin)
  , time_slice(time_slice)
  , holder(holder)
{
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();
   
}

void WireCell2dToy::LowmemTiling::DivideWires(int wire_limit, int min_wire){
  
  // loop over all the parent wires,  pick up one parent wire
  for (auto it = pwire_wires_map.begin(); it!=pwire_wires_map.end(); it++){
    MergeGeomWire *pwire = (MergeGeomWire*)it->first;
    GeomWireSelection wires = it->second;
    if (!wire_type_map[pwire]) continue; // if the parent wire is bad, no need to do anything ...
    if (wires.size()<=1) continue; // if the number of wires are zero or one, no need to do anything
    // ensure the sorting is working ... 
    MergeGeomWireSet ordered_wire_set;
    int saved_size = 0;
    for (int i=0;i!=wires.size();i++){
      MergeGeomWire *wire = (MergeGeomWire*)wires.at(i);
      ordered_wire_set.insert(wire);
      if (ordered_wire_set.size() == saved_size){
	for (auto it1 = ordered_wire_set.begin(); it1 != ordered_wire_set.end(); it1++){
	  if (wire->get_allwire().front()->index() == (*it1)->get_allwire().front()->index() &&
	      wire->get_allwire().back()->index() == (*it1)->get_allwire().back()->index()){
	    //std::cout << wire << " " << (*it1) << std::endl;
	    replace_wire(wire,*it1);
	    break;
	  }
	}
      }else{
	saved_size = ordered_wire_set.size();
      }
      //      std::cout << wire << " " << wire->get_allwire().front()->index() << " " << wire->get_allwire().back()->index() << std::endl;
    }
    //std::cout << ordered_wire_set.size() << " " << pwire_wires_map[pwire].size() << std::endl;
    if (ordered_wire_set.size()==1) continue; // if there is 
    //std::cout << wires.size() << " " << ordered_wire_set.size() << std::endl;

    
    
    
    // figure out a way to create small wires ...
    int start_wire_index;
    int end_wire_index;
    start_wire_index = (*ordered_wire_set.begin())->get_allwire().front()->index();
    std::vector<std::pair<int,int>> saved_results;
    
    while(start_wire_index <= (*(ordered_wire_set.rbegin()))->get_allwire().back()->index()){
      
      // std::cout << start_wire_index << " " << (*(ordered_wire_set.rbegin()))->get_allwire().back()->index() << std::endl;
      
      int flag = 0;
      // find the first end wire ... 
      for (auto it1 = ordered_wire_set.begin(); it1!= ordered_wire_set.end(); it1++){
	if (start_wire_index >= (*it1)->get_allwire().front()->index() && 
	    start_wire_index <= (*it1)->get_allwire().back()->index()){
	  //std::cout << start_wire_index << " " << (*it1)->get_allwire().front()->index() << " " << (*it1)->get_allwire().back()->index() << std::endl;
	  end_wire_index = (*it1)->get_allwire().back()->index();
	  break;
	}
      }
      
      
      // do the sorting
      std::set<int> sorted_index;
      for (auto it1 = ordered_wire_set.begin(); it1!= ordered_wire_set.end(); it1++){
	int temp_index =  (*it1)->get_allwire().front()->index() - 1;
	if (temp_index >=start_wire_index && temp_index <=end_wire_index){ 
	  end_wire_index = temp_index;
	  if (temp_index - start_wire_index >= min_wire)
	    sorted_index.insert(temp_index);
	}
	temp_index = (*it1)->get_allwire().back()->index();
	if (temp_index >=start_wire_index && temp_index <=end_wire_index){ 
	  end_wire_index = temp_index;
	  if (temp_index - start_wire_index >= min_wire)
	    sorted_index.insert(temp_index);
	}
      }
      

      // only do this when satisfy the condition ... 
      if (sorted_index.size()>0){
	for (auto it2 = sorted_index.begin(); it2!= sorted_index.end(); it2++){
	  end_wire_index = *it2;
	  int count = 0;
	  for (auto it1 = ordered_wire_set.begin(); it1!= ordered_wire_set.end(); it1++){
	    int temp_start_index = (*it1)->get_allwire().front()->index();
	    int temp_end_index = (*it1)->get_allwire().back()->index();
	    if (end_wire_index - temp_start_index >= min_wire && 
		temp_end_index - start_wire_index >= min_wire && 
		(temp_start_index != start_wire_index || 
		 temp_end_index != end_wire_index) &&
		(temp_start_index <= start_wire_index ||
		 temp_end_index >= end_wire_index)
		){
	      count ++;
	    }
	    
	    //    	  std::cout << temp_start_index << " " << temp_end_index << " " << start_wire_index << " " << end_wire_index << " " << count << std::endl;
	  }
	  if (count >0 && count < wire_limit){
	    saved_results.push_back(std::make_pair(start_wire_index,end_wire_index));
	    start_wire_index = end_wire_index +1;
	    flag = 1;
	    break;
	  }
	  //	std::cout << count << std::endl;
	} 
      }

      
      if (flag==0) {
	start_wire_index ++;
	if (start_wire_index < end_wire_index +1)
	  start_wire_index = end_wire_index + 1;
	end_wire_index = start_wire_index;
      }
      //std::cout << flag << " " << start_wire_index << " " << end_wire_index << " " << (*(ordered_wire_set.rbegin()))->get_allwire().back()->index() << std::endl;
    }
   
    if (saved_results.size()==0) continue;
    
    // create new cells, link them to the old cells
    // replace the old cell
    // replace the old wires 
    

    GeomCellSelection old_cells;
    GeomWireSelection old_wires;
    for (int i=0;i!=saved_results.size();i++){
      // find the old wire
      MergeGeomWire *old_wire;
      for (auto it1 = ordered_wire_set.begin(); it1!= ordered_wire_set.end(); it1++){
	if ( (*it1)->get_allwire().front()->index()<= saved_results.at(i).first && 
	     (*it1)->get_allwire().back()->index() >= saved_results.at(i).second && 
	     ( (*it1)->get_allwire().front()->index()!= saved_results.at(i).first ||
	       (*it1)->get_allwire().back()->index() != saved_results.at(i).second )){
	  old_wire = *it1;
	  
	  if (find(old_wires.begin(),old_wires.end(),old_wire) != old_wires.end()) continue;
	  old_wires.push_back(old_wire);
	  // create new wires
	  GeomWireSelection new_wires;
	  if (old_wire->get_allwire().front()->index() < saved_results.at(i).first){
	    GeomWireSelection nwires;
	     for (int j=old_wire->get_allwire().front()->index(); j!= saved_results.at(i).first; j++){
	       //std::cout << j << std::endl;
	       nwires.push_back(gds.by_planeindex(old_wire->get_allwire().front()->plane(),j));
	     }
	     MergeGeomWire *new_wire = new MergeGeomWire(0,nwires);
	     new_wires.push_back(new_wire);
	   }
	  if (old_wire->get_allwire().back()->index() > saved_results.at(i).second){
	    GeomWireSelection nwires;
	    for (int j=saved_results.at(i).second+1; j!= old_wire->get_allwire().back()->index() + 1; j++){
	      nwires.push_back(gds.by_planeindex(old_wire->get_allwire().front()->plane(),j));
	    }
	    MergeGeomWire *new_wire = new MergeGeomWire(0,nwires);
	    new_wires.push_back(new_wire);
	  }
	  {
	    GeomWireSelection nwires;
	    for (int j=saved_results.at(i).first; j!= saved_results.at(i).second+1; j++){
	      nwires.push_back(gds.by_planeindex(old_wire->get_allwire().front()->plane(),j));
	    }
	    MergeGeomWire *new_wire = new MergeGeomWire(0,nwires);
	    new_wires.push_back(new_wire);
	  }
	  
	  //std::cout << new_wires.size() << std::endl;

	  
	  for (auto it2 = wire_cells_map[old_wire].begin(); it2!= wire_cells_map[old_wire].end(); it2++){
	    // find the old cell
	    SlimMergeGeomCell *old_cell = (SlimMergeGeomCell*)(*it2);
	    // std::cout << "A: " <<  old_wire << " " << wire_cells_map[old_wire].size() << " " << old_cell << " " << create_single_cells(old_cell).size() << std::endl;
	    
	    // create new cells
	    MergeGeomWire *uwire=0;
	    MergeGeomWire *vwire=0;
	    MergeGeomWire *wwire=0;
	    GeomWireSelection temp_wires =  cell_wires_map[old_cell];
	    for (int j=0;j!=temp_wires.size();j++){
	      if (temp_wires.at(j) != old_wire){
		if ( ((MergeGeomWire*)temp_wires.at(j))->get_allwire().front()->plane() == WirePlaneType_t(0)){
		  uwire = (MergeGeomWire*)temp_wires.at(j);
		}else if ( ((MergeGeomWire*)temp_wires.at(j))->get_allwire().front()->plane() == WirePlaneType_t(1)){
		  vwire = (MergeGeomWire*)temp_wires.at(j);
		}else{
		  wwire = (MergeGeomWire*)temp_wires.at(j);
		}
	      }
	    }
	    int plane_flag;
	    if (uwire ==0 ){
	      plane_flag = 0;
	    }else if (vwire==0){
	      plane_flag = 1;
	    }else{
	      plane_flag = 2;
	    }
	    for (int j=0;j!=new_wires.size();j++){
	      SlimMergeGeomCell *new_cell;
	      if (plane_flag==0){
		uwire = (MergeGeomWire*)new_wires.at(j);
		new_cell = create_slim_merge_cell((MergeGeomWire*)new_wires.at(j),vwire,wwire);
	      }else if (plane_flag==1){
		vwire = (MergeGeomWire*)new_wires.at(j);
		new_cell = create_slim_merge_cell(uwire,(MergeGeomWire*)new_wires.at(j),wwire);
	      }else{
		wwire = (MergeGeomWire*)new_wires.at(j);
		new_cell = create_slim_merge_cell(uwire,vwire,(MergeGeomWire*)new_wires.at(j));
	      }
	      //std::cout << "B: " << plane_flag << " " << create_single_cells(new_cell).size() << std::endl;
	      // add in new cells and its associated wires into maps
	      MergeGeomWire *pwire = (MergeGeomWire*)wire_pwire_map[old_wire];
	      
	      // create new mwires 
	      GeomWireSelection uwires = new_cell->get_uwires();
	      GeomWireSelection vwires = new_cell->get_vwires();
	      GeomWireSelection wwires = new_cell->get_wwires();
	      MergeGeomWire *mwire_u = new MergeGeomWire(0,uwires);
	      MergeGeomWire *mwire_v = new MergeGeomWire(0,vwires);
	      MergeGeomWire *mwire_w = new MergeGeomWire(0,wwires);
	    
	      // create the map  
	      // cell to wires
	      GeomWireSelection wires;
	      wires.push_back(mwire_u);
	      wires.push_back(mwire_v);
	      wires.push_back(mwire_w);
	      cell_wires_map[new_cell] = wires;
	    
	      // wire to cells
	      GeomCellSelection cells;
	      cells.push_back(new_cell);
	      wire_cells_map[mwire_u] = cells;
	      wire_cells_map[mwire_v] = cells;
	      wire_cells_map[mwire_w] = cells;
	    
	      // wire to parent wire
	      // wire types
	      // parent wire to wires
	      if (plane_flag==0){
		wire_pwire_map[mwire_u] = pwire;
		wire_pwire_map[mwire_v] = wire_pwire_map[vwire];
		wire_pwire_map[mwire_w] = wire_pwire_map[wwire];

		pwire_wires_map[pwire].push_back(mwire_u);
		pwire_wires_map[wire_pwire_map[vwire]].push_back(mwire_v);
		pwire_wires_map[wire_pwire_map[wwire]].push_back(mwire_w);

		wire_type_map[mwire_u] = true;
		wire_type_map[mwire_v] = wire_type_map[vwire];     
		wire_type_map[mwire_w] = wire_type_map[wwire];
	      }else if (plane_flag==1){
		wire_pwire_map[mwire_u] = wire_pwire_map[uwire];
		wire_pwire_map[mwire_v] = pwire;
		wire_pwire_map[mwire_w] = wire_pwire_map[wwire];

		pwire_wires_map[wire_pwire_map[uwire]].push_back(mwire_u);
		pwire_wires_map[pwire].push_back(mwire_v);
		pwire_wires_map[wire_pwire_map[wwire]].push_back(mwire_w);

		wire_type_map[mwire_u] = wire_type_map[uwire];
		wire_type_map[mwire_v] = true;     
		wire_type_map[mwire_w] = wire_type_map[wwire];
	      }else{
		wire_pwire_map[mwire_u] = wire_pwire_map[uwire];
		wire_pwire_map[mwire_v] = wire_pwire_map[vwire];
		wire_pwire_map[mwire_w] = pwire;

		pwire_wires_map[wire_pwire_map[uwire]].push_back(mwire_u);
		pwire_wires_map[wire_pwire_map[vwire]].push_back(mwire_v);
		pwire_wires_map[pwire].push_back(mwire_w);

		wire_type_map[mwire_u] = wire_type_map[uwire];
		wire_type_map[mwire_v] = wire_type_map[vwire];     
		wire_type_map[mwire_w] = true;
	      }
	    
	    
	    }
	    // std::cout << uwire << " " << vwire << " " << wwire << std::endl;
	    
	    old_cells.push_back(old_cell);
	  }

	  // remove new wires ...
	  for (int j=0;j!=new_wires.size();j++){
	    delete new_wires.at(j);
	  }
	}
      }
    }

    // remove old cell
    for (int j=0;j!=old_cells.size();j++){
      remove_cell((SlimMergeGeomCell*)old_cells.at(j));
    }
    // remove old wires
    for (int j=0;j!=old_wires.size();j++){
      remove_wire((MergeGeomWire*)old_wires.at(j));
    }

    // std::cout << "abc " << " " << start_wire_index << " " << end_wire_index << std::endl;
    // for (auto it1 = ordered_wire_set.begin(); it1!= ordered_wire_set.end(); it1++){
    //   std::cout << (*it1)->get_allwire().front()->index() << " " << (*it1)->get_allwire().back()->index() << " " << (*it1)->get_allwire().front()->plane() << " " << wires.size() << std::endl;
    // }
   
  }

 
  
}

void WireCell2dToy::LowmemTiling::MergeWires(){
  // basically the original mergetoytiling algorithm ... 
  
  // loop over parent wire maps, 
  for (auto it = pwire_wires_map.begin(); it!= pwire_wires_map.end(); it++){
    MergeGeomWire *pwire = (MergeGeomWire*)it->first;
    GeomWireSelection allwires = it->second;
    if (wire_type_map[pwire]){
      //only do merge if the parent wire is good ...
      //while(further_mergewire(pwire_wires_map[pwire]));
      GeomWireSelection wires = pwire_wires_map[pwire];

      // std::cout << wire_cells_map.size() << std::endl;
      // std::cout << "a: " << pwire_wires_map[pwire].size() << " " << wires.size() << " " << further_mergewire(wires)  << std::endl;
      
      // wires = pwire_wires_map[pwire];
      // std::cout << "b: " << pwire_wires_map[pwire].size() << " " << wires.size() << " " << further_mergewire(wires)  << std::endl;
      // std::cout << wire_cells_map.size() << std::endl;
     
      while(further_mergewire(wires)){
	//	std::cout << pwire_wires_map.size() << " " << wires.size() << " " << wire_cells_map.size() << std::endl;
  	//pwire_wires_map[pwire] = wires;
       	wires = pwire_wires_map[pwire];
      }
    }
  }

  calculate_merged_wire_charge();

  for (auto it=cell_wires_map.begin(); it!= cell_wires_map.end(); it++){
    SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)it->first;
    mcell->OrderWires();

    // for (size_t i=0;i!=mcell->get_uwires().size();i++){
    //   std::cout << i << " " << mcell->get_uwires().at(i)->ident() << " " << mcell->get_uwires().at(i)->index() << std::endl;
    // }
  }
  
}

void WireCell2dToy::LowmemTiling::calculate_merged_wire_charge(){
  // set relative uncertainties to zero
  float relative_err_ind = 0.0; // 15% relative uncertainties for induction plane 
  float relative_err_col = 0.0; // 5% relative uncertainties for collection plane
  
  for (auto it = wire_cells_map.begin(); it!= wire_cells_map.end(); it++){
    MergeGeomWire *mwire = (MergeGeomWire*)it->first;
    if (wire_type_map[mwire]){
      // calculate merged good wire charge
      wirechargemap[mwire] = 0;
      wirecharge_errmap[mwire] = 0;
      GeomWireSelection wires = mwire->get_allwire();
      for (size_t i=0; i!=wires.size();i++){
	const GeomWire* swire = wires.at(i);
	// if (wirechargemap.find(swire)==wirechargemap.end() ||
	//     wirecharge_errmap.find(swire)==wirecharge_errmap.end()){
	//   std::cout << "Wrong! " << std::endl;
	// }
	wirechargemap[mwire] += wirechargemap[swire];
	wirecharge_errmap[mwire] += pow(wirecharge_errmap[swire],2);
      }
      if (wires.at(0)->plane()==WirePlaneType_t(2)){
	wirecharge_errmap[mwire] += pow(wirechargemap[mwire] * relative_err_col,2);
      }else{
	wirecharge_errmap[mwire] += pow(wirechargemap[mwire] * relative_err_ind,2);
      }
      wirecharge_errmap[mwire] = sqrt(wirecharge_errmap[mwire]);

      //std::cout << mwire << " " << wires.size()  << " " << wires.at(0)->plane() << " " << wirechargemap[mwire] << " " << wirecharge_errmap[mwire] << std::endl;
    }else{
      // calculate bad wire charge and uncertainties, charge = 0, error = -1
      wirechargemap[mwire] = 0;
      wirecharge_errmap[mwire] = -1;
    }
  }

}

int WireCell2dToy::LowmemTiling::further_mergewire(GeomWireSelection& allwire){
  
  WireCell::GeomWireSelection tempwire = allwire;
  allwire.clear();
  
  for (int i=0;i!=tempwire.size();i++){
    MergeGeomWire *wire = (MergeGeomWire*)tempwire[i];
    int flag=0;
    for (int k=0;k!=allwire.size();k++){
      if (((MergeGeomWire*)allwire[k])->AddWire(*wire)){
       	// now need to do something ... 
	replace_wire(wire,(MergeGeomWire*)allwire[k]);
	flag = 1;
    	break;
      }
    }
      
    if(flag==0){
      allwire.push_back(wire);
    }
  }
  
  int diff = tempwire.size() - allwire.size();  
  return diff;
}


bool WireCell2dToy::LowmemTiling::replace_wire(WireCell::MergeGeomWire *old_wire, WireCell::MergeGeomWire *wire){
  // replace the old_wire by the new wire
  // assume the new wire contains the old wire ... 
  // deal with the wire type map  // no point to do this if the wire is bad ... 
  if (wire_type_map.find(old_wire) != wire_type_map.end()){
    wire_type_map.erase(old_wire);
  }
  wire_type_map[wire] = true;

  // deal with the wire-wire maps
  if (wire_pwire_map.find(old_wire) != wire_pwire_map.end()){
    MergeGeomWire *pwire = (MergeGeomWire*)wire_pwire_map[old_wire];
  
    GeomWireSelection& wires = pwire_wires_map[pwire];
    auto it = find(wires.begin(),wires.end(),old_wire);
    if (it!=wires.end())
      wires.erase(it); // remove the old one
    if (find(wires.begin(),wires.end(),wire) == wires.end()){
      wires.push_back(wire); // push in the new one
    }
  
    wire_pwire_map.erase(old_wire);
    wire_pwire_map[wire] = pwire; // get the new one registered ...
  }
  

  // deal with the cell-wire maps
  if (wire_cells_map.find(old_wire)!=wire_cells_map.end()){
    // for each of these cell, remove the old wire and add in the new wire
    GeomCellSelection& cells = wire_cells_map[old_wire];
    
    if (wire_cells_map.find(wire)==wire_cells_map.end()){
      GeomCellSelection cells1;
      wire_cells_map[wire] = cells1;
    }

    for (int i=0;i!=cells.size();i++){
      const GeomCell *cell = cells.at(i);
      // replace the old wire with the new wire ... 
      GeomWireSelection& wires = cell_wires_map[cell];
      auto it = find(wires.begin(),wires.end(),old_wire);
      if (it!=wires.end())
	wires.erase(it);
      if (find(wires.begin(),wires.end(),wire)==wires.end())
	wires.push_back(wire);

      // add cells to this new wire cell ... 
      if (find(wire_cells_map[wire].begin(), wire_cells_map[wire].end(), cell) == wire_cells_map[wire].end()){
      	wire_cells_map[wire].push_back(cell);
      }
    }

    wire_cells_map.erase(old_wire);
  }
  
  
  delete old_wire;
}

bool WireCell2dToy::LowmemTiling::remove_wire_clear(MergeGeomWire *wire){
  // do the cell-wire maps
  if (wire_cells_map.find(wire)!=wire_cells_map.end()){
    wire_cells_map.erase(wire);
  }
  
  //do the wire-wire maps
  if (wire_pwire_map.find(wire) != wire_pwire_map.end()){
    MergeGeomWire *pwire = (MergeGeomWire*)wire_pwire_map[wire];
    GeomWireSelection& wires = pwire_wires_map[pwire];
    auto it = find(wires.begin(),wires.end(),wire);
    if (it!=wires.end())
      wires.erase(it);
    wire_pwire_map.erase(wire);
  }
  
  //do wire type map 
  if (wire_type_map.find(wire) != wire_type_map.end()){
    wire_type_map.erase(wire);
    delete wire;
  }
  return true;
}



bool WireCell2dToy::LowmemTiling::remove_wire(MergeGeomWire *wire){
  // do the cell-wire maps
  if (wire_cells_map.find(wire)!=wire_cells_map.end()){
    GeomCellSelection& cells = wire_cells_map[wire];
    if(cells.size()!=0){
      return false;
    }else{
      wire_cells_map.erase(wire);
    }
  }
  
  //do the wire-wire maps
  if (wire_pwire_map.find(wire) != wire_pwire_map.end()){
    MergeGeomWire *pwire = (MergeGeomWire*)wire_pwire_map[wire];
  
    GeomWireSelection& wires = pwire_wires_map[pwire];
    auto it = find(wires.begin(),wires.end(),wire);
    if (it!=wires.end())
      wires.erase(it);
  
    wire_pwire_map.erase(wire);
  }
  

  //do wire type map 
  if (wire_type_map.find(wire) != wire_type_map.end()){
    wire_type_map.erase(wire);
    delete wire;
  }
 

  return true;
}

void WireCell2dToy::LowmemTiling::Erase_Cell(SlimMergeGeomCell *cell){

  {
    // remove it from three_good_wire_cells
    auto it = find(three_good_wire_cells.begin(), three_good_wire_cells.end(), cell);
    if (it!= three_good_wire_cells.end())
      three_good_wire_cells.erase(it);
  }
  
  {
    // remove it from two_good_wire_cells
    auto it = find(two_good_wire_cells.begin(), two_good_wire_cells.end(), cell);
    if (it!= two_good_wire_cells.end())
      two_good_wire_cells.erase(it);
  }
  // remove it from one_good_wire_cells
  {
    // remove it from two_good_wire_cells
    auto it = find(one_good_wire_cells.begin(), one_good_wire_cells.end(), cell);
    if (it!= one_good_wire_cells.end())
      one_good_wire_cells.erase(it);
  }
  
  // remove it from map ...
  
  if (cell_wires_map.find(cell) != cell_wires_map.end()){
    // find all the wires connect to the cells
    GeomWireSelection& wires = cell_wires_map[cell];
    for (int i = 0; i!=wires.size();i++){
      MergeGeomWire *wire = (MergeGeomWire*)wires.at(i);
      // remove the cell from the wire map
      GeomCellSelection& cells = wire_cells_map[wire];
      auto it = find(cells.begin(),cells.end(),cell);
      if (it!=cells.end())
	cells.erase(it);

      // if the wire does not connect anything, remove it ... 
      if(cells.size()==0){
	remove_wire(wire); 
      }

      
      // remove the cell 
      cell_wires_map.erase(cell);
    }
  }
   
  // remove it from holder (no need ...)
  
}

bool WireCell2dToy::LowmemTiling::remove_cell(SlimMergeGeomCell *cell){
  if (cell_wires_map.find(cell) != cell_wires_map.end()){
    // find all the wires connect to the cells
    GeomWireSelection& wires = cell_wires_map[cell];
    for (int i = 0; i!=wires.size();i++){
      MergeGeomWire *wire = (MergeGeomWire*)wires.at(i);
      // remove the cell from the wire map
      GeomCellSelection& cells = wire_cells_map[wire];
      auto it = find(cells.begin(),cells.end(),cell);
      if (it!=cells.end())
	cells.erase(it);

      // if the wire does not connect anything, remove it ... 
      if(cells.size()==0){
	remove_wire(wire); 
      }

      // remove the cell 
      cell_wires_map.erase(cell);
    }
    delete cell;
  }
  
  return false;
}


// bool WireCell2dToy::LowmemTiling::check_crossing(float dis_u, float u_pitch, float dis_v, float v_pitch,
// 						 float min_w, float max_w){
//   std::vector<Vector> psave(5);
  
//   for (int i=0;i!=5;i++){
//     bool flag;
//     if (i==0){
//       flag = gds.crossing_point(dis_u-u_pitch/2.,dis_v-v_pitch/2.,kUwire,kVwire, psave[0]);
//     }else if (i==1){
//       flag = gds.crossing_point(dis_u+u_pitch/2.,dis_v+v_pitch/2.,kUwire,kVwire, psave[0]);
//     }else if (i==2){
//       flag = gds.crossing_point(dis_u,dis_v,kUwire,kVwire, psave[0]);
//     }else if (i==3){
//       flag = gds.crossing_point(dis_u+u_pitch/2.,dis_v-v_pitch/2.,kUwire,kVwire, psave[0]);
//     }else if (i==4){
//       flag = gds.crossing_point(dis_u-u_pitch/2.,dis_v+v_pitch/2.,kUwire,kVwire, psave[0]);
//     }
 
//     if (flag){
//       if (psave[0].z > max_w && psave[0].z < max_w)
// 	return true;
//     }
//   }
  
//   return false;
// }

bool WireCell2dToy::LowmemTiling::check_crossing(const WireCell::GeomWire* wire1, const WireCell::GeomWire* wire2, float pitch1, float pitch2, WirePlaneType_t plane, float min, float max, float tolerance){
  
  float dis1 = gds.wire_dist(*wire1);
  float dis2 = gds.wire_dist(*wire2);  
  
  bool flag;
  float dis;
  std::vector<Vector> puv_save(5);
  
  flag = gds.crossing_point(dis1-pitch1/4.,dis2-pitch2/4.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1-pitch1/4.,dis2+pitch2/4.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1+pitch1/4.,dis2-pitch2/4.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1+pitch1/4.,dis2+pitch2/4.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1-pitch1/2.,dis2-pitch2/2.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1-pitch1/2.,dis2+pitch2/2.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1+pitch1/2.,dis2-pitch2/2.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1+pitch1/2.,dis2+pitch2/2.,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  flag = gds.crossing_point(dis1,dis2,wire1->plane(),wire2->plane(), puv_save[0]);
  dis = gds.wire_dist(puv_save[0],plane);
  if (flag){
    if (dis >= min+tolerance && dis <= max-tolerance)
      return true;
  }
  
  return false;
}

bool WireCell2dToy::LowmemTiling::test_point(WireCell::PointVector& pcell, float dis_u, float dis_v,
					     float bmin_w, float bmax_w, float u_pitch, float v_pitch,
					     float& dis1){
  bool flag_qx = false;
  bool flag[10];
  std::vector<Vector> puv_save(9);
  flag[0] = gds.crossing_point(dis_u-u_pitch/2.,dis_v-v_pitch/2.,kUwire,kVwire, puv_save[0]);
  flag[1] = gds.crossing_point(dis_u+u_pitch/2.,dis_v-v_pitch/2.,kUwire,kVwire, puv_save[1]);
  flag[2] = gds.crossing_point(dis_u-u_pitch/2.,dis_v+v_pitch/2.,kUwire,kVwire, puv_save[2]);
  flag[3] = gds.crossing_point(dis_u+u_pitch/2.,dis_v+v_pitch/2.,kUwire,kVwire, puv_save[3]);
  flag[4] = gds.crossing_point(dis_u,dis_v,kUwire,kVwire, puv_save[4]);
  // flag[5] = gds.crossing_point(dis_u-u_pitch/4.,dis_v-v_pitch/4.,kUwire,kVwire, puv_save[5]);
  // flag[6] = gds.crossing_point(dis_u+u_pitch/4.,dis_v-v_pitch/4.,kUwire,kVwire, puv_save[6]);
  // flag[7] = gds.crossing_point(dis_u-u_pitch/4.,dis_v+v_pitch/4.,kUwire,kVwire, puv_save[7]);
  // flag[8] = gds.crossing_point(dis_u+u_pitch/4.,dis_v+v_pitch/4.,kUwire,kVwire, puv_save[8]);

  

  float dis[9];
  for (int i=0;i!=5;i++){
    dis[i] = gds.wire_dist(puv_save[i],kYwire);
    if (flag[i] && dis[i] > bmin_w && dis[i] < bmax_w){
      pcell.push_back(puv_save[i]);
      flag_qx = true;
    }
  }
  dis1 = dis[4];
  
  return flag_qx;

}

void WireCell2dToy::LowmemTiling::test_cross(WireCell::PointVector& pcell, float dis, WireCell::WirePlaneType_t plane, float low_limit, float high_limit, float pitch, float pitch1, float w_pitch, const WireCell::GeomWire* wire1, const WireCell::GeomWire *wire2, int dir, float bmin_w, float bmax_w){
  
  Vector pcross;
  float start_pos;
  float end_pos;
  float step_size;
  
  WireCell::WirePlaneType_t plane1;
  if (plane == kUwire){
    plane1 = kVwire;
  }else if (plane == kVwire){
    plane1 = kUwire;
  }
  
  if (dir==1){
    start_pos = low_limit - 2 * w_pitch;
    end_pos = high_limit + 2.5 * w_pitch;
    step_size = w_pitch;
  }else{
    start_pos = high_limit + 2 * w_pitch;
    end_pos = low_limit - 2.5 * w_pitch;
    step_size = -1 * w_pitch;
  }
  
  bool flag[10];
  float dis2[10];
  bool flag_qx;
  std::vector<Vector> puv_save(10);
  
  float step = 0;
  
  while ( fabs(step) < fabs(end_pos - start_pos)){
    //std::cout << step << " " << end_pos - start_pos <<  " " << step_size << std::endl;
    float pos = start_pos + step;
    
    if (gds.crossing_point(pos,dis,kYwire,plane,pcross)){
      const GeomWire *uwire_3 = gds.closest(pcross,plane1);
      if (uwire_3->index()>=wire1->index() && uwire_3->index()<= wire2->index()){
	float dis1 = gds.wire_dist(*uwire_3);
	flag[0] = gds.crossing_point(dis1-pitch1/2.,dis-pitch/2.,plane1,plane, puv_save[0]);
	flag[1] = gds.crossing_point(dis1+pitch1/2.,dis-pitch/2.,plane1,plane, puv_save[1]);
	flag[2] = gds.crossing_point(dis1-pitch1/2.,dis+pitch/2.,plane1,plane, puv_save[2]);
	flag[3] = gds.crossing_point(dis1+pitch1/2.,dis+pitch/2.,plane1,plane, puv_save[3]);
	flag[4] = gds.crossing_point(dis1,dis,plane1,plane, puv_save[4]);
	// flag[5] = gds.crossing_point(dis1-pitch1/4.,dis-pitch/4.,plane1,plane, puv_save[5]);
	// flag[6] = gds.crossing_point(dis1+pitch1/4.,dis-pitch/4.,plane1,plane, puv_save[6]);
	// flag[7] = gds.crossing_point(dis1-pitch1/4.,dis+pitch/4.,plane1,plane, puv_save[7]);
	// flag[8] = gds.crossing_point(dis1+pitch1/4.,dis+pitch/4.,plane1,plane, puv_save[8]);
	flag_qx = false;
	
	for (int i=0;i!=5;i++){
	  dis2[i] = gds.wire_dist(puv_save[i],kYwire);
	  if (flag[i] && dis2[i] > bmin_w && dis2[i] < bmax_w){
	    pcell.push_back(puv_save[i]);
	    flag_qx = true;
	  }
	}
	
	if (flag_qx) break;
      }
    }
    
    step += step_size;
  }

}



// void WireCell2dToy::LowmemTiling::test_crossing(PointVector& pcell, float dis_u, float dis_v, float bmin_w, float bmax_w, float u_pitch, float v_pitch, float w_pitch, const GeomWire* uwire_1, const GeomWire* uwire_2, const GeomWire *vwire_1, const GeomWire *vwire_2, int dir_u, int dir_v){
  

// 	//flag_qx = false;
// 	for (int i=0;i!=5;i++){
// 	  dis[i] = gds.wire_dist(puv_save[i],kYwire);
// 	  if (flag[i] && dis[i] > bmin_w && dis[i] < bmax_w){
// 	    pcell.push_back(puv_save[i]);
// 	    //  flag_qx = true;
// 	  }
// 	}
	
// 	// if (flag_qx){
// 	//   break;
// 	// }
//       }
//     }
      
//     if (gds.crossing_point(bmax_w,dis_v,kYwire,kVwire,pcross[0])){
//       // quick calculation
//       const GeomWire *uwire_4 = gds.closest(pcross[0],WirePlaneType_t(0));
//       for (int k=0;k!=5;k++){
// 	int index = uwire_4->index()+(-2+k);
// 	if (index < uwire_1->index() || 
// 	    index > uwire_2->index()) continue;
// 	const GeomWire *uwire_3 = gds.by_planeindex(WirePlaneType_t(0), index);
// 	dis_u1 = gds.wire_dist(*uwire_3);
// 	flag[0] = gds.crossing_point(dis_u1-u_pitch/2.,dis_v-v_pitch/2.,kUwire,kVwire, puv_save[0]);
// 	flag[1] = gds.crossing_point(dis_u1+u_pitch/2.,dis_v-v_pitch/2.,kUwire,kVwire, puv_save[1]);
// 	flag[2] = gds.crossing_point(dis_u1-u_pitch/2.,dis_v+v_pitch/2.,kUwire,kVwire, puv_save[2]);
// 	flag[3] = gds.crossing_point(dis_u1+u_pitch/2.,dis_v+v_pitch/2.,kUwire,kVwire, puv_save[3]);
// 	flag[4] = gds.crossing_point(dis_u1,dis_v,kUwire,kVwire, puv_save[4]);
// 	// flag[5] = gds.crossing_point(dis_u1-u_pitch/4.,dis_v-v_pitch/4.,kUwire,kVwire, puv_save[5]);
// 	// flag[6] = gds.crossing_point(dis_u1+u_pitch/4.,dis_v-v_pitch/4.,kUwire,kVwire, puv_save[6]);
// 	// flag[7] = gds.crossing_point(dis_u1-u_pitch/4.,dis_v+v_pitch/4.,kUwire,kVwire, puv_save[7]);
// 	// flag[8] = gds.crossing_point(dis_u1+u_pitch/4.,dis_v+v_pitch/4.,kUwire,kVwire, puv_save[8]);
// 	//  flag_qx = false;
// 	for (int i=0;i!=5;i++){
// 	  dis[i] = gds.wire_dist(puv_save[i],kYwire);
// 	  if (flag[i] && dis[i] > bmin_w && dis[i] < bmax_w){
// 	    pcell.push_back(puv_save[i]);
// 	    // flag_qx = true;
// 	  }
// 	}
	
// 	// if (flag_qx){
// 	//   break;
// 	// }
//       }
//     }
    
//     // look at the V wire
//     if (gds.crossing_point(bmin_w,dis_u,kYwire,kUwire,pcross[0])){
//       const GeomWire *vwire_4 = gds.closest(pcross[0],WirePlaneType_t(1));
//       for (int k=0;k!=5;k++){
// 	int index = vwire_4->index()+(-2+k);
// 	if (index < vwire_1->index() || 
// 	    index > vwire_2->index()) continue;
// 	const GeomWire *vwire_3 = gds.by_planeindex(WirePlaneType_t(1), index);
// 	dis_v1 = gds.wire_dist(*vwire_3);
// 	flag[0] = gds.crossing_point(dis_u-u_pitch/2.,dis_v1-v_pitch/2.,kUwire,kVwire, puv_save[0]);
// 	flag[1] = gds.crossing_point(dis_u+u_pitch/2.,dis_v1-v_pitch/2.,kUwire,kVwire, puv_save[1]);
// 	flag[2] = gds.crossing_point(dis_u-u_pitch/2.,dis_v1+v_pitch/2.,kUwire,kVwire, puv_save[2]);
// 	flag[3] = gds.crossing_point(dis_u+u_pitch/2.,dis_v1+v_pitch/2.,kUwire,kVwire, puv_save[3]);
// 	flag[4] = gds.crossing_point(dis_u,dis_v1,kUwire,kVwire, puv_save[4]);
// 	// flag[5] = gds.crossing_point(dis_u-u_pitch/4.,dis_v1-v_pitch/4.,kUwire,kVwire, puv_save[5]);
// 	// flag[6] = gds.crossing_point(dis_u+u_pitch/4.,dis_v1-v_pitch/4.,kUwire,kVwire, puv_save[6]);
// 	// flag[7] = gds.crossing_point(dis_u-u_pitch/4.,dis_v1+v_pitch/4.,kUwire,kVwire, puv_save[7]);
// 	// flag[8] = gds.crossing_point(dis_u+u_pitch/4.,dis_v1+v_pitch/4.,kUwire,kVwire, puv_save[8]);
// 	//flag_qx = false;
// 	for (int i=0;i!=5;i++){
// 	  dis[i] = gds.wire_dist(puv_save[i],kYwire);
// 	  if (flag[i] && dis[i] > bmin_w && dis[i] < bmax_w ){
// 	    pcell.push_back(puv_save[i]);
// 	    //  flag_qx = true;
// 	  }
// 	}
	
// 	// if (flag_qx){
// 	//   break;
// 	// }
//       }
//     }
    
//     if (gds.crossing_point(bmax_w,dis_u,kYwire,kUwire,pcross[0])){
//       const GeomWire *vwire_4 = gds.closest(pcross[0],WirePlaneType_t(1));
//       for (int k=0;k!=5;k++){
// 	int index = vwire_4->index()+(-2+k);
// 	if (index < vwire_1->index() || 
// 	    index > vwire_2->index()) continue;
// 	const GeomWire *vwire_3 = gds.by_planeindex(WirePlaneType_t(1), index);
// 	dis_v1 = gds.wire_dist(*vwire_3);
// 	flag[0] = gds.crossing_point(dis_u-u_pitch/2.,dis_v1-v_pitch/2.,kUwire,kVwire, puv_save[0]);
// 	flag[1] = gds.crossing_point(dis_u+u_pitch/2.,dis_v1-v_pitch/2.,kUwire,kVwire, puv_save[1]);
// 	flag[2] = gds.crossing_point(dis_u-u_pitch/2.,dis_v1+v_pitch/2.,kUwire,kVwire, puv_save[2]);
// 	flag[3] = gds.crossing_point(dis_u+u_pitch/2.,dis_v1+v_pitch/2.,kUwire,kVwire, puv_save[3]);
// 	flag[4] = gds.crossing_point(dis_u,dis_v1,kUwire,kVwire, puv_save[4]);
// 	// flag[5] = gds.crossing_point(dis_u-u_pitch/4.,dis_v1-v_pitch/4.,kUwire,kVwire, puv_save[5]);
// 	// flag[6] = gds.crossing_point(dis_u+u_pitch/4.,dis_v1-v_pitch/4.,kUwire,kVwire, puv_save[6]);
// 	// flag[7] = gds.crossing_point(dis_u-u_pitch/4.,dis_v1+v_pitch/4.,kUwire,kVwire, puv_save[7]);
// 	// flag[8] = gds.crossing_point(dis_u+u_pitch/4.,dis_v1+v_pitch/4.,kUwire,kVwire, puv_save[8]);
// 	//flag_qx = false;
// 	for (int i=0;i!=5;i++){
// 	  dis[i] = gds.wire_dist(puv_save[i],kYwire);
// 	  if (flag[i] && dis[i] > bmin_w && dis[i] < bmax_w){
// 	    pcell.push_back(puv_save[i]);
// 	    //flag_qx = true;
// 	  }
// 	}
	
	
// 	// if (flag_qx){
// 	//   break;
// 	// }
//       }
//     }
//   }
  
  
// }


WireCell::SlimMergeGeomCell* WireCell2dToy::LowmemTiling::create_slim_merge_cell(WireCell::MergeGeomWire *uwire, WireCell::MergeGeomWire *vwire, WireCell::MergeGeomWire *wwire){
  float u_pitch = gds.pitch(kUwire);
  float v_pitch = gds.pitch(kVwire);
  float w_pitch = gds.pitch(kYwire);
  
  float tolerance = 0.1 * units::mm / 2.;
  
  // find U plane wires
  float dis_u[3];
  const GeomWire *uwire_1 = uwire->get_allwire().front();
  const GeomWire *uwire_2 = uwire->get_allwire().back();
  dis_u[0] = gds.wire_dist(*uwire_1);
  dis_u[1] = gds.wire_dist(*uwire_2);
  // float bmin_u = dis_u[0] - u_pitch/2. - tolerance;
  // float bmax_u = dis_u[1] + u_pitch/2. + tolerance;

  // find V plane wires
  float dis_v[3];
  const GeomWire *vwire_1 = vwire->get_allwire().front();
  const GeomWire *vwire_2 = vwire->get_allwire().back();
  dis_v[0] = gds.wire_dist(*vwire_1);
  dis_v[1] = gds.wire_dist(*vwire_2);
  // float bmin_v = dis_v[0] - v_pitch/2. - tolerance;
  // float bmax_v = dis_v[1] + v_pitch/2. + tolerance;

  // define W plane range
  const GeomWire *wwire_1 = wwire->get_allwire().front();
  const GeomWire *wwire_2 = wwire->get_allwire().back();
  float bmin_w = gds.wire_dist(*wwire_1) - w_pitch/2. - tolerance;
  float bmax_w = gds.wire_dist(*wwire_2) + w_pitch/2. + tolerance;

  float dis_point[4]; //left, top, bottom, right
  bool flag_point[4];
  PointVector pcell;
  flag_point[0] = test_point(pcell, dis_u[0], dis_v[0], bmin_w, bmax_w, u_pitch, v_pitch, dis_point[0]);
  flag_point[1] = test_point(pcell, dis_u[0], dis_v[1], bmin_w, bmax_w, u_pitch, v_pitch, dis_point[1]);
  flag_point[2] = test_point(pcell, dis_u[1], dis_v[0], bmin_w, bmax_w, u_pitch, v_pitch, dis_point[2]);
  flag_point[3] = test_point(pcell, dis_u[1], dis_v[1], bmin_w, bmax_w, u_pitch, v_pitch, dis_point[3]);

  float low_limit, high_limit;
  
  if (!flag_point[0]){
    // Scan left + u 
    low_limit = bmin_w;
    if (dis_point[0] > low_limit) low_limit = dis_point[0];
    high_limit = bmax_w;
    if (dis_point[1] < high_limit) high_limit = dis_point[1];
    if (high_limit >= low_limit){
      test_cross(pcell,dis_u[0], WirePlaneType_t(0), low_limit, high_limit, u_pitch, v_pitch, w_pitch, vwire_1, vwire_2,1,bmin_w,bmax_w);
    }
    // Scan left + v
    high_limit = bmax_w;
    if (dis_point[2] < high_limit) high_limit = dis_point[2];
    if (high_limit >= low_limit){
      test_cross(pcell,dis_v[0], WirePlaneType_t(1), low_limit, high_limit, v_pitch, u_pitch, w_pitch, uwire_1, uwire_2,1,bmin_w,bmax_w);
    }
  }

  if (!flag_point[1]){
    // Scan top - u
    low_limit = bmin_w;
    if (dis_point[0] > low_limit) low_limit = dis_point[0];
    high_limit = bmax_w;
    if (dis_point[1] < high_limit) high_limit = dis_point[1];
    if (high_limit >= low_limit){
      test_cross(pcell,dis_u[0], WirePlaneType_t(0), low_limit, high_limit, u_pitch, v_pitch, w_pitch, vwire_1, vwire_2,-1,bmin_w,bmax_w);
    }
    // Scan top + v
    low_limit = bmin_w;
    if (dis_point[1] > low_limit) low_limit = dis_point[1];
    high_limit = bmax_w;
    if (dis_point[3] < high_limit) high_limit = dis_point[3];
    if (high_limit >= low_limit){
      test_cross(pcell,dis_v[1], WirePlaneType_t(1), low_limit, high_limit, v_pitch, u_pitch, w_pitch, uwire_1, uwire_2,1,bmin_w,bmax_w);
    }
  }
  
  if (!flag_point[2]){
    // Scan bottom +u 
    low_limit = bmin_w;
    if (dis_point[2] > low_limit) low_limit = dis_point[2];
    high_limit = bmax_w;
    if (dis_point[3] < high_limit) high_limit = dis_point[3];
    if (high_limit >= low_limit){
      test_cross(pcell,dis_u[1], WirePlaneType_t(0), low_limit, high_limit, u_pitch, v_pitch, w_pitch, vwire_1, vwire_2,1,bmin_w,bmax_w);
    }
    // Scan bottom -v
    low_limit = bmin_w;
    if (dis_point[0] > low_limit) low_limit = dis_point[0];
    high_limit = bmax_w;
    if (dis_point[2] < high_limit) high_limit = dis_point[2];
    if (high_limit >= low_limit){
      test_cross(pcell,dis_v[0], WirePlaneType_t(1), low_limit, high_limit, v_pitch, u_pitch, w_pitch, uwire_1, uwire_2,-1,bmin_w,bmax_w);
    }
  }
  
  if (!flag_point[3]){
    // Scan right -u
    low_limit = bmin_w;
    if (dis_point[2] > low_limit) low_limit = dis_point[2];
    high_limit = bmax_w;
    if (dis_point[3] < high_limit) high_limit = dis_point[3];
    if (high_limit >= low_limit){
      test_cross(pcell,dis_u[1], WirePlaneType_t(0), low_limit, high_limit, u_pitch, v_pitch, w_pitch, vwire_1, vwire_2,-1,bmin_w,bmax_w);
    }
    // Scan right -v
    low_limit = bmin_w;
    if (dis_point[1] > low_limit) low_limit = dis_point[1];
    high_limit = bmax_w;
    if (dis_point[3] < high_limit) high_limit = dis_point[3];
    if (high_limit >= low_limit){
      test_cross(pcell,dis_v[1], WirePlaneType_t(1), low_limit, high_limit, v_pitch, u_pitch, w_pitch, uwire_1, uwire_2,-1,bmin_w,bmax_w);
    }
  }

  // std::cout << flag_point[0] << " " << flag_point[1] << " " << flag_point[2] << " " << flag_point[3]
  // 	    << " " << dis_point[0] << " " << dis_point[1] << " " << dis_point[2] << " " << dis_point[3] << std::endl;

  //test_crossing(pcell, dis_u[0], dis_v[0], bmin_w, bmax_w, u_pitch, v_pitch, uwire_1, uwire_2, vwire_1, vwire_2,true,true);
  //test_crossing(pcell, dis_u[0], dis_v[1], bmin_w, bmax_w, u_pitch, v_pitch, uwire_1, uwire_2, vwire_1, vwire_2,true,true);
  //test_crossing(pcell, dis_u[1], dis_v[0], bmin_w, bmax_w, u_pitch, v_pitch, uwire_1, uwire_2, vwire_1, vwire_2,true,true);
  //test_crossing(pcell, dis_u[1], dis_v[1], bmin_w, bmax_w, u_pitch, v_pitch, uwire_1, uwire_2, vwire_1, vwire_2,true,true);
  
  
  
  // find the max and min
  if (pcell.size() >=1){
    //Creat a cell and then get all the wires in ... 
    double u_max = -1e9, u_min = 1e9;
    double v_max = -1e9, v_min = 1e9;
    double w_max = -1e9, w_min = 1e9;
    // tolerance = 0.;
    for (int k=0;k!=pcell.size();k++){
      double udist = gds.wire_dist(pcell.at(k),WirePlaneType_t(0));
      if (udist + tolerance> u_max) u_max = udist + tolerance;
      if (udist - tolerance< u_min) u_min = udist - tolerance;
      double vdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(1));
      if (vdist + tolerance> v_max) v_max = vdist + tolerance;
      if (vdist - tolerance< v_min) v_min = vdist + tolerance;
      double wdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(2));
      if (wdist + tolerance> w_max) w_max = wdist + tolerance;
      if (wdist - tolerance< w_min) w_min = wdist - tolerance;
    }
    
    int ident = holder.get_cell_no();
   
    

    const GeomWire* uwire_min = gds.closest(u_min,WirePlaneType_t(0));
    if (uwire_min->index() < uwire_1->index()){
      uwire_min = uwire_1;
    }else if (uwire_min->index() > uwire_2->index()){
      uwire_min = uwire_2;
    }
    const GeomWire* uwire_max = gds.closest(u_max,WirePlaneType_t(0));
    if (uwire_max->index() < uwire_1->index()){
      uwire_max = uwire_1;
    }else if (uwire_max->index() > uwire_2->index()){
      uwire_max = uwire_2;
    }
    const GeomWire* vwire_min = gds.closest(v_min,WirePlaneType_t(1));
    if (vwire_min->index() < vwire_1->index()){
      vwire_min = vwire_1;
    }else if (vwire_min->index() > vwire_2->index()){
      vwire_min = vwire_2;
    }
    const GeomWire* vwire_max = gds.closest(v_max,WirePlaneType_t(1));
    if (vwire_max->index() < vwire_1->index()){
      vwire_max = vwire_1;
    }else if (vwire_max->index() > vwire_2->index()){
      vwire_max = vwire_2;
    }
    const GeomWire* wwire_min = gds.closest(w_min,WirePlaneType_t(2));//.bounds(w_min,WirePlaneType_t(2)).second;
    if (wwire_min->index() < wwire_1->index()){
      wwire_min = wwire_1;
    }else if (wwire_min->index() > wwire_2->index()){
      wwire_min = wwire_2;
    }
    const GeomWire* wwire_max = gds.closest(w_max,WirePlaneType_t(2));//.bounds(w_max,WirePlaneType_t(2)).first;
    if (wwire_max->index() < wwire_1->index()){
      wwire_max = wwire_1;
    }else if (wwire_max->index() > wwire_2->index()){
      wwire_max = wwire_2;
    }
    u_min = gds.wire_dist(*uwire_min)-u_pitch/2.;
    u_max = gds.wire_dist(*uwire_max)+u_pitch/2.;
    v_min = gds.wire_dist(*vwire_min)-v_pitch/2.;
    v_max = gds.wire_dist(*vwire_max)+v_pitch/2.;
    w_min = gds.wire_dist(*wwire_min)-w_pitch/2.;
    w_max = gds.wire_dist(*wwire_max)+w_pitch/2.;
    
    // std::cout << uwire_min->index() << " " << uwire_max->index() << " " 
    // 	      << vwire_min->index() << " " << vwire_max->index() << " " 
    // 	      << wwire_min->index() << " " << wwire_max->index() << " " 
    // 	      << vwire_1->index() << " " << vwire_2->index() << " " 
    // 	      << w_min/units::m << " " << w_max / units::m << " " 
    // 	      << std::endl;
    
    
    GeomWireSelection ugroup,vgroup,wgroup;
    // find U group
    if (uwire_min!=0)
      ugroup.push_back(uwire_min);
    if (uwire_max!=uwire_min&&uwire_max!=0){
      ugroup.push_back(uwire_max);
      const GeomWire* uwire_p1 = gds.by_planeindex(WirePlaneType_t(0),uwire_min->index()+1);
      if (uwire_p1!=uwire_max&&uwire_p1!=0){
	ugroup.push_back(uwire_p1);
	const GeomWire* uwire_n1 = gds.by_planeindex(WirePlaneType_t(0),uwire_max->index()-1);
	if (uwire_n1!=uwire_p1&&uwire_n1!=0){
	  ugroup.push_back(uwire_n1);
	}
      }
    }
    // find V group
    if (vwire_min!=0)
      vgroup.push_back(vwire_min);
    if (vwire_max!=vwire_min&&vwire_max!=0){
      vgroup.push_back(vwire_max);
      const GeomWire* vwire_p1 = gds.by_planeindex(WirePlaneType_t(1),vwire_min->index()+1);
      if (vwire_p1!=vwire_max&&vwire_p1!=0){
	vgroup.push_back(vwire_p1);
	const GeomWire* vwire_n1 = gds.by_planeindex(WirePlaneType_t(1),vwire_max->index()-1);
	if (vwire_n1!=vwire_p1&&vwire_n1!=0){
	  vgroup.push_back(vwire_n1);
	}
      }
    }
    // find W group
    if (wwire_min!=0)
      wgroup.push_back(wwire_min);
    if (wwire_max!=wwire_min&&wwire_max!=0){
      wgroup.push_back(wwire_max);
      const GeomWire* wwire_p1 = gds.by_planeindex(WirePlaneType_t(2),wwire_min->index()+1);
      if (wwire_p1!=wwire_max&&wwire_p1!=0){
	wgroup.push_back(wwire_p1);
	const GeomWire* wwire_n1 = gds.by_planeindex(WirePlaneType_t(2),wwire_max->index()-1);
	if (wwire_n1!=wwire_p1&&wwire_n1!=0){
	  wgroup.push_back(wwire_n1);
	}
      }
    }
    // initialize map
    std::map<const GeomWire*, int> wiremap;
    for (auto it = ugroup.begin();it!=ugroup.end();it++){
      wiremap[*it]=0;
    }
    for (auto it = vgroup.begin();it!=vgroup.end();it++){
      wiremap[*it]=0;
    }
    for (auto it = wgroup.begin();it!=wgroup.end();it++){
      wiremap[*it]=0;
    }
    
    // std::cout << ugroup.size() << " " << vgroup.size() << " " << wgroup.size() << " " 
    // 	      << uwire_max->index()-uwire_min->index()+1 << " " 
    // 	      << vwire_max->index()-vwire_min->index()+1 << " " 
    // 	      << wwire_max->index()-wwire_min->index()+1 << " " 
    // 	      << std::endl;
    
    // check crossing 
    for (int i=0;i!=ugroup.size();i++){
      const GeomWire *uwire = ugroup.at(i);
      for (int j=0;j!=vgroup.size();j++){
	const GeomWire *vwire = vgroup.at(j);
	if (check_crossing(uwire, vwire, u_pitch, v_pitch, WirePlaneType_t(2), w_min, w_max, tolerance)){
	  wiremap[uwire]++;
	  wiremap[vwire]++;
	}
      }
    }
    for (int i=0;i!=wgroup.size();i++){
      const GeomWire *wwire = wgroup.at(i);
      for (int j=0;j!=vgroup.size();j++){
	const GeomWire *vwire = vgroup.at(j);
	if (check_crossing(wwire, vwire, w_pitch, v_pitch, WirePlaneType_t(0), u_min, u_max, tolerance)){
	  wiremap[wwire]++;
	  wiremap[vwire]++;
	}
      }
    }
    for (int i=0;i!=ugroup.size();i++){
      const GeomWire *uwire = ugroup.at(i);
      for (int j=0;j!=wgroup.size();j++){
	const GeomWire *wwire = wgroup.at(j);
	if (check_crossing(uwire, wwire, u_pitch, w_pitch, WirePlaneType_t(1), v_min, v_max, tolerance)){
	  wiremap[uwire]++;
	  wiremap[wwire]++;
	}
      }
    }

    
    // 	
    // 	for (int k=0;k!=wgroup.size();k++){
    // 	  const GeomWire *wwire = wgroup.at(k);
    // 	  if (check_crossing(uwire,vwire,wwire,u_pitch,v_pitch,w_pitch,tolerance)){
    // 	    wiremap[uwire]++;
    // 	    wiremap[vwire]++;
    // 	    wiremap[wwire]++;
    // 	  }
    // 	}
    //   }
    //   //      std::cout << wiremap[uwire] << " " << u_pitch << " " << tolerance << std::endl;
    // }

    int flag = 1;
    if (wiremap[uwire_min]==0 ){
      const GeomWire* wire_temp = gds.by_planeindex(WirePlaneType_t(0),uwire_min->index()+1);
      if (wiremap.find(wire_temp)!=wiremap.end() && wiremap[wire_temp]>0){ 
	uwire_min = wire_temp;
      }else{
	flag = 0;
      }
    }
    if (flag==1){
      if (wiremap[vwire_min]==0 ){
	const GeomWire* wire_temp = gds.by_planeindex(WirePlaneType_t(1),vwire_min->index()+1);
	if (wiremap.find(wire_temp)!=wiremap.end() && wiremap[wire_temp]>0){ 
	  vwire_min = wire_temp;
	}else{
	  flag = 0;
	}
      }
    }
    if (flag==1){
      if (wiremap[wwire_min]==0 ){
	const GeomWire* wire_temp = gds.by_planeindex(WirePlaneType_t(2),wwire_min->index()+1);
	if (wiremap.find(wire_temp)!=wiremap.end() && wiremap[wire_temp]>0){ 
	  wwire_min = wire_temp;
	}else{
	  flag = 0;
	}
      }
    }
    if (flag==1){
      if (wiremap[uwire_max]==0 ){
	const GeomWire* wire_temp = gds.by_planeindex(WirePlaneType_t(0),uwire_max->index()-1);
	if (wiremap.find(wire_temp)!=wiremap.end() && wiremap[wire_temp]>0){ 
	  uwire_max = wire_temp;
	}else{
	  flag = 0;
	}
      }
    }
    if (flag==1){
      if (wiremap[vwire_max]==0 ){
	const GeomWire* wire_temp = gds.by_planeindex(WirePlaneType_t(1),vwire_max->index()-1);
	if (wiremap.find(wire_temp)!=wiremap.end() && wiremap[wire_temp]>0){ 
	  vwire_max = wire_temp;
	}else{
	  flag = 0;
	}
      }
    }
    if (flag==1){
      if (wiremap[wwire_max]==0 ){
	const GeomWire* wire_temp = gds.by_planeindex(WirePlaneType_t(2),wwire_max->index()-1);
	if (wiremap.find(wire_temp)!=wiremap.end() && wiremap[wire_temp]>0){ 
	  wwire_max = wire_temp;
	}else{
	  flag = 0;
	}
      }
    }

   
   
    if (flag==1){
      SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
      mcell->SetTimeSlice(time_slice);
      holder.AddCell_No();
      // // Creat the Merge Wires
      // MergeGeomWire *mwire_u;
      // MergeGeomWire *mwire_v;
      // MergeGeomWire *mwire_w;


      // Inser U    
      for (int k=uwire_min->index();k!=uwire_max->index()+1;k++){
	const GeomWire *uwire1 = gds.by_planeindex(WirePlaneType_t(0),k);
	float charge = 0;
	if (wirechargemap.find(uwire1)!=wirechargemap.end())
	  charge = wirechargemap[uwire1];
	mcell->AddWire(uwire1,WirePlaneType_t(0),charge);
	// if (k==uwire_min->index()){
	//   mwire_u = new MergeGeomWire(ident,*uwire1);
	// }else{
	//   mwire_u->AddWire(*uwire1);
	// }
      }
      //Insert V
      for (int k=vwire_min->index();k!=vwire_max->index()+1;k++){
	const GeomWire *vwire1 = gds.by_planeindex(WirePlaneType_t(1),k);
	float charge = 0;
	if (wirechargemap.find(vwire1)!=wirechargemap.end())
	  charge = wirechargemap[vwire1];
	mcell->AddWire(vwire1,WirePlaneType_t(1),charge);
	// if (k==vwire_min->index()){
	//   mwire_v = new MergeGeomWire(ident,*vwire1);
	// }else{
	//   mwire_v->AddWire(*vwire1);
	// }
      }
      // Insert W
      for (int k=wwire_min->index();k!=wwire_max->index()+1;k++){
	const GeomWire *wwire1 = gds.by_planeindex(WirePlaneType_t(2),k);
	float charge = 0;
	if (wirechargemap.find(wwire1)!=wirechargemap.end())
	  charge = wirechargemap[wwire1];
	mcell->AddWire(wwire1,WirePlaneType_t(2),charge);
	// if (k==wwire_min->index()){
	//   mwire_w = new MergeGeomWire(ident,*wwire1);
	// }else{
	//   mwire_w->AddWire(*wwire1);
	// }
      }
      
     
      
      return mcell;
    }
  }
  
  return 0;
}
void WireCell2dToy::LowmemTiling::Print_maps(){
  std::cout << wire_type_map.size() << " " << cell_wires_map.size()
	    << " " << wire_cells_map.size() << " " << wire_pwire_map.size()
	    << " " << pwire_wires_map.size() << std::endl;
}

void WireCell2dToy::LowmemTiling::re_establish_maps(){
  std::map<const GeomCell*, const GeomWire*> cell_upwire_map;
  std::map<const GeomCell*, const GeomWire*> cell_vpwire_map;
  std::map<const GeomCell*, const GeomWire*> cell_wpwire_map;
  std::vector<const GeomCell*> all_cells;
  
  for (auto it = cell_wires_map.begin(); it!=cell_wires_map.end(); it++){
    const GeomCell *mcell = it->first;
    all_cells.push_back(mcell);
    GeomWireSelection wires = it->second;
    for (auto it1 = wires.begin(); it1!=wires.end(); it1++){
      MergeGeomWire *mwire = (MergeGeomWire*)(*it1);
      MergeGeomWire *p_mwire = (MergeGeomWire*)wire_pwire_map[mwire];
      if (mwire->get_allwire().at(0)->plane()==WirePlaneType_t(0)){
    	cell_upwire_map[mcell] = p_mwire;
      }else if (mwire->get_allwire().at(0)->plane()==WirePlaneType_t(1)){
    	cell_vpwire_map[mcell] = p_mwire;
      }else if (mwire->get_allwire().at(0)->plane()==WirePlaneType_t(2)){
    	cell_wpwire_map[mcell] = p_mwire;
      }
    }
  }

  //  std::cout << "T: " << cell_upwire_map.size() << " " << cell_vpwire_map.size() << " " << cell_wpwire_map.size() << " " << all_cells.size() << std::endl;
  
   // delete all wires ...
  GeomWireSelection wires;
  for (auto it = wire_pwire_map.begin(); it!=wire_pwire_map.end(); it++){
    MergeGeomWire *mwire = (MergeGeomWire*)it->first;
    wires.push_back(mwire);
  }
  for (auto it = wires.begin(); it!=wires.end(); it++){
    MergeGeomWire *mwire = (MergeGeomWire*)(*it);
    remove_wire_clear(mwire);
  }
  // clear all maps ...
  // wire_pwire_map.clear();
  // wire_cells_map.clear();
  cell_wires_map.clear();
  

  // fill in all maps again ... 
  for (auto it = all_cells.begin(); it!=all_cells.end(); it++){
    SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);

      // create new mwires 
    GeomWireSelection uwires = mcell->get_uwires();
    GeomWireSelection vwires = mcell->get_vwires();
    GeomWireSelection wwires = mcell->get_wwires();
    MergeGeomWire *mwire_u = new MergeGeomWire(0,uwires);
    MergeGeomWire *mwire_v = new MergeGeomWire(0,vwires);
    MergeGeomWire *mwire_w = new MergeGeomWire(0,wwires);

    // create the map  
    // cell to wires
    GeomWireSelection wires;
    wires.push_back(mwire_u);
    wires.push_back(mwire_v);
    wires.push_back(mwire_w);
    cell_wires_map[mcell] = wires;
    
    // wire to cells
    GeomCellSelection cells;
    cells.push_back(mcell);
    wire_cells_map[mwire_u] = cells;
    wire_cells_map[mwire_v] = cells;
    wire_cells_map[mwire_w] = cells;

    MergeGeomWire *uwire = (MergeGeomWire*)cell_upwire_map[mcell];
    MergeGeomWire *vwire = (MergeGeomWire*)cell_vpwire_map[mcell];
    MergeGeomWire *wwire = (MergeGeomWire*)cell_wpwire_map[mcell];

    // wire to parent wire
    wire_pwire_map[mwire_u] = uwire;
    wire_pwire_map[mwire_v] = vwire;
    wire_pwire_map[mwire_w] = wwire;
	  
    // wire types
    wire_type_map[mwire_u] = wire_type_map[uwire];
    wire_type_map[mwire_v] = wire_type_map[vwire];     
    wire_type_map[mwire_w] = wire_type_map[wwire];
	  

    //std::cout <<  wire_type_map[mwire_u] << " " << wire_type_map[mwire_v] << " " << wire_type_map[mwire_w]  << std::endl;

    //  std::cout << pwire_wires_map[uwire].size() << std::endl;
    
    //parent wire to wires
    if (pwire_wires_map.find(uwire)==pwire_wires_map.end()){
      GeomWireSelection wires;
      wires.push_back(mwire_u);
      pwire_wires_map[uwire] = wires;
    }else{
      pwire_wires_map[uwire].push_back(mwire_u);
    }
    
    if (pwire_wires_map.find(vwire)==pwire_wires_map.end()){
      GeomWireSelection wires;
      wires.push_back(mwire_v);
      pwire_wires_map[vwire] = wires;
    }else{
      pwire_wires_map[vwire].push_back(mwire_v);
    }
    
    if (pwire_wires_map.find(wwire)==pwire_wires_map.end()){
      GeomWireSelection wires;
      wires.push_back(mwire_w);
      pwire_wires_map[wwire] = wires;
    }else{
      pwire_wires_map[wwire].push_back(mwire_w);
    }

    // std::cout << pwire_wires_map[uwire].size() << std::endl;
    
  }
}

void WireCell2dToy::LowmemTiling::local_deghosting1(std::set<WireCell::SlimMergeGeomCell*>& good_mcells){
   
  std::map<const GeomWire*, float> wire_score_map;
  for (auto it = three_good_wire_cells.begin(); it!= three_good_wire_cells.end(); it++){
    SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);
    GeomWireSelection uwires = mcell->get_uwires();
    GeomWireSelection vwires = mcell->get_vwires();
    GeomWireSelection wwires = mcell->get_wwires();

    for(auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
      if (wire_score_map.find(*it1)==wire_score_map.end()){
	wire_score_map[*it1] = 1;
      }else{
	wire_score_map[*it1] = wire_score_map[*it1] + 1;
      }
    }
    for(auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
      if (wire_score_map.find(*it1)==wire_score_map.end()){
	wire_score_map[*it1] = 1;
      }else{
	wire_score_map[*it1] = wire_score_map[*it1] + 1;
      }
    }
    for(auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
      if (wire_score_map.find(*it1)==wire_score_map.end()){
	wire_score_map[*it1] = 1;
      }else{
	wire_score_map[*it1] = wire_score_map[*it1] + 1;
      }
    }
  }
  for (auto it = two_good_wire_cells.begin(); it!= two_good_wire_cells.end(); it++){
    SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);
    GeomWireSelection uwires = mcell->get_uwires();
    GeomWireSelection vwires = mcell->get_vwires();
    GeomWireSelection wwires = mcell->get_wwires();
    
    std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
    bool flag_u = true, flag_v = true, flag_w = true;
    for (auto it1 = bad_planes.begin(); it1!=bad_planes.end(); it1++){
      if (*it1 == WirePlaneType_t(0))
  	flag_u = false;
      if (*it1 == WirePlaneType_t(1))
  	flag_v = false;
      if (*it1 == WirePlaneType_t(2))
  	flag_w = false;
    }
    if (flag_u){
      for(auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
	if (wire_score_map.find(*it1)==wire_score_map.end()){
	  wire_score_map[*it1] = 1;
	}else{
	  wire_score_map[*it1] = wire_score_map[*it1] + 1;
	}
      }
    }
    if (flag_v){
      for(auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
	if (wire_score_map.find(*it1)==wire_score_map.end()){
	  wire_score_map[*it1] = 1;
	}else{
	  wire_score_map[*it1] = wire_score_map[*it1] + 1;
	}
      }
    }
    if (flag_w){
      for(auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
	if (wire_score_map.find(*it1)==wire_score_map.end()){
	  wire_score_map[*it1] = 1;
	}else{
	  wire_score_map[*it1] = wire_score_map[*it1] + 1;
	}
      }
    }
  }
  // finished filling fill in the wire score map ... 

  std::map<SlimMergeGeomCell*, float> mcell_high_score_map;
  
  for (auto it = two_good_wire_cells.begin(); it!= two_good_wire_cells.end(); it++){
    SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);
    
    GeomWireSelection uwires = mcell->get_uwires();
    GeomWireSelection vwires = mcell->get_vwires();
    GeomWireSelection wwires = mcell->get_wwires();
    
    std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
    bool flag_u = true, flag_v = true, flag_w = true;
    for (auto it1 = bad_planes.begin(); it1!=bad_planes.end(); it1++){
      if (*it1 == WirePlaneType_t(0))
  	flag_u = false;
      if (*it1 == WirePlaneType_t(1))
  	flag_v = false;
      if (*it1 == WirePlaneType_t(2))
  	flag_w = false;
    }

    float score_u=-1, score_v=-1, score_w=-1;
    
    if (flag_u){
      float sum1 = 0;
      float sum2 = 0;
      for(auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
	sum2 ++;
	sum1 += 1.0/wire_score_map[*it1];
      }
      score_u = sum1/sum2;
    }

    if (flag_v){
      float sum1 = 0;
      float sum2 = 0;
      for(auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
	sum2 ++;
	sum1 += 1.0/wire_score_map[*it1];
      }
      score_v = sum1/sum2;
    }

    if (flag_w){
      float sum1 = 0;
      float sum2 = 0;
      for(auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
	sum2 ++;
	sum1 += 1.0/wire_score_map[*it1];
      }
      score_w = sum1/sum2;
    }

    mcell_high_score_map[mcell] = std::max(std::max(score_u,score_v),score_w);
  }

  for (auto it = three_good_wire_cells.begin(); it!= three_good_wire_cells.end(); it++){
    SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);
    
    GeomWireSelection uwires = mcell->get_uwires();
    GeomWireSelection vwires = mcell->get_vwires();
    GeomWireSelection wwires = mcell->get_wwires();
    
    std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
    bool flag_u = true, flag_v = true, flag_w = true;
    for (auto it1 = bad_planes.begin(); it1!=bad_planes.end(); it1++){
      if (*it1 == WirePlaneType_t(0))
    	flag_u = false;
      if (*it1 == WirePlaneType_t(1))
    	flag_v = false;
      if (*it1 == WirePlaneType_t(2))
    	flag_w = false;
    }

    float score_u=-1, score_v=-1, score_w=-1;
    
    if (flag_u){
      float sum1 = 0;
      float sum2 = 0;
      for(auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
    	sum2 ++;
    	sum1 += 1.0/wire_score_map[*it1];
      }
      score_u = sum1/sum2;
    }

    if (flag_v){
      float sum1 = 0;
      float sum2 = 0;
      for(auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
    	sum2 ++;
    	sum1 += 1.0/wire_score_map[*it1];
      }
      score_v = sum1/sum2;
    }

    if (flag_w){
      float sum1 = 0;
      float sum2 = 0;
      for(auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
    	sum2 ++;
    	sum1 += 1.0/wire_score_map[*it1];
      }
      score_w = sum1/sum2;
    }

    mcell_high_score_map[mcell] = std::max(std::max(score_u,score_v),score_w);
    if (good_mcells.find(mcell)!=good_mcells.end())
      mcell_high_score_map[mcell] = 1;
  }

  
  GeomCellSelection to_be_removed;

  for (auto it = two_good_wire_cells.begin(); it!= two_good_wire_cells.end(); it++){
    SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);

    // const GeomWire *uwire = mcell->get_uwires().at(int(mcell->get_uwires().size()/2));
    // const GeomWire *vwire = mcell->get_vwires().at(int(mcell->get_vwires().size()/2));
    // Vector abc;
    // gds.crossing_point(*uwire, *vwire, abc);
    // std::cout << abc.x/units::cm << " " << abc.y/units::cm << " " << abc.z/units::cm << std::endl;
    
    GeomWireSelection uwires = mcell->get_uwires();
    GeomWireSelection vwires = mcell->get_vwires();
    GeomWireSelection wwires = mcell->get_wwires();
    
    GeomWireSelection mwires = cell_wires_map[mcell];
    std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
    int count = 0;
    for (auto it1=mwires.begin(); it1!=mwires.end(); it1++){
      MergeGeomWire *mwire = (MergeGeomWire*)(*it1);

      
      
      if (find(bad_planes.begin(),bad_planes.end(),mwire->get_allwire().front()->plane())==bad_planes.end()){
	int mcell_lwire, mcell_hwire;
	if (mwire->get_allwire().front()->plane()==WirePlaneType_t(0)){
	  mcell_lwire = uwires.front()->index();
	  mcell_hwire = uwires.back()->index();
	}else if (mwire->get_allwire().front()->plane()==WirePlaneType_t(1)){
	  mcell_lwire = vwires.front()->index();
	  mcell_hwire = vwires.back()->index();
	}else{
	  mcell_lwire = wwires.front()->index();
	  mcell_hwire = wwires.back()->index();
	}
	
	GeomCellSelection cells = wire_cells_map[mwire];
	for (auto it2 = cells.begin(); it2 != cells.end(); it2++){
	  SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*)(*it2);
	  if (mcell1==mcell) continue;
	  
	  // const GeomWire *uwire = mcell1->get_uwires().at(int(mcell1->get_uwires().size()/2));
	  // const GeomWire *vwire = mcell1->get_vwires().at(int(mcell1->get_vwires().size()/2));
	  // Vector abc;
	  // gds.crossing_point(*uwire, *vwire, abc);
	  // std::cout << abc.x/units::cm << " " << abc.y/units::cm << " " << abc.z/units::cm << std::endl;
	  

	  // std::cout << "a: " << mcell_high_score_map[mcell1] << " " << std::endl;
	  if (mcell_high_score_map[mcell1]>0.75){
	    int mcell1_lwire, mcell1_hwire;
	    if (mwire->get_allwire().front()->plane()==WirePlaneType_t(0)){
	      mcell1_lwire = mcell1->get_uwires().front()->index();
	      mcell1_hwire = mcell1->get_uwires().back()->index();
	    }else if (mwire->get_allwire().front()->plane()==WirePlaneType_t(1)){
	      mcell1_lwire = mcell1->get_vwires().front()->index();
	      mcell1_hwire = mcell1->get_vwires().back()->index();
	    }else{
	      mcell1_lwire = mcell1->get_wwires().front()->index();
	      mcell1_hwire = mcell1->get_wwires().back()->index();
	    }
	    int min_wire = std::max(mcell_lwire,mcell1_lwire);
	    int max_wire = std::min(mcell_hwire,mcell1_hwire);
	    //std::cout << "b: " << (max_wire - min_wire + 1.0)/(mcell_hwire - mcell_lwire + 1.0) << std::endl;
	    if ((max_wire - min_wire + 1.0)/(mcell_hwire - mcell_lwire + 1.0)>=0.75){
	      count ++;
	      break;
	    }
	  }
	}
	//	count ++;
      }
    }
    if (count == 2){
      to_be_removed.push_back(mcell);
    }
    // std::cout << mwires.size() << " " << count << std::endl;
  }
  for (auto it=to_be_removed.begin(); it!=to_be_removed.end(); it++){
    Erase_Cell((SlimMergeGeomCell*)(*it));
  }
  
}



GeomCellSelection WireCell2dToy::LowmemTiling::local_deghosting(std::set<SlimMergeGeomCell*>& potential_good_mcells, std::set<WireCell::SlimMergeGeomCell*>& good_mcells, bool flag_del){
  
  // do the local deghosting
  std::set<const GeomWire*> used_uwires;
  std::set<const GeomWire*> used_vwires;
  std::set<const GeomWire*> used_wwires;

  // find the wires for three-wire cells and also good cells
  for (auto it = three_good_wire_cells.begin(); it!= three_good_wire_cells.end(); it++){
    SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);
    if (good_mcells.find(mcell)!=good_mcells.end()){
      GeomWireSelection uwires = mcell->get_uwires();
      GeomWireSelection vwires = mcell->get_vwires();
      GeomWireSelection wwires = mcell->get_wwires();
      for (auto it1 = uwires.begin(); it1!= uwires.end(); it1++){
	used_uwires.insert(*it1);
      }
      for (auto it1 = vwires.begin(); it1!= vwires.end(); it1++){
	used_vwires.insert(*it1);
      }
      for (auto it1 = wwires.begin(); it1!= wwires.end(); it1++){
	used_wwires.insert(*it1);
      }
    }
  }

  GeomCellSelection to_be_removed;

  for (auto it = two_good_wire_cells.begin(); it!= two_good_wire_cells.end(); it++){
    SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);
    if (potential_good_mcells.find(mcell)==potential_good_mcells.end()){
      GeomWireSelection uwires = mcell->get_uwires();
      GeomWireSelection vwires = mcell->get_vwires();
      GeomWireSelection wwires = mcell->get_wwires();

      std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
      bool flag_u = true, flag_v = true, flag_w = true;
      for (auto it1 = bad_planes.begin(); it1!=bad_planes.end(); it1++){
  	if (*it1 == WirePlaneType_t(0))
  	  flag_u = false;
  	if (*it1 == WirePlaneType_t(1))
  	  flag_v = false;
  	if (*it1 == WirePlaneType_t(2))
  	  flag_w = false;
      }

      bool save_flag = false;

      if (flag_u && !save_flag){
  	int nwire = 0;
  	int nwire_common = 0;
  	for (auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
  	  nwire++;
  	  if (used_uwires.find((*it1))!=used_uwires.end())
  	    nwire_common++;
  	}
  	if (nwire_common != nwire) save_flag =true;
      }
      
      if (flag_v && !save_flag){
  	int nwire = 0;
  	int nwire_common = 0;
  	for (auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
  	  nwire++;
  	  if (used_vwires.find((*it1))!=used_vwires.end())
  	    nwire_common++;
  	}
  	if (nwire_common != nwire) save_flag =true;
      }
      
      if (flag_w && !save_flag){
  	int nwire = 0;
  	int nwire_common = 0;
  	for (auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
  	  nwire++;
  	  if (used_wwires.find((*it1))!=used_wwires.end())
  	    nwire_common++;
  	}
  	if (nwire_common != nwire) save_flag =true;
      }

      if (!save_flag) to_be_removed.push_back(mcell);
    }  
  }

  // std::cout << two_good_wire_cells.size() << " " << to_be_removed.size() << std::endl;
  //delete absolutely that can be deleted ... 
  for (auto it=to_be_removed.begin(); it!=to_be_removed.end(); it++){
    Erase_Cell((SlimMergeGeomCell*)(*it));
  }
  to_be_removed.clear();

  
  // std::cout << two_good_wire_cells.size() << std::endl;
  std::map<const GeomWire*, float> wire_score_map;
  for (auto it = two_good_wire_cells.begin(); it!= two_good_wire_cells.end(); it++){
    SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);
    GeomWireSelection uwires = mcell->get_uwires();
    GeomWireSelection vwires = mcell->get_vwires();
    GeomWireSelection wwires = mcell->get_wwires();
    
    std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
    bool flag_u = true, flag_v = true, flag_w = true;
    for (auto it1 = bad_planes.begin(); it1!=bad_planes.end(); it1++){
      if (*it1 == WirePlaneType_t(0))
  	flag_u = false;
      if (*it1 == WirePlaneType_t(1))
  	flag_v = false;
      if (*it1 == WirePlaneType_t(2))
  	flag_w = false;
    }
    if (flag_u){
      for(auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
	if (wire_score_map.find(*it1)==wire_score_map.end()){
	  wire_score_map[*it1] = 1;
	}else{
	  wire_score_map[*it1] = wire_score_map[*it1] + 1;
	}
      }
    }
    if (flag_v){
      for(auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
	if (wire_score_map.find(*it1)==wire_score_map.end()){
	  wire_score_map[*it1] = 1;
	}else{
	  wire_score_map[*it1] = wire_score_map[*it1] + 1;
	}
      }
    }
    if (flag_w){
      for(auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
	if (wire_score_map.find(*it1)==wire_score_map.end()){
	  wire_score_map[*it1] = 1;
	}else{
	  wire_score_map[*it1] = wire_score_map[*it1] + 1;
	}
      }
    }
  }

  
  for (auto it = two_good_wire_cells.begin(); it!= two_good_wire_cells.end(); it++){
    SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);
    
    GeomWireSelection uwires = mcell->get_uwires();
    GeomWireSelection vwires = mcell->get_vwires();
    GeomWireSelection wwires = mcell->get_wwires();
    
    std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
    bool flag_u = true, flag_v = true, flag_w = true;
    for (auto it1 = bad_planes.begin(); it1!=bad_planes.end(); it1++){
      if (*it1 == WirePlaneType_t(0))
  	flag_u = false;
      if (*it1 == WirePlaneType_t(1))
  	flag_v = false;
      if (*it1 == WirePlaneType_t(2))
  	flag_w = false;
    }

    float score_u=-1, score_v=-1, score_w=-1;
    
    if (flag_u){
      float sum1 = 0;
      float sum2 = 0;
      for(auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
	sum2 ++;
	if (used_uwires.find((*it1))==used_uwires.end()){
	  sum1 += 1.0/wire_score_map[*it1];
	}
      }
      score_u = sum1/sum2;
    }

    if (flag_v){
      float sum1 = 0;
      float sum2 = 0;
      for(auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
	sum2 ++;
	if (used_vwires.find((*it1))==used_vwires.end()){
	  sum1 += 1.0/wire_score_map[*it1];
	}
      }
      score_v = sum1/sum2;
    }

    if (flag_w){
      float sum1 = 0;
      float sum2 = 0;
      for(auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
	sum2 ++;
	if (used_wwires.find((*it1))==used_wwires.end()){
	  sum1 += 1.0/wire_score_map[*it1];
	}
      }
      score_w = sum1/sum2;
    }
    
    if (score_u <=0.5 && score_v <= 0.5 && score_w <=0.5 &&
	potential_good_mcells.find(mcell) == potential_good_mcells.end()){
      to_be_removed.push_back(mcell);
    }
    //    std::cout << "Xin: " << score_u << " " << score_v << " " << score_w << std::endl;
  }

  if (flag_del){
    for (auto it=to_be_removed.begin(); it!=to_be_removed.end(); it++){
      Erase_Cell((SlimMergeGeomCell*)(*it));
    }
    to_be_removed.clear();
  }
  
  
  return to_be_removed;
  
  // // loop through the two wire cell again
  // // find the ones that are absolute not matching with others ...
  // // for the rest, return it to be bad for consideration ... 
  // std::map<std::pair<int,int>,GeomCellSelection> proj_u;
  // std::map<std::pair<int,int>,GeomCellSelection> proj_v;
  // std::map<std::pair<int,int>,GeomCellSelection> proj_w;
  //

  //   if (flag_u){
  //     std::pair<int,int> u_1D = std::make_pair(uwires.front()->index(),uwires.back()->index());
  //     if (proj_u.find(u_1D)!=proj_u.end()){
  // 	GeomCellSelection& cells = proj_u[u_1D];
  // 	cells.push_back(mcell);
  //     }else{
  // 	bool flag_save = true;
  // 	std::vector<std::pair<int,int>> to_be_removed;
  // 	for (auto it1 = proj_u.begin(); it1!=proj_u.end(); it1++){
  // 	  std::pair<int,int> testing = it1->first;
  // 	  if (u_1D.first >= testing.first && u_1D.second <= testing.second){
  // 	    flag_save = false;
  // 	    break;
  // 	  }else if (u_1D.first <= testing.first && u_1D.second >= testing.second){
  // 	    to_be_removed.push_back(testing);
  // 	  }
  // 	}
  // 	if (flag_save){
  // 	  GeomCellSelection cells;
  // 	  cells.push_back(mcell);
  // 	  proj_u[u_1D] = cells;
  // 	}
  // 	for (auto it1 = to_be_removed.begin(); it1!=to_be_removed.end(); it1++){
  // 	  proj_u.erase((*it1));
  // 	}
  //     }
  //   }

    
  //   if (flag_v){
  //     std::pair<int,int> v_1D = std::make_pair(vwires.front()->index(),vwires.back()->index());
  //     if (proj_v.find(v_1D)!=proj_v.end()){
  // 	GeomCellSelection& cells = proj_v[v_1D];
  // 	cells.push_back(mcell);
  //     }else{
  // 	bool flag_save = true;
  // 	std::vector<std::pair<int,int>> to_be_removed;
  // 	for (auto it1 = proj_v.begin(); it1!=proj_v.end(); it1++){
  // 	  std::pair<int,int> testing = it1->first;
  // 	  if (v_1D.first >= testing.first && v_1D.second <= testing.second){
  // 	    flag_save = false;
  // 	    break;
  // 	  }else if (v_1D.first <= testing.first && v_1D.second >= testing.second){
  // 	    to_be_removed.push_back(testing);
  // 	  }
  // 	}
  // 	if (flag_save){
  // 	  GeomCellSelection cells;
  // 	  cells.push_back(mcell);
  // 	  proj_v[v_1D] = cells;
  // 	}
  // 	for (auto it1 = to_be_removed.begin(); it1!=to_be_removed.end(); it1++){
  // 	  proj_v.erase((*it1));
  // 	}
  //     }
  //   }

    
  //   if (flag_w){
  //     std::pair<int,int> w_1D = std::make_pair(wwires.front()->index(),wwires.back()->index());
  //     if (proj_w.find(w_1D)!=proj_w.end()){
  // 	GeomCellSelection& cells = proj_w[w_1D];
  // 	cells.push_back(mcell);
  //     }else{
  // 	bool flag_save = true;
  // 	std::vector<std::pair<int,int>> to_be_removed;
  // 	for (auto it1 = proj_w.begin(); it1!=proj_w.end(); it1++){
  // 	  std::pair<int,int> testing = it1->first;
  // 	  if (w_1D.first >= testing.first && w_1D.second <= testing.second){
  // 	    flag_save = false;
  // 	    break;
  // 	  }else if (w_1D.first <= testing.first && w_1D.second >= testing.second){
  // 	    to_be_removed.push_back(testing);
  // 	  }
  // 	}
  // 	if (flag_save){
  // 	  GeomCellSelection cells;
  // 	  cells.push_back(mcell);
  // 	  proj_w[w_1D] = cells;
  // 	}
  // 	for (auto it1 = to_be_removed.begin(); it1!=to_be_removed.end(); it1++){
  // 	  proj_w.erase((*it1));
  // 	}
  //     }
  //   }
  // }
  

  // std::set<const GeomCell*> saved;
  // // examine ... 
  // for (auto it = proj_u.begin(); it!=proj_u.end(); it++){
  //   GeomCellSelection& cells = it->second;
  //   if (cells.size()==1)
  //     saved.insert(cells.front());
  // }
  // for (auto it = proj_v.begin(); it!=proj_v.end(); it++){
  //   GeomCellSelection& cells = it->second;
  //   if (cells.size()==1)
  //     saved.insert(cells.front());
  // }
  // for (auto it = proj_w.begin(); it!=proj_w.end(); it++){
  //   GeomCellSelection& cells = it->second;
  //   if (cells.size()==1)
  //     saved.insert(cells.front());
  // }

  // for (auto it = saved.begin(); it!=saved.end(); it++){
  //   SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);
  //   GeomWireSelection uwires = mcell->get_uwires();
  //   GeomWireSelection vwires = mcell->get_vwires();
  //   GeomWireSelection wwires = mcell->get_wwires();
    
  //   std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
  //   bool flag_u = true, flag_v = true, flag_w = true;
  //   for (auto it1 = bad_planes.begin(); it1!=bad_planes.end(); it1++){
  //     if (*it1 == WirePlaneType_t(0))
  // 	flag_u = false;
  //     if (*it1 == WirePlaneType_t(1))
  // 	flag_v = false;
  //     if (*it1 == WirePlaneType_t(2))
  // 	flag_w = false;
  //   }
    
  //   if (flag_u){
  //     for (auto it1 = uwires.begin(); it1!= uwires.end(); it1++){
  // 	used_uwires.insert(*it1);
  //     }
  //   }
  //   if (flag_v){
  //     for (auto it1 = vwires.begin(); it1!= vwires.end(); it1++){
  // 	used_vwires.insert(*it1);
  //     }
  //   }
  //   if (flag_w){
  //     for (auto it1 = wwires.begin(); it1!= wwires.end(); it1++){
  // 	used_wwires.insert(*it1);
  //     }
  //   }
  // }

  
  // for (auto it = two_good_wire_cells.begin(); it!= two_good_wire_cells.end(); it++){
  //   SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);
  //   if (good_mcells.find(mcell)==good_mcells.end() && saved.find(mcell)==saved.end()){
  //     GeomWireSelection uwires = mcell->get_uwires();
  //     GeomWireSelection vwires = mcell->get_vwires();
  //     GeomWireSelection wwires = mcell->get_wwires();

  //     std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
  //     bool flag_u = true, flag_v = true, flag_w = true;
  //     for (auto it1 = bad_planes.begin(); it1!=bad_planes.end(); it1++){
  // 	if (*it1 == WirePlaneType_t(0))
  // 	  flag_u = false;
  // 	if (*it1 == WirePlaneType_t(1))
  // 	  flag_v = false;
  // 	if (*it1 == WirePlaneType_t(2))
  // 	  flag_w = false;
  //     }

  //     bool save_flag = false;

  //     if (flag_u && !save_flag){
  // 	int nwire = 0;
  // 	int nwire_common = 0;
  // 	for (auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
  // 	  nwire++;
  // 	  if (used_uwires.find((*it1))!=used_uwires.end())
  // 	    nwire_common++;
  // 	}
  // 	if (nwire_common <0.8*  nwire) save_flag =true;
  //     }
      
  //     if (flag_v && !save_flag){
  // 	int nwire = 0;
  // 	int nwire_common = 0;
  // 	for (auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
  // 	  nwire++;
  // 	  if (used_vwires.find((*it1))!=used_vwires.end())
  // 	    nwire_common++;
  // 	}
  // 	if (nwire_common <0.8* nwire) save_flag =true;
  //     }
      
  //     if (flag_w && !save_flag){
  // 	int nwire = 0;
  // 	int nwire_common = 0;
  // 	for (auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
  // 	  nwire++;
  // 	  if (used_wwires.find((*it1))!=used_wwires.end())
  // 	    nwire_common++;
  // 	}
  // 	if (nwire_common <0.8* nwire) save_flag =true;
  //     }

  //     if (!save_flag) to_be_removed.push_back(mcell);
  //   }  
  // }

  
  
 
  
  
  

  // GeomCellSelection to_be_saved;
  // for (auto it = proj_u.begin(); it!=proj_u.end(); it++){
  //   //std::cout << "U: " << it->first.first << " " << it->first.second << std::endl;
  //   GeomCellSelection& cells = it->second;
  //   if (cells.size()>1){
  //     bool flag_saved = true;
  //     for (auto it1 = cells.begin(); it1!=cells.end(); it1++){
  // 	if (saved.find((*it1))!=saved.end()){
  // 	  flag_saved = false;
  // 	  break;
  // 	}
  //     }
  //     if (flag_saved){
  // 	for (auto it1 = cells.begin(); it1!=cells.end(); it1++){
  // 	  to_be_saved.push_back(*it1);
  // 	}
  //     }
  //   }
  // }
  // for (auto it = proj_v.begin(); it!=proj_v.end(); it++){
  //   //std::cout << "V: " << it->first.first << " " << it->first.second << std::endl;
  //   GeomCellSelection& cells = it->second;
  //   if (cells.size()>1){
  //     bool flag_saved = true;
  //     for (auto it1 = cells.begin(); it1!=cells.end(); it1++){
  // 	if (saved.find((*it1))!=saved.end()){
  // 	  flag_saved = false;
  // 	  break;
  // 	}
  //     }
  //     if (flag_saved){
  // 	for (auto it1 = cells.begin(); it1!=cells.end(); it1++){
  // 	  to_be_saved.push_back(*it1);
  // 	}
  //     }
  //   }
  // }
  // for (auto it = proj_w.begin(); it!=proj_w.end(); it++){
  //   //std::cout << "W: " << it->first.first << " " << it->first.second << std::endl;
  //   GeomCellSelection& cells = it->second;
  //   if (cells.size()>1){
  //     bool flag_saved = true;
  //     for (auto it1 = cells.begin(); it1!=cells.end(); it1++){
  // 	if (saved.find((*it1))!=saved.end()){
  // 	  flag_saved = false;
  // 	  break;
  // 	}
  //     }
  //     if (flag_saved){
  // 	for (auto it1 = cells.begin(); it1!=cells.end(); it1++){
  // 	  to_be_saved.push_back(*it1);
  // 	}
  //     }
  //   }
  // }

  // //std::cout << proj_u.size() << " " << proj_v.size() << " " << proj_w.size() << " " << two_good_wire_cells.size() << " " << saved.size() << " " << to_be_saved.size() << std::endl;
  
  // for (auto it = to_be_saved.begin(); it!=to_be_saved.end(); it++){
  //   saved.insert((*it));
  // }

  
}

void WireCell2dToy::LowmemTiling::init_good_cells(const WireCell::Slice& slice, const WireCell::Slice& slice_err, std::vector<float>& uplane_rms, std::vector<float>& vplane_rms, std::vector<float>& wplane_rms){
  // form good wires group
  form_fired_merge_wires(slice, slice_err);
  
  bool flag_two_wire_cell = true;
  bool flag_one_wire_cell = false;

  // create three good wire cells & two good wire + one bad wire cells
  // U/V/W = 1/1/1
  for (int i=0;i!=fired_wire_u.size();i++){
    MergeGeomWire *uwire = (MergeGeomWire *)fired_wire_u.at(i);
    for (int j=0;j!=fired_wire_v.size();j++){
      MergeGeomWire *vwire = (MergeGeomWire *)fired_wire_v.at(j);
      for (int k=0;k!=fired_wire_w.size();k++){
	MergeGeomWire *wwire = (MergeGeomWire *)fired_wire_w.at(k);
	SlimMergeGeomCell *mcell = create_slim_merge_cell(uwire,vwire,wwire);
	
	if (mcell !=0) {
	  three_good_wire_cells.push_back(mcell);
	  
	  // create new mwires 
	  GeomWireSelection uwires = mcell->get_uwires();
	  GeomWireSelection vwires = mcell->get_vwires();
	  GeomWireSelection wwires = mcell->get_wwires();
	  MergeGeomWire *mwire_u = new MergeGeomWire(0,uwires);
	  MergeGeomWire *mwire_v = new MergeGeomWire(0,vwires);
	  MergeGeomWire *mwire_w = new MergeGeomWire(0,wwires);

	  // create the map  
	  // cell to wires
	  GeomWireSelection wires;
	  wires.push_back(mwire_u);
	  wires.push_back(mwire_v);
	  wires.push_back(mwire_w);
	  cell_wires_map[mcell] = wires;
      
	  // wire to cells
	  GeomCellSelection cells;
	  cells.push_back(mcell);
	  wire_cells_map[mwire_u] = cells;
	  wire_cells_map[mwire_v] = cells;
	  wire_cells_map[mwire_w] = cells;
      
	  // wire to parent wire
	  wire_pwire_map[mwire_u] = uwire;
	  wire_pwire_map[mwire_v] = vwire;
	  wire_pwire_map[mwire_w] = wwire;
	  
	  // wire types
	  wire_type_map[mwire_u] = wire_type_map[uwire];
	  wire_type_map[mwire_v] = wire_type_map[vwire];     
	  wire_type_map[mwire_w] = wire_type_map[wwire];
	  

	  //std::cout <<  wire_type_map[mwire_u] << " " << wire_type_map[mwire_v] << " " << wire_type_map[mwire_w]  << std::endl;
	  
	  //parent wire to wires
	  if (pwire_wires_map.find(uwire)==pwire_wires_map.end()){
	    GeomWireSelection wires;
	    wires.push_back(mwire_u);
	    pwire_wires_map[uwire] = wires;
	  }else{
	    pwire_wires_map[uwire].push_back(mwire_u);
	  }
	  
	  if (pwire_wires_map.find(vwire)==pwire_wires_map.end()){
	    GeomWireSelection wires;
	    wires.push_back(mwire_v);
	    pwire_wires_map[vwire] = wires;
	  }else{
	    pwire_wires_map[vwire].push_back(mwire_v);
	  }
	  
	  if (pwire_wires_map.find(wwire)==pwire_wires_map.end()){
	    GeomWireSelection wires;
	    wires.push_back(mwire_w);
	    pwire_wires_map[wwire] = wires;
	  }else{
	    pwire_wires_map[wwire].push_back(mwire_w);
	  }

	  // end fill all maps
	}
      }
    }
  }
  
  if (flag_two_wire_cell){
    //U/V/W = 1/1/0
    for (int i=0;i!=fired_wire_u.size();i++){
      MergeGeomWire *uwire = (MergeGeomWire *)fired_wire_u.at(i);
      for (int j=0;j!=fired_wire_v.size();j++){
	MergeGeomWire *vwire = (MergeGeomWire *)fired_wire_v.at(j);
	for (int k=0;k!=bad_wire_w.size();k++){
	  MergeGeomWire *wwire = (MergeGeomWire *)bad_wire_w.at(k);
	  
	  SlimMergeGeomCell *mcell = create_slim_merge_cell(uwire,vwire,wwire);
	  
	  if (mcell !=0) {
	    two_good_wire_cells.push_back(mcell);
	    mcell->add_bad_planes(WirePlaneType_t(2));
	    
  	    // create new mwires 
	    GeomWireSelection uwires = mcell->get_uwires();
	    GeomWireSelection vwires = mcell->get_vwires();
	    GeomWireSelection wwires = mcell->get_wwires();
	    MergeGeomWire *mwire_u = new MergeGeomWire(0,uwires);
	    MergeGeomWire *mwire_v = new MergeGeomWire(0,vwires);
	    MergeGeomWire *mwire_w = new MergeGeomWire(0,wwires);
	    
	    // create the map  
	    // cell to wires
	    GeomWireSelection wires;
	    wires.push_back(mwire_u);
	    wires.push_back(mwire_v);
	    wires.push_back(mwire_w);
	    cell_wires_map[mcell] = wires;
	    
	    // wire to cells
	    GeomCellSelection cells;
	    cells.push_back(mcell);
	    wire_cells_map[mwire_u] = cells;
	    wire_cells_map[mwire_v] = cells;
	    wire_cells_map[mwire_w] = cells;
	    
	    // wire to parent wire
	    wire_pwire_map[mwire_u] = uwire;
	    wire_pwire_map[mwire_v] = vwire;
	    wire_pwire_map[mwire_w] = wwire;
	    
	    // wire types
	    wire_type_map[mwire_u] = wire_type_map[uwire];
	    wire_type_map[mwire_v] = wire_type_map[vwire];     
	    wire_type_map[mwire_w] = wire_type_map[wwire];
	    
	    //parent wire to wires
	    if (pwire_wires_map.find(uwire)==pwire_wires_map.end()){
	      GeomWireSelection wires;
	      wires.push_back(mwire_u);
	      pwire_wires_map[uwire] = wires;
	    }else{
	      pwire_wires_map[uwire].push_back(mwire_u);
	    }
	    
	    if (pwire_wires_map.find(vwire)==pwire_wires_map.end()){
	      GeomWireSelection wires;
	      wires.push_back(mwire_v);
	      pwire_wires_map[vwire] = wires;
	    }else{
	      pwire_wires_map[vwire].push_back(mwire_v);
	    }
	    
	    if (pwire_wires_map.find(wwire)==pwire_wires_map.end()){
	      GeomWireSelection wires;
	      wires.push_back(mwire_w);
	      pwire_wires_map[wwire] = wires;
	    }else{
	      pwire_wires_map[wwire].push_back(mwire_w);
	    }
	    
	    // end fill all maps
	  }
	}
      }
    }
    
  
    //U/V/W = 1/0/1
    for (int i=0;i!=fired_wire_u.size();i++){
      MergeGeomWire *uwire = (MergeGeomWire *)fired_wire_u.at(i);
      for (int j=0;j!=bad_wire_v.size();j++){
	MergeGeomWire *vwire = (MergeGeomWire *)bad_wire_v.at(j);
	for (int k=0;k!=fired_wire_w.size();k++){
	  MergeGeomWire *wwire = (MergeGeomWire *)fired_wire_w.at(k);
	  
	  SlimMergeGeomCell *mcell = create_slim_merge_cell(uwire,vwire,wwire);
	  
	  if (mcell !=0) {
	    two_good_wire_cells.push_back(mcell);
	    mcell->add_bad_planes(WirePlaneType_t(1));
  	    // create new mwires 
	    GeomWireSelection uwires = mcell->get_uwires();
	    GeomWireSelection vwires = mcell->get_vwires();
	    GeomWireSelection wwires = mcell->get_wwires();
	    MergeGeomWire *mwire_u = new MergeGeomWire(0,uwires);
	    MergeGeomWire *mwire_v = new MergeGeomWire(0,vwires);
	    MergeGeomWire *mwire_w = new MergeGeomWire(0,wwires);
	    
	    // create the map  
	    // cell to wires
	    GeomWireSelection wires;
	    wires.push_back(mwire_u);
	    wires.push_back(mwire_v);
	    wires.push_back(mwire_w);
	    cell_wires_map[mcell] = wires;
	    
	    // wire to cells
	    GeomCellSelection cells;
	    cells.push_back(mcell);
	    wire_cells_map[mwire_u] = cells;
	    wire_cells_map[mwire_v] = cells;
	    wire_cells_map[mwire_w] = cells;
	    
	    // wire to parent wire
	    wire_pwire_map[mwire_u] = uwire;
	    wire_pwire_map[mwire_v] = vwire;
	    wire_pwire_map[mwire_w] = wwire;
	    
	    // wire types
	    wire_type_map[mwire_u] = wire_type_map[uwire];
	    wire_type_map[mwire_v] = wire_type_map[vwire];     
	    wire_type_map[mwire_w] = wire_type_map[wwire];
	    
	    //parent wire to wires
	    if (pwire_wires_map.find(uwire)==pwire_wires_map.end()){
	      GeomWireSelection wires;
	      wires.push_back(mwire_u);
	      pwire_wires_map[uwire] = wires;
	    }else{
	      pwire_wires_map[uwire].push_back(mwire_u);
	    }
	    
	    if (pwire_wires_map.find(vwire)==pwire_wires_map.end()){
	      GeomWireSelection wires;
	      wires.push_back(mwire_v);
	      pwire_wires_map[vwire] = wires;
	    }else{
	      pwire_wires_map[vwire].push_back(mwire_v);
	    }
	    
	    if (pwire_wires_map.find(wwire)==pwire_wires_map.end()){
	      GeomWireSelection wires;
	      wires.push_back(mwire_w);
	      pwire_wires_map[wwire] = wires;
	    }else{
	      pwire_wires_map[wwire].push_back(mwire_w);
	    }
	    
	    // end fill all maps
	  }
	}
      }
    }
    
    //U/V/W = 0/1/1
    for (int i=0;i!=bad_wire_u.size();i++){
      MergeGeomWire *uwire = (MergeGeomWire *)bad_wire_u.at(i);
      for (int j=0;j!=fired_wire_v.size();j++){
	MergeGeomWire *vwire = (MergeGeomWire *)fired_wire_v.at(j);
	for (int k=0;k!=fired_wire_w.size();k++){
	  MergeGeomWire *wwire = (MergeGeomWire *)fired_wire_w.at(k);
	  
	  SlimMergeGeomCell *mcell = create_slim_merge_cell(uwire,vwire,wwire);
	  
	  if (mcell !=0) {
	    two_good_wire_cells.push_back(mcell);
	    mcell->add_bad_planes(WirePlaneType_t(0));
	    
  	    // create new mwires 
	    GeomWireSelection uwires = mcell->get_uwires();
	    GeomWireSelection vwires = mcell->get_vwires();
	    GeomWireSelection wwires = mcell->get_wwires();
	    MergeGeomWire *mwire_u = new MergeGeomWire(0,uwires);
	    MergeGeomWire *mwire_v = new MergeGeomWire(0,vwires);
	    MergeGeomWire *mwire_w = new MergeGeomWire(0,wwires);
	    
	    // create the map  
	    // cell to wires
	    GeomWireSelection wires;
	    wires.push_back(mwire_u);
	    wires.push_back(mwire_v);
	    wires.push_back(mwire_w);
	    cell_wires_map[mcell] = wires;
	    
	    // wire to cells
	    GeomCellSelection cells;
	    cells.push_back(mcell);
	    wire_cells_map[mwire_u] = cells;
	    wire_cells_map[mwire_v] = cells;
	    wire_cells_map[mwire_w] = cells;
	    
	    // wire to parent wire
	    wire_pwire_map[mwire_u] = uwire;
	    wire_pwire_map[mwire_v] = vwire;
	    wire_pwire_map[mwire_w] = wwire;
	    
	    // wire types
	    wire_type_map[mwire_u] = wire_type_map[uwire];
	    wire_type_map[mwire_v] = wire_type_map[vwire];     
	    wire_type_map[mwire_w] = wire_type_map[wwire];
	    
	    //parent wire to wires
	    if (pwire_wires_map.find(uwire)==pwire_wires_map.end()){
	      GeomWireSelection wires;
	      wires.push_back(mwire_u);
	      pwire_wires_map[uwire] = wires;
	    }else{
	      pwire_wires_map[uwire].push_back(mwire_u);
	    }
	    
	    if (pwire_wires_map.find(vwire)==pwire_wires_map.end()){
	      GeomWireSelection wires;
	      wires.push_back(mwire_v);
	      pwire_wires_map[vwire] = wires;
	    }else{
	      pwire_wires_map[vwire].push_back(mwire_v);
	    }
	    
	    if (pwire_wires_map.find(wwire)==pwire_wires_map.end()){
	      GeomWireSelection wires;
	      wires.push_back(mwire_w);
	      pwire_wires_map[wwire] = wires;
	    }else{
	      pwire_wires_map[wwire].push_back(mwire_w);
	    }
	    
	    // end fill all maps
	  }
	}
      }
    }
    
  }

  if (flag_one_wire_cell)
    create_one_good_wire_cells();

  //create map

}

WireCell::PointVector WireCell2dToy::LowmemTiling::get_all_cell_centers(){
  // PointVector pcells;
  // for (auto it = cell_wires_map.begin(); it!= cell_wires_map.end(); it++){
  //   SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)it->first;
  //   const GeomWire *uwire = mcell->get_uwires().at(int(mcell->get_uwires().size()/2));
  //   const GeomWire *vwire = mcell->get_vwires().at(int(mcell->get_vwires().size()/2));
  //   Vector abc;
  //   gds.crossing_point(*uwire, *vwire, abc);
  //   pcells.push_back(abc);
  // }

  return points;
}

GeomWireSelection WireCell2dToy::LowmemTiling::get_all_good_wires(){
  GeomWireSelection wires;
  
  for (auto it = wire_cells_map.begin(); it!= wire_cells_map.end(); it++){
    MergeGeomWire *mwire = (MergeGeomWire*)it->first;
    if (wire_type_map[mwire]){
      wires.push_back(mwire);
    }
  }
  
  return wires;
}


GeomWireSelection WireCell2dToy::LowmemTiling::get_all_bad_wires(){
  GeomWireSelection wires;
  
  for (auto it = wire_cells_map.begin(); it!= wire_cells_map.end(); it++){
    MergeGeomWire *mwire = (MergeGeomWire*)it->first;
    if (wire_type_map[mwire]){
    }else{
      wires.push_back(mwire);
    }
  }
  
  return wires;
}




void WireCell2dToy::LowmemTiling::create_one_good_wire_cells(){
   // figure out the fired wires for each plane 
  std::set<const GeomWire*> fired_wires;
  for (int i = 0; i!= three_good_wire_cells.size(); i++){
    GeomWireSelection uwires = ((SlimMergeGeomCell*)three_good_wire_cells.at(i))->get_uwires();
    GeomWireSelection vwires = ((SlimMergeGeomCell*)three_good_wire_cells.at(i))->get_vwires();
    GeomWireSelection wwires = ((SlimMergeGeomCell*)three_good_wire_cells.at(i))->get_wwires();
    for (int j=0;j!=uwires.size();j++){
      fired_wires.insert(uwires.at(j));
    }
    for (int j=0;j!=vwires.size();j++){
      fired_wires.insert(vwires.at(j));
    }
    for (int j=0;j!=wwires.size();j++){
      fired_wires.insert(wwires.at(j));
    }
  }
  
  for (int i=0; i!=two_good_wire_cells.size();i++){
    GeomWireSelection uwires = ((SlimMergeGeomCell*)two_good_wire_cells.at(i))->get_uwires();
    GeomWireSelection vwires = ((SlimMergeGeomCell*)two_good_wire_cells.at(i))->get_vwires();
    GeomWireSelection wwires = ((SlimMergeGeomCell*)two_good_wire_cells.at(i))->get_wwires();

    GeomWireSelection wires;
    std::vector<WirePlaneType_t> bad_planes = ((SlimMergeGeomCell*)two_good_wire_cells.at(i))->get_bad_planes();
    
    if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(0)) == bad_planes.end()){
      wires.insert(wires.end(),uwires.begin(),uwires.end());
    }
    if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(1)) == bad_planes.end()){
      wires.insert(wires.end(),vwires.begin(),vwires.end());
    }
    if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(2)) == bad_planes.end()){
      wires.insert(wires.end(),wwires.begin(),wwires.end());
    }
    for (int j=0;j!=wires.size();j++){
      fired_wires.insert(wires.at(j));
    }
  }

  // int count = 0;
  // for (int i =0; i!= fired_wire_u.size();i++){
  //   count += ((MergeGeomWire*)fired_wire_u.at(i))->get_allwire().size();
  // }
  // for (int i =0; i!= fired_wire_v.size();i++){
  //   count += ((MergeGeomWire*)fired_wire_v.at(i))->get_allwire().size();
  // }
  // for (int i =0; i!= fired_wire_w.size();i++){
  //   count += ((MergeGeomWire*)fired_wire_w.at(i))->get_allwire().size();
  // }
  // std::cout << fired_wires.size() << " " << count << std::endl;

  WireCell::GeomWireWireMap temp_wire_map;
  
  std::set<const GeomWire*> leftover_wires;
  //  GeomWireSelection leftover_wires;
  for (int i =0; i!= fired_wire_u.size();i++){
    for (int j=0; j!= ((MergeGeomWire*)fired_wire_u.at(i))->get_allwire().size(); j++){
      const GeomWire *wire = ((MergeGeomWire*)fired_wire_u.at(i))->get_allwire().at(j);
      if (fired_wires.find(wire)==fired_wires.end()){
	leftover_wires.insert(wire);
	temp_wire_map[wire] = (MergeGeomWire*)fired_wire_u.at(i);
      }
    }
  }
  for (int i =0; i!= fired_wire_v.size();i++){
    for (int j=0; j!= ((MergeGeomWire*)fired_wire_v.at(i))->get_allwire().size(); j++){
      const GeomWire *wire = ((MergeGeomWire*)fired_wire_v.at(i))->get_allwire().at(j);
      if (fired_wires.find(wire)==fired_wires.end()){
	leftover_wires.insert(wire);
	temp_wire_map[wire] = (MergeGeomWire*)fired_wire_v.at(i);
      }
    }
  }
  for (int i =0; i!= fired_wire_w.size();i++){
    for (int j=0; j!= ((MergeGeomWire*)fired_wire_w.at(i))->get_allwire().size(); j++){
      const GeomWire *wire = ((MergeGeomWire*)fired_wire_w.at(i))->get_allwire().at(j);
      if (fired_wires.find(wire)==fired_wires.end()){
	leftover_wires.insert(wire);
	temp_wire_map[wire] = (MergeGeomWire*)fired_wire_w.at(i);
      }
    }
  }

  
  //std::cout << leftover_wires.size() << std::endl;

  GeomWireSelection remaining_fired_wire_u, remaining_fired_wire_v, remaining_fired_wire_w;
  // do U
  MergeGeomWire *mwire = 0;
  int ident = 0;
  int last_wire = -1;
  for (int i = 0; i!= nwire_u; i++){
    const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)0,i);
    if (leftover_wires.find(wire)!=leftover_wires.end()){
      if (mwire == 0){
   	mwire = new MergeGeomWire(ident,*wire);
	wire_type_map[mwire] = true;
	temp_wire_map[mwire] = temp_wire_map[wire];
	remaining_fired_wire_u.push_back(mwire);
   	ident ++;
      }else{
   	if (i==last_wire+1){
   	  mwire->AddWire(*wire);
   	}else{
   	  // create a new wire
   	  mwire = new MergeGeomWire(ident,*wire);
    	  wire_type_map[mwire] = true;
	  temp_wire_map[mwire] = temp_wire_map[wire];
	  remaining_fired_wire_u.push_back(mwire);
   	  ident ++;
   	}
      }
      last_wire = i;
    }
  }
  
  // V
  mwire = 0;
  last_wire = -1;
  for (int i = 0; i!= nwire_v; i++){
    const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)1,i);
    if (leftover_wires.find(wire)!=leftover_wires.end()){
      if (mwire == 0){
  	mwire = new MergeGeomWire(ident,*wire);
  	//mwire->SetTimeSlice(time_slice);
	wire_type_map[mwire] = true;
	temp_wire_map[mwire] = temp_wire_map[wire];
  	remaining_fired_wire_v.push_back(mwire);
  	ident ++;
      }else{
  	if (i==last_wire+1){
  	  mwire->AddWire(*wire);
  	}else{
  	  // create a new wire
  	  mwire = new MergeGeomWire(ident,*wire);
  	  //mwire->SetTimeSlice(time_slice);
	  wire_type_map[mwire] = true;
	  temp_wire_map[mwire] = temp_wire_map[wire];
  	  remaining_fired_wire_v.push_back(mwire);
  	  ident ++;
  	}
      }
      last_wire = i;
    }
  }
  
  // W
  mwire = 0;
  last_wire = -1;
  for (int i = 0; i!= nwire_w; i++){
    const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)2,i);
    if (leftover_wires.find(wire)!=leftover_wires.end()){
      if (mwire == 0){
  	mwire = new MergeGeomWire(ident,*wire);
  	//mwire->SetTimeSlice(time_slice);
	wire_type_map[mwire] = true;
	temp_wire_map[mwire] = temp_wire_map[wire];
  	remaining_fired_wire_w.push_back(mwire);
  	ident ++;
      }else{
  	if (i==last_wire+1){
  	  mwire->AddWire(*wire);
  	}else{
  	  // create a new wire
  	  mwire = new MergeGeomWire(ident,*wire);
  	  //mwire->SetTimeSlice(time_slice);
	  wire_type_map[mwire] = true;
	  temp_wire_map[mwire] = temp_wire_map[wire];
  	  remaining_fired_wire_w.push_back(mwire);
  	  ident ++;
  	}
      }
      last_wire = i;
    }
  }

  // int count = 0;
  // for (int i =0; i!= remaining_fired_wire_u.size();i++){
  //   count += ((MergeGeomWire*)remaining_fired_wire_u.at(i))->get_allwire().size();
  // }
  // for (int i =0; i!= remaining_fired_wire_v.size();i++){
  //   count += ((MergeGeomWire*)remaining_fired_wire_v.at(i))->get_allwire().size();
  // }
  // for (int i =0; i!= remaining_fired_wire_w.size();i++){
  //   count += ((MergeGeomWire*)remaining_fired_wire_w.at(i))->get_allwire().size();
  // }
  //std::cout << remaining_fired_wire_u.size() << " " << remaining_fired_wire_v.size() << " " << remaining_fired_wire_w.size() << " " << 0 << std::endl;

  
  GeomCellSelection temp_cells;
  //U/V/W = 1/0/0
  for (int i=0;i!=remaining_fired_wire_u.size();i++){
    MergeGeomWire *uwire = (MergeGeomWire *)remaining_fired_wire_u.at(i);
    for (int j=0;j!=bad_wire_v.size();j++){
      MergeGeomWire *vwire = (MergeGeomWire *)bad_wire_v.at(j);
      for (int k=0;k!=bad_wire_w.size();k++){
  	MergeGeomWire *wwire = (MergeGeomWire *)bad_wire_w.at(k);
	
  	SlimMergeGeomCell *mcell = create_slim_merge_cell(uwire,vwire,wwire);
	
  	if (mcell !=0) {
	  temp_cells.push_back(mcell);
	  mcell->add_bad_planes(WirePlaneType_t(1));
	  mcell->add_bad_planes(WirePlaneType_t(2));

	  // create new mwires 
	  GeomWireSelection uwires = mcell->get_uwires();
	  GeomWireSelection vwires = mcell->get_vwires();
	  GeomWireSelection wwires = mcell->get_wwires();
	  MergeGeomWire *mwire_u = new MergeGeomWire(0,uwires);
	  MergeGeomWire *mwire_v = new MergeGeomWire(0,vwires);
	  MergeGeomWire *mwire_w = new MergeGeomWire(0,wwires);

	  // create the map  
	  // cell to wires
	  GeomWireSelection wires;
	  wires.push_back(mwire_u);
	  wires.push_back(mwire_v);
	  wires.push_back(mwire_w);
	  cell_wires_map[mcell] = wires;
      
	  // wire to cells
	  GeomCellSelection cells;
	  cells.push_back(mcell);
	  wire_cells_map[mwire_u] = cells;
	  wire_cells_map[mwire_v] = cells;
	  wire_cells_map[mwire_w] = cells;
      
	  // wire to parent wire
	  wire_pwire_map[mwire_u] = temp_wire_map[uwire];
	  wire_pwire_map[mwire_v] = vwire;
	  wire_pwire_map[mwire_w] = wwire;
	  
	  // wire types
	  wire_type_map[mwire_u] = wire_type_map[temp_wire_map[uwire]];
	  wire_type_map[mwire_v] = wire_type_map[vwire];     
	  wire_type_map[mwire_w] = wire_type_map[wwire];
	  
	  //parent wire to wires
	  if (pwire_wires_map.find(temp_wire_map[uwire])==pwire_wires_map.end()){
	    GeomWireSelection wires;
	    wires.push_back(mwire_u);
	    pwire_wires_map[temp_wire_map[uwire]] = wires;
	  }else{
	    pwire_wires_map[temp_wire_map[uwire]].push_back(mwire_u);
	  }
	  
	  if (pwire_wires_map.find(vwire)==pwire_wires_map.end()){
	    GeomWireSelection wires;
	    wires.push_back(mwire_v);
	    pwire_wires_map[vwire] = wires;
	  }else{
	    pwire_wires_map[vwire].push_back(mwire_v);
	  }
	  
	  if (pwire_wires_map.find(wwire)==pwire_wires_map.end()){
	    GeomWireSelection wires;
	    wires.push_back(mwire_w);
	    pwire_wires_map[wwire] = wires;
	  }else{
	    pwire_wires_map[wwire].push_back(mwire_w);
	  }
	  

	}
      }
    }
  }
  
  //U/V/W = 0/1/0
  for (int i=0;i!=bad_wire_u.size();i++){
    MergeGeomWire *uwire = (MergeGeomWire *)bad_wire_u.at(i);
    for (int j=0;j!=remaining_fired_wire_v.size();j++){
      MergeGeomWire *vwire = (MergeGeomWire *)remaining_fired_wire_v.at(j);
      for (int k=0;k!=bad_wire_w.size();k++){
  	MergeGeomWire *wwire = (MergeGeomWire *)bad_wire_w.at(k);
	
  	SlimMergeGeomCell *mcell = create_slim_merge_cell(uwire,vwire,wwire);
	
  	if (mcell !=0) {
	  temp_cells.push_back(mcell);
	  mcell->add_bad_planes(WirePlaneType_t(0));
	  mcell->add_bad_planes(WirePlaneType_t(2));

	    // create new mwires 
	  GeomWireSelection uwires = mcell->get_uwires();
	  GeomWireSelection vwires = mcell->get_vwires();
	  GeomWireSelection wwires = mcell->get_wwires();
	  MergeGeomWire *mwire_u = new MergeGeomWire(0,uwires);
	  MergeGeomWire *mwire_v = new MergeGeomWire(0,vwires);
	  MergeGeomWire *mwire_w = new MergeGeomWire(0,wwires);

	  // create the map  
	  // cell to wires
	  GeomWireSelection wires;
	  wires.push_back(mwire_u);
	  wires.push_back(mwire_v);
	  wires.push_back(mwire_w);
	  cell_wires_map[mcell] = wires;
      
	  // wire to cells
	  GeomCellSelection cells;
	  cells.push_back(mcell);
	  wire_cells_map[mwire_u] = cells;
	  wire_cells_map[mwire_v] = cells;
	  wire_cells_map[mwire_w] = cells;
      
	  // wire to parent wire
	  wire_pwire_map[mwire_u] = uwire;
	  wire_pwire_map[mwire_v] = temp_wire_map[vwire];
	  wire_pwire_map[mwire_w] = wwire;
	  
	  // wire types
	  wire_type_map[mwire_u] = wire_type_map[uwire];
	  wire_type_map[mwire_v] = wire_type_map[temp_wire_map[vwire]];     
	  wire_type_map[mwire_w] = wire_type_map[wwire];
	  
	  //parent wire to wires
	  if (pwire_wires_map.find(uwire)==pwire_wires_map.end()){
	    GeomWireSelection wires;
	    wires.push_back(mwire_u);
	    pwire_wires_map[uwire] = wires;
	  }else{
	    pwire_wires_map[uwire].push_back(mwire_u);
	  }
	  
	  if (pwire_wires_map.find(temp_wire_map[vwire])==pwire_wires_map.end()){
	    GeomWireSelection wires;
	    wires.push_back(mwire_v);
	    pwire_wires_map[temp_wire_map[vwire]] = wires;
	  }else{
	    pwire_wires_map[temp_wire_map[vwire]].push_back(mwire_v);
	  }
	  
	  if (pwire_wires_map.find(wwire)==pwire_wires_map.end()){
	    GeomWireSelection wires;
	    wires.push_back(mwire_w);
	    pwire_wires_map[wwire] = wires;
	  }else{
	    pwire_wires_map[wwire].push_back(mwire_w);
	  }
	}
      }
    }
  }
  
  //U/V/W = 0/0/1
  for (int i=0;i!=bad_wire_u.size();i++){
    MergeGeomWire *uwire = (MergeGeomWire *)bad_wire_u.at(i);
    for (int j=0;j!=bad_wire_v.size();j++){
      MergeGeomWire *vwire = (MergeGeomWire *)bad_wire_v.at(j);
      for (int k=0;k!=remaining_fired_wire_w.size();k++){
  	MergeGeomWire *wwire = (MergeGeomWire *)remaining_fired_wire_w.at(k);
	
  	SlimMergeGeomCell *mcell = create_slim_merge_cell(uwire,vwire,wwire);
  	
	if (mcell !=0) {
	  temp_cells.push_back(mcell);
	  mcell->add_bad_planes(WirePlaneType_t(0));
	  mcell->add_bad_planes(WirePlaneType_t(1));

	    // create new mwires 
	  GeomWireSelection uwires = mcell->get_uwires();
	  GeomWireSelection vwires = mcell->get_vwires();
	  GeomWireSelection wwires = mcell->get_wwires();
	  MergeGeomWire *mwire_u = new MergeGeomWire(0,uwires);
	  MergeGeomWire *mwire_v = new MergeGeomWire(0,vwires);
	  MergeGeomWire *mwire_w = new MergeGeomWire(0,wwires);

	  // create the map  
	  // cell to wires
	  GeomWireSelection wires;
	  wires.push_back(mwire_u);
	  wires.push_back(mwire_v);
	  wires.push_back(mwire_w);
	  cell_wires_map[mcell] = wires;
      
	  // wire to cells
	  GeomCellSelection cells;
	  cells.push_back(mcell);
	  wire_cells_map[mwire_u] = cells;
	  wire_cells_map[mwire_v] = cells;
	  wire_cells_map[mwire_w] = cells;
      
	  // wire to parent wire
	  wire_pwire_map[mwire_u] = uwire;
	  wire_pwire_map[mwire_v] = vwire;
	  wire_pwire_map[mwire_w] = temp_wire_map[wwire];
	  
	  // wire types
	  wire_type_map[mwire_u] = wire_type_map[uwire];
	  wire_type_map[mwire_v] = wire_type_map[vwire];     
	  wire_type_map[mwire_w] = wire_type_map[temp_wire_map[wwire]];
	  
	  //parent wire to wires
	  if (pwire_wires_map.find(uwire)==pwire_wires_map.end()){
	    GeomWireSelection wires;
	    wires.push_back(mwire_u);
	    pwire_wires_map[uwire] = wires;
	  }else{
	    pwire_wires_map[uwire].push_back(mwire_u);
	  }
	  
	  if (pwire_wires_map.find(vwire)==pwire_wires_map.end()){
	    GeomWireSelection wires;
	    wires.push_back(mwire_v);
	    pwire_wires_map[vwire] = wires;
	  }else{
	    pwire_wires_map[vwire].push_back(mwire_v);
	  }
	  
	  if (pwire_wires_map.find(temp_wire_map[wwire])==pwire_wires_map.end()){
	    GeomWireSelection wires;
	    wires.push_back(mwire_w);
	    pwire_wires_map[temp_wire_map[wwire]] = wires;
	  }else{
	    pwire_wires_map[temp_wire_map[wwire]].push_back(mwire_w);
	  }

	}
      }
    }
  }

  
  // if connected to one of the existing cell, add them, or remove ... 
  for (int i=0;i!=temp_cells.size();i++){
    SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)temp_cells.at(i);
    int flag = 0;

    for (int j=0;j!=three_good_wire_cells.size();j++){
      SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*)three_good_wire_cells.at(j);
      if (mcell->Overlap(mcell1)){
	flag = 1;
	break;
      }
    }
    if (flag==0){
      for (int j=0;j!=two_good_wire_cells.size();j++){
	SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*)two_good_wire_cells.at(j);
	if (mcell->Overlap(mcell1)){
	  flag = 1;
	  break;
	}
      }
    }
   
    if (flag==1) {
      one_good_wire_cells.push_back(mcell);
    }else{
      remove_cell(mcell);
      //not_used_one_good_wire_cells.push_back(mcell);
    }
  }
  
  //for (int i=0;i!=not_used_one_good_wire_cells.size();i++){
    // remove_cell((SlimMergeGeomCell*)not_used_one_good_wire_cells.at(i));
  // delete not_used_one_good_wire_cells.at(i);
  // }
  
  //std::cout << not_used_one_good_wire_cells.size() << std::endl;
  // if (one_good_wire_cells.size()!=0){
  //   std::cout << one_good_wire_cells.size() << std::endl;
  // }
}


GeomCellSelection WireCell2dToy::LowmemTiling::create_single_cells(){
  GeomCellSelection cells;
  
  // // loop three wire cells
  // for (int i = 0; i!= three_good_wire_cells.size(); i++){
  //   GeomCellSelection temp_cells = create_single_cells((SlimMergeGeomCell*)three_good_wire_cells.at(i));
  //   cells.insert(cells.end(),temp_cells.begin(),temp_cells.end());
    
  //   points.push_back(temp_cells.at(int(temp_cells.size()/2))->center());
  // }

  // // loop two wire cells
  // for (int i=0; i!=two_good_wire_cells.size();i++){
  //   GeomCellSelection temp_cells = create_single_cells((SlimMergeGeomCell*)two_good_wire_cells.at(i));
  //   cells.insert(cells.end(),temp_cells.begin(),temp_cells.end());

  //   points.push_back(temp_cells.at(int(temp_cells.size()/2))->center());
  // }

  // // loop one wire cells
  // for (int i=0; i!=one_good_wire_cells.size();i++){
  //   GeomCellSelection temp_cells = create_single_cells((SlimMergeGeomCell*)one_good_wire_cells.at(i));
  //   cells.insert(cells.end(),temp_cells.begin(),temp_cells.end());

  //   points.push_back(temp_cells.at(int(temp_cells.size()/2))->center());
  // }

  for (auto it = cell_wires_map.begin(); it!= cell_wires_map.end(); it++){
    GeomCellSelection temp_cells = create_single_cells((SlimMergeGeomCell*)it->first);
    cells.insert(cells.end(),temp_cells.begin(),temp_cells.end());
    points.push_back(temp_cells.at(int(temp_cells.size()/2))->center());
  }
  
  

  return cells;
}


GeomCellSelection WireCell2dToy::LowmemTiling::create_single_cells(SlimMergeGeomCell *mcell){
  GeomCellSelection cells;
  GeomWireSelection wire_u = mcell->get_uwires();
  GeomWireSelection wire_v = mcell->get_vwires();
  GeomWireSelection wire_w = mcell->get_wwires();
  float tolerance = 0.1 * units::mm;
  float dis_u[3]={0.0},dis_v[3]={0.0},dis_w[3]={0.0},dis_puv[5]={0.0},dis_puw[5]={0.0},dis_pwv[5]={0.0};
  int ncell = 1;
  
  std::vector<float> udis,vdis,wdis;
  for (int i=0;i!=wire_u.size();i++){
    udis.push_back(gds.wire_dist(*wire_u[i]));
  }
  for (int j=0;j!=wire_v.size();j++){
    vdis.push_back(gds.wire_dist(*wire_v[j]));
  }
  for (int k=0;k!=wire_w.size();k++){
    wdis.push_back(gds.wire_dist(*wire_w[k]));
  }
					
  std::vector<Vector> puv_save(5), puw_save(5), pwv_save(5);
  float dis_puv_save[5]={0.0}, dis_puw_save[5]={0.0}, dis_pwv_save[5]={0.0};
  
  double u_pitch=0, v_pitch=0, w_pitch=0;
  u_pitch = gds.pitch(kUwire);
  v_pitch = gds.pitch(kVwire);
  w_pitch = gds.pitch(kYwire);

  if (wire_u.size()>=1&&wire_v.size()>=1){
    dis_u[0] = udis.at(0) - gds.pitch(kUwire)/2.;
    dis_u[1] = dis_u[0] + gds.pitch(kUwire);
    dis_u[2] = udis.at(0);

    dis_v[0] = vdis.at(0) - gds.pitch(kVwire)/2.;
    dis_v[1] = dis_v[0] + gds.pitch(kVwire);
    dis_v[2] = vdis.at(0);

    gds.crossing_point(dis_u[0],dis_v[0],kUwire,kVwire, puv_save[0]);
    gds.crossing_point(dis_u[0],dis_v[1],kUwire,kVwire, puv_save[1]);
    gds.crossing_point(dis_u[1],dis_v[1],kUwire,kVwire, puv_save[2]);
    gds.crossing_point(dis_u[1],dis_v[0],kUwire,kVwire, puv_save[3]);
    gds.crossing_point(dis_u[2],dis_v[2],kUwire,kVwire, puv_save[4]);

    for (int k=0;k!=5;k++){
      dis_puv_save[k] = gds.wire_dist(puv_save[k],kYwire);
    }
    
  }

  if (wire_u.size()>=1&&wire_w.size()>=1){
    dis_u[0] = udis.at(0) - gds.pitch(kUwire)/2.;
    dis_u[1] = dis_u[0] + gds.pitch(kUwire);
    dis_u[2] = udis.at(0);

    dis_w[0] = wdis.at(0) - w_pitch/2.;
    dis_w[1] = dis_w[0] + w_pitch;
    dis_w[2] = wdis.at(0);

    gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, puw_save[0]);
    gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, puw_save[1]);
    gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, puw_save[2]);
    gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, puw_save[3]);
    gds.crossing_point(dis_u[2],dis_w[2],kUwire,kYwire, puw_save[4]);

    for (int k=0;k!=5;k++){
      dis_puw_save[k] = gds.wire_dist(puw_save[k],kVwire);
    }

  }

  if (wire_v.size()>=1&&wire_w.size()>=1){

    dis_w[0] = wdis.at(0) - w_pitch/2.;
    dis_w[1] = dis_w[0] + w_pitch;
    dis_w[2] = wdis.at(0);

    dis_v[0] = vdis.at(0) - gds.pitch(kVwire)/2.;
    dis_v[1] = dis_v[0] + gds.pitch(kVwire);
    dis_v[2] = vdis.at(0);
    
    gds.crossing_point(dis_v[0],dis_w[0],kVwire,kYwire, pwv_save[0]);
    gds.crossing_point(dis_v[0],dis_w[1],kVwire,kYwire, pwv_save[1]);
    gds.crossing_point(dis_v[1],dis_w[1],kVwire,kYwire, pwv_save[2]);
    gds.crossing_point(dis_v[1],dis_w[0],kVwire,kYwire, pwv_save[3]);
    gds.crossing_point(dis_v[2],dis_w[2],kVwire,kYwire, pwv_save[4]);
    
    for (int k=0;k!=5;k++){
      dis_pwv_save[k] = gds.wire_dist(pwv_save[k],kUwire);
    }
  }


   for (int i=0;i!=wire_u.size();i++){
    dis_u[0] = udis.at(i) - u_pitch/2.;
    dis_u[1] = dis_u[0] + u_pitch;
    dis_u[2] = udis.at(i);

    for (int j=0;j!=wire_v.size();j++){
      //  if (wirechargemap[wire_u[i]] <100 && wirechargemap[wire_v[j]] <100 ) continue;
      dis_v[0] = vdis.at(j) - v_pitch/2.;
      dis_v[1] = dis_v[0] + v_pitch;
      dis_v[2] = vdis.at(j);
      
      //four vertices around
      //      std::vector<Vector> puv(4);
      std::vector<Vector> puv(5);
      
      //counter ++;

      if(!gds.crossing_point(dis_u[2],dis_v[2],kUwire,kVwire, puv[4])) continue;
      
      // Point pp;
      // pp.x = 0; pp.y = puv[4].y; pp.z = puv[4].z;
      // if (!gds.contained_yz(pp)) continue;

      //counter1 ++;
      dis_puv[4] = gds.wire_dist(puv[4],kYwire);
      for (int k=0;k!=4;k++){
	puv[k] = puv[4] + (puv_save[k]-puv_save[4]);
	dis_puv[k] = dis_puv[4] +(dis_puv_save[k]-dis_puv_save[4]);
      }

      for (int k=0;k!=wire_w.size();k++){

	int flag = 0;
	PointVector pcell;
	dis_w[0] = wdis.at(k) - w_pitch/2.;
  	dis_w[1] = dis_w[0] + w_pitch;//gds.wire_dist(*wire_w[k]) + gds.pitch(kYwire)/2.;	
	dis_w[2] = wdis.at(k);

	if (fabs(dis_w[0] - dis_puv[0])>3*units::cm) continue;
	
	//counter3++;
	
	for (int m = 0;m!=4;m++){
	  if (dis_puv[m] > dis_w[0]-tolerance/2. && dis_puv[m] < dis_w[1]+tolerance/2.){
	    flag = 1;
	    pcell.push_back(puv[m]);
	  }
	}


	if (flag==1 ) {
	  //counter ++;
	  std::vector<Vector> puw(5);
	  
	  gds.crossing_point(dis_u[2],dis_w[2],kUwire,kYwire, puw[4]);
	  dis_puw[4] = gds.wire_dist(puw[4],kVwire);

	  for (int k1=0;k1!=4;k1++){
	    puw[k1] = puw[4] + (puw_save[k1]-puw_save[4]);
	    dis_puw[k1] = dis_puw[4] +(dis_puw_save[k1]-dis_puw_save[4]);
	    if (dis_puw[k1] > dis_v[0]-tolerance/2. && dis_puw[k1] < dis_v[1]+tolerance/2.){
	      int flag_abc = 0;
	      for (int kk = 0; kk!=pcell.size();kk++){
	      	float dis = sqrt(pow(puw[k1].y-pcell.at(kk).y,2) + pow(puw[k1].z-pcell.at(kk).z,2));
	      	if (dis < tolerance) {
	      	  flag_abc = 1;
	      	  break;
	      	}
	      }
	      if (flag_abc == 0)
	    	pcell.push_back(puw[k1]);
	    }
	  }
	  std::vector<Vector> pwv(5);
	  gds.crossing_point(dis_v[2],dis_w[2],kVwire,kYwire, pwv[4]);
	  dis_pwv[4] = gds.wire_dist(pwv[4],kUwire);

	  for (int k1=0;k1!=4;k1++){
	    pwv[k1] = pwv[4] + (pwv_save[k1]-pwv_save[4]);
	    dis_pwv[k1] = dis_pwv[4] +(dis_pwv_save[k1]-dis_pwv_save[4]);
	    //	    dis_pwv[k] = gds.wire_dist(pwv[k],kUwire);
	    if (dis_pwv[k1] > dis_u[0]-tolerance/2. && dis_pwv[k1] < dis_u[1]+tolerance/2.){
	      int flag_abc = 0;
	      for (int kk = 0; kk!=pcell.size();kk++){
	      	float dis = sqrt(pow(pwv[k1].y-pcell.at(kk).y,2) + pow(pwv[k1].z-pcell.at(kk).z,2));
	      	if (dis < tolerance) {
	      	  flag_abc = 1;
	      	  break;
	      	}
	      }
	      if (flag_abc == 0)
	      	pcell.push_back(pwv[k1]);
	    }
	  }

	  
	  if (pcell.size()>=3){
	    //order all the points by phi angle
	    
	    const GeomCell *cell = 0;
	    
	    // old method
	    GeomCell *cell_t = new GeomCell(ncell,pcell);
	    cell_t->set_uwire(wire_u.at(i));
	    cell_t->set_vwire(wire_v.at(j));
	    cell_t->set_wwire(wire_w.at(k));
	    cell = cell_t;

	  

	    Point cell_center = cell->center();
	    if (gds.contained_yz(cell_center)){
	      
	      cells.push_back(cell);
	      ncell++;
	    }else{
	      delete cell;
	    }
	  }
	}
	
      } // W-loop
    } // V-loop
  } // U-loop
  


  
  return cells;
}

void WireCell2dToy::LowmemTiling::check_bad_cells(WireCell2dToy::LowmemTiling* tiling,WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map){
  
  int prev_time_slice = tiling->get_time_slice();

  int flag = 0;
  for (auto it = uplane_map.begin(); it!= uplane_map.end(); it++){
    int start_time_slice = it->second.first/ nrebin;
    int end_time_slice = it->second.second/ nrebin;
    if (prev_time_slice < start_time_slice && time_slice >= start_time_slice ||
	prev_time_slice <= end_time_slice && time_slice > end_time_slice )
      {
	flag = 1;
	break;
      }
    //    std::cout << start_time_slice << " " << end_time_slice << std::endl;
  }
  if (flag==0){
    for (auto it = vplane_map.begin(); it!= vplane_map.end(); it++){
      int start_time_slice = it->second.first/ nrebin;
      int end_time_slice = it->second.second/ nrebin;
      if (prev_time_slice < start_time_slice && time_slice >= start_time_slice ||
	prev_time_slice <= end_time_slice && time_slice > end_time_slice
	  ){
	flag = 1;
	break;
      }
      //    std::cout << start_time_slice << " " << end_time_slice << std::endl;
    }
  }
  if (flag==0){
    for (auto it = wplane_map.begin(); it!= wplane_map.end(); it++){
      int start_time_slice = it->second.first/ nrebin;
      int end_time_slice = it->second.second/ nrebin;
      if (prev_time_slice < start_time_slice && time_slice >= start_time_slice ||
	prev_time_slice <= end_time_slice && time_slice > end_time_slice
	  ){
	flag = 1;
	break;
      }
      //    std::cout << start_time_slice << " " << end_time_slice << std::endl;
    }
  }

  

  if (flag==1){
    std::cout << "Regenerate bad cells at time slice " << time_slice << std::endl;
    init_bad_cells(uplane_map,vplane_map,wplane_map);
  }else{
    //copy wires 
    for (int i=0;i!=tiling->get_bad_wire_u().size();i++){
      MergeGeomWire *mwire1 = (MergeGeomWire*)tiling->get_bad_wire_u().at(i);
      bad_wire_u.push_back(mwire1);
    }
    for (int i=0;i!=tiling->get_bad_wire_v().size();i++){
      MergeGeomWire *mwire1 = (MergeGeomWire*)tiling->get_bad_wire_v().at(i);
      bad_wire_v.push_back(mwire1);
    }
    for (int i=0;i!=tiling->get_bad_wire_w().size();i++){
      MergeGeomWire *mwire1 = (MergeGeomWire*)tiling->get_bad_wire_w().at(i);
      bad_wire_w.push_back(mwire1);
    }
    // copy cells;
    for (int i=0;i!=tiling->get_two_bad_wire_cells().size();i++){
      SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*)tiling->get_two_bad_wire_cells().at(i);
      two_bad_wire_cells.push_back(mcell1);
    }
  }
  
  // std::cout << bad_wire_u.size() << " " << bad_wire_v.size() << " " << bad_wire_w.size() << std::endl;

  
}



void WireCell2dToy::LowmemTiling::init_bad_cells(WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map){
  // form bad wires group
  form_bad_merge_wires(uplane_map, vplane_map, wplane_map);
  // create two bad wire cells // these are special ones ... 
  form_two_bad_cells();

  //  std::cout << two_bad_wire_cells.size() << std::endl; 
	    
}



void WireCell2dToy::LowmemTiling::form_two_bad_cells(){
  // form two bad wire cells ... taken from BadTiling ...  
  
  // U-V and insert Y ... 
  for (int i =0; i!=bad_wire_u.size();i++){
    const GeomWire *uwire_1 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().front();
    const GeomWire *uwire_2 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().back();
    float dis_u[3];
    float u_pitch = gds.pitch(kUwire);
    dis_u[0] = gds.wire_dist(*uwire_1) - u_pitch/2.;
    dis_u[1] = gds.wire_dist(*uwire_2) + u_pitch/2.;
    for (int j=0; j!=bad_wire_v.size();j++){
      const GeomWire *vwire_1 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().front();
      const GeomWire *vwire_2 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().back();
      float dis_v[3];
      float v_pitch = gds.pitch(kVwire);
      dis_v[0] = gds.wire_dist(*vwire_1) - v_pitch/2.;
      dis_v[1] = gds.wire_dist(*vwire_2) + v_pitch/2.;
      
      PointVector pcell;
      std::vector<Vector> pcross(4);

      bool flag1 = gds.crossing_point(dis_u[0],dis_v[0],kUwire,kVwire, pcross[0]); // check the inner point
      
      if (flag1){
	// fill the outer point
  	pcell.push_back(pcross[0]);
      }else{
  	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) - u_pitch/2.;

	  if (gds.crossing_point(dis_u[2],dis_v[0],kUwire,kVwire, pcross[0])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(k);
	  dis_v[2] = gds.wire_dist(*vwire_3) - v_pitch/2.;
	  if (gds.crossing_point(dis_u[0],dis_v[2],kUwire,kVwire, pcross[0])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
      }

      bool flag2 = gds.crossing_point(dis_u[0],dis_v[1],kUwire,kVwire, pcross[1]); 
      if (flag2){
  	pcell.push_back(pcross[1]);
      }else{
  	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) - u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_v[1],kUwire,kVwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size()-1-k);
	  dis_v[2] = gds.wire_dist(*vwire_3) + v_pitch/2.;
	  if (gds.crossing_point(dis_u[0],dis_v[2],kUwire,kVwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
      }
      
      bool flag3 = gds.crossing_point(dis_u[1],dis_v[0],kUwire,kVwire, pcross[2]); 
      if (flag3){
  	pcell.push_back(pcross[2]);
      }else{
  	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size()-1-k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) + u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_v[0],kUwire,kVwire, pcross[2])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(k);
	  dis_v[2] = gds.wire_dist(*vwire_3) - v_pitch/2.;
	  if (gds.crossing_point(dis_u[1],dis_v[2],kUwire,kVwire, pcross[2])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
      }
      
      bool flag4 = gds.crossing_point(dis_u[1],dis_v[1],kUwire,kVwire, pcross[3]);
      if (flag4){
  	pcell.push_back(pcross[3]);
      }else{
	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size()-1-k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) + u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_v[1],kUwire,kVwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[3]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size()-1-k);
	  dis_v[2] = gds.wire_dist(*vwire_3) + v_pitch/2.;
	  if (gds.crossing_point(dis_u[1],dis_v[2],kUwire,kVwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	    pcell.push_back(pcross[3]);
	    break;
	  }
	}
      }

      if (pcell.size() >=3){
	//Creat a cell and then get all the wires in ... 
	double u_max = -1e9, u_min = 1e9;
	double v_max = -1e9, v_min = 1e9;
	double w_max = -1e9, w_min = 1e9;
	for (int k=0;k!=pcell.size();k++){
	  double udist = gds.wire_dist(pcell.at(k),WirePlaneType_t(0));
	  if (udist > u_max) u_max = udist;
	  if (udist < u_min) u_min = udist;
	  double vdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(1));
	  if (vdist > v_max) v_max = vdist;
	  if (vdist < v_min) v_min = vdist;
	  double wdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(2));
	  if (wdist > w_max) w_max = wdist;
	  if (wdist < w_min) w_min = wdist;
	}
	//std::cout << u_max << " " << u_min << " " << v_max << " " << v_min << " " << w_max << " " << w_min << std::endl;
	//Create a cell
	int ident = holder.get_cell_no();
	SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
	mcell->SetTimeSlice(time_slice);
	mcell->AddBoundary(pcell);
	two_bad_wire_cells.push_back(mcell);
	holder.AddCell_No();
	
	//Insert U
	// for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	//   const GeomWire *uwire = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(k);
	//   double udist = gds.wire_dist(*uwire);
	//   if (udist>u_min && udist < u_max)
	//     mcell->AddWire(uwire,WirePlaneType_t(0));
	// }

	// for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	//   const GeomWire *vwire = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(k);
	//   double vdist = gds.wire_dist(*vwire);
	//   if (vdist>v_min && vdist < v_max)
	//     mcell->AddWire(vwire,WirePlaneType_t(1));
	// }
	// //Insert V
	// for (int k=0;k!=gds.wires_in_plane(WirePlaneType_t(2)).size();k++){
	//   const GeomWire *wwire = gds.wires_in_plane(WirePlaneType_t(2)).at(k);
	//   double wdist = gds.wire_dist(*wwire);
	//   if (wdist>w_min && wdist < w_max)
	//     mcell->AddWire(wwire,WirePlaneType_t(2));
	// }
	
	const GeomWire* uwire_min = gds.bounds(u_min,WirePlaneType_t(0)).second;
	const GeomWire* uwire_max = gds.bounds(u_max,WirePlaneType_t(0)).first;
	for (int k=uwire_min->index();k!=uwire_max->index()+1;k++){
	  const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),k);
	  float charge = 0;
	  if (wirechargemap.find(uwire)!=wirechargemap.end())
	    charge = wirechargemap[uwire];
	  mcell->AddWire(uwire,WirePlaneType_t(0),charge);
	}
	
	//Insert V
	const GeomWire* vwire_min = gds.bounds(v_min,WirePlaneType_t(1)).second;
	const GeomWire* vwire_max = gds.bounds(v_max,WirePlaneType_t(1)).first;
	// std::cout << vwire_min->index() << " " << vwire_max->index() << std::endl;
	for (int k=vwire_min->index();k!=vwire_max->index()+1;k++){
	  const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),k);
	  float charge = 0;
	  if (wirechargemap.find(vwire)!=wirechargemap.end())
	    charge = wirechargemap[vwire];
	  mcell->AddWire(vwire,WirePlaneType_t(1),charge);
	}

	//Insert W
	const GeomWire* wwire_min = gds.closest(w_min,WirePlaneType_t(2));//.bounds(w_min,WirePlaneType_t(2)).second;
	const GeomWire* wwire_max = gds.closest(w_max,WirePlaneType_t(2));//.bounds(w_max,WirePlaneType_t(2)).first;
	for (int k=wwire_min->index();k!=wwire_max->index()+1;k++){
	  const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),k);
	  float charge = 0;
	  if (wirechargemap.find(wwire)!=wirechargemap.end())
	    charge = wirechargemap[wwire];
	  mcell->AddWire(wwire,WirePlaneType_t(2),charge);
	}
	
	//	std::cout << mcell->get_uwires().size() << " " << mcell->get_vwires().size() << " " << mcell->get_wwires().size() << std::endl;
      }
    }
  }


  // U-W and insert V
  for ( int i = 0; i != bad_wire_u.size() ; i++){
    const GeomWire *uwire_1 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().front();
    const GeomWire *uwire_2 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().back();
    float dis_u[2];
    float u_pitch = gds.pitch(kUwire);
    dis_u[0] = gds.wire_dist(*uwire_1) - u_pitch/2.;
    dis_u[1] = gds.wire_dist(*uwire_2) + u_pitch/2.;
    

    for (int j = 0; j != bad_wire_w.size() ; j++){
      const GeomWire *wwire_1 = ((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().front();
      const GeomWire *wwire_2 = ((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().back();
      float dis_w[2];
      float w_pitch = gds.pitch(kYwire);
      dis_w[0] = gds.wire_dist(*wwire_1) - w_pitch/2.;
      dis_w[1] = gds.wire_dist(*wwire_2) + w_pitch/2.;

      //      std::cout << dis_u[0]/units::m << " " << dis_u[1]/units::m << " " << dis_w[0]/units::m << " " << dis_w[1]/units::m << std::endl;

      PointVector pcell;
      std::vector<Vector> pcross(4);
      
      
      bool flag1 = gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, pcross[0]); 
      if (flag1){
  	pcell.push_back(pcross[0]);
      }else{
  	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) - u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_w[0],kUwire,kYwire, pcross[0])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
	// scan w-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().at(k);
	  dis_w[2] = gds.wire_dist(*wwire_3) - w_pitch/2.;
	  if (gds.crossing_point(dis_u[0],dis_w[2],kUwire,kYwire, pcross[0])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
      }

      bool flag2 = gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, pcross[1]); 
      if (flag2){
  	pcell.push_back(pcross[1]);
      }else{
  	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) - u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_w[1],kUwire,kYwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
	// scan w-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().at(((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().size()-1-k);
	  dis_w[2] = gds.wire_dist(*wwire_3) + w_pitch/2.;
	  if (gds.crossing_point(dis_u[0],dis_w[2],kUwire,kYwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
      }
      
      bool flag3 = gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, pcross[2]); 
      if (flag3){
  	pcell.push_back(pcross[2]);
      }else{
  	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size()-1-k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) + u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_w[0],kUwire,kYwire, pcross[2])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
	// scan w-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().at(k);
	  dis_w[2] = gds.wire_dist(*wwire_3) - w_pitch/2.;
	  if (gds.crossing_point(dis_u[1],dis_w[2],kUwire,kYwire, pcross[2])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
      }
      
      bool flag4 = gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, pcross[3]);
      if (flag4){
  	pcell.push_back(pcross[3]);
      }else{
	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size();k++){
	  const GeomWire *uwire_3 = ((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().at(((MergeGeomWire*)bad_wire_u.at(i))->get_allwire().size()-1-k);
	  dis_u[2] =  gds.wire_dist(*uwire_3) + u_pitch/2.;
	  if (gds.crossing_point(dis_u[2],dis_w[1],kUwire,kYwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[3]);
	    break;
	  }
	}
	// scan w-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().at(((MergeGeomWire*)bad_wire_w.at(j))->get_allwire().size()-1-k);
	  dis_w[2] = gds.wire_dist(*wwire_3) + w_pitch/2.;
	  if (gds.crossing_point(dis_u[1],dis_w[2],kUwire,kYwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	    pcell.push_back(pcross[3]);
	    break;
	  }
	}
      }
     
      
    
      if(pcell.size() >=3){
	//Creat a cell and then get all the wires in ... 
	double u_max = -1e9, u_min = 1e9;
	double v_max = -1e9, v_min = 1e9;
	double w_max = -1e9, w_min = 1e9;
	for (int k=0;k!=pcell.size();k++){
	  double udist = gds.wire_dist(pcell.at(k),WirePlaneType_t(0));
	  if (udist > u_max) u_max = udist;
	  if (udist < u_min) u_min = udist;
	  double vdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(1));
	  if (vdist > v_max) v_max = vdist;
	  if (vdist < v_min) v_min = vdist;
	  double wdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(2));
	  if (wdist > w_max) w_max = wdist;
	  if (wdist < w_min) w_min = wdist;
	}
	//std::cout << u_max << " " << u_min << " " << v_max << " " << v_min << " " << w_max << " " << w_min << std::endl;
	//Create a cell
	int ident = holder.get_cell_no();
	SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
	mcell->SetTimeSlice(time_slice);
	mcell->AddBoundary(pcell);
	two_bad_wire_cells.push_back(mcell);
	holder.AddCell_No();
	//Insert U
	const GeomWire* uwire_min = gds.bounds(u_min,WirePlaneType_t(0)).second;
	const GeomWire* uwire_max = gds.bounds(u_max,WirePlaneType_t(0)).first;
	for (int k=uwire_min->index();k!=uwire_max->index()+1;k++){
	  const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),k);
	  float charge = 0;
	  if (wirechargemap.find(uwire)!=wirechargemap.end())
	    charge = wirechargemap[uwire];
	  mcell->AddWire(uwire,WirePlaneType_t(0),charge);
	}
	
	//Insert V
	const GeomWire* vwire_min = gds.closest(v_min,WirePlaneType_t(1));//.bounds(v_min,WirePlaneType_t(1)).second;
	const GeomWire* vwire_max = gds.closest(v_max,WirePlaneType_t(1));//.bounds(v_max,WirePlaneType_t(1)).first;
	// std::cout << vwire_min->index() << " " << vwire_max->index() << std::endl;
	for (int k=vwire_min->index();k!=vwire_max->index()+1;k++){
	  const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),k);
	  float charge = 0;
	  if (wirechargemap.find(vwire)!=wirechargemap.end())
	    charge = wirechargemap[vwire];
	  mcell->AddWire(vwire,WirePlaneType_t(1),charge);
	}

	//Insert W
	const GeomWire* wwire_min = gds.bounds(w_min,WirePlaneType_t(2)).second;
	const GeomWire* wwire_max = gds.bounds(w_max,WirePlaneType_t(2)).first;
	for (int k=wwire_min->index();k!=wwire_max->index()+1;k++){
	  const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),k);
	  float charge = 0;
	  if (wirechargemap.find(wwire)!=wirechargemap.end())
	    charge = wirechargemap[wwire];
	  mcell->AddWire(wwire,WirePlaneType_t(2),charge);
	}


	//	std::cout << mcell->get_uwires().size() << " " << mcell->get_vwires().size() << " " << mcell->get_wwires().size() << std::endl;
      }
    } 
  }


  // W-V and add U
   // deal with w-v
  for ( int i = 0; i != bad_wire_w.size() ; i++){
    const GeomWire *wwire_1 = ((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().front();
    const GeomWire *wwire_2 = ((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().back();
    float dis_w[2];
    float w_pitch = gds.pitch(kYwire);
    dis_w[0] = gds.wire_dist(*wwire_1) - w_pitch/2.;
    dis_w[1] = gds.wire_dist(*wwire_2) + w_pitch/2.;
    

    for (int j = 0; j != bad_wire_v.size() ; j++){
      const GeomWire *vwire_1 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().front();
      const GeomWire *vwire_2 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().back();
      float dis_v[2];
      float v_pitch = gds.pitch(kVwire);
      dis_v[0] = gds.wire_dist(*vwire_1) - v_pitch/2.;
      dis_v[1] = gds.wire_dist(*vwire_2) + v_pitch/2.;
      
      PointVector pcell;
      std::vector<Vector> pcross(4);
      
      bool flag1 = gds.crossing_point(dis_w[0],dis_v[0],kYwire,kVwire, pcross[0]); 
      if (flag1){
  	pcell.push_back(pcross[0]);
      }else{
  	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().at(k);
	  dis_w[2] =  gds.wire_dist(*wwire_3) - w_pitch/2.;
	  if (gds.crossing_point(dis_w[2],dis_v[0],kYwire,kVwire, pcross[0])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(k);
	  dis_v[2] = gds.wire_dist(*vwire_3) - v_pitch/2.;
	  if (gds.crossing_point(dis_w[0],dis_v[2],kYwire,kVwire, pcross[0])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[0].y,2) + pow(pcell.at(k1).z - pcross[0].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[0]);
	    break;
	  }
	}
      }

      bool flag2 = gds.crossing_point(dis_w[0],dis_v[1],kYwire,kVwire, pcross[1]); 
      if (flag2){
  	pcell.push_back(pcross[1]);
      }else{
  	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().at(k);
	  dis_w[2] =  gds.wire_dist(*wwire_3) - w_pitch/2.;
	  if (gds.crossing_point(dis_w[2],dis_v[1],kYwire,kVwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size()-1-k);
	  dis_v[2] = gds.wire_dist(*vwire_3) + v_pitch/2.;
	  if (gds.crossing_point(dis_w[0],dis_v[2],kYwire,kVwire, pcross[1])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[1].y,2) + pow(pcell.at(k1).z - pcross[1].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[1]);
	    break;
	  }
	}
      }
      
      bool flag3 = gds.crossing_point(dis_w[1],dis_v[0],kYwire,kVwire, pcross[2]); 
      if (flag3){
  	pcell.push_back(pcross[2]);
      }else{
  	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().at(((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().size()-1-k);
	  dis_w[2] =  gds.wire_dist(*wwire_3) + w_pitch/2.;
	  if (gds.crossing_point(dis_w[2],dis_v[0],kYwire,kVwire, pcross[2])){

	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(k);
	  dis_v[2] = gds.wire_dist(*vwire_3) - v_pitch/2.;
	  if (gds.crossing_point(dis_w[1],dis_v[2],kYwire,kVwire, pcross[2])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[2].y,2) + pow(pcell.at(k1).z - pcross[2].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[2]);
	    break;
	  }
	}
      }
      
      bool flag4 = gds.crossing_point(dis_w[1],dis_v[1],kYwire,kVwire, pcross[3]);
      if (flag4){
  	pcell.push_back(pcross[3]);
      }else{
	// scan u-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().size();k++){
	  const GeomWire *wwire_3 = ((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().at(((MergeGeomWire*)bad_wire_w.at(i))->get_allwire().size()-1-k);
	  dis_w[2] =  gds.wire_dist(*wwire_3) + w_pitch/2.;
	  if (gds.crossing_point(dis_w[2],dis_v[1],kYwire,kVwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	      pcell.push_back(pcross[3]);
	    break;
	  }
	}
	// scan v-wire
	for (int k=0;k!=((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size();k++){
	  const GeomWire *vwire_3 = ((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().at(((MergeGeomWire*)bad_wire_v.at(j))->get_allwire().size()-1-k);
	  dis_v[2] = gds.wire_dist(*vwire_3) + v_pitch/2.;
	  if (gds.crossing_point(dis_w[1],dis_v[2],kYwire,kVwire, pcross[3])){
	    int flag_qx = 1;
	    for (int k1 = 0;k1!=pcell.size();k1++){
	      if (sqrt(pow(pcell.at(k1).y - pcross[3].y,2) + pow(pcell.at(k1).z - pcross[3].z,2))<0.5*units::mm){
		flag_qx = 0;
		break;
	      }
	    }
	    if (flag_qx == 1)
	    pcell.push_back(pcross[3]);
	    break;
	  }
	}
      }
     
      
    
      if(pcell.size() >=3){
	//Creat a cell and then get all the wires in ... 
	double u_max = -1e9, u_min = 1e9;
	double v_max = -1e9, v_min = 1e9;
	double w_max = -1e9, w_min = 1e9;
	for (int k=0;k!=pcell.size();k++){
	  double udist = gds.wire_dist(pcell.at(k),WirePlaneType_t(0));
	  if (udist > u_max) u_max = udist;
	  if (udist < u_min) u_min = udist;
	  double vdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(1));
	  if (vdist > v_max) v_max = vdist;
	  if (vdist < v_min) v_min = vdist;
	  double wdist = gds.wire_dist(pcell.at(k),WirePlaneType_t(2));
	  if (wdist > w_max) w_max = wdist;
	  if (wdist < w_min) w_min = wdist;
	}
	//std::cout << u_max << " " << u_min << " " << v_max << " " << v_min << " " << w_max << " " << w_min << std::endl;
	//Create a cell
	int ident = holder.get_cell_no();
	SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
	mcell->SetTimeSlice(time_slice);
	mcell->AddBoundary(pcell);
	two_bad_wire_cells.push_back(mcell);
	holder.AddCell_No();

	//Insert U	
	const GeomWire* uwire_min = gds.closest(u_min,WirePlaneType_t(0));//.bounds(u_min,WirePlaneType_t(0)).second;
	const GeomWire* uwire_max = gds.closest(u_max,WirePlaneType_t(0));//.bounds(u_max,WirePlaneType_t(0)).first;
	for (int k=uwire_min->index();k!=uwire_max->index()+1;k++){
	  const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),k);
	  float charge = 0;
	  if (wirechargemap.find(uwire)!=wirechargemap.end())
	    charge = wirechargemap[uwire];
	  mcell->AddWire(uwire,WirePlaneType_t(0),charge);
	}
	
	//Insert V
	const GeomWire* vwire_min = gds.bounds(v_min,WirePlaneType_t(1)).second;
	const GeomWire* vwire_max = gds.bounds(v_max,WirePlaneType_t(1)).first;
	// std::cout << vwire_min->index() << " " << vwire_max->index() << std::endl;
	for (int k=vwire_min->index();k!=vwire_max->index()+1;k++){
	  const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),k);
	  float charge = 0;
	  if (wirechargemap.find(vwire)!=wirechargemap.end())
	    charge = wirechargemap[vwire];
	  mcell->AddWire(vwire,WirePlaneType_t(1),charge);
	}

	//Insert W
	const GeomWire* wwire_min = gds.bounds(w_min,WirePlaneType_t(2)).second;
	const GeomWire* wwire_max = gds.bounds(w_max,WirePlaneType_t(2)).first;
	for (int k=wwire_min->index();k!=wwire_max->index()+1;k++){
	  const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),k);
	  float charge = 0;
	  if (wirechargemap.find(wwire)!=wirechargemap.end())
	    charge = wirechargemap[wwire];
	  mcell->AddWire(wwire,WirePlaneType_t(2),charge);
	}


	// std::cout << mcell->get_uwires().size() << " " << mcell->get_vwires().size() << " " << mcell->get_wwires().size() << std::endl;
      }
    } 
  }

}


void WireCell2dToy::LowmemTiling::form_fired_merge_wires(const WireCell::Slice& slice, const WireCell::Slice& slice_err){
  WireCell::Channel::Group group = slice.group();
  WireCell::Channel::Group group_err = slice_err.group();
  
  
  //double sum = 0;
  for (int i=0;i!=group.size();i++){
    const WireCell::GeomWire *wire = gds.by_channel_segment(group.at(i).first,0);
    float charge = group.at(i).second;
    float charge_err = group_err.at(i).second;
    //  if (charge < 11) charge = 11;
    if (wirechargemap.find(wire) == wirechargemap.end()){
      //not found
      wirechargemap[wire] = charge;
      wirecharge_errmap[wire] = charge_err;
      //std::cout << wirechargemap[wire] << " " << wirecharge_errmap[wire] << std::endl;
    }else{
      //wirechargemap[wire] += charge;
    }
  }
  
  // do U
  MergeGeomWire *mwire = 0;
  int ident = 0;
  int last_wire = -1;
  for (int i = 0; i!= nwire_u; i++){
    const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)0,i);
    if (wirechargemap.find(wire)!=wirechargemap.end()){
      if (mwire == 0){
	mwire = new MergeGeomWire(ident,*wire);
	//mwire->SetTimeSlice(time_slice);
	wire_type_map[mwire] = true;
	fired_wire_u.push_back(mwire);
	ident ++;
      }else{
	if (i==last_wire+1){
	  mwire->AddWire(*wire);
	}else{
	  // create a new wire
	  mwire = new MergeGeomWire(ident,*wire);
	  // mwire->SetTimeSlice(time_slice);
	  wire_type_map[mwire] = true;
	  fired_wire_u.push_back(mwire);
	  ident ++;
	}
      }
      last_wire = i;
    }
  }
  
  // V
  mwire = 0;
  last_wire = -1;
  for (int i = 0; i!= nwire_v; i++){
    const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)1,i);
    if (wirechargemap.find(wire)!=wirechargemap.end()){
      if (mwire == 0){
	mwire = new MergeGeomWire(ident,*wire);
	//mwire->SetTimeSlice(time_slice);
	wire_type_map[mwire] = true;
	fired_wire_v.push_back(mwire);
	ident ++;
      }else{
	if (i==last_wire+1){
	  mwire->AddWire(*wire);
	}else{
	  // create a new wire
	  mwire = new MergeGeomWire(ident,*wire);
	  //mwire->SetTimeSlice(time_slice);
	  wire_type_map[mwire] = true;
	  fired_wire_v.push_back(mwire);
	  ident ++;
	}
      }
      last_wire = i;
    }
  }
  
  // W
  mwire = 0;
  last_wire = -1;
  for (int i = 0; i!= nwire_w; i++){
    const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)2,i);
    if (wirechargemap.find(wire)!=wirechargemap.end()){
      if (mwire == 0){
	mwire = new MergeGeomWire(ident,*wire);
	//mwire->SetTimeSlice(time_slice);
	wire_type_map[mwire] = true;
	fired_wire_w.push_back(mwire);
	ident ++;
      }else{
	if (i==last_wire+1){
	  mwire->AddWire(*wire);
	}else{
	  // create a new wire
	  mwire = new MergeGeomWire(ident,*wire);
	  //mwire->SetTimeSlice(time_slice);
	  wire_type_map[mwire] = true;
	  fired_wire_w.push_back(mwire);
	  ident ++;
	}
      }
      last_wire = i;
    }
  }

  
  // int sum = 0;
  // for (int i=0;i!=fired_wire_u.size();i++){
  //   MergeGeomWire *mwire = (MergeGeomWire*)fired_wire_u.at(i);
  //   sum += mwire->get_allwire().size();
  //   for (Int_t j=0;j!=mwire->get_allwire().size();j++){
  //     std::cout << i << " " << mwire->get_allwire().at(j)->index() << std::endl;
  //   }
  // }
  // for (int i=0;i!=fired_wire_v.size();i++){
  //   MergeGeomWire *mwire = (MergeGeomWire*)fired_wire_v.at(i);
  //   sum += mwire->get_allwire().size();
  // }
  // for (int i=0;i!=fired_wire_w.size();i++){
  //   MergeGeomWire *mwire = (MergeGeomWire*)fired_wire_w.at(i);
  //   sum += mwire->get_allwire().size();
  // }

  
  //std::cout << fired_wire_u.size() << " " << fired_wire_v.size() << " " << fired_wire_w.size() << " " << group.size() << " " << sum << std::endl;
  
}

void WireCell2dToy::LowmemTiling::form_bad_merge_wires(WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map){
  
  // do U
  MergeGeomWire *mwire = 0;
  int last_wire = -1;

  for (int i = 0; i!= nwire_u; i++){
    if (uplane_map.find(i)!=uplane_map.end()){
      if (time_slice >= uplane_map[i].first/nrebin && 
	  time_slice <= uplane_map[i].second/nrebin){
	const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)0,i);
	// first one 
	if (mwire == 0){
	  int ident = holder.get_wire_no();
	  mwire = new MergeGeomWire(ident,*wire);
	  wire_type_map[mwire] = false;
	  bad_wire_u.push_back(mwire);
	  holder.AddWire_No();
	}else{
	  if (i==last_wire+1){
	    mwire->AddWire(*wire);
	  }else{
	    // create a new wire
	    int ident = holder.get_wire_no();
	    mwire = new MergeGeomWire(ident,*wire);
	    wire_type_map[mwire] = false;
	    bad_wire_u.push_back(mwire);
	    holder.AddWire_No();
	  }
	}
	last_wire = i;
      }
    }
  }

  // now do V plane
  mwire = 0;
  last_wire = -1;

  for (int i = 0; i!= nwire_v; i++){
    if (vplane_map.find(i)!=vplane_map.end()){
      if (time_slice >= vplane_map[i].first/nrebin && 
	  time_slice <= vplane_map[i].second/nrebin){
	const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)1,i);
	// first one 
	if (mwire == 0){
	  int ident = holder.get_wire_no();
	  mwire = new MergeGeomWire(ident,*wire);
	  wire_type_map[mwire] = false;
	  bad_wire_v.push_back(mwire);
	  holder.AddWire_No();
	}else{
	  if (i==last_wire+1){
	    mwire->AddWire(*wire);
	  }else{
	    // create a new wire
	    int ident = holder.get_wire_no();
	    mwire = new MergeGeomWire(ident,*wire);
	    wire_type_map[mwire] = false;
	    bad_wire_v.push_back(mwire);
	    holder.AddWire_No();
	  }
	}
	last_wire = i;
      }
    }
  }
  
  // W plane
  mwire = 0;
  last_wire = -1;
  
  for (int i = 0; i!= nwire_w; i++){
    if (wplane_map.find(i)!=wplane_map.end()){
      if (time_slice >= wplane_map[i].first/nrebin && 
	  time_slice <= wplane_map[i].second/nrebin){
	const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)2,i);
	// first one 
	if (mwire == 0){
	  int ident = holder.get_wire_no();
	  mwire = new MergeGeomWire(ident,*wire);
	  wire_type_map[mwire] = false;
	  bad_wire_w.push_back(mwire);
	  holder.AddWire_No();
	}else{
	  if (i==last_wire+1){
	    mwire->AddWire(*wire);
	  }else{
	    // create a new wire
	    int ident = holder.get_wire_no();
	    mwire = new MergeGeomWire(ident,*wire);
	    wire_type_map[mwire] = false;
	    bad_wire_w.push_back(mwire);
	    holder.AddWire_No();
	  }
	}
	last_wire = i;
      }
    }
  }
  


  // int sum = 0;
  // for (int i=0;i!=bad_wire_w.size();i++){
  //   MergeGeomWire *mwire = (MergeGeomWire*)bad_wire_w.at(i);
  //   sum += mwire->get_allwire().size();
  //   for (int j=0;j!=mwire->get_allwire().size();j++){
  //     GeomWire *wire = (GeomWire*)mwire->get_allwire().at(j);
  //     std::cout << i << " " << wire->index() << std::endl;
  //   }
  // }
  // std::cout << uplane_map.size() << " " << bad_wire_u.size() << " " << sum << std::endl;
  //std::cout << vplane_map.size() << " " << bad_wire_v.size() << " " << sum << std::endl;
  //std::cout << wplane_map.size() << " " << bad_wire_w.size() << " " << sum << std::endl;

}

WireCell2dToy::LowmemTiling::~LowmemTiling(){
  // for (int i=0;i!=bad_wire_u.size();i++){
  //   delete bad_wire_u.at(i);
  // }
  // bad_wire_u.clear();
  
  // for (int i=0;i!=bad_wire_v.size();i++){
  //   delete bad_wire_v.at(i);
  // }
  // bad_wire_v.clear();

  // for (int i=0;i!=bad_wire_w.size();i++){
  //   delete bad_wire_w.at(i);
  // }
  // bad_wire_w.clear();
  
  // for (int i=0;i!=fired_wire_u.size();i++){
  //   delete fired_wire_u.at(i);
  // }
  // fired_wire_u.clear();
  
  // for (int i=0;i!=fired_wire_v.size();i++){
  //   delete fired_wire_v.at(i);
  // }
  // fired_wire_v.clear();

  // for (int i=0;i!=fired_wire_w.size();i++){
  //   delete fired_wire_w.at(i);
  // }
  // fired_wire_w.clear();
}

WireCell::GeomWireSelection WireCell2dToy::LowmemTiling::wires(const WireCell::GeomCell& cell) const{
  GeomWireSelection wires;
  return wires;
}
WireCell::GeomCellSelection WireCell2dToy::LowmemTiling::cells(const WireCell::GeomWire& wire) const{
  GeomCellSelection cells;
  return cells;
}

const WireCell::GeomCell* WireCell2dToy::LowmemTiling::cell(const WireCell::GeomWireSelection& wires) const{
  return 0;
}
