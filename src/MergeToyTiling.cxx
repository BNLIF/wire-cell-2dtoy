#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"

#include <cmath>

using namespace WireCell;




WireCell2dToy::MergeToyTiling::MergeToyTiling(WireCell2dToy::ToyTiling& tiling, int time_slice, int merge_strategy, int flag_remerge){
  IsRemerged = false;

  ncell = tiling.get_ncell();
  wire_u = tiling.get_wire_u();
  wire_v = tiling.get_wire_v();
  wire_w = tiling.get_wire_w();
  wirechargemap = tiling.wcmap();
  
  // goal is to create merged version of 
  // cell_all
  // cellmap
  // wiremap

  std::cout << "Merge Strategy " << merge_strategy << std::endl;

  if (merge_strategy == 1){
    
    // //start with wire_u
    for (int i =0;i!=wire_u.size();i++){
      //ntemp += tiling.cells(*tiling.get_wire_u()[i]).size();
      for (int j=0;j!=tiling.cells(*tiling.get_wire_u()[i]).size();j++){
	const GeomCell *cell = tiling.cells(*wire_u[i])[j];
	int flag=0;
	
	for (int k=0;k!=cell_all.size();k++){
	  if (((MergeGeomCell*)cell_all[k])->AddCell(*cell)){
	    flag = 1;
	    break;
	  }
	}
	
	if(flag==0){
	  MergeGeomCell *mcell = new MergeGeomCell(ncell,*cell);
	  mcell->SetTimeSlice(time_slice);
	  ncell++;
	  cell_all.push_back(mcell);
	}
      }
    }
    
    //int ntemp = 0;
    // for (int i=0;i!=cell_all.size();i++){
    //   MergeGeomCell *mcell = (MergeGeomCell*)cell_all.at(i);
    //   ntemp += mcell->get_allcell().size();
    // }
    // std::cout << ntemp << " " << cell_all.size() << std::endl;
    
    while(further_merge(cell_all,tiling.get_ncell(),time_slice));
    
  }else if (merge_strategy == 2){
    //try Brett's new merge strategy ... 
    //put cells into a list
    //std::cout << "Copy " << std::endl; 
    GeomCellSelection cells = tiling.get_allcell();
    GeomCellList cell_list(cells.begin(),cells.end());
    // // creat a list
    // for (int i = 0; i!= tiling.get_allcell().size();i++){
    //   cell_list.push_back(tiling.get_allcell().at(i));
    // }
    
    while(cell_list.size()!=0){
      
      // construct a new merged cell from the first element of the cell
      auto it_abc = cell_list.begin();
      const GeomCell *first_cell = *it_abc;
      MergeGeomCell *mcell = new MergeGeomCell(ncell,*first_cell);
      mcell->SetTimeSlice(time_slice);
      //std::cout << "Big: " << cell_list.size() << " " << cell_all.size() << first_cell->center().y << " " <<first_cell->center().z << std::endl;
      it_abc = cell_list.erase(it_abc);
      ncell++;
      cell_all.push_back(mcell);

      // loop through all the cells in the existing list 
      // and add them to this merge cell
      // need to keep a lits to save the ones that are successful
      GeomCellList suceed_list;

      // auto it =cell_list.begin();
      // while(it!=cell_list.end()){
      // 	for (auto it1 = it;it1!=cell_list.end();it1++){
      // 	  const GeomCell *current_cell = *it1;
      // 	  if (mcell->Connected(*first_cell,*current_cell)){
      // 	    mcell->AddNewCell(*current_cell);
      // 	    suceed_list.push_back(current_cell);
      // 	    //cell_list.erase(it);
      // 	    it = it1;
      // 	    break;
      // 	  }
      // 	}
      // 	it = cell_list.erase(it);
      // }
      
      auto it1 = it_abc;
      while(it1!=cell_list.end()){
	//      for (auto it = it_abc;it!=cell_list.end();it++){
	const GeomCell *current_cell = *it1;
	if (mcell->Connected(*first_cell,*current_cell)){
	  mcell->AddNewCell(*current_cell);
	  suceed_list.push_back(current_cell);
	  it1 = cell_list.erase(it1);
	}else{
	  it1++;
	}
      }
      //std::cout << suceed_list.size() << std::endl;

      //remove the succeed ones from the original list
      // for (auto it = suceed_list.begin(); it!=suceed_list.end(); it++){
      //  	const GeomCell *current_cell = *it;
      //  	cell_list.remove(current_cell);
      // }
      
      // Now, need to go through the existing
      while(suceed_list.size()!=0){
	
	const GeomCell *current_cell = *(suceed_list.begin()); // get first element
	suceed_list.erase(suceed_list.begin());
	GeomCellList temp_list;
	
	auto it = cell_list.begin();
	while(it!=cell_list.end()){
	  //	for (auto it = cell_list.begin();it!=cell_list.end();it++){
	  const GeomCell *current_cell1 = *it;
	  if (mcell->Connected(*current_cell,*current_cell1)){
	    mcell->AddNewCell(*current_cell1);
	    temp_list.push_back(current_cell1);
	    suceed_list.push_back(current_cell1);
	    it = cell_list.erase(it);
	  }else{
	    it ++;
	  }
	}
	
	//std::cout << "Small " << suceed_list.size() << " " << temp_list.size() << std::endl;
	//remove the succeed ones from the original list and add them into the suceed_list
	// for (auto it = temp_list.begin(); it!=temp_list.end(); it++){
	//   const GeomCell *current_cell1 = *it;
	//    cell_list.remove(current_cell1);
	  
	// }
      }
      

    }
  }else if (merge_strategy == 3){
    
    //create a lot of wire lists Map(wire) --> list of cells
    GeomWireSelection wires = tiling.get_allwire();
    GeomWireLMap wlmap;
    for (int i=0;i!=wires.size();i++){
      const GeomWire *wire = wires.at(i);
      GeomCellSelection cells = tiling.cells(*wire);
      GeomCellList abc(cells.begin(),cells.end());
      wlmap[wire] = abc;
    }
    

    //start from wire
    int i = 0;
    //    for (int i=0;i!=wires.size();i++){
    while(i!=wires.size()){
      
      const GeomWire *wire = wires.at(i);
      GeomCellList& abc = wlmap[wire];
      
      auto it = abc.begin();
      if (it!=abc.end()){
	const GeomCell *cell = *it;
	MergeGeomCell *mcell = new MergeGeomCell(ncell,*cell);
	mcell->SetTimeSlice(time_slice);
	//std::cout << "Big: " << i << " " << cell_all.size() << " " << abc.size()  << std::endl;
	ncell++;
	cell_all.push_back(mcell);
	it = abc.erase(it);
	//	std::cout << abc.size() << std::endl;
	//remove this cell from all its wires
	GeomWireSelection wires4 = tiling.wires(*cell);
	for (int k=0;k!=wires4.size();k++){
	  const GeomWire *wire4 = wires4.at(k);
	  if (wire4 != wire){
	    GeomCellList& abc2 = wlmap[wire4];
	    abc2.remove(cell);
	  }
	}

	
	GeomCellList suceed_list;
	//Now loop through all cells in the wires associated with this cell
	GeomWireSelection wires1 = tiling.wires(*cell);
	for (int j=0;j!=wires1.size();j++){
	  const GeomWire *wire1 = wires1.at(j);
	  GeomCellList& abc1 = wlmap[wire1];
	  auto it1 = abc1.begin();
	  while(it1!=abc1.end()){
	    //	  for (auto it1 = abc1.begin();it1!=abc1.end();it1++){
	    const GeomCell *cell1 = *it1;
	    //judge if the cell1 is connected
	    if (mcell->Connected(*cell,*cell1)){
	      mcell->AddNewCell(*cell1);
	      suceed_list.push_back(cell1);
	      it1 = abc1.erase(it1);

	      //remove the cell from the other wires
	      GeomWireSelection wires2 = tiling.wires(*cell1);
	      for (int k=0;k!=wires2.size();k++){
		const GeomWire *wire2 = wires2.at(k);
		if (wire2 != wire1){
		  GeomCellList& abc2 = wlmap[wire2];
		  abc2.remove(cell1);
		}
	      }
	    }else{
	      it1++;
	    }
	    
	  }
	}

	
	//Now deal with the suceed list
	while(suceed_list.size()!=0){
	  const GeomCell *cell1 = *(suceed_list.begin()); // get first element
	  suceed_list.erase(suceed_list.begin());
	  
	  GeomWireSelection wires2 = tiling.wires(*cell1);
	  for (int j=0;j!=wires2.size();j++){
	    const GeomWire *wire2 = wires2.at(j);
	    GeomCellList& abc2 = wlmap[wire2];
	    // now loop through abc2's cells

	    auto it1 = abc2.begin();
	    while(it1!=abc2.end()){
	      const GeomCell *cell2 = *it1;
	      //judge if the cell1 is connected
	      if (mcell->Connected(*cell1,*cell2)){
		mcell->AddNewCell(*cell2);
		suceed_list.push_back(cell2);
		it1 = abc2.erase(it1);

		//remove the cell from the other wires
		GeomWireSelection wires3 = tiling.wires(*cell2);
		for (int k=0;k!=wires3.size();k++){
		  const GeomWire *wire3 = wires3.at(k);
		  if (wire3 != wire2){
		    GeomCellList& abc3 = wlmap[wire3];
		    abc3.remove(cell2);
		  }
		}
	      }else{
		it1++;
	      }
	      
	     
	    }
	  }
	 }

 
      }else{
	i++;
      }
    }
    
  }
  // ntemp = 0;
  // for (int i=0;i!=cell_all.size();i++){
  //   MergeGeomCell *mcell = (MergeGeomCell*)cell_all.at(i);
  //   ntemp += mcell->get_allcell().size();
  // }
  // std::cout << ntemp << " " << cell_all.size()<< std::endl;

  MergeGeomCellSet mset;
  // rank the merged cell ... 
  for (int i=0;i!=cell_all.size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*)cell_all[i];
    //std::cout << mcell << " " << mcell->cross_section() << std::endl;
    mset.insert(mcell);
  }
  
  cell_all.clear();
  for (auto it = mset.begin();it!=mset.end();it++){
    cell_all.push_back(*it);
    //std::cout << (*it)->cross_section() << std::endl;
  }

  // ntemp = 0;
  // for (int i=0;i!=cell_all.size();i++){
  //   MergeGeomCell *mcell = (MergeGeomCell*)cell_all.at(i);
  //   ntemp += mcell->get_allcell().size();
  // }
  // std::cout << ntemp << " " << cell_all.size()<< std::endl;

  // Now merge all the wires   wire_all
  
  form_wiremap(tiling, time_slice);
  
  if (flag_remerge == 1){
    double dis = 0*units::mm;

    for (int i=0;i!=cell_all.size();i++){
      MergeGeomCell* mcell = (MergeGeomCell*)cell_all.at(i);
      //std::cout << "Before: " << mcell->boundary().size() << " " ; 
      mcell->Organize_edge_boundary();
      //std::cout << mcell->boundary().size() << std::endl;
    }

    int current_ncell = cell_all.size();
    int prev_ncell = -1;
    
    while (tiling.get_allcell().size()>10000 && cell_all.size() >0.45 * wire_all.size()){
      dis += sqrt(0.3*0.3/2.*tiling.get_allcell().size())/3.*units::cm;
      
      
      IsRemerged = true;
      while(further_merge(cell_all,tiling.get_ncell(),time_slice,dis));
      current_ncell = cell_all.size();
      
      if (current_ncell != prev_ncell){
	//clear all wire and all maps
	for (int i=0;i!=wire_all.size();i++){
	  delete wire_all[i];
	}
	wire_u.clear();
	wire_v.clear();
	wire_w.clear();
	wire_all.clear();
	cellmap.clear();
	wiremap.clear();
	cellmap1.clear();
	wiremap1.clear();
	wwmap.clear();
	wwsmap.clear();
	form_wiremap(tiling, time_slice);
      }
      
      prev_ncell = current_ncell;
    }
    
    current_ncell = cell_all.size();
    prev_ncell = -1;
    
    while (cell_all.size() > 2 * wire_all.size() && cell_all.size() - wire_all.size() > 50){
      dis += 6*units::mm;
      
      
      //std::cout << " Start to remerge blob " << cell_all.size() << " " << wire_all.size() << std::endl; 
      IsRemerged = true;
      
      
      while(further_merge(cell_all,tiling.get_ncell(),time_slice,dis));
      
      current_ncell = cell_all.size();
      if (current_ncell != prev_ncell){
	//clear all wire and all maps
	for (int i=0;i!=wire_all.size();i++){
	  delete wire_all[i];
	}
	wire_u.clear();
	wire_v.clear();
	wire_w.clear();
	wire_all.clear();
	cellmap.clear();
	wiremap.clear();
	cellmap1.clear();
	wiremap1.clear();
	wwmap.clear();
	wwsmap.clear();
	form_wiremap(tiling, time_slice);
      }
      prev_ncell = current_ncell;
    }

  }

}

void WireCell2dToy::MergeToyTiling::form_wiremap(WireCell2dToy::ToyTiling& tiling, int time_slice){
  int ident_wire = 50000;
  

  for (int i=0;i!=cell_all.size();i++){    
    GeomCellSelection call =  ((MergeGeomCell*)cell_all[i])->get_allcell();
    for (int k=0;k!=3;k++){
      WirePlaneType_t plane = (WirePlaneType_t)k;
      MergeGeomWire *mwire = 0;
      int flag = 0;

      for (int j=0;j!=call.size();j++){
   	GeomWireSelection wires = tiling.wires(*call[j]);
  	for (int nwire = 0; nwire!=wires.size();nwire++){

  	  // std::cout << wires[nwire]->plane() << " " << plane << std::endl;

   	  if (wires[nwire]->plane()==plane){
	    if (flag==0){
	      mwire = new MergeGeomWire(ident_wire,*wires[nwire]);
	      mwire->SetTimeSlice(time_slice);
	      
	      ident_wire++;
	      flag = 1;
	    }else {
	      mwire->AddWire(*wires[nwire]);
	    }
  	  }
  	}
      }
      
     wire_all.push_back(mwire);
      
     
    }
  }

  while(further_mergewire(wire_all,50000,time_slice));
  

  // Now construc the wire map;
  for (int i=0;i!=wire_all.size();i++){
    MergeGeomWire *mwire = (MergeGeomWire*)wire_all[i];
    GeomWireSelection wires = mwire->get_allwire();
    //std::cout << wires.size() << std::endl;
    wwsmap[mwire] = wires;
    for (int j=0;j!=wires.size();j++){
      wwmap[wires[j]] = mwire;
    }
  }


  // Now construct the map
  for (int i=0;i!=cell_all.size();i++){
    const MergeGeomCell *cell = (MergeGeomCell*)cell_all[i];
    GeomWireSelection wiresel;
        
    for (int j=0;j!=cell->get_allcell().size();j++){
      const GeomCell *scell = cell->get_allcell()[j];
      //std::cout << i << " " << scell->ident()<< " " << tiling.wires(*scell).size() << std::endl;
      for (int k=0;k!=tiling.wires(*scell).size();k++){
  	const GeomWire *wire = tiling.wires(*scell)[k];
  	wiresel.push_back(wire);

  	//also do the wiremap
  	if (wiremap1.find(wire) == wiremap1.end()){
  	  GeomCellSelection cellsel;
  	  cellsel.push_back(cell);
  	  wiremap1[wire]= cellsel;
  	}else{
  	  int flag = 0;
  	  for (int n=0;n!=wiremap1[wire].size();n++){
  	    if (cell == wiremap1[wire].at(n)){
  	      flag = 1;
  	      break;
  	    }
  	  }
  	  if(flag==0){
  	    wiremap1[wire].push_back(cell);
  	  }
  	}

      }
    }

    cellmap1[cell] = wiresel;
  }


  //Now construct the real map
  //loop through merged cells
  for (int i=0;i!=cell_all.size();i++){
    const MergeGeomCell *cell = (MergeGeomCell*)cell_all[i];
    for (int j=0;j!=cellmap1[cell].size();j++){
      const GeomWire *wire = (cellmap1[cell])[j];
      const MergeGeomWire *mwire = (MergeGeomWire*)wwmap[wire];

      if (cellmap.find(cell) == cellmap.end()){
	GeomWireSelection wiresel;
	wiresel.push_back(mwire);
	cellmap[cell] = wiresel;
      }else{
	GeomWireSelection wiresel = cellmap[cell];
	int flag = 0;
	for (int k=0; k!=wiresel.size();k++){
	  if (wiresel[k] == mwire){
	    flag = 1;
	    break;
	  }
	}
	if (flag==0){
	  cellmap[cell].push_back(mwire);
	}
      }
      
      
      if (wiremap.find(mwire) == wiremap.end()){
	GeomCellSelection cellsel;
	cellsel.push_back(cell);
	wiremap[mwire] = cellsel;
      }else{
	GeomCellSelection cellsel = wiremap[mwire];
	int flag = 0;
	for (int k=0; k!=cellsel.size();k++){
	  if (cellsel[k] == cell){
	    flag = 1;
	    break;
	  }
	}
	if (flag==0){
	  wiremap[mwire].push_back(cell);
	}
      }

    }
  }
}



WireCell2dToy::MergeToyTiling::~MergeToyTiling(){
  for (int i=0;i!=cell_all.size();i++){
    delete cell_all[i];
  }

  for (int i=0;i!=wire_all.size();i++){
    delete wire_all[i];
  }

  wire_u.clear();
  wire_v.clear();
  wire_w.clear();
  wire_all.clear();
  cell_all.clear();
  cellmap.clear();
  wiremap.clear();

  cellmap1.clear();
  wiremap1.clear();

  wirechargemap.clear();

  wwmap.clear();
  wwsmap.clear();
}


int WireCell2dToy::MergeToyTiling::further_mergewire(WireCell::GeomWireSelection &allwire, int nwire, int time_slice){

  WireCell::GeomWireSelection tempwire = allwire;
  allwire.clear();

  for (int i=0;i!=tempwire.size();i++){
    MergeGeomWire *wire = (MergeGeomWire*)tempwire[i];
    int flag=0;
    for (int k=0;k!=allwire.size();k++){
      if (((MergeGeomWire*)allwire[k])->AddWire(*wire)){
    	flag = 1;
    	break;
      }
    }
      
    if(flag==0){
      MergeGeomWire *mwire = new MergeGeomWire(nwire,*wire);
      mwire->SetTimeSlice(time_slice);
      nwire++;
      allwire.push_back(mwire);
    }
  }
  
  int diff = tempwire.size() - allwire.size();
  
  for (int i=0;i!=tempwire.size();i++){
    delete tempwire[i];
  }
  tempwire.clear();
  
  
  return diff;



  // if (allwire.size()>0){
  //   MergeGeomWire *wire = (MergeGeomWire*)allwire[0];
    
  //   std::cout << wire->get_allwire().size() << std::endl;
  //   //    MergeGeomWire *mwire = new MergeGeomWire(nwire,*wire);
  //   //mwire = 0;
  // }
  // return 0;
}




int WireCell2dToy::MergeToyTiling::further_merge(WireCell::GeomCellSelection &allcell, int ncell,int time_slice, double dis){
  WireCell::GeomCellSelection tempcell = allcell;
  allcell.clear();

  

  for (int i=0;i!=tempcell.size();i++){
    MergeGeomCell *cell = (MergeGeomCell*)tempcell[i];
      int flag=0;
      for (int k=0;k!=allcell.size();k++){
	


	if (((MergeGeomCell*)allcell[k])->AddCell(*cell,dis)){
	  flag = 1;
	  break;
	}
      }
      
      if(flag==0){
      	MergeGeomCell *mcell = new MergeGeomCell(ncell,*cell);
	mcell->SetTimeSlice(time_slice);
	ncell++;
	allcell.push_back(mcell);
      }
  }
  
  int diff = tempcell.size() - allcell.size();
  
  for (int i=0;i!=tempcell.size();i++){
    //tempcell[i] = 0;
    delete tempcell[i];
  }
  tempcell.clear();

 
  
  return diff;
}


const WireCell::GeomCell* WireCell2dToy::MergeToyTiling::cell(const WireCell::GeomWireSelection& wires) const
{
  return 0;
}

ClassImp(WireCell2dToy::MergeToyTiling);
