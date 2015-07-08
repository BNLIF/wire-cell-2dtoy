#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"

#include <cmath>

using namespace WireCell;




WireCell2dToy::MergeToyTiling::MergeToyTiling(WireCell2dToy::ToyTiling& tiling, int time_slice, int merge_strategy){
  

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
    std::cout << "Copy " << std::endl; 
    GeomCellList cell_list;
    // creat a list
    for (int i = 0; i!= tiling.get_allcell().size();i++){
      cell_list.push_back(tiling.get_allcell().at(i));
    }
    
    while(cell_list.size()!=0){
      std::cout << "Big: " << cell_list.size() << " " << cell_all.size() << std::endl;
      // construct a new merged cell from the first element of the cell
      const GeomCell *first_cell = *(cell_list.begin());
      MergeGeomCell *mcell = new MergeGeomCell(ncell,*first_cell);
      cell_list.remove(first_cell);
      ncell++;
      cell_all.push_back(mcell);

      // loop through all the cells in the existing list 
      // and add them to this merge cell
      // need to keep a lits to save the ones that are successful
      GeomCellList suceed_list;
      for (auto it = cell_list.begin();it!=cell_list.end();it++){
	const GeomCell *current_cell = *it;
	if (mcell->Connected(*first_cell,*current_cell)){
	  mcell->AddNewCell(*current_cell);
	  suceed_list.push_back(current_cell);
	}
      }
      //remove the succeed ones from the original list
      for (auto it = suceed_list.begin(); it!=suceed_list.end(); it++){
	const GeomCell *current_cell = *it;
	cell_list.remove(current_cell);
      }
      
      // Now, need to go through the existing
      while(suceed_list.size()!=0){
	
	const GeomCell *current_cell = *(suceed_list.begin()); // get first element
	suceed_list.remove(current_cell);
	GeomCellList temp_list;
	
	for (auto it = cell_list.begin();it!=cell_list.end();it++){
	  const GeomCell *current_cell1 = *it;
	  if (mcell->Connected(*current_cell,*current_cell1)){
	    mcell->AddNewCell(*current_cell1);
	    temp_list.push_back(current_cell1);
	  }
	}
	//std::cout << "Small " << suceed_list.size() << " " << temp_list.size() << std::endl;
	//remove the succeed ones from the original list and add them into the suceed_list
	for (auto it = temp_list.begin(); it!=temp_list.end(); it++){
	  const GeomCell *current_cell1 = *it;
	  cell_list.remove(current_cell1);
	  suceed_list.push_back(current_cell1);
	}
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
    tempwire[i] = 0;
  }
  
  return diff;



  // if (allwire.size()>0){
  //   MergeGeomWire *wire = (MergeGeomWire*)allwire[0];
    
  //   std::cout << wire->get_allwire().size() << std::endl;
  //   //    MergeGeomWire *mwire = new MergeGeomWire(nwire,*wire);
  //   //mwire = 0;
  // }
  // return 0;
}




int WireCell2dToy::MergeToyTiling::further_merge(WireCell::GeomCellSelection &allcell, int ncell,int time_slice){
  WireCell::GeomCellSelection tempcell = allcell;
  allcell.clear();
  
  for (int i=0;i!=tempcell.size();i++){
    MergeGeomCell *cell = (MergeGeomCell*)tempcell[i];
      int flag=0;
      for (int k=0;k!=allcell.size();k++){
	if (((MergeGeomCell*)allcell[k])->AddCell(*cell)){
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
    tempcell[i] = 0;
  }
  
  return diff;
}


const WireCell::GeomCell* WireCell2dToy::MergeToyTiling::cell(const WireCell::GeomWireSelection& wires) const
{
  return 0;
}

ClassImp(WireCell2dToy::MergeToyTiling);
