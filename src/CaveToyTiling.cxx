#include "WireCell2dToy/CaveToyTiling.h"

#include <cmath>

using namespace WireCell;

WireCell2dToy::CaveToyTiling::CaveToyTiling(WireCell2dToy::ToyTiling *toytiling1, WireCell2dToy::MergeToyTiling& mergetiling, WireCell2dToy::ToyMatrix& toymatrix){
  toytiling = toytiling1;
  // create cell_all_save, 
  GeomCellCellMap tmp_ccmap;
  GeomCellCellMap tmp_ccmap_inv;
  GeomCellSelection icellall = mergetiling.get_allcell();
  for (int i=0; i!= icellall.size();i++){
    const GeomCell* cell = icellall.at(i);
    double charge = toymatrix.Get_Cell_Charge(cell,1);
    double charge_err = toymatrix.Get_Cell_Charge(cell,2);
    
    if (charge_err !=0 && charge !=0){
      MergeGeomCell *mcell = new MergeGeomCell(10000,*((const MergeGeomCell*)cell));
      cell_all_save.push_back(mcell); //only save the cells that are good
      tmp_ccmap[mcell] = (GeomCell*) cell;
      tmp_ccmap_inv[cell] = mcell;
    }
  }
  //create wire_all_save, 
  GeomWireWireMap tmp_wwmap;
  GeomWireWireMap tmp_wwmap_inv;
  GeomWireSelection iwireall = mergetiling.get_allwire();
  for (int i=0;i!=iwireall.size();i++){
    const GeomWire* wire = iwireall.at(i);
    MergeGeomWire *mwire = new MergeGeomWire(10000,*((const MergeGeomWire*)wire));
    wire_all_save.push_back(mwire);
    tmp_wwmap[mwire] = (GeomWire*) wire;
    tmp_wwmap_inv[wire] = mwire;
  }
  // and do map of (need to have constructor in mergegeomcell and mergegeomwire)  
  // cellmap_save and wiremap_save  (difficult?)
  // need to have a mapping of wire to wire and cell to cell
  for (int i=0;i!=cell_all_save.size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*) cell_all_save.at(i);
    GeomCell *cell = tmp_ccmap[mcell];
    GeomWireSelection wires = mergetiling.wires(*cell);
    GeomWireSelection new_wires;
    for (int j=0;j!=wires.size();j++){
      const GeomWire* wire = wires.at(j);
      GeomWire *mwire = tmp_wwmap_inv[wire];
      new_wires.push_back(mwire);
    }
    cellmap_save[mcell] = new_wires;
  }
  for (int i=0;i!=wire_all_save.size();i++){
    MergeGeomWire *mwire = (MergeGeomWire*)wire_all_save.at(i);
    GeomWire *wire = tmp_wwmap[mwire];
    GeomCellSelection cells = mergetiling.cells(*wire);
    GeomCellSelection new_cells;
    for (int j=0;j!=cells.size();j++){
      const GeomCell* cell = cells.at(j);
      if (tmp_ccmap_inv.count(cell)==1){
	GeomCell *mcell = tmp_ccmap_inv[cell];
	new_cells.push_back(mcell);
      }
    }
    wiremap_save[mwire] = new_cells;
  }
  
  
  // for (int i=0;i!=wire_all_save.size();i++){
  //   GeomCellSelection cells = wiremap_save[wire_all_save.at(i)];
  //   if (cells.size()==0){ // add all the cells back
  //     GeomWire* wire = tmp_wwmap[wire_all_save.at(i)];
  //     GeomCellSelection cells1 = mergetiling.cells(*wire);
  //     std::cout << "abc" << cells1.size() << std::endl; 
  //     //Update cells
  //     for (int j=0;j!=cells1.size();j++){
  //      	const GeomCell* cell = cells1.at(j);
  //      	if (tmp_ccmap_inv.count(cell)==0){
  // 	  std::cout << cell->center().y << " " << cell->center().z << std::endl;
  //      	  //MergeGeomCell *mcell = new MergeGeomCell(10000,*((const MergeGeomCell*)cell));
  //     // 	  // update all cell
  //      	  //cell_all_save.push_back(mcell); //only save the cells that are good 
  //      	  //tmp_ccmap[mcell] = (GeomCell*) cell;
  //      	  //tmp_ccmap_inv[cell] = mcell;
  //     // 	  cells.push_back(mcell);
  //     // 	  //update cell map;
  //     // 	  GeomWireSelection wires = mergetiling.wires(*cell);
  //     // 	  GeomWireSelection new_wires;
  //     // 	  for (int k=0;k!=wires.size();k++){
  //     // 	    const GeomWire* wire1 = wires.at(k);
  //     // 	    GeomWire *mwire = tmp_wwmap_inv[wire1];
  //     // 	    new_wires.push_back(mwire);
  //     // 	  }
  //     // 	  cellmap_save[mcell] = new_wires;
  //      	}
  //     }
  //   }
  // }
  // check
  // std::cout << cell_all_save.size() << " " << wire_all_save.size() << std::endl;
  // for (int i=0;i!=cell_all_save.size();i++){
  //   GeomWireSelection wires = cellmap_save[cell_all_save.at(i)];
  //   std::cout << wires.size() << " " << wires.at(0)->ident() << std::endl;
  // }
  // for (int i=0;i!=wire_all_save.size();i++){
  //   GeomCellSelection cells = wiremap_save[wire_all_save.at(i)];
  //   std::cout << cells.size() << " " << std::endl;
  // }
  //std::cout << cells.size() << " " <<  std::endl;

  
  
  // and then create cell_to_remove all the cells at the boundary
  // how???
  for (int i=0;i!=cell_all_save.size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*)cell_all_save[i];
    mcell->FindEdges(); 
    //    std::cout << mcell->cross_section() << std::endl; //area ranked from small to big
    // rank from big cell to smaller cell 
    // then rank smaller area to larger area
  }
  

  
  // create a function to move current to save 
  // create a function to construct current from save
  // key is to create new wires from merged cells and do all the mapping

  //when the cell_to_remove is empty, create a new one
  //move current_cell to somewhere ... 
}

WireCell2dToy::CaveToyTiling::~CaveToyTiling(){
  for (int i=0;i!=cell_all.size();i++){
    delete cell_all[i];
  }

  for (int i=0;i!=wire_all.size();i++){
    delete wire_all[i];
  }

  for (int i=0;i!=cell_all_save.size();i++){
    delete cell_all_save[i];
  }

  for (int i=0;i!=wire_all_save.size();i++){
    delete wire_all_save[i];
  }
}

const WireCell::GeomCell* WireCell2dToy::CaveToyTiling::cell(const WireCell::GeomWireSelection& wires) const
{
  return 0;
}

