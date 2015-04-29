#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/MergeGeomCell.h"
#include <cmath>

using namespace WireCell;

WireCell2dToy::MergeToyTiling::MergeToyTiling(WireCell2dToy::ToyTiling tiling){
  ncell = tiling.get_ncell();
  wire_u = tiling.get_wire_u();
  wire_v = tiling.get_wire_v();
  wire_w = tiling.get_wire_w();
  wire_all = tiling.get_wire_all();

  // goal is to create merged version of 
  // cell_all
  // cellmap
  // wiremap

  //start with wire_u
  for (int i =0;i!=wire_u.size();i++){
    for (int j=0;j!=tiling.cells(*wire_u[i]).size();j++){
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
	ncell++;
	cell_all.push_back(mcell);
      }
    }
  }
  
  while(further_merge(cell_all,tiling.get_ncell()));
    
}

int WireCell2dToy::MergeToyTiling::further_merge(WireCell::GeomCellSelection &allcell, int ncell){
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
