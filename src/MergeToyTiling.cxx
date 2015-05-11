#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"

#include <cmath>

using namespace WireCell;

WireCell2dToy::MergeToyTiling::MergeToyTiling(WireCell2dToy::ToyTiling& tiling){
  ncell = tiling.get_ncell();
  wire_u = tiling.get_wire_u();
  wire_v = tiling.get_wire_v();
  wire_w = tiling.get_wire_w();
  
  // goal is to create merged version of 
  // cell_all
  // cellmap
  // wiremap

  // //start with wire_u
  for (int i =0;i!=wire_u.size();i++){
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
      	ncell++;
	cell_all.push_back(mcell);
      }
    }
  }
  
  while(further_merge(cell_all,tiling.get_ncell()));
    


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
	      ident_wire++;
	      flag = 1;
	    }else {
	      mwire->AddWire(*wires[nwire]);
	    }
  	  }
  	}
      }
      
      // if ( flag!= call.size()){
      // 	for (int j=0;j!=call.size();j++){
      // 	  GeomWireSelection wires = tiling.wires(*call[j]);
      // 	  for (int nwire = 0; nwire!=wires.size();nwire++){
      // 	    std::cout << wires[nwire]->plane() << " " << wires[nwire]->ident() << " " << plane << std::endl;
      // 	  }
      // 	}
      // }
      
      // if (flag==1)
      //  	std::cout << mwire->get_allwire().size() << " " << call.size() << std::endl;
      
      wire_all.push_back(mwire);
    }
  }

 

  

  while(further_mergewire(wire_all,50000));
  


  // // // Now construct the map
  // // for (int i=0;i!=cell_all.size();i++){
  // //   const MergeGeomCell *cell = (MergeGeomCell*)cell_all[i];
  // //   GeomWireSelection wiresel;
        
  // //   for (int j=0;j!=cell->get_allcell().size();j++){
  // //     const GeomCell *scell = cell->get_allcell()[j];
  // //     //std::cout << i << " " << scell->ident()<< " " << tiling.wires(*scell).size() << std::endl;
  // //     for (int k=0;k!=tiling.wires(*scell).size();k++){
  // // 	const GeomWire *wire = tiling.wires(*scell)[k];
  // // 	wiresel.push_back(wire);

  // // 	//also do the wiremap
  // // 	if (wiremap1.find(wire) == wiremap1.end()){
  // // 	  GeomCellSelection cellsel;
  // // 	  cellsel.push_back(cell);
  // // 	  wiremap1[wire]= cellsel;
  // // 	}else{
  // // 	  int flag = 0;
  // // 	  for (int n=0;n!=wiremap1[wire].size();n++){
  // // 	    if (cell == wiremap1[wire].at(n)){
  // // 	      flag = 1;
  // // 	      break;
  // // 	    }
  // // 	  }
  // // 	  if(flag==0){
  // // 	    wiremap1[wire].push_back(cell);
  // // 	  }
  // // 	}

  // //     }
  // //   }

  // //   cellmap1[cell] = wiresel;
  // // }


  // //Now construct the real map
  
  

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
  wirechargemap.clear();
}


int WireCell2dToy::MergeToyTiling::further_mergewire(WireCell::GeomWireSelection &allwire, int nwire){

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
