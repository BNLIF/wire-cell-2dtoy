#include "WireCell2dToy/BlobToyTiling.h"


#include <cmath>

using namespace WireCell;

WireCell2dToy::BlobToyTiling::BlobToyTiling(WireCell2dToy::ToyTiling& toytiling, WireCell2dToy::MergeToyTiling& mergetiling, WireCell2dToy::ToyMatrix& toymatrix, int time_slice, int num_merge_wire){
  num_wire = num_merge_wire;
  tiling = &toytiling;
  //all the single wires
  wirechargemap = mergetiling.wcmap();

  GeomWireSet1 tmp_set;
  for (int i=0;i!= mergetiling.get_wire_u().size();i++){
    tmp_set.insert(mergetiling.get_wire_u().at(i));
  }
  for (auto it = tmp_set.begin();it!=tmp_set.end();it++){
    wire_u.push_back(*it) ;
  }
  tmp_set.clear();
  
  for (int i=0;i!= mergetiling.get_wire_v().size();i++){
    tmp_set.insert(mergetiling.get_wire_v().at(i));
  }
  for (auto it = tmp_set.begin();it!=tmp_set.end();it++){
    wire_v.push_back(*it) ;
  }
  tmp_set.clear();
  
  for (int i=0;i!= mergetiling.get_wire_w().size();i++){
    tmp_set.insert(mergetiling.get_wire_w().at(i));
  }
  for (auto it = tmp_set.begin();it!=tmp_set.end();it++){
    wire_w.push_back(*it) ;
  }
  tmp_set.clear();
  


  //initialize ncell;
  ncell = 0;

  GeomCellSelection cellall = mergetiling.get_allcell();
  for (int i=0;i!=cellall.size();i++){
    MergeGeomCell *mcell = (WireCell::MergeGeomCell*)cellall[i];
    if (toymatrix.Get_Cell_Charge(mcell)+toymatrix.Get_Cell_Charge(mcell,2)>0){
      GeomCellSelection acell = mcell->get_allcell();
      //put all the single cells in here
      remain_cells.insert(remain_cells.end(),acell.begin(),acell.end());
      ncell += acell.size();
    }
  }
  
  
  //std::cout << ncell << std::endl;
  //std::cout << "U " << std::endl; 
  int ident_wire = 100000;
  int flag = 0;
  MergeGeomWire *mwire = 0;
  
  for (int i=0;i!=wire_u.size();i++){
    if (flag == 0 ){
      mwire = new MergeGeomWire(ident_wire,*wire_u[i]);
      ident_wire ++;
      flag = wire_u[i]->index();
    }else{
      int abc = wire_u.size()/num_wire/2.;
      if (abc < 1) abc = 1;
      if (mwire->get_allwire().size()<abc && fabs(wire_u[i]->index()-flag)==1){
	mwire->AddWire(*wire_u[i]);
	flag = wire_u[i]->index();
      }else{
	wire_u1.push_back(mwire);
	wire_all.push_back(mwire);
	mwire = new MergeGeomWire(ident_wire,*wire_u[i]);
	ident_wire ++;
	flag = wire_u[i]->index();
      }
    }
    //std::cout << wire_u[i]->index() << std::endl;
  }
  wire_u1.push_back(mwire);
  wire_all.push_back(mwire);

  //  std::cout << wire_all.size() << std::endl;
 
  
  for (int i=0;i!=wire_v.size();i++){
    if (flag == 0 ){
      mwire = new MergeGeomWire(ident_wire,*wire_v[i]);
      ident_wire ++;
      flag = wire_v[i]->index();
    }else{
      int abc = wire_v.size()/num_wire/2.;
      if (abc < 1) abc = 1;
      if (mwire->get_allwire().size()<abc && fabs(wire_v[i]->index()-flag)==1){
	mwire->AddWire(*wire_v[i]);
	flag = wire_v[i]->index();
      }else{
	wire_v1.push_back(mwire);
	wire_all.push_back(mwire);
	mwire = new MergeGeomWire(ident_wire,*wire_v[i]);
	ident_wire ++;
	flag = wire_v[i]->index();
      }
    }
    //std::cout << wire_v[i]->index() << std::endl;
  }
  wire_v1.push_back(mwire);
  wire_all.push_back(mwire);
  
  for (int i=0;i!=wire_w.size();i++){
    if (flag == 0 ){
      mwire = new MergeGeomWire(ident_wire,*wire_w[i]);
      ident_wire ++;
      flag = wire_w[i]->index();
    }else{
      int abc = wire_w.size()/num_wire;
      if (abc < 1) abc = 1;
      if (mwire->get_allwire().size()< abc&& fabs(wire_w[i]->index()-flag)==1){
	mwire->AddWire(*wire_w[i]);
	flag = wire_w[i]->index();
      }else{
	wire_w1.push_back(mwire);
	wire_all.push_back(mwire);
	mwire = new MergeGeomWire(ident_wire,*wire_w[i]);
	ident_wire ++;
	flag = wire_w[i]->index();
      }
    }
    //std::cout << wire_w[i]->index() << std::endl;
  }
  wire_w1.push_back(mwire);
  wire_all.push_back(mwire);
 
  //given a set of merged wire, find the single cells that are contained in remain_cells and form a merged cell
  int ident_cell = 1000;
  for (int i =0; i!=wire_u1.size();i++){
    for (int j=0;j!=wire_v1.size();j++){
      for (int k=0;k!=wire_w1.size();k++){
	MergeGeomCell *mcell = FormMergeCell((MergeGeomWire*)wire_u1[i],(MergeGeomWire*)wire_v1[j], (MergeGeomWire*)wire_w1[k],ident_cell,time_slice);
	if (mcell!=0){
	  cell_all.push_back(mcell);
	  GeomWireSelection wiresel;
	  wiresel.push_back(wire_u1[i]);
	  wiresel.push_back(wire_v1[j]);
	  wiresel.push_back(wire_w1[k]);	
	  cellmap[mcell]=wiresel;

	  //fill wiremap
	  if (wiremap.find(wire_u1[i]) == wiremap.end()){
	    //not found
	    GeomCellSelection cellsel;
	    cellsel.push_back(mcell);
	    wiremap[wire_u1[i]]=cellsel;
	  }else{
	    //found
	    wiremap[wire_u1[i]].push_back(mcell);
	  }

	  if (wiremap.find(wire_v1[j]) == wiremap.end()){
	    //not found
	    GeomCellSelection cellsel;
	    cellsel.push_back(mcell);
	    wiremap[wire_v1[j]]=cellsel;
	  }else{
	    //found
	    wiremap[wire_v1[j]].push_back(mcell);
	  }

	  if (wiremap.find(wire_w1[k]) == wiremap.end()){
	    //not found
	    GeomCellSelection cellsel;
	    cellsel.push_back(mcell);
	    wiremap[wire_w1[k]]=cellsel;
	  }else{
	    //found
	    wiremap[wire_w1[k]].push_back(mcell);
	  }

	}
      }
    }
  }

  //std::cout << cell_all.size() << " " << wire_all.size() << std::endl;

}


MergeGeomCell* WireCell2dToy::BlobToyTiling::FormMergeCell(MergeGeomWire* mwireu, MergeGeomWire* mwirev, MergeGeomWire* mwirew,int ident_cell,int time_slice) {
  MergeGeomCell *mcell = 0;

  for (int i =0; i!= mwireu->get_allwire().size();i++){
    for (int j =0; j!= mwirev->get_allwire().size();j++){
      for (int k =0; k!= mwirew->get_allwire().size();k++){
	GeomWireSelection wires;
	wires.push_back(mwireu->get_allwire().at(i));
	wires.push_back(mwirev->get_allwire().at(j));
	wires.push_back(mwirew->get_allwire().at(k));

	const GeomCell* cell = tiling->cell(wires);
	
	if (cell!=0){
	  auto it = find(remain_cells.begin(),remain_cells.end(),cell);
	  if (it!=remain_cells.end()){
	    if (mcell == 0){
	      mcell = new MergeGeomCell(ident_cell,*cell);
	      mcell->SetTimeSlice(time_slice);
	    }else{
	      mcell->AddCell(*cell);
	    }
	  }
	}
      }
    }
  }

  

  return mcell;
}


const WireCell::GeomCell* WireCell2dToy::BlobToyTiling::cell(const WireCell::GeomWireSelection& wires) const
{
  return 0;
}


WireCell2dToy::BlobToyTiling::~BlobToyTiling(){
  for (int i=0;i!=cell_all.size();i++){
    delete cell_all[i];
  }

  for (int i=0;i!=wire_all.size();i++){
    delete wire_all[i];
  }

  wire_u.clear();
  wire_v.clear();
  wire_w.clear();
  
  wire_u1.clear();
  wire_v1.clear();
  wire_w1.clear();
  
  wire_all.clear();
  cell_all.clear();
  cellmap.clear();
  wiremap.clear();

  remain_cells.clear();
  wirechargemap.clear();

}

ClassImp(WireCell2dToy::BlobToyTiling);
