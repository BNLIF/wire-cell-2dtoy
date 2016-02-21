#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"
#include "WireCellData/Singleton.h"
#include "WireCellData/TPCParams.h"

#include <cmath>

using namespace WireCell;

double WireCell2dToy::MergeToyTiling::time_convert = 0.32;
double WireCell2dToy::MergeToyTiling::dis_offset = -256;


WireCell2dToy::MergeToyTiling::MergeToyTiling(WireCell2dToy::ToyTiling& tiling, int time_slice, int merge_strategy, int flag_remerge)
{
  IsRemerged = false;

  ncell = tiling.get_ncell();
  wire_u = tiling.get_wire_u();
  wire_v = tiling.get_wire_v();
  wire_w = tiling.get_wire_w();
  //wirechargemap = tiling.wcmap();
  
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
    
    while ((cell_all.size() > 2 * wire_all.size() && cell_all.size() - wire_all.size() > 50) || 
	   (cell_all.size()>80 && cell_all.size() > 0.7 * wire_all.size()) ||
	   (cell_all.size()>100)){
      
      // two wire pitch
      dis += 2 * Singleton<TPCParams>::Instance().get_pitch();//6*units::mm;
      
      
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

  }else if (flag_remerge == 2){
    // merge things that are different by 1 wire gap ... 
    
    //    double dis = (3.0 + 1.5)*units::mm; //slightly bigger than 1 wire pitch
    double dis = Singleton<TPCParams>::Instance().get_pitch() * 1.5;
    for (int i=0;i!=cell_all.size();i++){
      MergeGeomCell* mcell = (MergeGeomCell*)cell_all.at(i);
      //std::cout << "Before: " << mcell->boundary().size() << " " ; 
      mcell->Organize_edge_boundary();
      //std::cout << mcell->boundary().size() << std::endl;
    }

    int current_ncell = cell_all.size();
    int prev_ncell = -1;
          
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
    
  }


  //find the edge cells ... 
  
  for (int i=0;i!=cell_all.size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*)cell_all.at(i);
    mcell->FindEdges(); // find ege wires
    for (int j=0;j!=mcell->get_edge_wires().size();j++){
      const GeomWire* edge_wire = mcell->get_edge_wires().at(j);
      for (int k=0;k!=tiling.cells(*edge_wire).size();k++){
	const GeomCell *tmp_cell = tiling.cells(*edge_wire).at(k);
	GeomCellSelection allcells = mcell->get_allcell();
	GeomCellSelection edgecells = mcell->get_edgecells();
	auto it1 = find(allcells.begin(),allcells.end(), tmp_cell);
	auto it2 = find(edgecells.begin(),edgecells.end(),tmp_cell);
	if (it1 != allcells.end() && it2 == edgecells.end())
	  mcell->get_edge_cells().push_back(tmp_cell);
      }
    }
    //std::cout <<"xin: " << mcell->get_allcell().size() << " " << mcell->get_edge_cells().size() << " " << mcell->get_edgecells().size() << " " << mcell->get_edge_wires().size() << std::endl;
    
  }
  //  deghost();
  
}


WireCell2dToy::MergeToyTiling::MergeToyTiling(const DetectorGDS& gds, WireCell2dToy::ToyTiling& tiling, int time_slice){
  
  //check cells
  // for (int i = 0; i!=tiling.get_allcell().size(); i++){
  //   const GeomCell *cell = tiling.get_allcell().at(i);
  //   std::cout << i << " " << tiling.wires(*cell).size() << std::endl;
  // }
  //check wires;
  // for (int i=0;i!=tiling.get_allwire().size();i++){
  //   const GeomWire *wire = tiling.get_allwire().at(i);
  //   std::cout << i << " " << tiling.cells(*wire).size() << std::endl;
  // }
  
  GeomWireSelection all_wires = tiling.get_allwire();
  GeomWireLMap wlmap;
  for (int i=0;i!=all_wires.size();i++){
    const GeomWire *wire = all_wires.at(i);
    GeomCellSelection cells = tiling.cells(*wire);
    if (cells.size() >0){
      GeomCellList abc(cells.begin(),cells.end());
      wlmap[wire] = abc;
    }
  }

  
  // do merge cells apa by apa
  int n_cryos = gds.ncryos();
  for (int ii=0;ii!=n_cryos;ii++){
    int n_apa = gds.napa(ii);
    for (int jj=0;jj!=n_apa;jj++){
      
      GeomWireSelection wires;
      for (int kk=0;kk!=all_wires.size();kk++){
	const GeomWire *wire = all_wires.at(kk);
	if (wire->cryo() == ii && wire->apa() == jj && 
	    tiling.cells(*wire).size()>0){
	  wires.push_back(wire);
	}
      }

      // std::cout << ii << " " << jj << " " << wires.size() << std::endl;
      
      int i = 0;
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
      } // end while
      
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
      
      //      std::cout << ii << " " << jj << " " << cell_all.size() << std::endl;
    
    }
  }
  
  
  form_wiremap(gds,tiling,time_slice);
  

  // form a new mapping for merged cell to merged cell 
  for (int i=0;i!=cell_all.size();i++){
    const MergeGeomCell *mcell = (MergeGeomCell*)cell_all[i];
    GeomCellSelection cells;
    mcmcsmap[mcell] = cells;

    for (int j=0;j!=cellmap[mcell].size();j++){
      const MergeGeomWire *mwire = (MergeGeomWire*)cellmap[mcell].at(j);
      for (int k=0;k!=wiremap[mwire].size();k++){
	const MergeGeomCell *mcell1 = (MergeGeomCell*)wiremap[mwire].at(k);
	if (mcell1 == mcell) continue;
	
	// judge if two cells share any wires and distance apart is 1 m
	double dis = sqrt(pow(mcell->center().y - mcell1->center().y,2)
			  + pow(mcell->center().z - mcell1->center().z,2));
	
	if (dis > 1 *units::m){
	  int flag_common_wire = 0;
	  for (int i1 = 0; i1!= cellmap1[mcell].size(); i1++){
	    auto it = find( cellmap1[mcell1].begin(), cellmap1[mcell1].end(),cellmap1[mcell].at(i1));
	    if (it != cellmap1[mcell1].end()){
	      flag_common_wire = 1;
	      break;
	    }
	  }
	  if (flag_common_wire == 0){
	    auto it = find(mcmcsmap[mcell].begin(),mcmcsmap[mcell].end(),mcell1);
	    if (it == mcmcsmap[mcell].end()){
	      mcmcsmap[mcell].push_back(mcell1);
	      // std::cout << "Check " 
	      // 		<< " " << mcell->center().y - mcell1->center().y
	      // 		<< " " << mcell->center().z - mcell1->center().z
	      // 		<< std::endl;
	    }
	  }
	}
      }
    }
    // std::cout << "Check " <<  mcmcsmap[mcell].size() << std::endl;
  }


  for (int i=0;i!=cell_all.size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*)cell_all.at(i);
    mcell->FindEdges(); // find ege wires
    // std::cout << "xin1: " << mcell->get_allcell().size() << " " << mcell->get_edge_cells().size() << " " << mcell->get_edgecells().size() << " " << mcell->get_edge_wires().size() << std::endl;

    for (int j=0;j!=mcell->get_edge_wires().size();j++){
      const GeomWire* edge_wire = mcell->get_edge_wires().at(j);
      //std::cout << "xin1: " << j << " " << tiling.cells(*edge_wire).size() << std::endl;
      for (int k=0;k!=tiling.cells(*edge_wire).size();k++){
	const GeomCell *tmp_cell = tiling.cells(*edge_wire).at(k);
	GeomCellSelection allcells = mcell->get_allcell();
	GeomCellSelection edgecells = mcell->get_edgecells();
	auto it1 = find(allcells.begin(),allcells.end(), tmp_cell);
	auto it2 = find(edgecells.begin(),edgecells.end(),tmp_cell);
	if (it1 != allcells.end() && it2 == edgecells.end())
	  mcell->get_edge_cells().push_back(tmp_cell);
      }
      
    }
    //std::cout <<"xin: " << mcell->get_allcell().size() << " " << mcell->get_edge_cells().size() << " " << mcell->get_edgecells().size() << " " << mcell->get_edge_wires().size() << std::endl;
    

  }

}


void WireCell2dToy::MergeToyTiling::form_wiremap(const DetectorGDS& gds, WireCell2dToy::ToyTiling& tiling, int time_slice){
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

  // std::cout <<cell_all.size() << " " << wire_all.size() << std::endl;
  while(further_mergewire(wire_all,50000,time_slice));
  //std::cout <<cell_all.size() << " " << wire_all.size() << std::endl;
  

  //deal with associations ... 
  
  
  


  
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
  
  // fill in charge
  std::map<int,float> ccmap = tiling.ccmap();
  for (int i=0;i!=wire_all.size();i++){
    MergeGeomWire *mwire = (MergeGeomWire*)wire_all[i];
    float charge = 0;
    GeomWireSelection wires = mwire->get_allwire();
    std::set<int> channels;
    for (int j=0;j!=wires.size();j++){
      int channel = wires.at(j)->channel();
      if (channels.find(channel) != channels.end()){
      }else{
	channels.insert(channel);
	charge += ccmap[channel];
      }
    }
    wirechargemap[mwire] = charge;
  }

  // 


  //check ... 
  // for (int i=0;i!=cell_all.size();i++){
  //   std::cout << i << " " << cellmap[cell_all.at(i)].size() << std::endl;
  // }
  // for (int i=0;i!=wire_all.size();i++){
  //   std::cout << i << " " << wiremap[wire_all.at(i)].size() << " " << wirechargemap[wire_all.at(i)]<< std::endl;
  // }

}




void WireCell2dToy::MergeToyTiling::deghost(GeomCellSelection &good_mcells){
  
  // fill the three-wires cells and two-wires cells 
  two_wires_cells.clear();
  three_wires_cells.clear();
  for (int i=0;i!=cell_all.size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*) cell_all.at(i);
    GeomWireSelection wires = cellmap[mcell];
    int nwire = 0;
    for (int j=0;j!=wires.size();j++){
      MergeGeomWire *mwire = (MergeGeomWire*)wires.at(j);
      if (wirechargemap[mwire] > 10){
	nwire ++;
      }else{
	//std::cout << "Xin1: " << wirechargemap[mwire] <<std::endl;
      }
    }
    //std::cout << "Xin: " << wires.size() << " " << nwire << std::endl;
    auto it = find(good_mcells.begin(),good_mcells.end(),mcell);

    // if (it != good_mcells.end()){
    //   std::cout << "Xin2: " << std::endl;
    // }

    if (nwire == 3 || it != good_mcells.end()){
      three_wires_cells.push_back(mcell);
    }else{
      two_wires_cells.push_back(mcell);
    }
  }
  // remove two_wires cells if they share anything with three_wires_cells
  GeomWireSelection used_wires;
  for (int i=0;i!=three_wires_cells.size();i++){
    GeomWireSelection wires = cellmap[three_wires_cells.at(i)];
    for (int j=0;j!=wires.size();j++){
      auto it = find(used_wires.begin(),used_wires.end(),wires.at(j));
      if (it == used_wires.end())
	used_wires.push_back(wires.at(j));
    }
  }
  GeomCellSelection to_be_removed;
  for (int i=0;i!=two_wires_cells.size();i++){
    GeomWireSelection wires = cellmap[two_wires_cells.at(i)];
    for (int j=0;j!=wires.size();j++){
      auto it = find(used_wires.begin(),used_wires.end(),wires.at(j));
      if (it != used_wires.end()){
	to_be_removed.push_back(two_wires_cells.at(i));
	break;
      }
    }
  }
  
  for (int i=0;i!=to_be_removed.size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*)to_be_removed.at(i);
    // remove from cell_all
    auto it1 = find(cell_all.begin(),cell_all.end(),mcell);
    cell_all.erase(it1);
    // remove from two_wire_cells
    auto it2 = find(two_wires_cells.begin(),two_wires_cells.end(),mcell);
    two_wires_cells.erase(it2);
    // remove from cellmap
    cellmap.erase(mcell);
    // remove from wiremap
    for (auto it3 = wiremap.begin(); it3 != wiremap.end(); it3++){
      auto it4 = find(it3->second.begin(),it3->second.end(),mcell);
      if (it4 != it3->second.end())
	it3->second.erase(it4);
    }


    // delete the cell
    delete mcell;
  }
  
    

  // first figure out how many cells are single wire cell
  int flag1 = 1;
  while(flag1==1){
    flag1 = 0;
    single_wire_cells.clear();
    GeomCellSelection to_be_removed_cells;
    for (int i=0;i!=cell_all.size();i++){
      MergeGeomCell *mcell = (MergeGeomCell*) cell_all.at(i);
      GeomWireSelection wires = cellmap[mcell];
      int flag = 0;
      for (int j=0;j!=wires.size();j++){
  	MergeGeomWire *mwire = (MergeGeomWire*)wires.at(j);
  	if (wirechargemap[mwire] > 10){
  	  if (wiremap[mwire].size()==1){
  	    flag = 1;
  	    break;
  	  }
  	}
      }
      if (flag==1)
  	single_wire_cells.push_back(mcell);
    }
  
    //figure out to be removed cells
    for (int i=0;i!=single_wire_cells.size();i++){
      MergeGeomCell *mcell = (MergeGeomCell*)single_wire_cells.at(i);
      GeomWireSelection wires = cellmap[mcell];
      for (int j=0;j!=wires.size();j++){
  	MergeGeomWire *mwire = (MergeGeomWire*)wires.at(j);
  	GeomCellSelection cells = wiremap[mwire];
  	for (int k=0;k!=cells.size();k++){
  	  auto it = find(single_wire_cells.begin(),single_wire_cells.end(),
  			 cells.at(k));
	  auto it2 = find(two_wires_cells.begin(),two_wires_cells.end(),
			  cells.at(k));
  	  if (it == single_wire_cells.end() && it2!=two_wires_cells.end()){
  	    auto it1 = find(to_be_removed_cells.begin(), to_be_removed_cells.end(), cells.at(k));
	    
  	    if (it1 == to_be_removed_cells.end())
  	      to_be_removed_cells.push_back(cells.at(k));
  	  }
  	}
      }
    }
  
    if (to_be_removed_cells.size() !=0 ) flag1 = 1; 
    
    //Now remove all cells and then reconstruct the two maps;
    for (int i=0;i!=to_be_removed_cells.size();i++){
      MergeGeomCell *mcell = (MergeGeomCell*)to_be_removed_cells.at(i);
      auto it = find(cell_all.begin(),cell_all.end(),mcell);
      cell_all.erase(it);
      delete mcell;
    }
    GeomCellMap cellmap_save = cellmap;
    cellmap.clear();
    wiremap.clear();
    for (int i=0;i!=cell_all.size();i++){
      cellmap[cell_all.at(i)] = cellmap_save[cell_all.at(i)];
      for (int j=0;j!=cellmap[cell_all.at(i)].size();j++){
  	if (wiremap.find(cellmap[cell_all.at(i)].at(j)) == wiremap.end()){
  	  GeomCellSelection cells;
  	  cells.push_back(cell_all.at(i));
  	  wiremap[cellmap[cell_all.at(i)].at(j)] = cells;
  	}else{
  	  wiremap[cellmap[cell_all.at(i)].at(j)].push_back(cell_all.at(i));
  	}
      }
    }

    if (flag1 == 0){
      //  std::cout << "All Cells " << cell_all.size() << " " << single_wire_cells.size() << " " << to_be_removed_cells.size() << std::endl;
    }

  }

  // find the remaining cells ...
  single_wire_cells.clear();
  for (int i=0;i!=cell_all.size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*) cell_all.at(i);
    GeomWireSelection wires = cellmap[mcell];
    int flag = 0;
    for (int j=0;j!=wires.size();j++){
      MergeGeomWire *mwire = (MergeGeomWire*)wires.at(j);
      if (wirechargemap[mwire] > 10){
	if (wiremap[mwire].size()==1){
	  flag = 1;
	  break;
	}
      }
    }
    if (flag==1)
      single_wire_cells.push_back(mcell);
  }

  GeomCellSelection remaining_cells;
  for (int i=0;i!=cell_all.size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*) cell_all.at(i);
    auto it1 = find(single_wire_cells.begin(),single_wire_cells.end(),mcell);
    auto it2 = find(three_wires_cells.begin(),three_wires_cells.end(),mcell);
    if (it1==single_wire_cells.end() && it2 == three_wires_cells.end()){
      remaining_cells.push_back(mcell);
    }
  }

  // std::cout << "Xin1: " << remaining_cells.size() << std::endl;

  // remove some of the remaining cells? 
  while(remaining_cells.size()>0){
    // std::cout << "Xin2: " << remaining_cells.size() << std::endl;
    // remove one cell
    MergeGeomCell *mcell = (MergeGeomCell*)remaining_cells.back();
    remaining_cells.pop_back();
    
    auto it = find(cell_all.begin(),cell_all.end(),mcell);
    cell_all.erase(it);
    delete mcell;

    // reconstruct map
    GeomCellMap cellmap_save = cellmap;
    cellmap.clear();
    wiremap.clear();
    for (int i=0;i!=cell_all.size();i++){
      cellmap[cell_all.at(i)] = cellmap_save[cell_all.at(i)];
      for (int j=0;j!=cellmap[cell_all.at(i)].size();j++){
  	if (wiremap.find(cellmap[cell_all.at(i)].at(j)) == wiremap.end()){
  	  GeomCellSelection cells;
  	  cells.push_back(cell_all.at(i));
  	  wiremap[cellmap[cell_all.at(i)].at(j)] = cells;
  	}else{
  	  wiremap[cellmap[cell_all.at(i)].at(j)].push_back(cell_all.at(i));
  	}
      }
    }
    
    GeomCellSelection to_be_removed_mcells;
    // re-calculate single wire and move them out of remaining_cells
    for (int i=0;i!=remaining_cells.size();i++){
      MergeGeomCell *mcell1 = (MergeGeomCell*) remaining_cells.at(i);
      GeomWireSelection wires = cellmap[mcell1];
      int flag = 0;
      for (int j=0;j!=wires.size();j++){
  	MergeGeomWire *mwire = (MergeGeomWire*)wires.at(j);
  	if (wirechargemap[mwire] > 10){
  	  if (wiremap[mwire].size()==1){
  	    flag = 1;
  	    break;
  	  }
  	}
      }
      if (flag==1){
	to_be_removed_mcells.push_back(mcell1); 
      }
    }
    
    for (int i=0;i!=to_be_removed_mcells.size();i++){
      auto it3 = find(remaining_cells.begin(),remaining_cells.end(),to_be_removed_mcells.at(i));
      remaining_cells.erase(it3);
    }

  }
  
  //  std::cout << "Check: " << cell_all.size() << " " << three_wires_cells.size() << " " << two_wires_cells.size() << std::endl;
  
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
  

  //std::cout << wire_all.size() << std::endl;
  WireChargeMap wcmap = tiling.wcmap();
  for (int i=0;i!=wire_all.size();i++){
    MergeGeomWire *mwire = (MergeGeomWire*)wire_all[i];
    GeomWireSelection wires = mwire->get_allwire();
    float charge = 0;
    for (int j=0;j!=wires.size();j++){
      charge += wcmap[wires.at(j)];
    }
    // std::cout << i << " " << charge << std::endl;
    wirechargemap[mwire] = charge;
  }
  //


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
