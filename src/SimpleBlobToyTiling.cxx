#include "WireCell2dToy/SimpleBlobToyTiling.h"
using namespace WireCell;
WireCell2dToy::SimpleBlobToyTiling::SimpleBlobToyTiling(WireCell2dToy::ToyTiling& toytiling1, WireCell2dToy::MergeToyTiling& mergetiling1, WireCell2dToy::ToyMatrix& toymatrix1,
							WireCell2dToy::MergeToyTiling& prev_mergetiling, WireCell2dToy::ToyMatrix& prev_toymatrix,
							WireCell2dToy::MergeToyTiling& next_mergetiling, WireCell2dToy::ToyMatrix& next_toymatrix){
  toytiling = &toytiling1;
  mergetiling = &mergetiling1;
  toymatrix = &toymatrix1;
  

  //save the first pass results
  WireCell::GeomCellSelectionV corner_mcells;
  WireCell::GeomCellSelectionV corner_smcells;
  

  nsimple_blob = 0;
  MergeGeomCell *mcorner_cell;
  if (toymatrix->GetSimpleBlobReduction()){
    GeomCellSelection mcells = mergetiling->get_allcell();
    for (int i=0;i!=mcells.size();i++){
      MergeGeomCell *mcell = (MergeGeomCell*)mcells.at(i);
      if (mcell->IsSimpleBlob() && mcell->IsBlob()){
	CellIndexMap indexmap = mcell->get_cornercells_index();
	
	for (int j=0;j!=12;j++){
	  int index1 = mcell->index1(j);
	  int index2 = mcell->index2(j);
	  GeomCellSelection cells = mcell->get_cornercells(index1,index2);
	  if (cells.size()!=0){
	    //std::cout << index1 << " " << index2 << " " << cells.size() << std::endl;
	    int flag = -1;
	    for (int k=0;k!=cells.size();k++){
	      if (k==0){
		mcorner_cell = new MergeGeomCell(10000,*cells.at(k));
		cell_rank[mcorner_cell] = 0;
		mcorner_cell->ewires.insert(index1);
		mcorner_cell->ewires.insert(index2);
	      }else{
		mcorner_cell->AddCell(*cells.at(k));
	      }
	      if (indexmap[cells.at(k)] == 3) {
		flag=k;
		cell_rank[mcorner_cell] += 5; // should be adjusted later?
	      }
	    }
	    int flag1 = 0;
	    
	    // if (flag!=-1){
	    //   for (int k=0;k!=corner_smcells[nsimple_blob].size();k++){
	    // 	// see if any previous merged cell contain this cell
	    // 	// if so, add all cells to this merged cell 
	    // 	// and delete this cell
	    // 	// change flag1 value;
	    // 	GeomCellSelection tcells = ((MergeGeomCell*)corner_smcells[nsimple_blob].at(k))->get_allcell();
	    // 	auto it = find(tcells.begin(),tcells.end(),cells.at(flag));
	    // 	if (it!=tcells.end()){
	    // 	  delete mcorner_cell;
	    // 	  flag1 = 1;
	    // 	  ((MergeGeomCell*)corner_smcells[nsimple_blob].at(k))->ewires.insert(index1);
	    // 	  ((MergeGeomCell*)corner_smcells[nsimple_blob].at(k))->ewires.insert(index2);
	    // 	  for (int kk = 0;kk!=cells.size();kk++){
	    // 	    auto it1 = find(tcells.begin(),tcells.end(),cells.at(kk));
	    // 	    if (it1==tcells.end()){
	    // 	      ((MergeGeomCell*)corner_smcells[nsimple_blob].at(k))->AddCell(*cells.at(kk));
	    // 	    }
	    // 	  }
	    // 	}
	    //   }
	    // }
	    //save the special cells ...  //need to add time later ... 
	    if (flag1==0){
	      if (flag!=-1){
		if (corner_smcells.size() < nsimple_blob + 1){
		  GeomCellSelection cellss;
		  corner_smcells.push_back(cellss);
		}
		corner_smcells[nsimple_blob].push_back(mcorner_cell);
		
	      }else{
		if (corner_mcells.size() < nsimple_blob + 1){
		  GeomCellSelection cellss;
		  corner_mcells.push_back(cellss);
		}
		corner_mcells[nsimple_blob].push_back(mcorner_cell);
		  
	      }
	    }
	    
	  }

	}

	// going through smcells array and properly calculate index
	for (auto it = corner_smcells[nsimple_blob].begin(); it!=corner_smcells[nsimple_blob].end();it++){
	  MergeGeomCell *mcell = (MergeGeomCell*) (*it);
	  for (int j=0;j!=prev_mergetiling.get_allcell().size();j++){
	    MergeGeomCell *prev_mcell = (MergeGeomCell*) prev_mergetiling.get_allcell().at(j);
	    if (prev_toymatrix.Get_Cell_Charge(prev_mcell)>2000){
	      int temp_val = mcell->Overlap1(*prev_mcell);
	      if (temp_val){
		cell_rank[mcell] += temp_val;
		break;
	      }
	    }
	  }
	  //see next time slice
	  for (int j=0;j!=next_mergetiling.get_allcell().size();j++){
	    MergeGeomCell *next_mcell = (MergeGeomCell*) next_mergetiling.get_allcell().at(j);
	    if (next_toymatrix.Get_Cell_Charge(next_mcell)>2000){
	      int temp_val = mcell->Overlap1(*next_mcell);
	      if (temp_val){
		cell_rank[mcell] += temp_val;
		break;
	      }
	    }
	  }
	}
	
	// going through the mcells array and judge if any of them are special
	// if so move to the smcells
	for (auto it = corner_mcells[nsimple_blob].begin(); it!=corner_mcells[nsimple_blob].end();it++){
	  MergeGeomCell *mcell = (MergeGeomCell*) (*it);
	  int flag = 0;
	  //see previous time slice
	  for (int j=0;j!=prev_mergetiling.get_allcell().size();j++){
	    MergeGeomCell *prev_mcell = (MergeGeomCell*) prev_mergetiling.get_allcell().at(j);
	    if (prev_toymatrix.Get_Cell_Charge(prev_mcell)>2000){
	      int temp_val = mcell->Overlap1(*prev_mcell);
	      if (temp_val){
		cell_rank[mcell] += temp_val;
		flag = 1;
	      	break;
	      }
	    }
	  }
	  //see next time slice
	  for (int j=0;j!=next_mergetiling.get_allcell().size();j++){
	    MergeGeomCell *next_mcell = (MergeGeomCell*) next_mergetiling.get_allcell().at(j);
	    if (next_toymatrix.Get_Cell_Charge(next_mcell)>2000){
	      int temp_val = mcell->Overlap1(*next_mcell);
	      if (temp_val){
		cell_rank[mcell] += temp_val;
		flag = 1;
		break;
	      }
	    }
	  }
	  if (flag==1){
	    corner_smcells[nsimple_blob].push_back(mcell);
	  }
	}
	for (auto it = corner_smcells[nsimple_blob].begin(); it!=corner_smcells[nsimple_blob].end();it++){
	  MergeGeomCell *mcell = (MergeGeomCell*) (*it);
	  auto itq = find( corner_mcells[nsimple_blob].begin(),  corner_mcells[nsimple_blob].end(), mcell);
	  if (itq != corner_mcells[nsimple_blob].end()){
	    corner_mcells[nsimple_blob].erase(itq);
	  }
	}

	



	//now put everything into hypo_ccells
	if (hypo_ccells.size() < nsimple_blob + 1){
	  GeomCellSelectionV qcv; // vector of vector
	  
	  // first deal with smcells
	  // do a loop, and insert the first element, 
	  // do a loop, and merge anything that can be merged, erase, if none
	  // insert the second element until nothing in the smcells;
	  int flag_new = 1; //insert?
	  int flag_merge = 0; //merge?
	  while(corner_smcells[nsimple_blob].size()){
	    if (flag_merge == 1 && corner_smcells[nsimple_blob].size()!=0){
	      int temp_flag = 0;
	      //do the merge part
	      for (int j=0;j!=corner_smcells[nsimple_blob].size();j++){
	        MergeGeomCell* mcell1 = (MergeGeomCell*) corner_smcells[nsimple_blob].at(j);
		for (int k = 0;k!=qcv.size();k++){
		  for (int kk=0;kk!=qcv.at(k).size();kk++){
		    MergeGeomCell* mcell2 = (MergeGeomCell*) qcv.at(k).at(kk);
		    if (mcell1->Overlap1(*mcell2)){
		      qcv.at(k).push_back(mcell1);
		      corner_smcells[nsimple_blob].erase(corner_smcells[nsimple_blob].begin() + j);
		      temp_flag = 1;
		      break;
		    }
		  }
		  if (temp_flag==1) break;
		}
		if (temp_flag==1) break;
	      }
	      if (temp_flag == 0){
		flag_merge = 0;
		flag_new = 1;
	      }
	    }
	    
	    
	    if (flag_new == 1 && corner_smcells[nsimple_blob].size()!=0){
	      //insert an element and delete one ... 
	      GeomCellSelection qc;
	      qc.push_back(corner_smcells[nsimple_blob].at(corner_smcells[nsimple_blob].size()-1));
	      qcv.push_back(qc);
	      corner_smcells[nsimple_blob].pop_back();
	      flag_new = 0;
	      flag_merge = 1;
	    }
	  }
	  
	  
	  // then deal with mcells
	  if (corner_mcells[nsimple_blob].size()!=0)
	    qcv.push_back(corner_mcells[nsimple_blob]);

	  // flag_new = 1; //insert?
	  // flag_merge = 0; //merge?
	  // while(corner_mcells[nsimple_blob].size()){
	  //   if (flag_merge == 1 && corner_mcells[nsimple_blob].size()!=0){
	  //     int temp_flag = 0;
	  //     //do the merge part
	  //     for (int j=0;j!=corner_mcells[nsimple_blob].size();j++){
	  //       MergeGeomCell* mcell1 = (MergeGeomCell*) corner_mcells[nsimple_blob].at(j);
	  // 	for (int k = 0;k!=qcv.size();k++){
	  // 	  for (int kk=0;kk!=qcv.at(k).size();kk++){
	  // 	    MergeGeomCell* mcell2 = (MergeGeomCell*) qcv.at(k).at(kk);
	  // 	    if (mcell1->Overlap1(*mcell2)){
	  // 	      qcv.at(k).push_back(mcell1);
	  // 	      corner_mcells[nsimple_blob].erase(corner_mcells[nsimple_blob].begin() + j);
	  // 	      temp_flag = 1;
	  // 	      break;
	  // 	    }
	  // 	  }
	  // 	  if (temp_flag==1) break;
	  // 	}
	  // 	if (temp_flag==1) break;
	  //     }
	  //     if (temp_flag == 0){
	  // 	flag_merge = 0;
	  // 	flag_new = 1;
	  //     }
	  //   }
	  //   if (flag_new == 1 && corner_mcells[nsimple_blob].size()!=0){
	  //     //insert an element and delete one ... 
	  //     GeomCellSelection qc;
	  //     qc.push_back(corner_mcells[nsimple_blob].at(corner_mcells[nsimple_blob].size()-1));
	  //     qcv.push_back(qc);
	  //     corner_mcells[nsimple_blob].pop_back();
	  //     flag_new = 0;
	  //     flag_merge = 1;
	  //   }
	  // }

	  
	  hypo_ccells.push_back(qcv); //save things in ... 
	}


	nsimple_blob ++;
	// if (nsimple_blob >= 10) {
	//   break;
	// }
	
      }
    }

    
    
    std::cout << "SimpleBlobTiling: "<< nsimple_blob << " " << corner_smcells[0].size() << " " << corner_mcells[0].size() << " " << hypo_ccells.at(0).size() << std::endl;
    for (int j=0;j!=hypo_ccells.at(0).size();j++){
      std::cout << hypo_ccells.at(0).at(j).size() << " ";
      for (int k = 0;k!= hypo_ccells.at(0).at(j).size();k++){
	std::cout << cell_rank[hypo_ccells.at(0).at(j).at(k)] << " ";
      }
      std::cout << std::endl;
    }
    // for (int j=0;j!=corner_smcells[0].size();j++){
    //   std::cout << ((MergeGeomCell*)corner_smcells[0].at(j))->get_allcell().size() << " ";
    //   for (auto it = ((MergeGeomCell*)corner_smcells[0].at(j))->ewires.begin(); it!= ((MergeGeomCell*)corner_smcells[0].at(j))->ewires.end(); it++){
    // 	std::cout << *it << " ";
    //   }
    //   std::cout << std::endl;
    // }

  }
}

WireCell2dToy::SimpleBlobToyTiling::~SimpleBlobToyTiling(){
  // for (int i=0;i!=nsimple_blob;i++){
  //   for (int j=0;j!=corner_mcells[i].size();j++){
  //     delete corner_mcells[i].at(j);
  //   }
  //   for (int j=0;j!=corner_smcells[i].size();j++){
  //     delete corner_smcells[i].at(j);
  //   }
  // }

  for (int i = 0 ;i!=hypo_ccells.size(); i++){
    for (int j = 0; j!=hypo_ccells.at(i).size();j++){
      for (int k=0;k!=hypo_ccells.at(i).at(j).size();k++){
	delete hypo_ccells.at(i).at(j).at(k);
      }
    }
  }
}
