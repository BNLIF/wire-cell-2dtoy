#include "WireCell2dToy/SimpleBlobToyTiling.h"
using namespace WireCell;
WireCell2dToy::SimpleBlobToyTiling::SimpleBlobToyTiling(WireCell2dToy::ToyTiling& toytiling1, WireCell2dToy::MergeToyTiling& mergetiling1, WireCell2dToy::ToyMatrix& toymatrix1,
							WireCell2dToy::MergeToyTiling& prev_mergetiling, WireCell2dToy::ToyMatrix& prev_toymatrix,
							WireCell2dToy::MergeToyTiling& next_mergetiling, WireCell2dToy::ToyMatrix& next_toymatrix){
  toytiling = &toytiling1;
  mergetiling = &mergetiling1;
  toymatrix = &toymatrix1;
  
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
		mcorner_cell->ewires.insert(index1);
		mcorner_cell->ewires.insert(index2);
	      }else{
		mcorner_cell->AddCell(*cells.at(k));
	      }
	      if (indexmap[cells.at(k)] == 3) flag=k;
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
		// if (corner_smcells.at(nsimple_blob) == corner_smcells.end()){
		//   GeomCellSelection cellss;
		//   corner_smcells.push_back(cellss);
		// }
		corner_smcells[nsimple_blob].push_back(mcorner_cell);
		
	      }else{
		// if (corner_mcells.at(nsimple_blob) == corner_mcells.end()){
		//   GeomCellSelection cellss;
		//   corner_mcells.push_back(cellss);
		// }
		corner_mcells[nsimple_blob].push_back(mcorner_cell);
		  
	      }
	    }
	    
	  }

	}

	// going through the mcells array and judge if any of them are special
	// if so move to the smcells
	for (auto it = corner_mcells[nsimple_blob].end()-1; it>=corner_mcells[nsimple_blob].begin();it--){
	  MergeGeomCell *mcell = (MergeGeomCell*) (*it);
	  
	  // for (int j=0;j!=mcell->get_allcell().size();j++){
	  //   std::cout << j << " " << mcell->get_allcell().at(j)->center().y << " " << mcell->get_allcell().at(j)->center().z << " " << std::endl;
	  // }

	  int flag = 0;
	  //see previous time slice
	  for (int j=0;j!=prev_mergetiling.get_allcell().size();j++){
	    MergeGeomCell *prev_mcell = (MergeGeomCell*) prev_mergetiling.get_allcell().at(j);
	    if (prev_toymatrix.Get_Cell_Charge(prev_mcell)>2000){
	      
	      // for (int k=0;k!=prev_mcell->get_allcell().size();k++){
	      // 	std::cout << k << " " << prev_mcell->get_allcell().at(k)->center().y << " " << prev_mcell->get_allcell().at(k)->center().z << " " << std::endl;
	      // }
	      // std::cout << it - corner_mcells[nsimple_blob].begin() << " " << mcell->get_allcell().size() << " " << prev_mcell->get_allcell().size() << " " << mcell->Overlap(*prev_mcell) << std::endl;

	      if (mcell->Overlap1(*prev_mcell)){
	      	std::cout << mcell->center().y << " " << mcell->center().z << " " <<
	      	  prev_mcell->center().y << " " << prev_mcell->center().z << std::endl;
	      	//std::cout << it - corner_mcells[nsimple_blob].begin();
	      	flag = 1;
	      	break;
	      }
	    }
	  }
	  
	  //see next time slice
	  if (flag==0){
	    for (int j=0;j!=next_mergetiling.get_allcell().size();j++){
	      MergeGeomCell *next_mcell = (MergeGeomCell*) next_mergetiling.get_allcell().at(j);
	      if (next_toymatrix.Get_Cell_Charge(next_mcell)>2000){
		if (mcell->Overlap1(*next_mcell)){
		  flag = 1;
		  break;
		}
	      }
	    }
	  }
	  
	  
	  if (flag==1){
	    corner_smcells[nsimple_blob].push_back(mcell);
	    corner_mcells[nsimple_blob].pop_back();
	  }
	  
	}



	nsimple_blob ++;
	// if (nsimple_blob >= 10) {
	//   break;
	// }
	
      }
    }

    
    
    std::cout << "SimpleBlobTiling: "<< nsimple_blob << " " << corner_smcells[0].size() << " " << corner_mcells[0].size() << std::endl;
    for (int j=0;j!=corner_smcells[0].size();j++){
      std::cout << ((MergeGeomCell*)corner_smcells[0].at(j))->get_allcell().size() << " ";
      for (auto it = ((MergeGeomCell*)corner_smcells[0].at(j))->ewires.begin(); it!= ((MergeGeomCell*)corner_smcells[0].at(j))->ewires.end(); it++){
	std::cout << *it << " ";
      }
      std::cout << std::endl;
    }

  }
}

WireCell2dToy::SimpleBlobToyTiling::~SimpleBlobToyTiling(){
  for (int i=0;i!=nsimple_blob;i++){
    for (int j=0;j!=corner_mcells[i].size();j++){
      delete corner_mcells[i].at(j);
    }
    for (int j=0;j!=corner_smcells[i].size();j++){
      delete corner_smcells[i].at(j);
    }
  }
}
