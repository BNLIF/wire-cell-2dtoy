#include "WireCell2dToy/SimpleBlobToyTiling.h"
#include "TRandom.h"
#include "TDecompSVD.h"

// #include <lsp/nnls.h>

//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>

//#include <sstream>


using namespace WireCell;

double WireCell2dToy::SimpleBlobToyTiling::Get_Cell_Charge( const WireCell::GeomCell *cell, int flag )  {
  // flag == 1 charge
  // flag == 2 uncertainty
  int index = scimap_save[cell];
  if (flag==1){
    return Cx_save.at(index);
  }else if (flag==2){
    return dCx_save.at(index);
  }
}


WireCell2dToy::SimpleBlobToyTiling::SimpleBlobToyTiling(WireCell2dToy::ToyTiling& toytiling1, WireCell2dToy::MergeToyTiling& mergetiling1, WireCell2dToy::ToyMatrix& toymatrix1,
							WireCell2dToy::MergeToyTiling& prev_mergetiling, WireCell2dToy::ToyMatrix& prev_toymatrix,
							WireCell2dToy::MergeToyTiling& next_mergetiling, WireCell2dToy::ToyMatrix& next_toymatrix, int recon_t){
  toytiling = &toytiling1;
  mergetiling = &mergetiling1;
  toymatrix = &toymatrix1;
  recon_threshold = recon_t;
  

  nsimple_blob = 0;
  MergeGeomCell *mcorner_cell;
  if (toymatrix->GetSimpleBlobReduction()){
    GeomCellSelection mcells = mergetiling->get_allcell();
    for (int i=0;i!=mcells.size();i++){
      MergeGeomCell *mcell = (MergeGeomCell*)mcells.at(i);
      if (mcell->IsSimpleBlob() && mcell->IsBlob()){
	//save all the simple blob
	sbcells.push_back(mcell);
	if (corner_smcells.size() < nsimple_blob + 1){
	  GeomCellSelection cellss;
	  corner_smcells.push_back(cellss);
	}
	if (corner_mcells.size() < nsimple_blob + 1){
	  GeomCellSelection cellss;
	  corner_mcells.push_back(cellss);
	}
	
	CellIndexMap indexmap = mcell->get_cornercells_index();
	
	for (int j=0;j!=12;j++){
	  int index1 = mcell->index1(j);
	  int index2 = mcell->index2(j);
	  GeomCellSelection cells = mcell->get_cornercells(index1,index2);
	  if (cells.size()!=0){
	    //std::cout << index1 << " " << index2 << " " << cells.size() << std::endl;
	    int flag = -1;
	    //form merged cell and start to fill rank ... 
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
		
		corner_smcells[nsimple_blob].push_back(mcorner_cell);
		
	      }else{
		
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
	    if (prev_toymatrix.Get_Cell_Charge(prev_mcell)>recon_threshold){
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
	    if (next_toymatrix.Get_Cell_Charge(next_mcell)>recon_threshold){
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
	    if (prev_toymatrix.Get_Cell_Charge(prev_mcell)>recon_threshold){
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
	    if (next_toymatrix.Get_Cell_Charge(next_mcell)>recon_threshold){
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
	
	//save stuff into the first, second and third ... 
	Organize(nsimple_blob);
	
	nsimple_blob ++;
	// if (nsimple_blob >= 10) {
	//   break;
	// }
      }
    }

   

    //form the hypothesis and do the test ... 
    ncount = 0;

    //std::cout << "xin1" << std::endl;
    FormHypo();
    //std::cout << "xin1" << std::endl;
    DoTiling();
    //std::cout << "xin1" << std::endl;
    cur_chi2 = CalChi2();
    //std::cout << "xin1" << std::endl;
    SaveResult();
    //std::cout << "xin1" << std::endl;
    while(cur_chi2 > 5 && ncount < 100){
      //std::cout << cur_chi2 << std::endl;
      
      FormHypo();
      //std::cout << "xin1" << std::endl;
      DoTiling();
      //std::cout << "xin1" << std::endl;
      cur_chi2 = CalChi2();
      //std::cout << "xin2" << std::endl;
      //std::cout << cur_chi2 << std::endl;
      if (cur_chi2 < chi2_save)
    	SaveResult();
    }

    std::cout << "DeBlob Chi2: " << cur_chi2 << " " << cell_all.size() << std::endl;
    
    
    
    // std::cout << "SimpleBlobTiling: "<< nsimple_blob << " " << corner_smcells[0].size() << " " << corner_mcells[0].size() << " " << hypo_ccells.at(0).size() << std::endl;
    // for (int j=0;j!=hypo_ccells.at(0).size();j++){
    //   std::cout << hypo_ccells.at(0).at(j).size() << " ";
    //   for (int k = 0;k!= hypo_ccells.at(0).at(j).size();k++){
    // 	std::cout << cell_rank[hypo_ccells.at(0).at(j).at(k)] << " ";
    //   }
    //   std::cout << std::endl;
    // }
   

  }
}

double WireCell2dToy::SimpleBlobToyTiling::CalChi2(){
  //using namespace boost::numeric::ublas;

  double chi2;
  Buildup_index();
  TMatrixD MA(swindex,scindex);
  TMatrixD MAT(scindex,swindex);
 
  TMatrixD Vy(swindex,swindex);
  TMatrixD Vy_inv(swindex,swindex);
 
  TVectorD Wy(swindex);
  TVectorD Wy_pred(swindex);
  
  TMatrixD MC(scindex,scindex);
  TMatrixD MC_inv(scindex,scindex);

  TMatrixD Vx(scindex,scindex);
  TMatrixD Vx_inv(scindex,scindex);


  TVectorD VCx(scindex);
  
  // fill wire charge
  WireChargeMap wcmap = toytiling->wcmap();
  WireChargeMap wcemap = toytiling->wcemap();
  for (int i=0; i!=wire_all.size();i++){
    int index = swimap[wire_all.at(i)];
    float charge = wcmap[wire_all.at(i)];
    float charge_err = wcemap[wire_all.at(i)];
    Wy[index] =charge;
    Vy(index,index) = charge_err*charge_err/1e6;
    //WirePlaneType_t plane = wire_all.at(i)->plane();
    //double charge_noise;
    // if (plane == 0){
    //   charge_noise = 14000*0.05;
    // }else if (plane==1){
    //   charge_noise = 14000*0.03;
    // }else if (plane==2){
    //   charge_noise = 14000*0.02;
    // }
    // Vy(index,index) = (charge_noise*charge_noise + 0.05*0.05 * charge*charge)/1e6;
  }
  
  //fill MA
  for (int i=0;i!=wire_all.size();i++){
    int index = swimap[wire_all.at(i)];
    for (int j=0;j!=wiremap[wire_all.at(i)].size();j++){
      int index1 = scimap[wiremap[wire_all.at(i)].at(j)];
      MA(index,index1) = 1;
    }
  }
  
  MAT.Transpose(MA);
  Vy_inv = Vy;
  Vy_inv.Invert();

  MC = MAT * Vy_inv * MA;


  // //Now try NNLS  need to improve
  // int m = scindex, n = swindex;

  // std::cout << m << " " << n << std::endl;

  // matrix < double > NNLS_A(n,m);
  // vector < double > NNLS_b(n);
    
  // //fill value ...
  // for (int i = 0;i!=n;i++){
  //   NNLS_b(i) = Wy[i];
  //   //std::cout << i << " " << Wy[i] << std::endl;
  //   for (int j=0;j!=m;j++){
  //     NNLS_A(i,j) = MA(i,j);
  //     //std::cout << i << " " << j << " " << MA(i,j) << std::endl;
  //   }
  // }

  // lsp::nnls< matrix< double >, vector< double >  > nnls( NNLS_A, NNLS_b );
  // vector< double > NNLS_x(m);
  // matrix< double > NNLS_cov(m,m);
  // nnls.solve( NNLS_x , NNLS_cov );
  
  // // save results and calculate chi2 
  // for (int i=0;i!=scindex;i++){
  //   (*Cx)[i] = NNLS_x(i);
  //   (*dCx)[i] = sqrt(NNLS_cov(i,i));
  // }


  //SVD solution
  TDecompSVD svd(MC);
  MC_inv = svd.Invert();
  
  VCx = MC_inv * MAT * Vy_inv * Wy;
  Vx_inv = MAT * Vy_inv * MA;
  TDecompSVD svd1(Vx_inv);
  Vx = svd1.Invert();

  Cx.clear();
  dCx.clear();
  for (int i = 0; i!=scindex; i++){
    if (Vx(i,i)>=0){
      dCx.push_back(sqrt( Vx(i,i)) * 1000.);
    }else{
      dCx.push_back(0);
    }
    Cx.push_back(VCx[i]);
    //std::cout << i << " " << (*Cx)[i] << " " << (*dCx)[i] << std::endl;
  }
  //  std::cout << std::endl;


  Wy_pred = MA * VCx;
  TVectorD r1 = Vy_inv * (Wy - Wy_pred);
  TVectorD r2 = Wy - Wy_pred;
  chi2 = r1 * r2 / 1e6;


  return chi2;
}

void WireCell2dToy::SimpleBlobToyTiling::Buildup_index(){
  scindex = 0;
  swindex = 0;
  
  scimap.clear();
  swimap.clear();
  

  
  for (int j=0;j!=cell_all.size();j++){
    //construct merged cell index
    const GeomCell *cell = cell_all[j];
    if (scimap.find(cell) == scimap.end()){
      scimap[cell] = scindex;
      scindex ++;
      
      const GeomWireSelection wires = cellmap[cell];
      for (int k=0;k!=wires.size();k++){
	//construct merged wire index
	const GeomWire *wire = wires[k];
	if (swimap.find(wire) == swimap.end()){
	  swimap[wire] = swindex;
	  swindex ++;
	}
      }
    }
  }
  
}


void WireCell2dToy::SimpleBlobToyTiling::SaveResult(){
  chi2_save = cur_chi2;
  //  std::cout << "abc1 " <<std::endl;
  // if (hypo_save.size() !=0){
  //   for (int i = 0;i!=hypo_save.size();i++){
  //     for (int j=0;j!=hypo_save.at(i).size();j++){
  // 	delete hypo_save.at(i).at(j);
  //     }
  //     hypo_save.at(i).clear();
  //   }
  //   hypo_save.clear();
  // }
  // hypo_save = cur_hypo;
  //std::cout << "abc2 " <<std::endl;

  wire_u_save.clear();
  wire_v_save.clear();
  wire_w_save.clear();
  wire_all_save.clear();
  cell_all_save.clear();

  wire_u_save = wire_u;
  wire_v_save = wire_v;
  wire_w_save = wire_w;
  wire_all_save = wire_all;
  cell_all_save = cell_all;

   Cx_save.clear();
  dCx_save.clear();
  Cx_save = Cx;
  dCx_save = dCx;

  //  std::cout << "abc3 " <<std::endl;

  cellmap_save.clear();
  wiremap_save.clear();

  cellmap_save = cellmap;
  wiremap_save = wiremap;
  
  //  std::cout << "abc4 " <<std::endl;

  scimap_save.clear();
  swimap_save.clear();
  
  scindex_save = scindex;
  swindex_save = swindex;
  swimap_save = swimap;
  scimap_save = scimap;

  //  std::cout << "abc6 " <<std::endl;

 
}


void WireCell2dToy::SimpleBlobToyTiling::DoTiling(){
  
  GeomWireSelection mwires; //save all the merged wires
  GeomCellSelection mcells; // save other cells, and any associated cells
  GeomCellSelection mcells_no;
  //use merge tiling to do stuff

  //save the first pass merge wires
  for (int i=0;i!=sbcells.size();i++){
    MergeGeomCell* mcell = (MergeGeomCell*) sbcells.at(i);
    if (flag_cell.at(i)!=0){
      GeomWireSelection wires = mergetiling->wires(*mcell);
      for (int j=0;j!=wires.size();j++){
	mwires.push_back(wires.at(j));
      }
    }
  }

  //std::cout << "xin1: " << mcells.size() << " " << mwires.size() << std::endl; 

  //save the first pass merged wires;
  for (int i=0;i!=mwires.size();i++){
    MergeGeomWire* mwire = (MergeGeomWire*) mwires.at(i);
    GeomCellSelection cells = mergetiling->cells(*mwire);
    for (int j=0;j!=cells.size();j++){
      if (toymatrix->Get_Cell_Charge(cells.at(j))>recon_threshold){
	auto it = find(sbcells.begin(),sbcells.end(),cells.at(j));
	if (it==sbcells.end()){
	  mcells.push_back(cells.at(j));
	}
      }else{
	mcells_no.push_back(cells.at(j));
      }
    }
  }
  // std::cout << "xin2: " << mcells.size() << " " << mwires.size() << std::endl; 

  
  //save all ... 
  int flag = 1;
  while(flag){
    flag = 0;
    //do wires
    for (int i=0;i!=mcells.size();i++){
      MergeGeomCell* mcell = (MergeGeomCell*)mcells.at(i);
      GeomWireSelection wires = mergetiling->wires(*mcell);
      for (int j=0;j!=wires.size();j++){
	MergeGeomWire* mwire = (MergeGeomWire*)wires.at(j);
	auto it = find(mwires.begin(),mwires.end(),mwire);
	if (it == mwires.end()){
	  flag = 1;
	  mwires.push_back(mwire);
	}
      }
    }
    //do cells
    for (int i=0;i!=mwires.size();i++){
      MergeGeomWire* mwire = (MergeGeomWire*)mwires.at(i);
      GeomCellSelection cells = mergetiling->cells(*mwire);
      for (int j=0;j!=cells.size();j++){
	MergeGeomCell* mcell = (MergeGeomCell*)cells.at(j);
	auto it = find(mcells.begin(),mcells.end(),mcell);
	auto it1 = find(sbcells.begin(),sbcells.end(),mcell);
	auto it2 = find(mcells_no.begin(),mcells_no.end(),mcell);
	if (it==mcells.end() && it1 == sbcells.end() && it2 == mcells_no.end()){
	  flag = 1;
	  mcells.push_back(mcell);
	}
      }
    }
  }
  //std::cout << "xin3: " << mcells.size() << " " << mwires.size() << std::endl; 



  //use hypothesis to save stuff
  //save all the cells;
  cell_all.clear();
  for (int i = 0; i!=sbcells.size();i++){
    WireCell2dToy::HypoSelection hypos = cur_hypo.at(i);
    MergeGeomCell* mcell = (MergeGeomCell*) sbcells.at(i);
    for (int j=0;j!=mcell->get_allcell().size();j++){
      const GeomCell* cell = mcell->get_allcell().at(j);
      for (int k=0;k!=hypos.size();k++){
	WireCell2dToy::ToyHypothesis *hypo = (WireCell2dToy::ToyHypothesis*) hypos.at(k);

	//std::cout << i << " " << j << " " << k << " " << hypo->IsInside(*cell) << std::endl;

	if (hypo->IsInside(*cell)){
	  cell_all.push_back(cell);
	  break;
	}
      }
    }
  }
  
  for (int i=0;i!=mcells.size();i++){
    MergeGeomCell* mcell = (MergeGeomCell*) mcells.at(i);
    for (int j=0;j!=mcell->get_allcell().size();j++){
      const GeomCell* cell = mcell->get_allcell().at(j);
      cell_all.push_back(cell);
    }
  }
  
  wire_all.clear();
  wire_u.clear();
  wire_v.clear();
  wire_w.clear();

  for (int i=0;i!=mwires.size();i++){
    MergeGeomWire* mwire = (MergeGeomWire*)mwires.at(i);
    for (int j=0;j!=mwire->get_allwire().size();j++){
      const GeomWire* wire = mwire->get_allwire().at(j);
      wire_all.push_back(wire);
      if (wire->plane() == ((WirePlaneType_t)0)){
  	wire_u.push_back(wire);
      }else if (wire->plane() == ((WirePlaneType_t)1)){
  	wire_v.push_back(wire);
      }else{
  	wire_w.push_back(wire);
      }
    }
  }
    
  //use toy tiling to do association
  cellmap.clear();
  wiremap.clear();

  GeomCellMap cmap = toytiling->cmap();
  GeomWireMap wmap = toytiling->wmap();
  
  for (int i=0;i!=cell_all.size();i++){
    const GeomCell* cell = cell_all.at(i);
    cellmap[cell] = cmap[cell];
  }
  for (int i=0;i!=wire_all.size();i++){
    const GeomWire* wire = wire_all.at(i);
    GeomCellSelection cells = wmap[wire];
    GeomCellSelection ncells;
    for (int j=0;j!=cells.size();j++){
      const GeomCell* cell = cells.at(j);
      auto it = find(cell_all.begin(),cell_all.end(),cell);
      if (it != cell_all.end()){
  	ncells.push_back(cell);
      }
    }
    wiremap[wire] = ncells;
  }
  
  
  
}


void WireCell2dToy::SimpleBlobToyTiling::FormHypo(){
  // std::cout << ncount << std::endl;
   ClearHypo();

   //std::cout << ncount << std::endl;
  
  if (ncount == 0){
    for (int i=0;i!=nsimple_blob;i++){
      WireCell2dToy::HypoSelection hypos;
      if (flag_cell.at(i) == 1){
	// only one hypothsis to connect highest and second highest
	MergeGeomCell *mcell1 = (MergeGeomCell*) first_cell.at(i).at(0);
	MergeGeomCell *mcell2 = (MergeGeomCell*) second_cell.at(i).at(0);
	WireCell2dToy::ToyHypothesis *hypo = new WireCell2dToy::ToyHypothesis(*mcell1,*mcell2);
	hypos.push_back(hypo);
      }else if (flag_cell.at(i) == 2){
	// two hypothesis, highest and second highest, and then highest to the third
	MergeGeomCell *mcell1 = (MergeGeomCell*) first_cell.at(i).at(0);
	MergeGeomCell *mcell2 = (MergeGeomCell*) second_cell.at(i).at(0);
	MergeGeomCell *mcell3 = (MergeGeomCell*) other_cell.at(i).at(0);
	WireCell2dToy::ToyHypothesis *hypo = new WireCell2dToy::ToyHypothesis(*mcell1,*mcell2);
	hypos.push_back(hypo);
	
	//do the smallest distance
	double dis1 = pow(mcell1->center().y-mcell3->center().y,2) + pow(mcell1->center().z-mcell3->center().z,2);
	double dis2 = pow(mcell2->center().y-mcell3->center().y,2) + pow(mcell2->center().z-mcell3->center().z,2);
	
	if (dis1 < dis2){
	  hypo = new WireCell2dToy::ToyHypothesis(*mcell1,*mcell3);
	}else{
	  hypo = new WireCell2dToy::ToyHypothesis(*mcell2,*mcell3);
	}
	hypos.push_back(hypo);
      }else if (flag_cell.at(i) == 3){
	// two hypothesis, highest to third, and second highest to third
	MergeGeomCell *mcell1 = (MergeGeomCell*) first_cell.at(i).at(0);
	MergeGeomCell *mcell2 = (MergeGeomCell*) second_cell.at(i).at(0);
	MergeGeomCell *mcell3 = (MergeGeomCell*) other_cell.at(i).at(0);

	MergeGeomCell *mcell4;
	std::set<int> edge_wires;
      	for (auto it = mcell1->ewires.begin(); it!=mcell1->ewires.end(); it++){
	  edge_wires.insert(*it);
	}
	for (auto it = mcell2->ewires.begin(); it!=mcell2->ewires.end(); it++){
	  edge_wires.insert(*it);
	}
	for (auto it = mcell3->ewires.begin(); it!=mcell3->ewires.end(); it++){
	  edge_wires.insert(*it);
	}
	std::vector<int> not_found;
	for (int j=0;j!=6;j++){
	  if (edge_wires.find(j)==edge_wires.end()){
	    not_found.push_back(j);
	  }
	}
	for (int j=1;j<other_cell.at(i).size();j++){
	  mcell4 = (MergeGeomCell*) other_cell.at(i).at(j);
	  int flag = 1;
	  for (int k=0;k!=not_found.size();k++){
	    auto it = mcell4->ewires.find(not_found.at(k));
	    if (it == mcell4->ewires.end()){
	      flag = 0;
	      break;
	    }
	  }
	  if (flag==1) break;
	}
	WireCell2dToy::ToyHypothesis *hypo = new WireCell2dToy::ToyHypothesis(*mcell1,*mcell3);
	hypos.push_back(hypo);
	hypo = new WireCell2dToy::ToyHypothesis(*mcell2,*mcell4);
	hypos.push_back(hypo);
      }
      cur_hypo.push_back(hypos);
    }
  }else{
    //start to randomize ... 
    for (int i=0;i!=nsimple_blob;i++){
      WireCell2dToy::HypoSelection hypos;
      if (flag_cell.at(i) == 1){
	// single track (change to two tracks)
	int n = gRandom->Uniform(0,first_cell.at(i).size());
	MergeGeomCell *mcell1 = (MergeGeomCell*) first_cell.at(i).at(n);
	n = gRandom->Uniform(0,second_cell.at(i).size());
	MergeGeomCell *mcell2 = (MergeGeomCell*) second_cell.at(i).at(n);

	if (other_cell.at(i).size()!=0){
	  n = gRandom->Uniform(0,other_cell.at(i).size());
	  MergeGeomCell *mcell3 = (MergeGeomCell*) other_cell.at(i).at(n);	
	  WireCell2dToy::ToyHypothesis *hypo = new WireCell2dToy::ToyHypothesis(*mcell1,*mcell2);
	  hypos.push_back(hypo);
	  n = gRandom->Uniform(0.1,1.9);
	  if (n==0){
	    hypo = new WireCell2dToy::ToyHypothesis(*mcell1,*mcell3);
	    hypos.push_back(hypo);
	  }else{
	    hypo = new WireCell2dToy::ToyHypothesis(*mcell2,*mcell3);
	    hypos.push_back(hypo);
	  }
	}else{
	  WireCell2dToy::ToyHypothesis *hypo = new WireCell2dToy::ToyHypothesis(*mcell1,*mcell2);
	  hypos.push_back(hypo);
	}
      }else if (flag_cell.at(i) == 2){
	// two hypothesis, highest and second highest, and then highest to the third
	int n = gRandom->Uniform(0,first_cell.at(i).size());
	MergeGeomCell *mcell1 = (MergeGeomCell*) first_cell.at(i).at(n);
	n = gRandom->Uniform(0,second_cell.at(i).size());
	MergeGeomCell *mcell2 = (MergeGeomCell*) second_cell.at(i).at(n);
	
	if (other_cell.at(i).size()!=0){
	  n = gRandom->Uniform(0,other_cell.at(i).size());
	  MergeGeomCell *mcell3 = (MergeGeomCell*) other_cell.at(i).at(n);
	  int nnn = 0;
	  while(mcell3 == mcell2 && nnn<20){
	    n = gRandom->Uniform(0,second_cell.at(i).size());
	    MergeGeomCell *mcell2 = (MergeGeomCell*) second_cell.at(i).at(n);
	    n = gRandom->Uniform(0,other_cell.at(i).size());
	    MergeGeomCell *mcell3 = (MergeGeomCell*) other_cell.at(i).at(n);
	    nnn ++;
	  }
	  
	  WireCell2dToy::ToyHypothesis *hypo = new WireCell2dToy::ToyHypothesis(*mcell1,*mcell2);
	  hypos.push_back(hypo);
	  n = gRandom->Uniform(0.1,1.9);
	  if (n==0){
	    hypo = new WireCell2dToy::ToyHypothesis(*mcell1,*mcell3);
	    hypos.push_back(hypo);
	  }else{
	    hypo = new WireCell2dToy::ToyHypothesis(*mcell2,*mcell3);
	    hypos.push_back(hypo);
	  }
	}else{
	  WireCell2dToy::ToyHypothesis *hypo = new WireCell2dToy::ToyHypothesis(*mcell1,*mcell2);
	  hypos.push_back(hypo);
	}
      }else if (flag_cell.at(i) == 3){
	// two hypothesis, highest to third, and second highest to third
	int n = gRandom->Uniform(0,first_cell.at(i).size());
	MergeGeomCell *mcell1 = (MergeGeomCell*) first_cell.at(i).at(n);
	n = gRandom->Uniform(0,second_cell.at(i).size());
	MergeGeomCell *mcell2 = (MergeGeomCell*) second_cell.at(i).at(n);

	if (other_cell.at(i).size()!=0){
	  n = gRandom->Uniform(0,other_cell.at(i).size());
	  MergeGeomCell *mcell3 = (MergeGeomCell*) other_cell.at(i).at(n);
	  
	  int nnn = 0;
	  while(mcell3 == mcell2 && nnn<20){
	    n = gRandom->Uniform(0,second_cell.at(i).size());
	    MergeGeomCell *mcell2 = (MergeGeomCell*) second_cell.at(i).at(n);
	    n = gRandom->Uniform(0,other_cell.at(i).size());
	    MergeGeomCell *mcell3 = (MergeGeomCell*) other_cell.at(i).at(n);
	    nnn ++;
	  }
	
	  MergeGeomCell *mcell4;
	  std::set<int> edge_wires;
	  for (auto it = mcell1->ewires.begin(); it!=mcell1->ewires.end(); it++){
	    edge_wires.insert(*it);
	  }
	  for (auto it = mcell2->ewires.begin(); it!=mcell2->ewires.end(); it++){
	    edge_wires.insert(*it);
	  }
	  for (auto it = mcell3->ewires.begin(); it!=mcell3->ewires.end(); it++){
	    edge_wires.insert(*it);
	  }
	  std::vector<int> not_found;
	  for (int j=0;j!=6;j++){
	    if (edge_wires.find(j)==edge_wires.end()){
	      not_found.push_back(j);
	    }
	  }
	  
	  //need to improve later .. no random here
	  for (int j=1;j<other_cell.at(i).size();j++){
	    mcell4 = (MergeGeomCell*) other_cell.at(i).at(j);
	    int flag = 1;
	    for (int k=0;k!=not_found.size();k++){
	      auto it = mcell4->ewires.find(not_found.at(k));
	      if (it == mcell4->ewires.end()){
		flag = 0;
		break;
	      }
	    }
	    if (flag==1) break;
	  }
	  
	  n = gRandom->Uniform(0.1,2.9);
	  if (n==0){
	    WireCell2dToy::ToyHypothesis *hypo = new WireCell2dToy::ToyHypothesis(*mcell1,*mcell2);
	    hypos.push_back(hypo);
	    hypo = new WireCell2dToy::ToyHypothesis(*mcell3,*mcell4);
	    hypos.push_back(hypo);
	  }else if (n==1){
	    WireCell2dToy::ToyHypothesis *hypo = new WireCell2dToy::ToyHypothesis(*mcell1,*mcell3);
	    hypos.push_back(hypo);
	    hypo = new WireCell2dToy::ToyHypothesis(*mcell2,*mcell4);
	    hypos.push_back(hypo);
	  }else if (n==2){
	    WireCell2dToy::ToyHypothesis *hypo = new WireCell2dToy::ToyHypothesis(*mcell1,*mcell4);
	    hypos.push_back(hypo);
	    hypo = new WireCell2dToy::ToyHypothesis(*mcell2,*mcell3);
	    hypos.push_back(hypo);
	  }
	}else{
	  WireCell2dToy::ToyHypothesis *hypo = new WireCell2dToy::ToyHypothesis(*mcell1,*mcell2);
	  hypos.push_back(hypo);
	}
	
      }
      cur_hypo.push_back(hypos);
    }
  }
  ncount ++;
}

void WireCell2dToy::SimpleBlobToyTiling::ClearHypo(){
  if (cur_hypo.size() !=0){
    for (int i = 0;i!=cur_hypo.size();i++){
      for (int j=0;j!=cur_hypo.at(i).size();j++){
	delete cur_hypo.at(i).at(j);
      }
      cur_hypo.at(i).clear();
    }
    cur_hypo.clear();
  }
}


void WireCell2dToy::SimpleBlobToyTiling::Organize(int nsimple_blob){
  //fill the first rank, second rank and third rank and flag ... 
  // loop through the smcells and find the highest rank, fill in first rank, 
  int max = 0 ;
  int max_bin;

  flag_cell.push_back(0);
  
  //find the maximum ranked cell
  for (int i=0;i!=corner_smcells.at(nsimple_blob).size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*) corner_smcells.at(nsimple_blob).at(i);
    if (cell_rank[mcell] > max){
      max = cell_rank[mcell];
      max_bin = i;
    }
  }

  if (max !=0){
    //if found
    GeomCellSelection cells;
    cells.push_back(corner_smcells.at(nsimple_blob).at(max_bin));


    // loop through everything and fill in everything connected
    int flag = 1;
    while(flag){
      flag = 0;
      for (int i=0;i!=corner_smcells.at(nsimple_blob).size();i++){
	MergeGeomCell *mcell = (MergeGeomCell*) corner_smcells.at(nsimple_blob).at(i);
	auto it = find(cells.begin(),cells.end(),mcell);
	if (it==cells.end()){
	  for (int j=0;j!=cells.size();j++){
	    MergeGeomCell *mcell1 = (MergeGeomCell*) cells.at(j);
	    if (mcell1->Overlap1(*mcell) ){
	      flag = 1;
	      cells.push_back(mcell);
	      break;
	    }
	  }
	}
      }
    }
    first_cell.push_back(cells);
    
    // loop through the smcells and find the second highest rank, fill in second rank
    max = 0 ;
    for (int i=0;i!=corner_smcells.at(nsimple_blob).size();i++){
      MergeGeomCell *mcell = (MergeGeomCell*) corner_smcells.at(nsimple_blob).at(i);
      auto it = find(first_cell.at(nsimple_blob).begin(),first_cell.at(nsimple_blob).end(),mcell);
      if (it == first_cell.at(nsimple_blob).end() && cell_rank[mcell] > max){
	max = cell_rank[mcell];
	max_bin = i;
      }
    }

    if (max!=0){ //if found
      GeomCellSelection cells1;
      cells1.push_back(corner_smcells.at(nsimple_blob).at(max_bin));
      
  

      // loop through everything and fill in everything connected
      flag = 1;
      while(flag){
	flag = 0;
	for (int i=0;i!=corner_smcells.at(nsimple_blob).size();i++){
	  MergeGeomCell *mcell = (MergeGeomCell*) corner_smcells.at(nsimple_blob).at(i);
	
	  auto it = find(cells1.begin(),cells1.end(),mcell);
	  auto it1 = find(first_cell.at(nsimple_blob).begin(),first_cell.at(nsimple_blob).end(),mcell);
	  if (it==cells1.end() && it1 == first_cell.at(nsimple_blob).end()){
	    for (int j=0;j!=cells1.size();j++){
	      MergeGeomCell *mcell1 = (MergeGeomCell*) cells1.at(j);
	      if (mcell1->Overlap1(*mcell)){
		flag = 1;
		cells1.push_back(mcell);
		break;
	      }
	    }
	  }
	}
      }
      second_cell.push_back(cells1);
      // do the third group and fill the flag ... 
      
      std::set<int> edge_wires;
      for (int i=0;i!=first_cell.at(nsimple_blob).size();i++){
	MergeGeomCell* mcell = (MergeGeomCell*) first_cell.at(nsimple_blob).at(i);
	for (auto it = mcell->ewires.begin(); it!=mcell->ewires.end(); it++){
	  edge_wires.insert(*it);
	}
      }
      for (int i=0;i!=second_cell.at(nsimple_blob).size();i++){
	MergeGeomCell* mcell = (MergeGeomCell*) second_cell.at(nsimple_blob).at(i);
	for (auto it = mcell->ewires.begin(); it!=mcell->ewires.end(); it++){
	  edge_wires.insert(*it);
	}
      }

      // for (auto it = edge_wires.begin();it!=edge_wires.end();it++){
      // 	std::cout << *it << " ";
      // }
      // std::cout << std::endl;
      if (edge_wires.size()==6){
	//put everything into the third group
	flag_cell.at(nsimple_blob) = 1;
	GeomCellSelection cells;
	for (int i=0;i!=corner_smcells.at(nsimple_blob).size();i++){
	  MergeGeomCell *mcell = (MergeGeomCell*) corner_smcells.at(nsimple_blob).at(i);
	  auto it = find(second_cell.at(nsimple_blob).begin(),second_cell.at(nsimple_blob).end(),mcell);
	  auto it1 = find(first_cell.at(nsimple_blob).begin(),first_cell.at(nsimple_blob).end(),mcell);
	  if (it==second_cell.at(nsimple_blob).end() && it1 == first_cell.at(nsimple_blob).end()){
	    cells.push_back(mcell);
	  }
	}
	for (int i=0;i!=corner_mcells.at(nsimple_blob).size();i++){
	  MergeGeomCell *mcell = (MergeGeomCell*) corner_mcells.at(nsimple_blob).at(i);
	  cells.push_back(mcell);
	}
	other_cell.push_back(cells);
      }else{
	//two cases
	GeomCellSelection cells;
	std::vector<int> not_found;
	for (int i=0;i!=6;i++){
	  if (edge_wires.find(i)==edge_wires.end()){
	    not_found.push_back(i);
	  }
	}
	int flag = 0;
	for (int i=0;i!=corner_smcells.at(nsimple_blob).size();i++){
	  MergeGeomCell *mcell = (MergeGeomCell*) corner_smcells.at(nsimple_blob).at(i);
	  int flag1 = 1;
	  for (int j=0;j!=not_found.size();j++){
	    if (mcell->ewires.find(not_found.at(j))==mcell->ewires.end()){
	      flag1 = 0;
	    }
	  }
	  if (flag1==1){
	    cells.push_back(mcell);
	    flag = 1;
	  }
	}
	for (int i=0;i!=corner_mcells.at(nsimple_blob).size();i++){
	  MergeGeomCell *mcell = (MergeGeomCell*) corner_mcells.at(nsimple_blob).at(i);
	  int flag1 = 1;
	  for (int j=0;j!=not_found.size();j++){
	    if (mcell->ewires.find(not_found.at(j))==mcell->ewires.end()){
	      flag1 = 0;
	    }
	  }
	  if (flag1==1){
	    cells.push_back(mcell);
	    flag = 1;
	  }
	}
		
	if (flag==0){
	  //find the ones containing any of them ... 
	  for (int i=0;i!=corner_smcells.at(nsimple_blob).size();i++){
	    MergeGeomCell *mcell = (MergeGeomCell*) corner_smcells.at(nsimple_blob).at(i);
	    int flag1 = 0;
	    for (int j=0;j!=not_found.size();j++){
	      if (mcell->ewires.find(not_found.at(j))!=mcell->ewires.end()){
		flag1 = 1;
	      }
	    }
	    if (flag1==1){
	      cells.push_back(mcell);
	    }
	  }
	  for (int i=0;i!=corner_mcells.at(nsimple_blob).size();i++){
	    MergeGeomCell *mcell = (MergeGeomCell*) corner_mcells.at(nsimple_blob).at(i);
	    int flag1 = 0;
	    for (int j=0;j!=not_found.size();j++){
	      if (mcell->ewires.find(not_found.at(j))!=mcell->ewires.end()){
		flag1 = 1;
	      }
	    }
	    if (flag1==1){
	      cells.push_back(mcell);
	    }
	  }
	  flag_cell.at(nsimple_blob) = 3;
	}else{
	  flag_cell.at(nsimple_blob) = 2;
	}
	other_cell.push_back(cells);

	// fill in all the stuff with rank > 8
	for (int i = 0; i!=corner_smcells.at(nsimple_blob).size();i++){
	  MergeGeomCell *mcell = (MergeGeomCell*) corner_smcells.at(nsimple_blob).at(i);
	  auto it = find(cells1.begin(),cells1.end(),mcell);
      	  auto it1 = find(first_cell.at(nsimple_blob).begin(),first_cell.at(nsimple_blob).end(),mcell);
      	  if (it==cells1.end() && it1 == first_cell.at(nsimple_blob).end()){
      	    //std::cout << cell_rank[mcell] << std::endl;
      	    if (cell_rank[mcell] >=3){
      	      second_cell.at(nsimple_blob).push_back(mcell);
      	    }
      	  }
	}
	
      }
    }else{

      // only find the first ... 
      GeomCellSelection cells;
      std::set<int> edge_wires;
      for (int i=0;i!=first_cell.at(nsimple_blob).size();i++){
	MergeGeomCell* mcell = (MergeGeomCell*) first_cell.at(nsimple_blob).at(i);
	for (auto it = mcell->ewires.begin(); it!=mcell->ewires.end(); it++){
	  edge_wires.insert(*it);
	}
      }
      std::vector<int> found;
	for (int i=0;i!=6;i++){
	  if (edge_wires.find(i)!=edge_wires.end()){
	    found.push_back(i);
	  }
	}
      

      //second first try to find the ones does not contain the same as the existing one
      for (int i=0;i!=corner_smcells.at(nsimple_blob).size();i++){
	MergeGeomCell *mcell = (MergeGeomCell*) corner_smcells.at(nsimple_blob).at(i);
	int flag1 = 1;
	for (int j=0;j!=found.size();j++){
	  if (mcell->ewires.find(found.at(j))!=mcell->ewires.end()){
	    flag1 = 0;
	  }
	}
	if (flag1==1){
	  cells.push_back(mcell);
	}
      }
      for (int i=0;i!=corner_mcells.at(nsimple_blob).size();i++){
	MergeGeomCell *mcell = (MergeGeomCell*) corner_mcells.at(nsimple_blob).at(i);
	int flag1 = 1;
	for (int j=0;j!=found.size();j++){
	  if (mcell->ewires.find(found.at(j))!=mcell->ewires.end()){
	    flag1 = 0;
	  }
	}
	if (flag1==1){
	  cells.push_back(mcell);
	}
      }
      second_cell.push_back(cells);
      //put the rest to other 
      cells.clear();
      for (int i=0;i!=corner_smcells.at(nsimple_blob).size();i++){
	MergeGeomCell *mcell = (MergeGeomCell*) corner_smcells.at(nsimple_blob).at(i);
	auto it = find(second_cell.at(nsimple_blob).begin(),second_cell.at(nsimple_blob).end(),mcell);
	auto it1 = find(first_cell.at(nsimple_blob).begin(),first_cell.at(nsimple_blob).end(),mcell);
	if (it==second_cell.at(nsimple_blob).end() && it1 == first_cell.at(nsimple_blob).end()){
	  cells.push_back(mcell);
	}
      }
      for (int i=0;i!=corner_mcells.at(nsimple_blob).size();i++){
	MergeGeomCell *mcell = (MergeGeomCell*) corner_mcells.at(nsimple_blob).at(i);
	auto it = find(second_cell.at(nsimple_blob).begin(),second_cell.at(nsimple_blob).end(),mcell);
	auto it1 = find(first_cell.at(nsimple_blob).begin(),first_cell.at(nsimple_blob).end(),mcell);
	if (it==second_cell.at(nsimple_blob).end() && it1 == first_cell.at(nsimple_blob).end()){
	  cells.push_back(mcell);
	}
      }
      other_cell.push_back(cells);
      flag_cell.at(nsimple_blob) = 1;
    }
  }else{
    GeomCellSelection cells;
    first_cell.push_back(cells);
    second_cell.push_back(cells);
    other_cell.push_back(cells);
  }
  

  std::cout << first_cell.at(0).size() << " " << second_cell.at(0).size() << " " << other_cell.at(0).size() << " " << flag_cell.at(0) << std::endl;

  // //now put everything into hypo_ccells
  // if (hypo_ccells.size() < nsimple_blob + 1){
  //   GeomCellSelectionV qcv; // vector of vector
  
  //   // first deal with smcells
  //   // do a loop, and insert the first element, 
  //   // do a loop, and merge anything that can be merged, erase, if none
  //   // insert the second element until nothing in the smcells;
  //   int flag_new = 1; //insert?
  //   int flag_merge = 0; //merge?
  //   while(corner_smcells[nsimple_blob].size()){
  //     if (flag_merge == 1 && corner_smcells[nsimple_blob].size()!=0){
  //       int temp_flag = 0;
  //       //do the merge part
  //       for (int j=0;j!=corner_smcells[nsimple_blob].size();j++){
  //         MergeGeomCell* mcell1 = (MergeGeomCell*) corner_smcells[nsimple_blob].at(j);
  // 	for (int k = 0;k!=qcv.size();k++){
  // 	  for (int kk=0;kk!=qcv.at(k).size();kk++){
  // 	    MergeGeomCell* mcell2 = (MergeGeomCell*) qcv.at(k).at(kk);
  // 	    if (mcell1->Overlap1(*mcell2)){
  // 	      qcv.at(k).push_back(mcell1);
  // 	      corner_smcells[nsimple_blob].erase(corner_smcells[nsimple_blob].begin() + j);
  // 	      temp_flag = 1;
  // 	      break;
  // 	    }
  // 	  }
  // 	  if (temp_flag==1) break;
  // 	}
  // 	if (temp_flag==1) break;
  //       }
  //       if (temp_flag == 0){
  // 	flag_merge = 0;
  // 	flag_new = 1;
  //       }
  //     }
  
  
  //     if (flag_new == 1 && corner_smcells[nsimple_blob].size()!=0){
  //       //insert an element and delete one ... 
  //       GeomCellSelection qc;
  //       qc.push_back(corner_smcells[nsimple_blob].at(corner_smcells[nsimple_blob].size()-1));
  //       qcv.push_back(qc);
  //       corner_smcells[nsimple_blob].pop_back();
  //       flag_new = 0;
  //       flag_merge = 1;
  //     }
  //   }
  
  //   // then deal with mcells
  //   if (corner_mcells[nsimple_blob].size()!=0)
  //     qcv.push_back(corner_mcells[nsimple_blob]);
  
  //   hypo_ccells.push_back(qcv); //save things in ... 
  // }
  
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

  // for (int i = 0 ;i!=hypo_ccells.size(); i++){
  //   for (int j = 0; j!=hypo_ccells.at(i).size();j++){
  //     for (int k=0;k!=hypo_ccells.at(i).at(j).size();k++){
  // 	delete hypo_ccells.at(i).at(j).at(k);
  //     }
  //   }
  // }
}

ClassImp(WireCell2dToy::SimpleBlobToyTiling);
