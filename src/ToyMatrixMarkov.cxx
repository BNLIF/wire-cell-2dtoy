#include "WireCell2dToy/ToyMatrixMarkov.h"

using namespace WireCell;
#include "TMath.h"
#include "TRandom.h"

void WireCell2dToy::ToyMatrixMarkov::find_subset(WireCell2dToy::ToyMatrixKalman &toymatrix,WireCell2dToy::ToyMatrix &toymatrix1, std::vector<int>& vec){
  //std::cout << toymatrix.Get_already_removed().size() << " " << toymatrix.Get_no_need_remove().size() << " " << vec.size() << std::endl;

  if (toymatrix.Get_numz()!=0){
    for (int i=0;i!=toymatrix.Get_mcindex();i++){
      auto it1 = find(toymatrix.Get_already_removed().begin(),toymatrix.Get_already_removed().end(),i);
      auto it2 = find(toymatrix.Get_no_need_remove().begin(),toymatrix.Get_no_need_remove().end(),i);
      if (it1 == toymatrix.Get_already_removed().end() && it2 == toymatrix.Get_no_need_remove().end()){
	std::vector<int> already_removed = toymatrix.Get_already_removed();
	already_removed.push_back(i);
	WireCell2dToy::ToyMatrixKalman kalman(already_removed,toymatrix.Get_no_need_remove(),toymatrix1,0,0);

	//	std::cout << "abc: " << i << " " <<kalman.Get_numz() << " " << toymatrix.Get_numz() << std::endl;
	if (kalman.Get_numz()==toymatrix.Get_numz()){
	}else{
	  find_subset(kalman,toymatrix1,vec);
	  //move to "no need to remove"??
	  toymatrix.Get_no_need_remove().push_back(i);
	  if (vec.size()!=0) break;
	}
      }
    }

    if (vec.size() == 0){
      for (int i=0;i!=toymatrix.Get_mcindex();i++){
	auto it1 = find(toymatrix.Get_already_removed().begin(),toymatrix.Get_already_removed().end(),i);
	auto it2 = find(toymatrix.Get_no_need_remove().begin(),toymatrix.Get_no_need_remove().end(),i);
	if (it1 == toymatrix.Get_already_removed().end() && it2 == toymatrix.Get_no_need_remove().end()){
	  for (int j=0;j!=toymatrix.Get_mcindex();j++){
	    if (i==j) continue;
	    auto it3 = find(toymatrix.Get_already_removed().begin(),toymatrix.Get_already_removed().end(),j);
	    auto it4 = find(toymatrix.Get_no_need_remove().begin(),toymatrix.Get_no_need_remove().end(),j);
	    if (it3 == toymatrix.Get_already_removed().end() && it4 == toymatrix.Get_no_need_remove().end()){
	      std::vector<int> already_removed = toymatrix.Get_already_removed();
	      already_removed.push_back(i);
	      already_removed.push_back(j);
	      WireCell2dToy::ToyMatrixKalman kalman(already_removed,toymatrix.Get_no_need_remove(),toymatrix1,0,0);
	      
	      // std::cout << "abc: " << i << " " << j << " " << kalman.Get_numz() << " " << toymatrix.Get_numz() << std::endl;

	      if (kalman.Get_numz()==toymatrix.Get_numz()){
	      }else{
		find_subset(kalman,toymatrix1,vec);
		//move to "no need to remove"??
		// toymatrix.Get_no_need_remove().push_back(i);
		// toymatrix.Get_no_need_remove().push_back(j);
		if (vec.size()!=0) break;
	      }
	    }
	  }
	  if (vec.size()!=0) break;
	}
      }
    }
    

  }else{
    for (int i=0;i!=toymatrix.Get_mcindex();i++){
      auto it1 = find(toymatrix.Get_already_removed().begin(),toymatrix.Get_already_removed().end(),i);
      if (it1 == toymatrix.Get_already_removed().end()){
	vec.push_back(i);
      }
    }
  }
}

WireCell2dToy::ToyMatrixMarkov::ToyMatrixMarkov(WireCell2dToy::ToyMatrix &toycur,WireCell2dToy::MergeToyTiling &mergecur, WireCell::GeomCellSelection *allmcell1, WireCell::GeomCellSelection &cells, int recon_t1, int recon_t2)
  : penalty_ncpt(0)
{
  ncount = 0;
  first_flag = 0;
  toymatrix = &toycur;
  allmcell = allmcell1;
  mcindex = toymatrix->Get_mcindex();
  recon_threshold1 = recon_t1;
  recon_threshold2 = recon_t2;
  use_time_threshold = 0.75;

  GeomCellSelection allmcell_c = mergecur.get_allcell();

  //form two vectors 
  std::vector<int> already_removed; //dummy
  
  for (int i=0;i!=allmcell_c.size();i++){
    MergeGeomCell *mcell_c = (MergeGeomCell*)allmcell_c[i];
    int index_c = toycur.Get_mcindex(mcell_c);
    
    cell_penal[index_c] = 0;

    auto it = find(cells.begin(),cells.end(),mcell_c);
    if (it !=cells.end()){
      use_time.push_back(index_c);
    }
    
  }

  for (int i=0;i!=allmcell_c.size();i++){
    auto it = find(use_time.begin(),use_time.end(),i);
    if (it==use_time.end())
      already_removed.push_back(i);
  }
  use_time.clear();
  
  //initialize
  //toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(*toymatrix);  // hold the current results 

  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, use_time, *toymatrix,1,0);
  std::cout << "With Time: " << toymatrixkalman->Get_numz() << " " << allmcell_c.size() << " " <<  already_removed.size() << std::endl;
  // Find a sub-set that is not degenerated
  // put things into use_time
  find_subset(*toymatrixkalman,toycur,use_time);
  already_removed.clear();
  
  std::cout << use_time.size() << std::endl;
  // recalculate
  delete toymatrixkalman;
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, use_time, toycur,1,0);
  toymatrix->Set_Solve_Flag(0);
  toymatrix->Set_chi2(-1);

  if (use_time.size() !=0){
    std::cout << "Start making guesses " << std::endl;
    if (allmcell_c.size()<50){
      while (ncount < 1e4 
	     && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	     && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	     && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	     && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800)
	     && (cur_chi2 > 1*(cur_dof+0.1) || ncount < 10)
	     ){
	if (ncount == 1 || (ncount % 1000 ==0&&ncount >0) )
	  std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
	
	make_guess();
	next_chi2 = toymatrix->Get_Chi2();
	next_dof = toymatrix->Get_ndf();
      }
    }else if (allmcell_c.size()<75){
      while (ncount < 5000 
	     && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	     && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	     && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	     && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800)
	     && (cur_chi2 > 1*(cur_dof+0.1) || ncount < 10)
	     ){
	if (ncount == 1 || (ncount % 1000 ==0&&ncount >0) )
	  std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
	
	make_guess();
	next_chi2 = toymatrix->Get_Chi2();
	next_dof = toymatrix->Get_ndf();
      }
    }else if (allmcell_c.size()<110){
      while (ncount < 2000
	     && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	     && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	     && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	     && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800)
	     && (cur_chi2 > 1*(cur_dof+0.1) || ncount < 10)
	     ){
	if (ncount == 1 || (ncount % 1000 ==0&&ncount >0) )
	  std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
	
	make_guess();
	next_chi2 = toymatrix->Get_Chi2();
	next_dof = toymatrix->Get_ndf();
      }
      //    }else if (allmcell_c.size()<150){
    }else{
      for (int i=0;i!=200;i++){
	make_guess();
	next_chi2 = toymatrix->Get_Chi2();
	next_dof = toymatrix->Get_ndf();
      }
    }
  
  }else{
    if (allmcell_c.size()<30){
      while (ncount < 1e4 
	     && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	     && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	     && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	     && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800)
	     && (cur_chi2 > 1*(cur_dof+0.1) || ncount < 10)
	     ){
	if (ncount == 1 || (ncount % 1000 ==0&&ncount >0) )
	  std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
	
	make_guess();
	next_chi2 = toymatrix->Get_Chi2();
	next_dof = toymatrix->Get_ndf();
      }
    }else if (allmcell_c.size()<60){
      while (ncount < 5000 
	     && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	     && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	     && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	     && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800) 
	     && (cur_chi2 > 1*(cur_dof+0.1) || ncount < 10) 
	     ){
	if (ncount == 1 || (ncount % 1000 ==0&&ncount >0) )
	  std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
	
	make_guess();
	next_chi2 = toymatrix->Get_Chi2();
	next_dof = toymatrix->Get_ndf();
      }
    }
  }
}

WireCell2dToy::ToyMatrixMarkov::ToyMatrixMarkov(WireCell2dToy::ToyMatrix &toybefore, WireCell2dToy::ToyMatrix &toycur, WireCell2dToy::ToyMatrix &toyafter, WireCell2dToy::MergeToyTiling &mergebefore, WireCell2dToy::MergeToyTiling &mergecur, WireCell2dToy::MergeToyTiling &mergeafter, WireCell::GeomCellSelection *allmcell1,int recon_t1, int recon_t2, double penalty, double penalty_ncpt)
  : penalty_ncpt(penalty_ncpt)
{
  ToyMatrixMarkov(&toybefore, &toycur, &toyafter, &mergebefore, &mergecur, &mergeafter, allmcell, recon_t1, recon_t2, penalty, penalty_ncpt);
}


WireCell2dToy::ToyMatrixMarkov::ToyMatrixMarkov(WireCell2dToy::ToyMatrix *toybefore, WireCell2dToy::ToyMatrix *toycur, WireCell2dToy::ToyMatrix *toyafter, WireCell2dToy::MergeToyTiling *mergebefore, WireCell2dToy::MergeToyTiling *mergecur, WireCell2dToy::MergeToyTiling *mergeafter, WireCell::GeomCellSelection *allmcell1,int recon_t1, int recon_t2, double penalty, double penalty_ncpt)
  : penalty_ncpt(penalty_ncpt)
{
  ncount = 0;
  first_flag = 0;
  toymatrix = toycur;
  allmcell = allmcell1;
  mcindex = toymatrix->Get_mcindex();
  recon_threshold1 = recon_t1;
  recon_threshold2 = recon_t2;
  use_time_threshold = 0.75;
  
   //find good cells with time information and then use them ... 
  
  GeomCellSelection allmcell_p;
  if (mergebefore != 0)
    allmcell_p = mergebefore->get_allcell();
  GeomCellSelection allmcell_c = mergecur->get_allcell();
  GeomCellSelection allmcell_n;
  if (mergeafter != 0)
    allmcell_n = mergeafter->get_allcell();
  
  //form two vectors 
  std::vector<int> already_removed; //dummy
 
  for (int i=0;i!=allmcell_c.size();i++){
    MergeGeomCell *mcell_c = (MergeGeomCell*)allmcell_c[i];
    int index_c = toycur->Get_mcindex(mcell_c);
    
    int flag_before = 0;
    if (toybefore != 0){
      if (toybefore->Get_Solve_Flag()!=0 ){
	for (int j=0;j!=allmcell_p.size();j++){
	  MergeGeomCell *mcell_p = (MergeGeomCell*)allmcell_p[j];
	  int index_p = toybefore->Get_mcindex(mcell_p);
	  double charge = toybefore->Get_Cell_Charge(mcell_p,1);
	  if ( charge > recon_threshold2 && mcell_c->Overlap(*mcell_p)){
	    auto it = find(use_time.begin(),use_time.end(),index_c);
	    if (it == use_time.end()){
	      use_time.push_back(index_c);
	    }
	    flag_before = 1;
	  }
	}
      }
    }

    if (flag_before == 1){
      cell_penal[index_c] += penalty;
    }

    int flag_after = 0;
    if (toyafter !=0){
      if (toyafter->Get_Solve_Flag()!=0){
	for (int j=0;j!=allmcell_n.size();j++){
	  MergeGeomCell *mcell_n = (MergeGeomCell*)allmcell_n[j];
	  int index_n = toyafter->Get_mcindex(mcell_n);
	  double charge = toyafter->Get_Cell_Charge(mcell_n,1);
	  if ( charge > recon_threshold2 && mcell_c->Overlap(*mcell_n)){
	    auto it = find(use_time.begin(),use_time.end(),index_c);
	    if (it == use_time.end()){
	      use_time.push_back(index_c);
	    }
	    flag_after = 1;
	  }
	}
      }
    }

    if (flag_after == 1){
      cell_penal[index_c] += penalty;
    }
  }
  
  // fill in the cell_connect map;
  GeomCellCellsMap& cells_map = mergecur->get_not_compatible_cells_map();
  GeomCellSelection cells = mergecur->get_allcell();
  for (int i=0;i!=cells.size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*)cells.at(i);
    int index = toycur->Get_mcindex(mcell);
    if (cells_map[mcell].size() > 0){
      std::vector<int> cells_indices;
      for (int j=0;j!=cells_map[mcell].size();j++){
	int index_c = toycur->Get_mcindex(cells_map[mcell].at(j));
	cells_indices.push_back(index_c);
      }
      cells_ncpt[index] = cells_indices;
    }
  }



  for (int i=0;i!=allmcell_c.size();i++){
    auto it = find(use_time.begin(),use_time.end(),i);
    if (it==use_time.end())
      already_removed.push_back(i);
  }
  use_time.clear();
  
  //initialize
  //toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(*toymatrix);  // hold the current results 
  
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, use_time, *toymatrix,1,0);
  std::cout << "With Time: " << toymatrixkalman->Get_numz() << " " << allmcell_c.size() << " " <<  already_removed.size() << std::endl;
  // Find a sub-set that is not degenerated
  // put things into use_time
  find_subset(*toymatrixkalman,*toycur,use_time);
  already_removed.clear();
  
  //std::cout << use_time.size() << std::endl;
  // recalculate
  delete toymatrixkalman;
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, use_time, *toycur,1,0);
  toymatrix->Set_Solve_Flag(0);
  toymatrix->Set_chi2(-1);

  if (use_time.size() != 0 ){
    
    std::cout << "Start making guesses " << std::endl;
    if (allmcell_c.size()<50){
      while (ncount < 1e4 
	     && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	     && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	     && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	     && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800) 
	     && (cur_chi2 > 1*(cur_dof+0.1) || ncount < 10) 
	     ){
	if (ncount == 1 || (ncount % 1000 ==0&&ncount >0) )
	  std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
	
	make_guess();
	next_chi2 = toymatrix->Get_Chi2();
	next_dof = toymatrix->Get_ndf();
      }
    }else if (allmcell_c.size()<75){
      while (ncount < 5000 
	     && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	     && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	     && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	     && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800) 
	     && (cur_chi2 > 1*(cur_dof+0.1) || ncount < 10) 
	     ){
	if (ncount == 1 || (ncount % 1000 ==0&&ncount >0) )
	  std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
	
	make_guess();
	next_chi2 = toymatrix->Get_Chi2();
	next_dof = toymatrix->Get_ndf();
      }
    }else if (allmcell_c.size()<110){
      while (ncount < 2000 
	     && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	     && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	     && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	     && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800) 
	     && (cur_chi2 > 1*(cur_dof+0.1) || ncount < 10) 
	     ){
	if (ncount == 1 || (ncount % 1000 ==0&&ncount >0) )
	  std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
	
	make_guess();
	next_chi2 = toymatrix->Get_Chi2();
	next_dof = toymatrix->Get_ndf();
      }
    }else{
      //    }else if (allmcell_c.size()<150){
      for (int i=0;i!=200;i++){
	make_guess();
	next_chi2 = toymatrix->Get_Chi2();
	next_dof = toymatrix->Get_ndf();
      }
    }
  }else{ // test use ... 
    if (allmcell_c.size()<30){
      while (ncount < 1e4 
	     && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	     && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	     && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	     && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800) 
	     && (cur_chi2 > 1*(cur_dof+0.1) || ncount < 10) 
	     ){
	if (ncount == 1 || (ncount % 1000 ==0&&ncount >0) )
	  std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
	
	make_guess();
	next_chi2 = toymatrix->Get_Chi2();
	next_dof = toymatrix->Get_ndf();
      }
    }else if (allmcell_c.size()<60){
      while (ncount < 5000 
	     && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	     && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	     && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	     && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800) 
	     && (cur_chi2 > 1*(cur_dof+0.1) || ncount < 10) 
	     ){
	if (ncount == 1 || (ncount % 1000 ==0&&ncount >0) )
	  std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
	
	make_guess();
	next_chi2 = toymatrix->Get_Chi2();
	next_dof = toymatrix->Get_ndf();
      }
    }
  }
}



WireCell2dToy::ToyMatrixMarkov::ToyMatrixMarkov(WireCell2dToy::ToyMatrix *toymatrix1,WireCell::GeomCellSelection *allmcell1, int recon_t1, int recon_t2)
  : penalty_ncpt(0)
{
  ncount = 0;
  first_flag = 0;
  toymatrix = toymatrix1; //save the matrix in here
  allmcell = allmcell1;
  mcindex = toymatrix->Get_mcindex();
  recon_threshold1 = recon_t1;
  recon_threshold2 = recon_t2;
  use_time_threshold = 0.75;

  for (int i=0;i!=mcindex;i++){
    cell_penal[i] = 0;
  }


  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(*toymatrix,0);  // hold the current results 
  
  if (allmcell->size() < 50){
    while (ncount < 1e4 
	   && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	   && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	   && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	   && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800) 
	   && (cur_chi2 > (cur_dof+0.1) || ncount < 10) 
	   ){
      if (ncount == 1 || (ncount % 1000 ==0&&ncount >0) )
	std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
      
      make_guess();
      next_chi2 = toymatrix->Get_Chi2();
      next_dof = toymatrix->Get_ndf();
    }
  }else if (allmcell->size()<75){
    while (ncount < 5000
	   && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	   && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	   && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	   && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800) 
	   && (cur_chi2 > (cur_dof+0.1) || ncount < 10) 
	   ){
      if (ncount == 1 || (ncount % 1000 ==0&&ncount >0) )
	std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
      
      make_guess();
      next_chi2 = toymatrix->Get_Chi2();
      next_dof = toymatrix->Get_ndf();
    }
  }else if (allmcell->size()<110){
    while (ncount < 2000 
	   && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	   && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	   && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	   && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800) 
	   && (cur_chi2 > (cur_dof+0.1) || ncount < 10) 
	   ){
      if (ncount ==1 || (ncount % 1000 ==0&&ncount >0) )
	std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
      
      make_guess();
      next_chi2 = toymatrix->Get_Chi2();
      next_dof = toymatrix->Get_ndf();
    }
  }else{
    //  }else if (allmcell->size()<150){
    for (int i=0;i!=200;i++){
      make_guess();
      next_chi2 = toymatrix->Get_Chi2();
      next_dof = toymatrix->Get_ndf();
    }
  }
  
  //  std::cout << cur_chi2 << " " << ncount << std::endl;
  

  

  
}

WireCell2dToy::ToyMatrixMarkov::~ToyMatrixMarkov(){
  delete toymatrixkalman;
}

void WireCell2dToy::ToyMatrixMarkov::make_guess(){
  
  ncount ++; 
  
  if (first_flag ==0){
    nlevel = 0;
    //    std::cout << toymatrixkalman->Get_mcindex() << " " << toymatrixkalman->Get_numz() << std::endl;
    Iterate(*toymatrixkalman);
    first_flag = 1;
  }
  //need to save the results???
  
  cur_chi2 = toymatrix->Get_Chi2();
 

  cur_dof = toymatrix->Get_ndf();

  toymatrix->Update_pred();

  double_t max_res = 0;

  for (int i=0;i!= (*allmcell).size();i++){
    MergeGeomCell *mcell = (MergeGeomCell*)(*allmcell).at(i);
    double charge = toymatrix->Get_Cell_Charge(mcell,1);
    double charge_err = toymatrix->Get_Cell_Charge(mcell,2);
    
    cur_cell_index.push_back(toymatrix->Get_mcindex(mcell));

    if (charge ==0 && charge_err ==0 ){
      cur_cell_status.push_back(0); // removed for matrix calculation
    }else{
      cur_cell_status.push_back(1); // not removed
    }
    
    if (charge > recon_threshold1){ // hard coded to be fixed later
      cur_cell_pol.push_back(1); // on
    }else{
      cur_cell_pol.push_back(0); //off
    }

    //calculate residual for each cell, part of toymatrix
    // with residual and polarity, one can calculate probability
    cell_res.push_back(toymatrix->Get_residual(mcell));

    if (toymatrix->Get_residual(mcell) > max_res) 
      max_res = toymatrix->Get_residual(mcell);

    cell_charge.push_back(charge);

    if (cur_cell_status.at(i)==1){
      double ratio;
      
      if (cell_res.at(i)!=0){
	ratio = cell_charge.at(i)/cell_res.at(i);
      }else{
	ratio = cell_charge.at(i)/10.;
      }
      
      if (ratio > 10.){
	cell_prob.push_back(0.8); 
      }else if (ratio <=10 && ratio>=-10.){
	cell_prob.push_back(0.6);
      }else{
	cell_prob.push_back(0.2);
      }
    }else{
      if (cell_res.at(i) > 1000){
	cell_prob.push_back(0.6);
      }else{
	cell_prob.push_back(0.4);
      }
    }

    auto it = find(use_time.begin(),use_time.end(),toymatrix->Get_mcindex(mcell));
    if (it != use_time.end()){

      if (cell_prob.at(i) < use_time_threshold)
	cell_prob.at(i) = use_time_threshold;
      
    }

    //   std::cout << cur_cell_status.at(i) << " " << cur_cell_pol.at(i) <<
    //   " " << cell_res.at(i) << " " << std::endl;
  }


  cell_set.clear();
  //rank stuff ... 
  // first probability, second bias
  for (int i=0;i!= (*allmcell).size();i++){
    CellRankPair a(i,cell_prob.at(i)+gRandom->Uniform(0,1)/100.);
    cell_set.insert(a);
  }
  
  CellRankSet xp,middle,bad_one;
  
  

  
  toymatrixkalman->Get_no_need_remove().clear();
  toymatrixkalman->Get_already_removed().clear();
  
  //fill three cases
  for (auto it= cell_set.begin(); it!=cell_set.end();it++){
    //std::cout << (*it).first << " " << (*it).second << std::endl;
    int index = (*it).first ;
    double score = (*it).second;
    CellRankPair a(index,score);
    double r = gRandom->Uniform(0,1);
    if (r < cell_prob.at(index)){ 
      // start to insert
      if (xp.size()==0){
	
	xp.insert(a);
	toymatrixkalman->Get_no_need_remove().push_back(index);
      }else{
	//first push in and judge 
	toymatrixkalman->Get_no_need_remove().push_back(index);
	if (toymatrixkalman->Cal_numz(*toymatrix)==0){
	  xp.insert(a);
	}else{
	  bad_one.insert(a);
	  toymatrixkalman->Get_no_need_remove().pop_back();
	}
      }
    }else{
      middle.insert(a);
    } 
  }

  //std::cout << xp.size() << " " << middle.size() << " " << bad_one.size() << std::endl;

  while (middle.size()!=0){
    std::vector<int> del_list;
    for (auto it = middle.begin();it!=middle.end();it++){
      double r = gRandom->Uniform(0,1);
      int index = (*it).first ;
      double score = (*it).second;
      CellRankPair a(index,score);
      
      if (r < cell_prob.at(index)){ 
	del_list.push_back(index);
	//  	middle.erase(it);
	toymatrixkalman->Get_no_need_remove().push_back(index);
	if (toymatrixkalman->Cal_numz(*toymatrix)==0){
	  xp.insert(a);
	}else{
	  bad_one.insert(a);
	  toymatrixkalman->Get_no_need_remove().pop_back();
	}
      }
    }
    
    for (int i = 0; i!=del_list.size();i++){
      for (auto it = middle.begin();it!=middle.end();it++){
	if ((*it).first == del_list.at(i)){
	  middle.erase(it);
	  break;
	}
      }
    }
    
  }

  //std::cout << xp.size() << " " << middle.size() << " " << bad_one.size() << std::endl;
 
  for (int i=0;i!=mcindex;i++){
    auto it = find(toymatrixkalman->Get_no_need_remove().begin(),toymatrixkalman->Get_no_need_remove().end(),i);
    if (it==toymatrixkalman->Get_no_need_remove().end())
      toymatrixkalman->Get_already_removed().push_back(i);
  }
  toymatrixkalman->Get_no_need_remove().clear();
  
  double chi2_p = 0;
  for (int j = 0; j!=toymatrixkalman->Get_already_removed().size(); j++){
    chi2_p += cell_penal[toymatrixkalman->Get_already_removed().at(j)];
  }

  if (penalty_ncpt>0){
    //add the penalty due to not compatible cells (ncpt)
    //int flag_test = 0;
    for (auto it = cells_ncpt.begin(); it!= cells_ncpt.end(); it++){
      int index1 = it->first;
      if (find(toymatrixkalman->Get_already_removed().begin(),toymatrixkalman->Get_already_removed().end(),index1) != toymatrixkalman->Get_already_removed().end()) continue;
      std::vector<int> indices = it->second;
      
      for (int j=0;j!=indices.size();j++){
	int index2 = indices.at(j);
	if (find(toymatrixkalman->Get_already_removed().begin(),toymatrixkalman->Get_already_removed().end(),index2) == toymatrixkalman->Get_already_removed().end()){
	  chi2_p += penalty_ncpt;
	  // flag_test = 1;
	  // break;
	}
      }
      // if (flag_test ==1) break;
    }
    // if (flag_test == 1){
    //   chi2_p += penalty_ncpt;
    // }
  }


  toymatrixkalman->Set_penalty(chi2_p);

  //why initiate again???  //calculate chi2 ... 
  toymatrixkalman->init(*toymatrix);
  //std::cout << toymatrixkalman->Get_already_removed().size() << std::endl;
}

void WireCell2dToy::ToyMatrixMarkov::Iterate(WireCell2dToy::ToyMatrixKalman &toykalman){
  nlevel++;
  if (toykalman.Get_numz()!=0&&toymatrix->Get_Chi2()<0){
    for (int i=0;i!=toykalman.Get_mcindex();i++){
      auto it1 = find(toykalman.Get_already_removed().begin(),toykalman.Get_already_removed().end(),i);
      auto it2 = find(toykalman.Get_no_need_remove().begin(),toykalman.Get_no_need_remove().end(),i);
      if (it1 == toykalman.Get_already_removed().end() && it2 == toykalman.Get_no_need_remove().end()){
	std::vector<int> already_removed = toykalman.Get_already_removed();
	already_removed.push_back(i);

	// test ... 
	if (toykalman.Get_no_need_remove().size() + already_removed.size() >= toykalman.Get_mcindex()){
	  //std::cout << already_removed.size() << " " << toykalman.Get_no_need_remove().size() << std::endl;

	  //  for (int i=0;i!=allmcell_c.size();i++){
	  //   auto it = find(use_time.begin(),use_time.end(),i);
	  //   if (it==use_time.end())
	  //     already_removed.push_back(i);
	  // }
	  // use_time.clear();
	  
	  // //initialize
	  // //toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(*toymatrix);  // hold the current results 
	  
	  // toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, use_time, *toymatrix,1,0);
	  // std::cout << "With Time: " << toymatrixkalman->Get_numz() << " " << allmcell_c.size() << " " <<  already_removed.size() << std::endl;
	  // // Find a sub-set that is not degenerated
	  // // put things into use_time
	  // find_subset(*toymatrixkalman,toycur,use_time);

	  std::vector<int> temp_already_removed;
	  std::vector<int> temp_no_need_remove;
	  for (int k=0;k!=toykalman.Get_mcindex();k++){
	    auto it4 = find(toykalman.Get_no_need_remove().begin(),toykalman.Get_no_need_remove().end(),k);
	    if (it4 == toykalman.Get_no_need_remove().end())
	      temp_already_removed.push_back(k);
	  }
	  

	  //temp_no_need_remove.clear();
	  //initialize
	  WireCell2dToy::ToyMatrixKalman temp_tk(temp_already_removed, temp_no_need_remove, *toymatrix,0,0);
	  //	  std::cout << temp_tk.Get_numz() << " " << temp_already_removed.size() << " " << temp_no_need_remove.size() << std::endl;
	  // Find a sub-set that is not degenerated
	  find_subset(temp_tk,*toymatrix,temp_no_need_remove);

	  if (temp_no_need_remove.size()!=0){
	    toykalman.Get_no_need_remove().clear();
	    for (int k=0;k!=temp_no_need_remove.size();k++){
	      toykalman.Get_no_need_remove().push_back(temp_no_need_remove.at(k));
	    }
	    already_removed.clear();
	  }
	  //std::cout << already_removed.size() << " " << toykalman.Get_no_need_remove().size() << " " << temp_no_need_remove.size() << std::endl;
	}
	

	double chi2_p = 0;
	for (int j = 0; j!=already_removed.size(); j++){
	  chi2_p += cell_penal[already_removed.at(j)];
	}

	if (penalty_ncpt>0){
	  //add the penalty due to not compatible cells (ncpt)
	  //int flag_test = 0;
	  for (auto it = cells_ncpt.begin(); it!= cells_ncpt.end(); it++){
	    int index1 = it->first;
	    if (find(already_removed.begin(),already_removed.end(),index1) != already_removed.end()) continue;
	    std::vector<int> indices = it->second;
	    
	    for (int j=0;j!=indices.size();j++){
	      int index2 = indices.at(j);
	      if (find(already_removed.begin(),already_removed.end(),index2) == already_removed.end()){
		chi2_p += penalty_ncpt;
		// flag_test = 1;
		// break;
	      }
	    }
	    // if (flag_test ==1) break;
	  }
	  // if (flag_test == 1){
	  //   chi2_p += penalty_ncpt;
	  // }
	}




	//std::cout << nlevel << " " << already_removed.size() << " " << toykalman.Get_no_need_remove().size() << " " << toykalman.Get_numz() << std::endl;

	WireCell2dToy::ToyMatrixKalman kalman(already_removed,toykalman.Get_no_need_remove(),*toymatrix,0,1,chi2_p);
	//std::cout << kalman.Get_numz() << std::endl;
	// this part seems to be able to improve
	// to get the first solution, one should not take too much time ...
	// consider to improve in the future ...  i.e. add cut on the number of get_already_removed vs. what's left ... 
	Iterate(kalman);
	nlevel--;
	
	toykalman.Get_no_need_remove().push_back(i);
	if (toymatrix->Get_Chi2()>0) break;
      }
    }
  }
}
