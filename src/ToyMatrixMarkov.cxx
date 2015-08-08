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
	
	if (kalman.Get_numz()==toymatrix.Get_numz()){
	}else{
	  find_subset(kalman,toymatrix1,vec);
	  //move to "no need to remove"??
	  toymatrix.Get_no_need_remove().push_back(i);
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

WireCell2dToy::ToyMatrixMarkov::ToyMatrixMarkov(WireCell2dToy::ToyMatrix &toycur,WireCell2dToy::MergeToyTiling &mergecur, WireCell::GeomCellSelection *allmcell1, WireCell::GeomCellSelection &cells, int recon_t1, int recon_t2){
  ncount = 0;
  first_flag = 0;
  toymatrix = &toycur;
  allmcell = allmcell1;
  mcindex = toymatrix->Get_mcindex();
  recon_threshold1 = recon_t1;
  recon_threshold2 = recon_t2;
  use_time_threshold = 0.95;

  GeomCellSelection allmcell_c = mergecur.get_allcell();

  //form two vectors 
  std::vector<int> already_removed; //dummy
  
  for (int i=0;i!=allmcell_c.size();i++){
    MergeGeomCell *mcell_c = (MergeGeomCell*)allmcell_c[i];
    int index_c = toycur.Get_mcindex(mcell_c);
    
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
  
  // recalculate
  delete toymatrixkalman;
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, use_time, toycur,1,0);
  toymatrix->Set_Solve_Flag(0);
  toymatrix->Set_chi2(-1);


  
  while (ncount < 1e4 
	 && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	 && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	 && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	 && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800) 
	 ){
    if (ncount == 1 || (ncount % 1000 ==0&&ncount >0) )
      std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
    
    make_guess();
    next_chi2 = toymatrix->Get_Chi2();
    next_dof = toymatrix->Get_ndf();
  }


}


WireCell2dToy::ToyMatrixMarkov::ToyMatrixMarkov(WireCell2dToy::ToyMatrix &toybefore, WireCell2dToy::ToyMatrix &toycur, WireCell2dToy::ToyMatrix &toyafter, WireCell2dToy::MergeToyTiling &mergebefore, WireCell2dToy::MergeToyTiling &mergecur, WireCell2dToy::MergeToyTiling &mergeafter, WireCell::GeomCellSelection *allmcell1,int recon_t1, int recon_t2){
  ncount = 0;
  first_flag = 0;
  toymatrix = &toycur;
  allmcell = allmcell1;
  mcindex = toymatrix->Get_mcindex();
  recon_threshold1 = recon_t1;
  recon_threshold2 = recon_t2;
  use_time_threshold = 0.75;
  
   //find good cells with time information and then use them ... 
  
   GeomCellSelection allmcell_p = mergebefore.get_allcell();
  GeomCellSelection allmcell_c = mergecur.get_allcell();
  GeomCellSelection allmcell_n = mergeafter.get_allcell();
  
  //form two vectors 
  std::vector<int> already_removed; //dummy
 
  for (int i=0;i!=allmcell_c.size();i++){
    MergeGeomCell *mcell_c = (MergeGeomCell*)allmcell_c[i];
    int index_c = toycur.Get_mcindex(mcell_c);
    
    if (toybefore.Get_Solve_Flag()!=0 ){
      for (int j=0;j!=allmcell_p.size();j++){
	MergeGeomCell *mcell_p = (MergeGeomCell*)allmcell_p[j];
	int index_p = toybefore.Get_mcindex(mcell_p);
	double charge = toybefore.Get_Cell_Charge(mcell_p,1);
	if ( charge > recon_threshold2 && mcell_c->Overlap(*mcell_p)){
	  auto it = find(use_time.begin(),use_time.end(),index_c);
	  if (it == use_time.end()){
	    use_time.push_back(index_c);
	  }
	}
      }
    }
    if (toyafter.Get_Solve_Flag()!=0){
      for (int j=0;j!=allmcell_n.size();j++){
	MergeGeomCell *mcell_n = (MergeGeomCell*)allmcell_n[j];
	int index_n = toyafter.Get_mcindex(mcell_n);
	double charge = toyafter.Get_Cell_Charge(mcell_n,1);
	if ( charge > recon_threshold2 && mcell_c->Overlap(*mcell_n)){
	  auto it = find(use_time.begin(),use_time.end(),index_c);
	  if (it == use_time.end()){
	    use_time.push_back(index_c);
	  }
	}
      }
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
  
  // recalculate
  delete toymatrixkalman;
  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(already_removed, use_time, toycur,1,0);
  toymatrix->Set_Solve_Flag(0);
  toymatrix->Set_chi2(-1);


  
  while (ncount < 1e4 
	 && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	 && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	 && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	 && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800) 
	 ){
    if (ncount == 1 || (ncount % 1000 ==0&&ncount >0) )
      std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
    
    make_guess();
    next_chi2 = toymatrix->Get_Chi2();
    next_dof = toymatrix->Get_ndf();
  }
}





WireCell2dToy::ToyMatrixMarkov::ToyMatrixMarkov(WireCell2dToy::ToyMatrix *toymatrix1,WireCell::GeomCellSelection *allmcell1, int recon_t1, int recon_t2){
  ncount = 0;
  first_flag = 0;
  toymatrix = toymatrix1; //save the matrix in here
  allmcell = allmcell1;
  mcindex = toymatrix->Get_mcindex();
  recon_threshold1 = recon_t1;
  recon_threshold2 = recon_t2;
  use_time_threshold = 0.75;

  toymatrixkalman = new WireCell2dToy::ToyMatrixKalman(*toymatrix,0);  // hold the current results 
  
  while (ncount < 1e4 
	 && (cur_chi2 > 5*(cur_dof+0.1) || ncount < 6000) 
	 && (cur_chi2 > 4*(cur_dof+0.1) || ncount < 3000) 
	 && (cur_chi2 > 3*(cur_dof+0.1) || ncount < 1500) 
	 && (cur_chi2 > 2*(cur_dof+0.1) || ncount < 800) 
	 ){
    if (ncount == 1 || (ncount % 1000 ==0&&ncount >0) )
      std::cout << "MCMC: " << ncount << " " << cur_chi2 << " " << cur_dof << std::endl;
    
    make_guess();
    next_chi2 = toymatrix->Get_Chi2();
    next_dof = toymatrix->Get_ndf();
  }
  
  //  std::cout << cur_chi2 << " " << ncount << std::endl;
  

  

  
}

WireCell2dToy::ToyMatrixMarkov::~ToyMatrixMarkov(){
  delete toymatrixkalman;
}

void WireCell2dToy::ToyMatrixMarkov::make_guess(){
  
  ncount ++; 
  
  if (first_flag ==0){
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
  

  //why initiate again???  //calculate chi2 ... 
  toymatrixkalman->init(*toymatrix);
  //std::cout << toymatrixkalman->Get_already_removed().size() << std::endl;
}

void WireCell2dToy::ToyMatrixMarkov::Iterate(WireCell2dToy::ToyMatrixKalman &toykalman){
  if (toykalman.Get_numz()!=0&&toymatrix->Get_Chi2()<0){
    for (int i=0;i!=toykalman.Get_mcindex();i++){
      auto it1 = find(toykalman.Get_already_removed().begin(),toykalman.Get_already_removed().end(),i);
      auto it2 = find(toykalman.Get_no_need_remove().begin(),toykalman.Get_no_need_remove().end(),i);
      if (it1 == toykalman.Get_already_removed().end() && it2 == toykalman.Get_no_need_remove().end()){
	std::vector<int> already_removed = toykalman.Get_already_removed();
	already_removed.push_back(i);
	WireCell2dToy::ToyMatrixKalman kalman(already_removed,toykalman.Get_no_need_remove(),*toymatrix,0);
	
	Iterate(kalman);
	
	toykalman.Get_no_need_remove().push_back(i);
	if (toymatrix->Get_Chi2()>0) break;
      }
    }
  }
}
