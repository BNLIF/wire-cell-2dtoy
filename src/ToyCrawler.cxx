#include "WireCell2dToy/ToyCrawler.h"

using namespace WireCell;

WireCell2dToy::ToyCrawler::ToyCrawler(MergeSpaceCellSelection& mcells, int flag){

  CreateClusterTrack(mcells);
  FormGraph(); 

  //std::cout << "Merge Clusters " << std::endl; 
  // Merge Cluster ...
  MergeCTrack();  
  // Further merge trying to extend into other tracks ... 
  FurtherExtendCTrack();
  
  PurgeMergeCTrack();
  

  //this is not a good assumption
  CleanUpCTTrack(flag);
  
  
  // for (int i=0;i!=all_mergeclustertrack.size();i++){
  //   MergeClusterTrack *mct = all_mergeclustertrack.at(i);
  //   std::cout << i << " " << mct->Get_allmcells().size() << std::endl;
  // }
}

void WireCell2dToy::ToyCrawler::CleanUpCTTrack(int flag){
  //Sort the track first
  MergeClusterTrackSet MCT_set;
  for (int i =0;i!=all_mergeclustertrack.size();i++){
    MergeClusterTrack *mct = all_mergeclustertrack.at(i);
    MCT_set.insert(mct);
  }
  all_mergeclustertrack.clear();
  for (auto it = MCT_set.begin();it!=MCT_set.end();it++){
    MergeClusterTrack *mct = *it;
    all_mergeclustertrack.push_back(mct);
    //std::cout << mct->Get_allmcells().size() << std::endl;
  }

  MergeClusterTrackSelection to_be_removed;
  MergeClusterTrackSelection examined;
  
  for (int i=0;i!=all_mergeclustertrack.size();i++){
    MergeClusterTrack *mct = all_mergeclustertrack.at(i);
    examined.push_back(mct);
    
    for (int j = 0;j!=mct->Get_allmcells().size();j++){
      MergeSpaceCell *mcell = mct->Get_allmcells().at(j);
      for (int k=0; k!= mcells_mct_map[mcell].size();k++){
	MergeClusterTrack *mct1 =  mcells_mct_map[mcell].at(k);
	auto it1 = find(to_be_removed.begin(),to_be_removed.end(),mct1);
	auto it2 = find(examined.begin(),examined.end(),mct1);
	if (it1==to_be_removed.end() && it2 == examined.end()){
	  
	  int n_common = 0;
	  int n_diff = 0;
	  float ncells_common = 0;
	  float ncells_diff = 0;
 
	  for (int i1 = 0; i1!=mct1->Get_allmcells().size();i1++){
	    MergeSpaceCell *mcell1 = mct1->Get_allmcells().at(i1);
	    auto it3 = find(mct->Get_allmcells().begin(),mct->Get_allmcells().end(),mcell1);
	    if (it3 == mct->Get_allmcells().end()){
	      n_diff ++;
	      if (flag == 1){
		ncells_diff += mcell1->Get_all_spacecell().size();
	      }else{
		ncells_diff += mcell1->Get_Charge();
	      }
	    }else{
	      n_common ++;
	      if (flag == 1){
		ncells_common += mcell1->Get_all_spacecell().size();
	      }else{
		ncells_common += mcell1->Get_Charge();
	      }
	    }
	  }

	  if ((flag==1)&&(
			  (n_common > n_diff && n_diff < 2 && ncells_common >= 20*ncells_diff) )){
	    //merge
	    to_be_removed.push_back(mct1);
	    mct->AddTrack(mct1);
	  }else if ((flag==2)&&
		    (n_common >= n_diff && n_diff < 2 && ncells_common >= 20*ncells_diff) ){
	    //merge
	    to_be_removed.push_back(mct1);
	    mct->AddTrack(mct1);
	  }else{
	    //not merge
	    examined.push_back(mct1);
	  }
	}
      }
    }
  }

  
  for (int i=0;i!=to_be_removed.size();i++){
    auto it = find(all_mergeclustertrack.begin(),all_mergeclustertrack.end(),to_be_removed.at(i));
    all_mergeclustertrack.erase(it);
  }

  for (int i=0;i!=all_mergeclustertrack.size();i++){
    MergeClusterTrack *mct = all_mergeclustertrack.at(i);
    mct->Organize();
  }
  
  

  // Update map again ... 
  mcells_mct_map.clear();
  for (int i=0;i!=all_mergeclustertrack.size();i++){
    MergeClusterTrack *mct = all_mergeclustertrack.at(i);
    for (int j=0; j!=mct->Get_allmcells().size() ;j++){
      MergeSpaceCell *mcell = mct->Get_allmcells().at(j);
      
      if (mcells_mct_map.find(mcell) == mcells_mct_map.end()){
	MergeClusterTrackSelection temp1;
	temp1.push_back(mct);
	mcells_mct_map[mcell] = temp1;
      }else{
	mcells_mct_map[mcell].push_back(mct);
      }
    }
  }



}


void WireCell2dToy::ToyCrawler::PurgeMergeCTrack(){
  //remove merged cluster in which all the merge space cell are accounted for ... 
  MergeClusterTrackSelection temp1;
  for (int i=0;i!=all_mergeclustertrack.size();i++){
    MergeClusterTrack *mct = all_mergeclustertrack.at(i);
    int flag = 0;
    for (int j=0; j!=mct->Get_allmcells().size() ;j++){
      MergeSpaceCell *mcell = mct->Get_allmcells().at(j);
      
      if (mcells_mct_map[mcell].size()==1){
	flag = 1;
	break;
      }
    }
    if (flag == 0){
      temp1.push_back(mct);
    }

  }

  for (int i=0;i!=temp1.size();i++){
    auto it = find(all_mergeclustertrack.begin(),all_mergeclustertrack.end(),temp1.at(i));
    all_mergeclustertrack.erase(it);
  }
  temp1.clear();



  // Update map again ... 
  mcells_mct_map.clear();
  for (int i=0;i!=all_mergeclustertrack.size();i++){
    MergeClusterTrack *mct = all_mergeclustertrack.at(i);
    for (int j=0; j!=mct->Get_allmcells().size() ;j++){
      MergeSpaceCell *mcell = mct->Get_allmcells().at(j);
      
      if (mcells_mct_map.find(mcell) == mcells_mct_map.end()){
	MergeClusterTrackSelection temp1;
	temp1.push_back(mct);
	mcells_mct_map[mcell] = temp1;
      }else{
	mcells_mct_map[mcell].push_back(mct);
      }
    }
  }

}



void WireCell2dToy::ToyCrawler::FurtherExtendCTrack(){
  
  
  MergeClusterTrackSelection temp1;
  
  int flag3 = 1;
  flag3 = 0;
  float cut_angle = 15;
  
  for (int i=0;i!=all_mergeclustertrack.size();i++){
    MergeClusterTrack *mct = all_mergeclustertrack.at(i);
    MergeSpaceCellSelection mcells = mct->Get_allmcells();
    
    if (mcells.size()>=2){
      MergeSpaceCell *fvertex = mcells.front();
      MergeSpaceCell *fvertex_next = mcells.at(1);
      //determine direction
      int flag1 = 0;
      if (fvertex->Get_Center().x - fvertex_next->Get_Center().x>0){
	flag1 = 1;
      }else if (fvertex->Get_Center().x - fvertex_next->Get_Center().x<0){
	flag1 = -1;
      }
      
      if (flag1!=0){
	mct->SC_Hough(fvertex->Get_Center());
	double theta = mct->Get_Theta();
	double phi = mct->Get_Phi();
	
	mct->SC_Hough(fvertex->Get_Center(),-1,2);
	double theta_m = mct->Get_Theta();
	double phi_m = mct->Get_Phi();
	
	//save the one that satisfy the requirement
	MergeSpaceCellSelection temp;
	MergeSpaceCell *cur_msc = fvertex;
	MergeSpaceCellSelection cur_cells = mcells_map[cur_msc];
	
	int flag2 = 1;
	while(flag2 == 1){
	  //insert temp
	  flag2 = 0;
	  for (int j=0;j!=cur_cells.size();j++){
	    MergeSpaceCell *next_msc = cur_cells.at(j);
	    if ((next_msc->Get_Center().x-cur_msc->Get_Center().x)*flag1>0){
	      //judge overlap?
	      if (next_msc->CrossCell(fvertex->Get_Center(),theta_m,phi_m)){
		temp.push_back(next_msc);
		cur_msc = next_msc;
		cur_cells = mcells_map[cur_msc];
		flag2 = 1;
		break;
	      }
	    }
	  }
	}
	
	//Now temp contains all the cells to be merged
	
	// add one cell in, and then figure out whether to add its track
	flag2 = 1;
	while(flag2){
	  flag2 = 0;
	  for (int j = 0;j!=temp.size();j++){
	    MergeSpaceCell *mcell = temp.at(j);
	    auto it = find(mct->Get_allmcells().begin(),mct->Get_allmcells().end(),mcell);
	    if (it == mct->Get_allmcells().end()){
	      MergeClusterTrackSelection MCTS = mcells_mct_map[mcell];
	      // insert it in
	      int insert_flag = 1;
	      // for (int k=0;k!=MCTS.size();k++){
	      // 	if ( MCTS.at(k)->CheckCell(mcell)){
	      // 	  insert_flag = 1;
	      // 	  break;
	      // 	}
	      // }
	      // std::cout << "abc: " << insert_flag << std::endl;
	      
	      if (insert_flag ==1){
		mct->Add(mcell,-1);
		
		// Now find the associated clusterTrack
		// judge whether to put in temp1 to be deleted
		for (int k=0;k!=MCTS.size();k++){
		  MCTS.at(k)->SC_Hough(mcell->Get_Center());
		  float theta1 = MCTS.at(k)->Get_Theta();
		  float phi1 = MCTS.at(k)->Get_Phi();

		  MCTS.at(k)->SC_Hough(mcell->Get_Center(),-1,2);
		  float theta1_m = MCTS.at(k)->Get_Theta();
		  float phi1_m = MCTS.at(k)->Get_Phi();

		  if (MCTS.at(k)->CheckCell(mcell) ){
		    if ((fabs(theta1+theta-3.1415926)<cut_angle/180.*3.1415926 // 5 degrees
			 && fabs(fabs(phi1-phi)-3.1415926)<cut_angle/180.*3.1415926) ||
			(fabs(theta1_m+theta_m-3.1415926)<cut_angle/180.*3.1415926 // 5 degrees
			 && fabs(fabs(phi1_m-phi_m)-3.1415926)<cut_angle/180.*3.1415926)){
		      //add this Merged Track
		      mct->Add(MCTS.at(k),mcell,-1);
		      temp1.push_back(MCTS.at(k));
		      flag2 = 1;
		    }
		  }
		}
	      
	      }
	      
	      
	    }
	  }
	} //while(flag2)
	  //std::cout << i << " " << temp.size() << " " << temp1.size() << std::endl;
	if (temp1.size()>0) {
	  flag3 = 1;
	  goto abc;
	}
      }
      
      
      //deal with the front
      MergeSpaceCell *bvertex = mcells.back();
      MergeSpaceCell *bvertex_next = mcells.at(mcells.size()-2);
      
      
      flag1 = 0;
      if (bvertex->Get_Center().x - bvertex_next->Get_Center().x>0){
      	flag1 = 1;
      }else if (bvertex->Get_Center().x - bvertex_next->Get_Center().x<0){
      	flag1 = -1;
      }
      
      if (flag1!=0){
      	mct->SC_Hough(bvertex->Get_Center());
      	double theta = mct->Get_Theta();
      	double phi = mct->Get_Phi();
      	mct->SC_Hough(bvertex->Get_Center(),-1,2);
      	double theta_m = mct->Get_Theta();
      	double phi_m = mct->Get_Phi();
	
      	//save the one that satisfy the requirement
      	MergeSpaceCellSelection temp;
      	MergeSpaceCell *cur_msc = bvertex;
      	MergeSpaceCellSelection cur_cells = mcells_map[cur_msc];
	
      	int flag2 = 1;
      	while(flag2 == 1){
      	  //insert temp
      	  flag2 = 0;
      	  for (int j=0;j!=cur_cells.size();j++){
      	    MergeSpaceCell *next_msc = cur_cells.at(j);
      	    if ((next_msc->Get_Center().x-cur_msc->Get_Center().x)*flag1>0){
      	      //judge overlap?
      	      if (next_msc->CrossCell(bvertex->Get_Center(),theta,phi)){
      		temp.push_back(next_msc);
      		cur_msc = next_msc;
      		cur_cells = mcells_map[cur_msc];
      		flag2 = 1;
      		break;
      	      }
      	    }
      	  }
      	}
	
      	//Now temp contains all the cells to be merged
	
      	// add one cell in, and then figure out whether to add its track
      	flag2 = 1;
      	while(flag2){
      	  flag2 = 0;
      	  for (int j = 0;j!=temp.size();j++){
      	    MergeSpaceCell *mcell = temp.at(j);
      	    auto it = find(mct->Get_allmcells().begin(),mct->Get_allmcells().end(),mcell);
      	    if (it == mct->Get_allmcells().end()){
      	      MergeClusterTrackSelection MCTS = mcells_mct_map[mcell];
      	      int insert_flag = 1;
	      
      	      // for (int k=0;k!=MCTS.size();k++){
      	      // 	if ( MCTS.at(k)->CheckCell(mcell) ){
      	      // 	  insert_flag = 1;
      	      // 	  break;
      	      // 	}
      	      // }
	      
      	      if (insert_flag == 1){
      		// insert it in
      		mct->Add(mcell,1);
		
      		// Now find the associated clusterTrack
      		// judge whether to put in temp1 to be deleted
      		for (int k=0;k!=MCTS.size();k++){
      		  MCTS.at(k)->SC_Hough(mcell->Get_Center());
      		  float theta1 = MCTS.at(k)->Get_Theta();
      		  float phi1 = MCTS.at(k)->Get_Phi();
      		  MCTS.at(k)->SC_Hough(mcell->Get_Center(),-1,2);
      		  float theta1_m = MCTS.at(k)->Get_Theta();
      		  float phi1_m = MCTS.at(k)->Get_Phi();

      		  if (MCTS.at(k)->CheckCell(mcell) ){
      		    if ((fabs(theta1+theta-3.1415926)<cut_angle/180.*3.1415926 // 5 degrees
      			 && fabs(fabs(phi1-phi)-3.1415926)<cut_angle/180.*3.1415926)
      			||(fabs(theta1_m+theta_m-3.1415926)<cut_angle/180.*3.1415926 // 5 degrees
      			 && fabs(fabs(phi1_m-phi_m)-3.1415926)<cut_angle/180.*3.1415926)
      			){
		      
      		      //std::cout << theta << " " << theta1 << " " << phi << " " << phi1 << std::endl;
      		      //add this Merged Track
      		      mct->Add(MCTS.at(k),mcell,1);
      		      temp1.push_back(MCTS.at(k));
      		      flag2 = 1;
      		    }
      		  }
      		}
      	      }
	      
      	    }
      	  }
      	} //while(flag2)
      	  // std::cout << i << " " << temp.size() << " " << temp1.size() << std::endl;
	
      	if (temp1.size()>0) {
      	  flag3 = 1;
      	  goto abc;
      	}
      }


      
      
    }
  } // remove the merge cluster after this one
  
 abc:
  
  //std::cout << all_mergeclustertrack.size() << std::endl;
  for (int i=0;i!=temp1.size();i++){
    auto it = find(all_mergeclustertrack.begin(),all_mergeclustertrack.end(),temp1.at(i));
    all_mergeclustertrack.erase(it);
  }
  temp1.clear();
  //std::cout << all_mergeclustertrack.size() << std::endl;
  // }
  
  //Update Map 
  
  mcells_mct_map.clear();
  for (int i=0;i!=all_mergeclustertrack.size();i++){
    MergeClusterTrack *mct = all_mergeclustertrack.at(i);
    for (int j=0; j!=mct->Get_allmcells().size() ;j++){
      MergeSpaceCell *mcell = mct->Get_allmcells().at(j);
      
      if (mcells_mct_map.find(mcell) == mcells_mct_map.end()){
	MergeClusterTrackSelection temp1;
	temp1.push_back(mct);
	mcells_mct_map[mcell] = temp1;
      }else{
	mcells_mct_map[mcell].push_back(mct);
      }
    }
  }
    
    
    
}


void WireCell2dToy::ToyCrawler::MergeCTrack(){

  for (int i = 0;i!=all_clustertrack.size();i++){
    auto it = find(used_clustertrack.begin(),used_clustertrack.end(),all_clustertrack.at(i));
    if (it == used_clustertrack.end()){
      MergeClusterTrack *mct = new MergeClusterTrack(all_clustertrack.at(i));
      
      MergeSpaceCellSelection vertices; // save existing vertices
      all_mergeclustertrack.push_back(mct);
      used_clustertrack.push_back(all_clustertrack.at(i));
      vertices.push_back(all_clustertrack.at(i)->Get_FirstMSCell());
      vertices.push_back(all_clustertrack.at(i)->Get_LastMSCell());

      // start to add stuff ... 
      while (vertices.size()!=0){
	MergeSpaceCell *vertex = vertices.back();
	vertices.pop_back();

	// find the track inside MergeClusterTrack which contain this vertex
	ClusterTrack* old_cct = mct->GetClusterTrack(vertex);
	mct->SC_Hough(vertex->Get_Center());
	float theta1 = mct->Get_Theta();
	float phi1 = mct->Get_Phi();

	mct->SC_Hough(vertex->Get_Center(),-1,2);
	float theta1_m = mct->Get_Theta();
	float phi1_m = mct->Get_Phi();


	
	ClusterTrackSelection cts = ms_ct_map[vertex];

	//std::cout << theta1 << " " << phi1 << " " << cts.size() << std::endl;

	int flag_special = 0;
	for (int j=0;j!=cts.size();j++){
	  if (cts.at(j)->Get_allmcells().size()>3)
	    flag_special ++;
	}


	for (int j=0;j!=cts.size();j++){
	  ClusterTrack *cct = cts.at(j);
	  auto it1 = find(used_clustertrack.begin(),used_clustertrack.end(), cct);
	  if (it1==used_clustertrack.end()){
	    //Judge angle matching ... 
	    int flag = 0;
	    // write the main alg. ...
	    cct->SC_Hough(vertex->Get_Center());
	    float theta2 = cct->Get_Theta();
	    float phi2 = cct->Get_Phi();

	    cct->SC_Hough(vertex->Get_Center(),-1,2);
	    float theta2_m = cct->Get_Theta();
	    float phi2_m = cct->Get_Phi();
	    
	    float cut_angle;
	    if (flag_special > 2){
	      cut_angle = 10;
	    }else{
	      cut_angle = 25;
	    }
	    float cut_angle1 = 15;

	    float shift_angle = 100;

	    if ((fabs(theta1+theta2-3.1415926)<cut_angle/180.*3.1415926 // 5 degrees
	    	 && fabs(fabs(phi1-phi2)-3.1415926)<cut_angle/180.*3.1415926)
	    	){
	      flag = 1;
	    }

	    if ((fabs(theta1_m+theta2_m-3.1415926)<cut_angle1/180.*3.1415926 // 5 degrees
	    	 && fabs(fabs(phi1_m-phi2_m)-3.1415926)<cut_angle1/180.*3.1415926)
	    	){
	      flag = 1;
	    }

	    
	    // std::cout << theta1_m/3.1415926*180. << " " << theta2_m/3.1415926*180.
	    // 	      << " " << phi1_m/3.1415926*180. << " " << phi2_m/3.1415926*180. 
	    // 	      << std::endl;

	    

	    if (flag == 0 ){
	      int cross_num = cct->CrossNum(vertex->Get_Center(), theta1_m,phi1_m);
	      if ( cross_num == cct->Get_allmcells().size()){
	      	flag = 1;
	      }
	    }
	    
	    // if (flag==0){
	    //   theta2 = cct->Get_CTheta(vertex->Get_Center());
	    //   phi2 = cct->Get_CPhi(vertex->Get_Center());
	    //   if ((fabs(theta1+theta2-3.1415926)<cut_angle/180.*3.1415926 // 5 degrees
	    // 	   && fabs(fabs(phi1-phi2)-3.1415926)<cut_angle/180.*3.1415926)
	    // 	  ||(fabs(theta1-theta2-shift_angle)<cut_angle1/180.*3.1415926 // 5 degrees
	    // 	     && fabs(phi1-phi2)<cut_angle1/180.*3.1415926)
	    // 	  ){
	    // 	flag = 1;
	    //   }
	    // }

	    // if (flag == 0 ){
	    //   mct->SC_Hough(vertex->Get_Center(), 3 * units::cm);
	    //   cct->SC_Hough(vertex->Get_Center(), 3 * units::cm);
	      
	    //   theta1 = mct->Get_Theta();
	    //   phi1 = mct->Get_Phi();

	    //   theta2 = cct->Get_Theta();
	    //   phi2 = cct->Get_Phi();

	    //   if ((fabs(theta1+theta2-3.1415926)<cut_angle/180.*3.1415926 // 5 degrees
	    // 	   && fabs(fabs(phi1-phi2)-3.1415926)<cut_angle/180.*3.1415926)
	    // 	  ||(fabs(theta1-theta2-shift_angle)<cut_angle1/180.*3.1415926 // 5 degrees
	    // 	     && fabs(phi1-phi2)<cut_angle1/180.*3.1415926)
	    // 	  ){
	    // 	flag = 1;
	    //   }
	    // }

	    // std::cout <<  i << " " << find(all_clustertrack.begin(),all_clustertrack.end(),old_cct)-all_clustertrack.begin()
	    // 	      << " " << find(all_clustertrack.begin(),all_clustertrack.end(),cct)-all_clustertrack.begin() << " "
	    // 	      << cts.size() << " " 
	    // 	      << flag << " " << theta1/3.1415926*180. << " " << phi1/3.1415926*180. << " "
	    // 	      << theta2/3.1415926*180. << " " 
	    // 	      <<  phi2/3.1415926*180. << " "
	    // 	      << cct->Get_CTheta(vertex->Get_Center())/3.1415926*180. << " " 
	    // 	      << cct->Get_CPhi(vertex->Get_Center())/3.1415926*180. << " " 
	    // 	      << old_cct->Get_CTheta(vertex->Get_Center())/3.1415926*180. << " " 
	    // 	      << old_cct->Get_CPhi(vertex->Get_Center())/3.1415926*180. << " " 
	    // 	      << vertex->Get_Center().x/units::cm << " " 
	    // 	      << vertex->Get_Center().y/units::cm << " " 
	    // 	      << vertex->Get_Center().z/units::cm << " " 
	    // 	      << std::endl;

	    //flag = 0;

	    // if good, save it
	    if (flag==1){
	      mct->Add(cct,vertex);
	      used_clustertrack.push_back(cct);
	      if (cct->Get_FirstMSCell()!=vertex) vertices.push_back(cct->Get_FirstMSCell());
	      if (cct->Get_LastMSCell()!=vertex) vertices.push_back(cct->Get_LastMSCell());
	    }
	  }
	}
      }
    }
  }

  //Form map ...

  for (int i=0;i!=all_mergeclustertrack.size();i++){
    MergeClusterTrack *mct = all_mergeclustertrack.at(i);
    for (int j=0; j!=mct->Get_allmcells().size() ;j++){
      MergeSpaceCell *mcell = mct->Get_allmcells().at(j);
      
      if (mcells_mct_map.find(mcell) == mcells_mct_map.end()){
  	MergeClusterTrackSelection temp1;
  	temp1.push_back(mct);
  	mcells_mct_map[mcell] = temp1;
      }else{
  	mcells_mct_map[mcell].push_back(mct);
      }

    }
  }

}


void WireCell2dToy::ToyCrawler::CreateClusterTrack(MergeSpaceCellSelection& mcells){

  // form associations ...
  for (int i = 0; i!=mcells.size();i++){
    MergeSpaceCell *mcell1 = mcells.at(i);
    Point m1center = mcell1->Get_Center();
    double thickness = mcell1->thickness();
    
    MergeSpaceCellSelection mcell1_sel;

    for (int j=0;j!=mcells.size();j++){
      MergeSpaceCell *mcell2 = mcells.at(j);
      Point m2center = mcell2->Get_Center();

      if (fabs(fabs(m1center.x - m2center.x)-thickness) < 0.1*thickness ){
	if (mcell1->Overlap(*mcell2))
	    mcell1_sel.push_back(mcell2);
      }
      
      
    }
    MergeSpaceCellSelection mcell2_sel;
    mcells_save[mcell1] = mcell2_sel;
    mcells_map[mcell1] = mcell1_sel; // form structure   adjacent ones ...
    //    mcells_counter[mcell1] = 0;  // initiliaization
  }




  //std::cout << mcells_map.size() << " " << mcells_counter.size() << std::endl;

  //for (int qx = 0; qx!=14;qx++){
  //
  while(used_mcells.size()!=mcells_map.size()){
    //  std::cout << "Xin: " << qx << " " << used_mcells.size() << " " << mcells_map.size() << std::endl;
    //start to construct ClusterTrack ... first one 
    ClusterTrack *ctrack;
    // find the start point
    int flag1 = 0;
    for (int i=0;i!=mcells.size();i++){
      MergeSpaceCell *mcell = mcells.at(i);
      auto it1 = find(used_mcells.begin(),used_mcells.end(),mcell);
      auto it2 = find(end_mcells.begin(),end_mcells.end(),mcell);
            
      if (it1 == used_mcells.end()  // not insided the used ones
	  || (it2 != end_mcells.end()// inside the end ones
	      && mcells_map[mcell].size()> mcells_save[mcell].size() //not used all    
	      )) {
	
	// std::cout << it1 - used_mcells.end() << " " << 
	// it2 - end_mcells.end() << " " << 
	// mcells_map[mcell].size() << " " << mcells_save[mcell].size()
	// 	      <<std::endl;
	
	ctrack = new ClusterTrack(mcell);
	if (it2==end_mcells.end())
	  end_mcells.push_back(mcell);
	if (it1 == used_mcells.end())
	  used_mcells.push_back(mcell);
	//mcells_counter[mcell] ++;
	flag1 = 1;
	break;
      }
    }
    
    
    if (flag1 == 1){
      int flag = 0;
      while(flag==0){
	MergeSpaceCell *cur_cell = ctrack->Get_LastMSCell();
	int num_daughter = mcells_map[cur_cell].size();
	//      std::cout << "nd: " << num_daughter << std::endl;
	if (num_daughter == 0 ){ // no daugher, end
	  flag = 1;
	}else if (num_daughter == 1){
	  MergeSpaceCell *dcell = mcells_map[cur_cell].at(0);
	  auto it1 = find(ctrack->Get_allmcells().begin(),ctrack->Get_allmcells().end(),dcell);
	  if (it1 != ctrack->Get_allmcells().end()){
	    flag = 1; // already existed ... end it 
	  }else{
	    auto it2 = find(used_mcells.begin(),used_mcells.end(),dcell);
	    if (it2 == used_mcells.end()){
	      //not used yet? add it
	      
	      //Now need to judge angle ... 
	      	      
	      if (ctrack->AddMSCell(dcell)){
		used_mcells.push_back(dcell); // continue
		flag = 0;
	      }else{
		flag = 1;
	      }
	      	    

	    }else{
	      auto it3 = find(end_mcells.begin(),end_mcells.end(),dcell);
	      flag = 1;
	      if (it3 != end_mcells.end()){
		// already used, add it, and end it
		ctrack->AddMSCell(dcell);
	      }
	    }
	  }
	} // # of daughter  == 1
	else if (num_daughter == 2){
	  // two cases,    just one element
	  if (ctrack->Get_allmcells().size()==1){
	    for (int i=0; i!=num_daughter;i++){
	      MergeSpaceCell *dcell = mcells_map[cur_cell].at(i);
	      auto it1 = find(ctrack->Get_allmcells().begin(),ctrack->Get_allmcells().end(),dcell);
	      //std::cout << it1 - ctrack->Get_allmcells().end() << std::endl;
	      if (it1 == ctrack->Get_allmcells().end()){
		// not used yet
		auto it2 = find(used_mcells.begin(),used_mcells.end(),dcell);
		if (it2 == used_mcells.end()){
		  //not used yet? add it
		  if (ctrack->AddMSCell(dcell)){
		    used_mcells.push_back(dcell); // continue
		    flag = 0;
		    break;
		  }else{
		    flag = 1;
		  }
		}else{
		  auto it3 = find(end_mcells.begin(),end_mcells.end(),dcell);
		  flag = 1;
		  if (it3 != end_mcells.end()){
		    // already used, add it, and end it
		    auto it4 = find(mcells_save[cur_cell].begin(),mcells_save[cur_cell].end(),dcell);
		    if (it4 == mcells_save[cur_cell].end()){
		      if (ctrack->AddMSCell(dcell)){
			break;
		      }
		    }
		  }
		}
	      }else{
		flag = 1;
	      }
	    }
	  }else{
	    // more than one element
	    MergeSpaceCell *dcell1 = mcells_map[cur_cell].at(0);
	    MergeSpaceCell *dcell2 = mcells_map[cur_cell].at(1);
	    auto it1 = find(ctrack->Get_allmcells().begin(),ctrack->Get_allmcells().end(),dcell1);
	    auto it2 = find(ctrack->Get_allmcells().begin(),ctrack->Get_allmcells().end(),dcell2);
	    
	    if (it1 == ctrack->Get_allmcells().end() && it2 == ctrack->Get_allmcells().end()){
	      flag = 1; // end it
	    }else if (it1 == ctrack->Get_allmcells().end() && it2 != ctrack->Get_allmcells().end()){
	      // look at it1
	      auto it3 = find(used_mcells.begin(),used_mcells.end(),dcell1);
	      if (it3 == used_mcells.end()){
		//not used yet? add it
		if (ctrack->AddMSCell(dcell1)){
		  used_mcells.push_back(dcell1); // continue
		  flag = 0;
		}else{
		  flag = 1;
		}
	      }else{
		auto it4 = find(end_mcells.begin(),end_mcells.end(),dcell1);
		flag = 1;
		if (it4 != end_mcells.end()){
		  // already used, add it, and end it
		  //if (mcells_counter[dcell1]+1 < mcells_map[dcell1].size()){
		  ctrack->AddMSCell(dcell1);
		  //}
		}
	      }
	      
	    }else if (it1 != ctrack->Get_allmcells().end() && it2 == ctrack->Get_allmcells().end()){
	      // look at it2
	      auto it3 = find(used_mcells.begin(),used_mcells.end(),dcell2);
	      if (it3 == used_mcells.end()){
		//not used yet? add it
		if (ctrack->AddMSCell(dcell2)){
		  used_mcells.push_back(dcell2); // continue
		  flag = 0;
		}else{
		  flag = 1;
		}
	      }else{
		auto it4 = find(end_mcells.begin(),end_mcells.end(),dcell2);
		flag = 1;
		if (it4 != end_mcells.end()){
		  // already used, add it, and end it
		  //if (mcells_counter[dcell2]+1 < mcells_map[dcell2].size()){
		  ctrack->AddMSCell(dcell2);
		  //}
		}
	      }
	    }else if (it1 != ctrack->Get_allmcells().end() && it2 != ctrack->Get_allmcells().end()){
	      flag = 1;
	    }
	  }
	}
	
	if (num_daughter > 2){
	  // two cases,    just one element
	  if (ctrack->Get_allmcells().size()==1){
	    for (int i=0; i!=num_daughter;i++){
	      MergeSpaceCell *dcell = mcells_map[cur_cell].at(i);
	      auto it1 = find(ctrack->Get_allmcells().begin(),ctrack->Get_allmcells().end(),dcell);
	      
	      //std::cout << "Xin2: " <<it1 - ctrack->Get_allmcells().end() << std::endl;
	      
	      if (it1 == ctrack->Get_allmcells().end()){
		// not used yet
		auto it2 = find(used_mcells.begin(),used_mcells.end(),dcell);
		//std::cout << "Xin3: " << it2 - used_mcells.end() << std::endl;
		if (it2 == used_mcells.end()){
		  //not used yet? add it
		  if (ctrack->AddMSCell(dcell)){
		    used_mcells.push_back(dcell); // continue
		    flag = 0;
		    break;
		  }else{
		    flag = 1;
		  }
		}else{
		  auto it3 = find(end_mcells.begin(),end_mcells.end(),dcell);
		  flag = 1;
		  if (it3 != end_mcells.end()){
		    // already used, add it, and end it
		    //		  std::cout << mcells_counter[dcell] << " " << mcells_map[dcell].size() << std::endl;
		    auto it4 = find(mcells_save[cur_cell].begin(),mcells_save[cur_cell].end(),dcell);
		    if (it4 == mcells_save[cur_cell].end()){
		      //if (mcells_counter[dcell]+1 < mcells_map[dcell].size()){
		      if (ctrack->AddMSCell(dcell)){
			break;
		      }
		    }
		  }
		}
	      }else{
		flag = 1;
	      }
	      
	      //std::cout << i << " " << flag << std::endl;
	    }
	  }else{
	    flag = 1; // end
	  }
	}
      } // while(flag==0)
      
      // deal with the last element
      if (ctrack->Get_allmcells().size()>1){
	auto it = find(end_mcells.begin(),end_mcells.end(),ctrack->Get_LastMSCell());
	if (it == end_mcells.end())
	  end_mcells.push_back(ctrack->Get_LastMSCell());
      }
      
      if (ctrack->Get_allmcells().size()!=1){
	MergeSpaceCell *fcell = ctrack->Get_FirstMSCell();
	MergeSpaceCell *fcell_next = ctrack->Get_allmcells().at(1);
	
	MergeSpaceCell *lcell = ctrack->Get_LastMSCell();
	MergeSpaceCell *lcell_next = ctrack->Get_allmcells().at(ctrack->Get_allmcells().size()-2);
	
	mcells_save[fcell].push_back(fcell_next);
	mcells_save[lcell].push_back(lcell_next);
      }else{
	MergeSpaceCell *fcell = ctrack->Get_FirstMSCell();
	for (int qq = 0;qq!=mcells_map[fcell].size();qq++){
	  mcells_save[fcell].push_back(mcells_map[fcell].at(qq));
	}
      }
      
      
      all_clustertrack.push_back(ctrack);
    } 
  } // while loop
}

void WireCell2dToy::ToyCrawler::FormGraph(){
  ct_ms_map.clear();
  ms_ct_map.clear();

  for (int i=0;i!=all_clustertrack.size();i++){
    ClusterTrack *ctrack = all_clustertrack.at(i);
    MergeSpaceCell *FMSCell = ctrack->Get_FirstMSCell();
    MergeSpaceCell *LMSCell = ctrack->Get_LastMSCell();
    
    MergeSpaceCellSelection temp;
    temp.push_back(FMSCell);
    temp.push_back(LMSCell);
    ct_ms_map[ctrack] = temp;
    
    if (ms_ct_map.find(FMSCell) == ms_ct_map.end()){
      ClusterTrackSelection temp1;
      temp1.push_back(ctrack);
      ms_ct_map[FMSCell] = temp1;
    }else{
      ms_ct_map[FMSCell].push_back(ctrack);
    }
    
    if (ms_ct_map.find(LMSCell) == ms_ct_map.end()){
      ClusterTrackSelection temp1;
      temp1.push_back(ctrack);
      ms_ct_map[LMSCell] = temp1;
    }else{
      ms_ct_map[LMSCell].push_back(ctrack);
    }
        
  }

  // need to expand the distance cut a little bit ... 

}


WireCell2dToy::ToyCrawler::~ToyCrawler(){
  for (int i=0;i!=all_clustertrack.size();i++){
    delete all_clustertrack.at(i);
  }
  for (int i=0;i!=all_mergeclustertrack.size();i++){
    delete all_mergeclustertrack.at(i);
  }
  all_mergeclustertrack.clear();
  all_clustertrack.clear();
}
