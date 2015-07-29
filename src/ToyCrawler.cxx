#include "WireCell2dToy/ToyCrawler.h"

using namespace WireCell;

WireCell2dToy::ToyCrawler::ToyCrawler(MergeSpaceCellSelection& mcells){

  CreateClusterTrack(mcells);
  FormGraph(); 
  //std::cout << "Merge Clusters " << std::endl; 
  MergeCTrack();
   
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
	old_cct->SC_Hough(vertex->Get_Center());
	float theta1 = old_cct->Get_Theta();
	float phi1 = old_cct->Get_Phi();

	//
	

	ClusterTrackSelection cts = ms_ct_map[vertex];

	//std::cout << theta1 << " " << phi1 << " " << cts.size() << std::endl;

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
	    
	    float cut_angle = 15;
	    float cut_angle1 = 5;


	    if ((fabs(theta1+theta2-3.1415926)<cut_angle/180.*3.1415926 // 5 degrees
		 && fabs(fabs(phi1-phi2)-3.1415926)<cut_angle/180.*3.1415926)
		||(fabs(theta1-theta2-100)<cut_angle1/180.*3.1415926 // 5 degrees
		 && fabs(phi1-phi2)<cut_angle1/180.*3.1415926)
		){
	      flag = 1;
	    }

	    if (flag == 0 ){
	      
	    }

	    if (flag == 0 ){
	    
	    }

	    std::cout <<  i << " " << find(all_clustertrack.begin(),all_clustertrack.end(),old_cct)-all_clustertrack.begin()
		      << " " << find(all_clustertrack.begin(),all_clustertrack.end(),cct)-all_clustertrack.begin() << " "
		      << cts.size() << " " 
		      << flag << " " << (theta1 + theta2-3.1415926)/3.1415926*180. << " " << (fabs(phi1- phi2)-3.1415926)/3.1415926*180. << std::endl;


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
    mcells_map[mcell1] = mcell1_sel; // form structure 
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
	    ctrack->AddMSCell(dcell);
	    used_mcells.push_back(dcell); // continue
	    flag = 0;
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
		ctrack->AddMSCell(dcell);
		used_mcells.push_back(dcell); // continue
		flag = 0;
		break;
	      }else{
		auto it3 = find(end_mcells.begin(),end_mcells.end(),dcell);
		flag = 1;
		if (it3 != end_mcells.end()){
		  // already used, add it, and end it
		  auto it4 = find(mcells_save[cur_cell].begin(),mcells_save[cur_cell].end(),dcell);
		  if (it4 == mcells_save[cur_cell].end()){
		    ctrack->AddMSCell(dcell);
		    break;
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
	      ctrack->AddMSCell(dcell1);
	      used_mcells.push_back(dcell1); // continue
	      flag = 0;
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
	      ctrack->AddMSCell(dcell2);
	      used_mcells.push_back(dcell2); // continue
	      flag = 0;
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
		ctrack->AddMSCell(dcell);
		used_mcells.push_back(dcell); // continue
		flag = 0;
		break;
	      }else{
		auto it3 = find(end_mcells.begin(),end_mcells.end(),dcell);
		flag = 1;
		if (it3 != end_mcells.end()){
		  // already used, add it, and end it

		  //		  std::cout << mcells_counter[dcell] << " " << mcells_map[dcell].size() << std::endl;
		  auto it4 = find(mcells_save[cur_cell].begin(),mcells_save[cur_cell].end(),dcell);
		  if (it4 == mcells_save[cur_cell].end()){
		  //if (mcells_counter[dcell]+1 < mcells_map[dcell].size()){
		    ctrack->AddMSCell(dcell);
		    break;
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


  }
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
}


WireCell2dToy::ToyCrawler::~ToyCrawler(){
  for (int i=0;i!=all_clustertrack.size();i++){
    delete all_clustertrack.at(i);
  }
}
