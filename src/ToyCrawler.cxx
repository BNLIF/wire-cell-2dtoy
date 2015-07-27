#include "WireCell2dToy/ToyCrawler.h"

using namespace WireCell;

WireCell2dToy::ToyCrawler::ToyCrawler(MergeSpaceCellSelection& mcells){
  // form associations ...
  for (int i = 0; i!=mcells.size();i++){
    MergeSpaceCell *mcell1 = mcells.at(i);
    Point m1center = mcell1->Get_Center();
    double thickness = mcell1->thickness();
    
    MergeSpaceCellSelection mcell1_sel;
    for (int j=0;j!=mcells.size();j++){
      MergeSpaceCell *mcell2 = mcells.at(j);
      Point m2center = mcell2->Get_Center();

      if (fabs(m1center.x - m2center.x)==thickness){
	if (mcell1->Overlap(*mcell2))
	    mcell1_sel.push_back(mcell2);
      }
      
      
    }
    mcells_map[mcell1] = mcell1_sel;
    mcells_counter[mcell1] = 0;
  }
  
  //start to construct ClusterTrack ... 
  ClusterTrack *ctrack;

  // find the start point
  for (int i=0;i!=mcells.size();i++){
    MergeSpaceCell *mcell = mcells.at(i);
    auto it1 = find(used_mcells.begin(),used_mcells.end(),mcell);
    auto it2 = find(end_mcells.begin(),end_mcells.end(),mcell);
    
    if (it1 == used_mcells.end()  // not insided the used ones
	|| (it2 != end_mcells.end()// inside the end ones
	    && mcells_map[mcell].size()!= mcells_counter[mcell] //not used all    
	    )) {
      ctrack = new ClusterTrack(mcell);
      
      if (it2==end_mcells.end())
	end_mcells.push_back(mcell);


      used_mcells.push_back(mcell);
      mcells_counter[mcell] ++;
      
      break;
    }
  }

  while(mcells_counter[ctrack->Get_LastMSCell()]+1 != mcells_map[ctrack->Get_LastMSCell()].size()){
    // start to look at the next ones
    for (int i=0;i!=mcells_map[ctrack->Get_LastMSCell()].size();i++){
      MergeSpaceCell *mcell = mcells_map[ctrack->Get_LastMSCell()].at(i);
      auto it1 = find(used_mcells.begin(),used_mcells.end(),mcell);
      if (it1 == used_mcells.end()){
	ctrack->AddMSCell(mcell);
	used_mcells.push_back(mcell);
	break;
      }
    }
  }
  
  //add the last element
  auto it = find(end_mcells.begin(),end_mcells.end(),ctrack->Get_LastMSCell());
  if (it == end_mcells.end())
    end_mcells.push_back(ctrack->Get_LastMSCell());
  mcells_counter[ctrack->Get_LastMSCell()] ++;

  all_clustertrack.push_back(ctrack);
  

}

WireCell2dToy::ToyCrawler::~ToyCrawler(){
  for (int i=0;i!=all_clustertrack.size();i++){
    delete all_clustertrack.at(i);
  }
}
