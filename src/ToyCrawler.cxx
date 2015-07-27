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
	mcell1_sel.push_back(mcell2);
      }
      
      
    }
    mcells_map[mcell1] = mcell1_sel;
  }
  
}

WireCell2dToy::ToyCrawler::~ToyCrawler(){
  for (int i=0;i!=all_clustertrack.size();i++){
    delete all_clustertrack.at(i);
  }
}
