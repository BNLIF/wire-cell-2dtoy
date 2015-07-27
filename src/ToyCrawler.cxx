#include "WireCell2dToy/ToyCrawler.h"

using namespace WireCell;

WireCell2dToy::ToyCrawler::ToyCrawler(MergeSpaceCellSelection& mcells){
}

WireCell2dToy::ToyCrawler::~ToyCrawler(){
  for (int i=0;i!=all_clustertrack.size();i++){
    delete all_clustertrack.at(i);
  }
}
