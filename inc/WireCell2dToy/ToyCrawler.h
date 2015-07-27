#ifndef WIRECELL_TOYCRAWLER_H
#define WIRECELL_TOYCRAWLER_H

#include "WireCellData/ClusterTrack.h"

namespace WireCell2dToy{
  
  class ToyCrawler {
  public:
    ToyCrawler(WireCell::MergeSpaceCellSelection& mcells);
    ~ToyCrawler();

    WireCell::ClusterTrackSelection& Get_allCT(){return all_clustertrack;};
    
  protected:
    WireCell::ClusterTrackSelection all_clustertrack;
    
  };
  
}


#endif
