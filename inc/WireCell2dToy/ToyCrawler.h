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

    WireCell::MergeSpaceCellSelection used_mcells; // save used mergespacecell
    WireCell::MergeSpaceCellSelection end_mcells; // save all the merge spacecell at the beginning/end of a clustertrack

    WireCell::MergeSpaceCellMap mcells_map;
    WireCell::MergeSpaceCellMap mcells_save;

    WireCell::MergeSpaceCellCounter mcells_counter;
  };
  
}


#endif
