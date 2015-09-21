#ifndef WIRECELL_TOYCRAWLER_H
#define WIRECELL_TOYCRAWLER_H

#include "WireCellData/ClusterTrack.h"
#include "WireCellData/MergeClusterTrack.h"
#include "WireCell2dToy/MergeToyTiling.h"

namespace WireCell2dToy{
  
  class ToyCrawler {
  public:
    ToyCrawler(WireCell::MergeSpaceCellSelection& mcells, int flag = 1, int flag1 = 1);
    ~ToyCrawler();

    void CreateClusterTrack(WireCell::MergeSpaceCellSelection& mcells);
    void FormGraph();
    void MergeCTrack(int flag = 1);
    void FurtherExtendCTrack(int flag_qx = 0);
    void PurgeMergeCTrack();

    void CleanUpCTTrack(int flag=1);

    void PrepareTracking();

    void UpdateMap();

    WireCell::ClusterTrackSelection& Get_allCT(){return all_clustertrack;};
    WireCell::MergeClusterTrackSelection& Get_allMCT(){return all_mergeclustertrack;};
    WireCell::MSpaceCellClusterMap& Get_ms_ct_map(){return ms_ct_map;};
    WireCell::MergeSpaceCellCounter& Get_mcells_counter(){return mcells_counter;};
    WireCell::MergeSpaceCellMap& Get_mcells_map(){return mcells_map;};

    WireCell::MergeSpaceCell* GetClosestMSC(WireCell::Point p, WireCell::MergeSpaceCellSelection& cells);

   

  protected:
    WireCell::ClusterTrackSelection all_clustertrack;
    //WireCell::ClusterTrackSelection used_clustertrack; // hold the used tracks

    WireCell::MergeClusterTrackSelection all_mergeclustertrack;
    WireCell::MSC_MCT_Map  mcells_mct_map;


    WireCell::MergeSpaceCellSelection used_mcells; // save used mergespacecell
    WireCell::MergeSpaceCellSelection end_mcells; // save all the merge spacecell at the beginning/end of a clustertrack

    WireCell::MergeSpaceCellMap mcells_map;
    WireCell::MergeSpaceCellMap mcells_save;

    
    WireCell::MergeSpaceCellCounter mcells_counter;
    WireCell::ClusterMSpaceCellMap ct_ms_map;
    WireCell::MSpaceCellClusterMap ms_ct_map;

  };
  
}


#endif
