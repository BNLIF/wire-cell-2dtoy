#ifndef WIRECELL_TOYCRAWLER_H
#define WIRECELL_TOYCRAWLER_H

#include "WCPData/ClusterTrack.h"
#include "WCPData/MergeClusterTrack.h"
#include "WCP2dToy/MergeToyTiling.h"

namespace WCP2dToy{
  
  class ToyCrawler {
  public:
    ToyCrawler(WCP::MergeSpaceCellSelection& mcells, int flag = 1, int flag1 = 1);
    ~ToyCrawler();

    void CreateClusterTrack(WCP::MergeSpaceCellSelection& mcells);
    void FormGraph();
    void MergeCTrack(int flag = 1);
    void FurtherExtendCTrack(int flag_qx = 0);
    void PurgeMergeCTrack();

    void CleanUpCTTrack(int flag=1);

    void PrepareTracking();

    void UpdateMap();

    WCP::ClusterTrackSelection& Get_allCT(){return all_clustertrack;};
    WCP::MergeClusterTrackSelection& Get_allMCT(){return all_mergeclustertrack;};
    WCP::MSpaceCellClusterMap& Get_ms_ct_map(){return ms_ct_map;};
    WCP::MergeSpaceCellCounter& Get_mcells_counter(){return mcells_counter;};
    WCP::MergeSpaceCellMap& Get_mcells_map(){return mcells_map;};

    WCP::MergeSpaceCell* GetClosestMSC(WCP::Point p, WCP::MergeSpaceCellSelection& cells);

   

  protected:
    WCP::ClusterTrackSelection all_clustertrack;
    //WCP::ClusterTrackSelection used_clustertrack; // hold the used tracks

    WCP::MergeClusterTrackSelection all_mergeclustertrack;
    WCP::MSC_MCT_Map  mcells_mct_map;


    WCP::MergeSpaceCellSelection used_mcells; // save used mergespacecell
    WCP::MergeSpaceCellSelection end_mcells; // save all the merge spacecell at the beginning/end of a clustertrack

    WCP::MergeSpaceCellMap mcells_map;
    WCP::MergeSpaceCellMap mcells_save;

    
    WCP::MergeSpaceCellCounter mcells_counter;
    WCP::ClusterMSpaceCellMap ct_ms_map;
    WCP::MSpaceCellClusterMap ms_ct_map;

  };
  
}


#endif
