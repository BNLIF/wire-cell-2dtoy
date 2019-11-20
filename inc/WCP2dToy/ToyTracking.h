#ifndef WIRECELL_TOYTRACKING_H
#define WIRECELL_TOYTRACKING_H

#include "WCPData/WCTrack.h"
#include "WCPData/WCVertex.h"
#include "WCPData/WCShower.h"
#include "WCP2dToy/ToyCrawler.h"


namespace WCP2dToy{
  typedef std::map<WCP::WCTrack*, WCP::WCVertexSelection> WCT_WCV_Map;
  typedef std::map<WCP::WCVertex*, WCP::WCTrackSelection> WCV_WCT_Map;

  class ToyTracking{
  public:
    ToyTracking(ToyCrawler& toycrawler, int tracking_type = 0);
    ~ToyTracking();
    
    WCP::WCTrackSelection& get_tracks(){return tracks;};

    WCP::WCTrackSelection& get_good_tracks(){return good_tracks;};
    WCP::WCTrackSelection& get_bad_tracks(){return bad_tracks;};
    WCP::WCTrackSelection& get_parallel_tracks(){return parallel_tracks;};
    WCP::WCTrackSelection& get_short_tracks(){return short_tracks;};
    
    WCP::WCVertexSelection& get_vertices(){return vertices;};
    WCP::WCVertexSelection& get_good_vertices(){return good_vertices;};  
    WCP::WCVertexSelection& get_bad_vertices(){return bad_vertices;}; 
    
    WCP::WCShowerSelection& get_showers(){return showers;};

    void IterateMergeTracks(WCP::MergeSpaceCellMap& mcells_map);
    void MergeTracks(WCP::MergeSpaceCellMap& mcells_map);
    void MergeTracks_no_shared_vertex(WCP::MergeSpaceCellMap& mcells_map, int type = 0);

    void create_new_tracks_from_leftover(ToyCrawler& toycrawler);

    void CreateVertices(ToyCrawler& toycrawler);
    void RemoveSame();
    void RemoveSameTrack();

    void MergeVertices(int flag  = 1);
    void BreakTracks();
    void OrganizeTracks(int flag = 1);
    void Associate();
    void CleanUpVertex();
    void Crawl();
   
    bool IsContained();
    bool IsThisShower(WCP2dToy::ToyCrawler& toycrawler);
    

    void CheckVertices(WCP2dToy::ToyCrawler& toycrawler);

    bool ExamineVertex(WCP::WCVertex* vertex, WCP2dToy::ToyCrawler& toycrawler);

    void update_maps(int flag = 0);
    void update_maps1();
    void fine_tracking(int flag = 0);
    void cleanup_bad_tracks();

    void deal_wiggle_tracks();
    bool grow_track_fill_gap(WCP2dToy::ToyCrawler& toycrawler);

    void form_parallel_tiny_tracks(WCP2dToy::ToyCrawler& toycrawler);
    void form_parallel_tiny_tracks_wovertex(WCP2dToy::ToyCrawler& toycrawler);
    

    void parallel_tracking(WCP::WCVertex *vertex, WCP::MergeSpaceCellSelection &mcells, WCP2dToy::ToyCrawler& toycrawler);

    bool track_shower_reco(WCP2dToy::ToyCrawler& toycrawler);
    void Cleanup_showers();

    void single_shower_reco(WCP2dToy::ToyCrawler& toycrawler);
    
    void cosmic_finder_all(WCP2dToy::ToyCrawler& toycrawler);
    void cosmic_finder_part(WCP2dToy::ToyCrawler& toycrawler);
    std::map<WCP::WCTrack*, std::vector<float>>& get_gt_angle_map(){return tracks_angle_map;};
    std::map<WCP::WCTrack*, std::vector<float>>& get_gt_pos_map(){return tracks_pos_map;};
    void fill_maps();

    

  protected: 
    std::map<WCP::WCTrack*, std::vector<float>> tracks_angle_map;
    std::map<WCP::WCTrack*, std::vector<float>> tracks_pos_map;

    WCP::WCTrackSelection tracks;
    WCP::WCVertexSelection vertices;
    WCP::WCShowerSelection showers;
    
    //WCP::WCVertexSelection wiggle_vertices;
    
    WCP::WCVertexSelection good_vertices;  
    WCP::WCVertexSelection bad_vertices;  
    
    WCP::WCTrackSelection good_tracks;
    WCP::WCTrackSelection bad_tracks;
    WCP::WCTrackCounter type3_tracks;

    
    WCP::WCTrackSelection parallel_tracks;
    WCP::WCTrackSelection short_tracks;
    
    //    WCP::MergeSpaceCellSelection orig_mcells;
    WCP::MergeSpaceCellSelection new_mcells;
    WCP::MergeSpaceCellMap1 new_mcells_map;  //map the new stuff ... 

    /* WCP::MSC_WCV_Map msc_wcv_map; */
    /* WCP::MCT_WCT_Map mct_wct_map; */

    WCT_WCV_Map wct_wcv_map;
    WCV_WCT_Map wcv_wct_map; 
    
  };
  
  typedef std::vector<ToyTracking*> ToyTrackingSelection;

}

#endif
