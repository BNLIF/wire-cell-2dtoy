#ifndef WIRECELL_TOYTRACKING_H
#define WIRECELL_TOYTRACKING_H

#include "WireCellData/WCTrack.h"
#include "WireCellData/WCVertex.h"
#include "WireCellData/WCShower.h"
#include "WireCell2dToy/ToyCrawler.h"

namespace WireCell2dToy{
  typedef std::map<WireCell::WCTrack*, WireCell::WCVertexSelection> WCT_WCV_Map;
  typedef std::map<WireCell::WCVertex*, WireCell::WCTrackSelection> WCV_WCT_Map;

  class ToyTracking{
  public:
    ToyTracking(ToyCrawler& toycrawler);
    ~ToyTracking();
    
    WireCell::WCTrackSelection& get_tracks(){return tracks;};

    WireCell::WCTrackSelection& get_good_tracks(){return good_tracks;};
    WireCell::WCTrackSelection& get_bad_tracks(){return bad_tracks;};
    WireCell::WCTrackSelection& get_parallel_tracks(){return parallel_tracks;};
    WireCell::WCTrackSelection& get_short_tracks(){return short_tracks;};
    
    WireCell::WCVertexSelection& get_vertices(){return vertices;};
    /* WireCell::WCVertexSelection& get_good_vertices(){return good_vertices;}; */
    /* WireCell::WCVertexSelection& get_bad_vertices(){return bad_vertices;}; */
    
    void CreateVertices(ToyCrawler& toycrawler);
    void RemoveSame();
    void MergeVertices(int flag  = 1);
    void BreakTracks();
    void OrganizeTracks(int flag = 1);
    void Associate();
    void CleanUpVertex();
    void Crawl();
   
    
    bool IsThisShower(WireCell2dToy::ToyCrawler& toycrawler);
    

    void CheckVertices(WireCell2dToy::ToyCrawler& toycrawler);

    bool ExamineVertex(WireCell::WCVertex* vertex, WireCell2dToy::ToyCrawler& toycrawler);

    void update_maps();
    void update_maps1();
    void fine_tracking(int flag = 0);
    void cleanup_bad_tracks();

    void deal_wiggle_tracks();
    bool grow_track_fill_gap(WireCell2dToy::ToyCrawler& toycrawler);

    void form_parallel_tiny_tracks(WireCell2dToy::ToyCrawler& toycrawler);
    void parallel_tracking(WireCell::WCVertex *vertex, WireCell::MergeSpaceCellSelection &mcells, WireCell2dToy::ToyCrawler& toycrawler);

    bool track_shower_reco(WireCell2dToy::ToyCrawler& toycrawler);
    void single_shower_reco(WireCell2dToy::ToyCrawler& toycrawler);
    

  protected: 
    WireCell::WCTrackSelection tracks;
    WireCell::WCVertexSelection vertices;
    WireCell::WCShowerSelection showers;
    
    /* WireCell::WCVertexSelection good_vertices; */
    /* WireCell::WCVertexSelection bad_vertices; */
    
    WireCell::WCTrackSelection good_tracks;
    WireCell::WCTrackSelection bad_tracks;
    WireCell::WCTrackCounter type3_tracks;

    
    WireCell::WCTrackSelection parallel_tracks;
    WireCell::WCTrackSelection short_tracks;
    
    //    WireCell::MergeSpaceCellSelection orig_mcells;
    WireCell::MergeSpaceCellSelection new_mcells;
    WireCell::MergeSpaceCellMap1 new_mcells_map; 

    /* WireCell::MSC_WCV_Map msc_wcv_map; */
    /* WireCell::MCT_WCT_Map mct_wct_map; */

    WCT_WCV_Map wct_wcv_map;
    WCV_WCT_Map wcv_wct_map; 
    
  };
}

#endif
