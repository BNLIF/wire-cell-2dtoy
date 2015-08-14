#ifndef WIRECELL_TOYTRACKING_H
#define WIRECELL_TOYTRACKING_H

#include "WireCellData/WCTrack.h"
#include "WireCellData/WCVertex.h"
#include "WireCell2dToy/ToyCrawler.h"

namespace WireCell2dToy{
  typedef std::map<WireCell::WCTrack*, WireCell::WCVertexSelection> WCT_WCV_Map;
  typedef std::map<WireCell::WCVertex*, WireCell::WCTrackSelection> WCV_WCT_Map;

  class ToyTracking{
  public:
    ToyTracking(ToyCrawler& toycrawler);
    ~ToyTracking();
    
    WireCell::WCTrackSelection& get_tracks(){return tracks;};
    WireCell::WCVertexSelection& get_vertices(){return vertices;};
    
    void CreateVertices(ToyCrawler& toycrawler);
    void RemoveSame();
    void MergeVertices();
    void BreakTracks();
    void OrganizeTracks();


  protected: 
    WireCell::WCTrackSelection tracks;
    WireCell::WCVertexSelection vertices;
    
    WireCell::MSC_WCV_Map msc_wcv_map;
    WireCell::MCT_WCT_Map mct_wct_map;

    WCT_WCV_Map wct_wcv_map;
    WCV_WCT_Map wcv_wct_map;
    
  };
}

#endif
