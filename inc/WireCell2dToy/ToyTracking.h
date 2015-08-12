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
  protected: 
    WCT_WCV_Map wct_wcv_map;
    WCV_WCT_Map wcv_wct_map;
  };
}

#endif
