#ifndef WIRECELL_WCCOSMIC_H
#define WIRECELL_WCCOSMIC_H

#include "WireCell2dToy/ToyTracking.h"

namespace WireCell2dToy{
  class WCCosmic{
  public:
    WCCosmic(ToyTrackingSelection& toytrackings);
    ~WCCosmic();

    void Add(WCCosmic *cosmic);
    ToyTrackingSelection& get_trackings(){return toytrackings;};
    WireCell::MergeSpaceCellSelection& get_mcells(){return mcells;};
    float get_theta(){return theta;};
    float get_phi(){return phi;};
    float cal_costh(WCCosmic *cosmic);
    float cal_dist(WCCosmic *cosmic);
    WireCell::Point get_center(){return center;};
    WireCell::PointVector& get_points(){return points;};
    float cal_pos(WireCell::MergeSpaceCell *mcell1);
    float cal_pos(WireCell::SpaceCell *cell);
    void Sort();
    void fill_points();
    void judge_cosmic();

    bool IsCosmic(){return cosmic_flag;};
    bool IsNearBy(ToyTracking *toytracking);

  protected:
    ToyTrackingSelection& toytrackings;
    WireCell::MergeSpaceCellSelection mcells;
    WireCell::ClusterTrack *ct;
    WireCell::Point center;
    WireCell::PointVector points;
    float theta;
    float phi;
    bool cosmic_flag;

  };
  
  struct MSC_Struct
  {
    float key;
    WireCell::MergeSpaceCell *mcell;
  MSC_Struct(float key, WireCell::MergeSpaceCell *mcell) : key(key), mcell(mcell) {}
  };

  struct less_than_key{
    inline bool operator() (const MSC_Struct& s1, const MSC_Struct& s2)
    {
      return (s1.key < s2.key);
    }
  };

  typedef std::vector<WCCosmic*> WCCosmicSelection;
}

#endif
