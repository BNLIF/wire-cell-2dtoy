#ifndef WIRECELL_WCCOSMIC_H
#define WIRECELL_WCCOSMIC_H

#include "WCP2dToy/ToyTracking.h"

namespace WCP2dToy{
  class WCCosmic{
  public:
    WCCosmic(ToyTrackingSelection& toytrackings);
    ~WCCosmic();

    void Add(WCCosmic *cosmic);
    ToyTrackingSelection& get_trackings(){return toytrackings;};
    WCP::MergeSpaceCellSelection& get_mcells(){return mcells;};
    float get_theta(){return theta;};
    float get_phi(){return phi;};
    float cal_costh(WCCosmic *cosmic);
    float cal_dist(WCCosmic *cosmic);
    WCP::Point get_center(){return center;};
    WCP::PointVector& get_points(){return points;};
    float cal_pos(WCP::MergeSpaceCell *mcell1);
    float cal_pos(WCP::SpaceCell *cell);
    void Sort();
    void fill_points();
    void judge_cosmic();

    bool IsCosmic(){return cosmic_flag;};
    bool IsNearBy(ToyTracking *toytracking);

  protected:
    ToyTrackingSelection& toytrackings;
    WCP::MergeSpaceCellSelection mcells;
    WCP::ClusterTrack *ct;
    WCP::Point center;
    WCP::PointVector points;
    float theta;
    float phi;
    bool cosmic_flag;

  };
  
  struct MSC_Struct
  {
    float key;
    WCP::MergeSpaceCell *mcell;
  MSC_Struct(float key, WCP::MergeSpaceCell *mcell) : key(key), mcell(mcell) {}
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
