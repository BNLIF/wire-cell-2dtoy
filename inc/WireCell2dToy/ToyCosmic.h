#ifndef WIRECELL_TOYCOSMIC_H
#define WIRECELL_TOYCOSMIC_H

#include "WireCell2dToy/ToyTracking.h"

namespace WireCell2dToy{
  class ToyCosmic{
  public:
    ToyCosmic(ToyTrackingSelection& trackings);
    ~ToyCosmic();
    
    std::vector<ToyTrackingSelection>& get_raw_candidates(){return cosmic_candidates;};

    bool IsConnected(ToyTracking *tracking1, ToyTracking *tracking2);
    bool IsConnected(WireCell::MergeSpaceCell *mcell1, WireCell::MergeSpaceCell *mcell2, float dis_cut);
    
  protected:
    ToyTrackingSelection& trackings;
    std::vector<ToyTrackingSelection> cosmic_candidates;
  };
}

#endif
