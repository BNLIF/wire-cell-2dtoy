#ifndef WIRECELL_WCCOSMIC_H
#define WIRECELL_WCCOSMIC_H

#include "WireCell2dToy/ToyTracking.h"

namespace WireCell2dToy{
  class WCCosmic{
  public:
    WCCosmic(ToyTrackingSelection& toytrackings);
    ~WCCosmic();
  protected:
    ToyTrackingSelection& toytrackings;
  };
  
  typedef std::vector<WCCosmic*> WCCosmicSelection;
}

#endif
