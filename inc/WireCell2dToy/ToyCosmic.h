#ifndef WIRECELL_TOYCOSMIC_H
#define WIRECELL_TOYCOSMIC_H

#include "WireCell2dToy/ToyTracking.h"

namespace WireCell2dToy{
  class ToyCosmic{
  public:
    ToyCosmic(ToyTrackingSelection& trackings);
    ~ToyCosmic();
  protected:
    ToyTrackingSelection& trackings;
  };
}

#endif
