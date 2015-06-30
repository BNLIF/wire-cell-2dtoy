#ifndef WIRECELL2dToy_TOYSIGNALSIMU_H
#define WIRECELL2dToy_TOYSIGNALSIMU_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"

namespace WireCell2dToy {
  class ToySignalSimuFDS : public WireCell::FrameDataSource
  {
  public:
    ToySignalSimuFDS(const WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, int bins_per_frame = 9600, int nframes_total = -1);
  private:
    const WireCell::FrameDataSource& fds;
    const WireCell::GeomDataSource& gds;
    int bins_per_frame, max_frames;
    
    
  };
}

#endif
