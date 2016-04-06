#ifndef WIRECELL2dToy_DATASIGNALGAUSROI_H
#define WIRECELL2dToy_DATASIGNALGAUSROI_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"
#include "WireCell2dToy/DataSignalWien_ROI.h"
#include "TH2I.h"

namespace WireCell2dToy {
  class DataSignalGausROIFDS : public WireCell::FrameDataSource
  {
  public:
    DataSignalGausROIFDS(DataSignalWienROIFDS& fds, int nframes_total);
    ~DataSignalGausROIFDS();

    virtual int size() const;
    virtual int jump(int frame_number);

  private:
    int  max_frames;
    DataSignalWienROIFDS& fds;
    
  };
}

#endif
