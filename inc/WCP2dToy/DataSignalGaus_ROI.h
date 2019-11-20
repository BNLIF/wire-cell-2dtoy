#ifndef WIRECELL2dToy_DATASIGNALGAUSROI_H
#define WIRECELL2dToy_DATASIGNALGAUSROI_H

#include "WCPNav/FrameDataSource.h"
#include "WCPNav/GeomDataSource.h"
#include "WCP2dToy/DataSignalWien_ROI.h"
#include "TH2I.h"

namespace WCP2dToy {
  class DataSignalGausROIFDS : public WCP::FrameDataSource
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
