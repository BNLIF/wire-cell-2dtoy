#ifndef WIRECELL2dToy_TOYSIGNALSIMU_H
#define WIRECELL2dToy_TOYSIGNALSIMU_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"
#include "TH1F.h"
#include "TGraph.h"

namespace WireCell2dToy {
  class ToySignalSimuFDS : public WireCell::FrameDataSource
  {
  public:
    ToySignalSimuFDS(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, int bins_per_frame = 9600, int nframes_total = -1);
    ~ToySignalSimuFDS();

    virtual int size() const;
    virtual int jump(int frame_number);

    void Save();
    
  private:
    WireCell::FrameDataSource& fds;
    const WireCell::GeomDataSource& gds;
    int bins_per_frame, max_frames;

    int nwire_u, nwire_v, nwire_w;

    TH1F **hu;
    TH1F **hv;
    TH1F **hw;

    TGraph *gu, *gv, *gw;

    TH1F *hur, *hvr, *hwr;
    
  };
}

#endif
