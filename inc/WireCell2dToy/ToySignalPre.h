#ifndef WIRECELL2dToy_TOYSIGNALPRE_H
#define WIRECELL2dToy_TOYSIGNALPRE_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"
#include "TH1F.h"

namespace WireCell2dToy {
  class ToySignalPreFDS : public WireCell::FrameDataSource
  {
  public:
    ToySignalPreFDS(WireCell::FrameDataSource& fds,  const WireCell::GeomDataSource& gds, int bins_per_frame = 9600, int nframes_total = -1);
    ~ToySignalPreFDS();

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

    TH1F *hfilter_u;
    TH1F *hfilter_v;
    TH1F *hfilter_w;
    
    TH1F *hfilter_time_u;
    TH1F *hfilter_time_v;
    TH1F *hfilter_time_w;
  };

}

#endif
