#ifndef WIRECELL2dToy_TOYSIGNALGAUS_H
#define WIRECELL2dToy_TOYSIGNALGAUS_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"
#include "TH1F.h"
#include "TGraph.h"

namespace WireCell2dToy {
  class ToySignalGausFDS : public WireCell::FrameDataSource
  {
  public:
    ToySignalGausFDS(WireCell::FrameDataSource& fds,  const WireCell::GeomDataSource& gds, int bins_per_frame1 = 9600, int nframes_total = -1);
    ~ToySignalGausFDS();

    virtual int size() const;
    virtual int jump(int frame_number);

    void Save();

  private:
    WireCell::FrameDataSource& fds;
    const WireCell::GeomDataSource& gds;
    int max_frames;
    
    int nwire_u, nwire_v, nwire_w;

    TH1F **hu;
    TH1F **hv;
    TH1F **hw;

    TH1F *hfilter_u;
    TH1F *hfilter_v;
    TH1F *hfilter_w;

    
  };

}

#endif
