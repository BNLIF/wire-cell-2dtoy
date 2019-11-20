#ifndef WIRECELL2dToy_TOYSIGNALPRE_H
#define WIRECELL2dToy_TOYSIGNALPRE_H

#include "WCPNav/FrameDataSource.h"
#include "WCPNav/GeomDataSource.h"
#include "TH1F.h"
#include "TGraph.h"

namespace WCP2dToy {
  class ToySignalPreFDS : public WCP::FrameDataSource
  {
  public:
    ToySignalPreFDS(WCP::FrameDataSource& fds,  const WCP::GeomDataSource& gds, int bins_per_frame1 = 9600, int nframes_total = -1);
    ~ToySignalPreFDS();

    virtual int size() const;
    virtual int jump(int frame_number);

    void Save();

  private:
    WCP::FrameDataSource& fds;
    const WCP::GeomDataSource& gds;
    int max_frames;
    
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

    TGraph *gu;
    TGraph *gv;
    TGraph *gw;
  };

}

#endif
