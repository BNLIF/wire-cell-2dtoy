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
    ToySignalSimuFDS(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, int bins_per_frame1 = 9600, int nframes_total = -1,float time_offset_uv = 0, float time_offset_uw = 0, int flag_random = 1);
    ~ToySignalSimuFDS();

    virtual int size() const;
    virtual int jump(int frame_number);
    //fixed it ...
    

    void Save();
    
  private:
    WireCell::FrameDataSource& fds;
    const WireCell::GeomDataSource& gds;
    int  max_frames;

    float time_offset_uv, time_offset_uw;
    int flag_random;
    
    int nwire_u, nwire_v, nwire_w;

    /* TH1F **hu; */
    /* TH1F **hv; */
    /* TH1F **hw; */

    TH1F *hu;
    TH1F *hv;
    TH1F *hw;
    


    TGraph *gu, *gv, *gw;

    TH1F *hur, *hvr, *hwr;
    
  };
}

#endif
