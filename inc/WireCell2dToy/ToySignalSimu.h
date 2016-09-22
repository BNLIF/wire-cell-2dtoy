#ifndef WIRECELL2dToy_TOYSIGNALSIMU_H
#define WIRECELL2dToy_TOYSIGNALSIMU_H


#include "WireCellSignal/ElectronicsConfig.h"
#include "WireCellSignal/ConvolutedResponse.h"
#include "WireCellSignal/DetGenerativeFDS.h"
#include "WireCellSignal/GenerativeFDS.h"
#include "WireCellSignal/GenNoise.h"

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"
#include "WireCellNav/DetectorGDS.h"
#include "WireCellSst/MCTruth.h"

#include "TH1F.h"
#include "TGraph.h"

namespace WireCell2dToy {
  class ToySignalSimuFDS : public WireCell::FrameDataSource
  {
  public:
    ToySignalSimuFDS(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, int bins_per_frame1 = 9600, int nframes_total = -1,float time_offset_uv = 0, float time_offset_uw = 0, int flag_random = 1, float overall_time_offset = 0, int overall_time_shift = 0);
    ToySignalSimuFDS(WireCell::FrameDataSource& fds, const WireCell::DetectorGDS& gds, int bins_per_frame1 = 9600, int nframes_total = -1,float time_offset_uv = 0, float time_offset_uw = 0, int flag_random = 1, float overall_time_offset = 0, int overall_time_shift = 0);
    ToySignalSimuFDS(WireCell::FrameDataSource& fds, const WireCell::DetectorGDS& gds, WireCellSignal::ElectronicsConfig& conf, int nframes_total = -1,float time_offset_uv = 0, float time_offset_uw = 0, int flag_random = 1, float overall_time_offset = 0, int overall_time_shift = 0);
    ToySignalSimuFDS(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, WireCellSignal::ElectronicsConfig& conf, int nframes_total = -1,float time_offset_uv = 0, float time_offset_uw = 0, int flag_random = 1, float overall_time_offset = 0, int overall_time_shift = 0);

    ~ToySignalSimuFDS();

    virtual int size() const;
    virtual int jump(int frame_number);
    //fixed it ...
   
    void Save();
    
  private:
    WireCell::FrameDataSource& fds;

    int gds_flag;
    const WireCell::GeomDataSource* gds;
    const WireCell::DetectorGDS* dgds;
    int  max_frames;
    int overall_time_shift;

    float time_offset_uv, time_offset_uw;
    int flag_random;
    
    int nwire_u, nwire_v, nwire_w;

    //test save ... 
    TH1F **hu1;
    TH1F **hv1;
    TH1F **hw1;
    //test save

    TH1F *hu;
    TH1F *hv;
    TH1F *hw;
    

    TGraph *gu, *gv, *gw;

    TH1F *hur, *hvr, *hwr;

  public:
    WireCellSignal::GenNoise *noise;
    WireCellSignal::ElectronicsConfig* config;
    
  };
}

#endif
