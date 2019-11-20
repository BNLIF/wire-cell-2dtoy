#ifndef WIRECELL2dToy_TOYSIGNALSIMU_H
#define WIRECELL2dToy_TOYSIGNALSIMU_H

#include "WCPNav/FrameDataSource.h"
#include "WCPNav/GeomDataSource.h"
#include "WCPNav/DetectorGDS.h"
#include "WCPSst/MCTruth.h"

#include "WCPSignal/ElectronicsConfig.h"
#include "WCPSignal/ConvolutedResponse.h"
#include "WCPSignal/DetGenerativeFDS.h"
//#include "WCPSignal/GenerativeFDS.h"
#include "WCPSignal/GenNoise.h"



#include "TH1F.h"
#include "TGraph.h"

namespace WCP2dToy {
  class ToySignalSimuFDS : public WCP::FrameDataSource
  {
  public:
    ToySignalSimuFDS(WCP::FrameDataSource& fds, const WCP::GeomDataSource& gds, int bins_per_frame1 = 9600, int nframes_total = -1,float time_offset_uv = 0, float time_offset_uw = 0, int flag_random = 1, float overall_time_offset = 0, int overall_time_shift = 0);
    ToySignalSimuFDS(WCP::FrameDataSource& fds, const WCP::DetectorGDS& gds, int bins_per_frame1 = 9600, int nframes_total = -1,float time_offset_uv = 0, float time_offset_uw = 0, int flag_random = 1, float overall_time_offset = 0, int overall_time_shift = 0);
    ToySignalSimuFDS(WCP::FrameDataSource& fds, const WCP::DetectorGDS& gds, WCPSignal::ElectronicsConfig& conf, int nframes_total = -1,float time_offset_uv = 0, float time_offset_uw = 0, int flag_random = 1, float overall_time_offset = 0, int overall_time_shift = 0);
    ToySignalSimuFDS(WCP::FrameDataSource& fds, const WCP::GeomDataSource& gds, WCPSignal::ElectronicsConfig& conf, int nframes_total = -1,float time_offset_uv = 0, float time_offset_uw = 0, int flag_random = 1, float overall_time_offset = 0, int overall_time_shift = 0);

    ~ToySignalSimuFDS();

    virtual int size() const;
    virtual int jump(int frame_number);
    //fixed it ...
   
    void Save();
    
  private:
    WCP::FrameDataSource& fds;

    int simulation_type; 
    // 1 for the simple simulation (1D the oldest one)
    // 2 for Xiaoyue's simulation (2D)

    int gds_flag;
    const WCP::GeomDataSource* gds;
    const WCP::DetectorGDS* dgds;
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
    WCPSignal::GenNoise *noise;
    WCPSignal::ElectronicsConfig* config;
    
  };
}

#endif
