#ifndef WIRECELL2dToy_TOYSIGNALWIEN_H
#define WIRECELL2dToy_TOYSIGNALWIEN_H

#include "WCPNav/FrameDataSource.h"
#include "WCPNav/GeomDataSource.h"
#include "WCPNav/DetectorGDS.h"

#include "TH1F.h"
#include "TGraph.h"

namespace WCP2dToy {
  class ToySignalWienFDS : public WCP::FrameDataSource
  {
  public:
    ToySignalWienFDS(WCP::FrameDataSource& fds, const WCP::GeomDataSource& gds, int bins_per_frame1 = 9600, int nframes_total = -1, float time_offset_uv = 0, float time_offset_uw = 0, float overall_time_offset = 0);
    ToySignalWienFDS(WCP::FrameDataSource& fds, const WCP::DetectorGDS& gds, int bins_per_frame1 = 9600, int nframes_total = -1, float time_offset_uv = 0, float time_offset_uw = 0, float overall_time_offset = 0);
    ~ToySignalWienFDS();

    virtual int size() const;
    virtual int jump(int frame_number);
    //fixed it ...
    
    std::vector <float>& get_uplane_rms(){return uplane_rms;};
    std::vector <float>& get_vplane_rms(){return vplane_rms;};
    std::vector <float>& get_wplane_rms(){return wplane_rms;};


    void Save();
    
  private:
    WCP::FrameDataSource& fds;
    int gds_flag;
    const WCP::GeomDataSource* gds;
    const WCP::DetectorGDS* dgds;
    int  max_frames;
    int nbin;

    float time_offset_uv;
    float time_offset_uw;
    float overall_time_offset;
    
    int nwire_u, nwire_v, nwire_w;

    /* TH1F **hu; */
    /* TH1F **hv; */
    /* TH1F **hw; */

    TH1F *hu;
    TH1F *hv;
    TH1F *hw;

    TGraph *gu, *gv, *gw;

    TH1F *hur, *hvr, *hwr;
    
    TH1 *hmr_u, *hpr_u;
    TH1 *hmr_v, *hpr_v;
    TH1 *hmr_w, *hpr_w;

    std::vector<float> uplane_rms;
    std::vector<float> vplane_rms;
    std::vector<float> wplane_rms;

    
  };
}

#endif
