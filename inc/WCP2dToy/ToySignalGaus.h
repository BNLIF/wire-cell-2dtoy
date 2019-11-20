#ifndef WIRECELL2dToy_TOYSIGNALGAUS_H
#define WIRECELL2dToy_TOYSIGNALGAUS_H

#include "WCPNav/FrameDataSource.h"
#include "WCPNav/GeomDataSource.h"
#include "WCPNav/DetectorGDS.h"

#include "TH1F.h"
#include "TGraph.h"

namespace WCP2dToy {
  class ToySignalGausFDS : public WCP::FrameDataSource
  {
  public:
    ToySignalGausFDS(WCP::FrameDataSource& fds,  const WCP::GeomDataSource& gds, int bins_per_frame1 = 9600, int nframes_total = -1, float time_offset_uv = 0, float time_offset_uw = 0, float overall_time_offset = 0);
    ToySignalGausFDS(WCP::FrameDataSource& fds,  const WCP::DetectorGDS& gds, int bins_per_frame1 = 9600, int nframes_total = -1, float time_offset_uv = 0, float time_offset_uw = 0, float overall_time_offset = 0);
    ~ToySignalGausFDS();

    virtual int size() const;
    virtual int jump(int frame_number);

    void Save();

  private:
    WCP::FrameDataSource& fds;
    int gds_flag;
    const WCP::GeomDataSource* gds;
    const WCP::DetectorGDS* dgds;
    int max_frames;
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

    TH1F *hfilter_time_gaus;
    TH1 *hfilter_gaus;
    TF1 *filter_g;

    TH1F *hur, *hvr, *hwr;

    TGraph *gu, *gv, *gw;

    TH1 *hmr_u, *hpr_u;
    TH1 *hmr_v, *hpr_v;
    TH1 *hmr_w, *hpr_w;
    
  };

}

#endif
