#ifndef WIRECELL2dToy_DATASIGNALWIENROI_H
#define WIRECELL2dToy_DATASIGNALWIENROI_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"
#include "WireCellData/GeomWire.h"
#include "TH1F.h"
#include "TGraph.h"

namespace WireCell2dToy {
  class DataSignalWienROIFDS : public WireCell::FrameDataSource
  {
  public:
    DataSignalWienROIFDS(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, WireCell::ChirpMap& umap, WireCell::ChirpMap& vmap, WireCell::ChirpMap& wmap, int bins_per_frame1 = 9600, int nframes_total = -1, float time_offset_uv = 0, float time_offset_uw = 0, float overall_time_offset = 0);
    ~DataSignalWienROIFDS();

    virtual int size() const;
    virtual int jump(int frame_number);
    //fixed it ...
    std::vector <float>& get_uplane_rms(){return uplane_rms;};
    std::vector <float>& get_vplane_rms(){return vplane_rms;};
    std::vector <float>& get_wplane_rms(){return wplane_rms;};

    void Save();
    
  private:
    WireCell::FrameDataSource& fds;
    const WireCell::GeomDataSource& gds;
    int  max_frames;
    int nbin;

    WireCell::ChirpMap& umap;
    WireCell::ChirpMap& vmap;
    WireCell::ChirpMap& wmap;

    float time_offset_uv;
    float time_offset_uw;
    float overall_time_offset;
    
    int nwire_u, nwire_v, nwire_w;

    /* TH1F **hu; */
    /* TH1F **hv; */
    /* TH1F **hw; */

    TH1F *hu; //good 
    TGraph **gu; //Now an array 
    TH1F *hur; // good
    TH1 *hmr_u, *hpr_u; //good?

    TH1F *hv;
    TH1F *hw;
    TGraph *gv, *gw;    
    TH1F *hvr, *hwr;

    /* TGraph *wiener_filter_u; */
    /* TGraph *wiener_filter_v; */
    
    
    
    TH1 *hmr_v, *hpr_v;
    TH1 *hmr_w, *hpr_w;
    
    std::vector<float> uplane_rms;
    std::vector<float> vplane_rms;
    std::vector<float> wplane_rms;

  };
}

#endif
