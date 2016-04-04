#ifndef WIRECELL2dToy_DATASIGNALWIENROI_H
#define WIRECELL2dToy_DATASIGNALWIENROI_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"
#include "WireCellData/GeomWire.h"
#include "TH1F.h"
#include "TH2F.h"
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

    std::vector<float> uplane_rms; // calibrated field response
    std::vector<float> vplane_rms; // calibrated field response
    std::vector<float> wplane_rms; // calibrated field response

    std::vector<float> uplane_rms_g;
    std::vector<float> vplane_rms_g;
    std::vector<float> wplane_rms_g;

    void Deconvolute_U_1D_c();
    void Deconvolute_V_1D_c();
    void Deconvolute_W_1D_g();

    void Deconvolute_U_2D_g(); // both with filter and without filter 
    void Deconvolute_V_1D_g(); // both with filter and without filter 

    // U plane  1D-U-c
    TH2F *hu_1D_c;
    TGraph *gu_1D_c;

    // 2D-U-g-f, 2D-U-g
    TH2F *hu_2D_g_f;
    TH2F *hu_2D_g;
    TGraph **gu_2D_g;
    
          
    //V plane 1-D-V-c
    TH2F *hv_1D_c;
    TGraph *gv_1D_c;

    //1-D-V-g-f, 1-D-V-g
    TH2F *hv_1D_g_f;
    TH2F *hv_1D_g;
    TGraph *gv_1D_g;

    //W-plane  1-D-W-g
    TH2F *hw_1D_g;
    TGraph *gw_1D_g;
    TGraph *gw_1D_c;

    
 

  };
}

#endif
