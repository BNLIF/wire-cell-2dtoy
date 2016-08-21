#ifndef WIRECELL2dToy_uBooNE_Data_2D_Deconvolution_H
#define WIRECELL2dToy_uBooNE_Data_2D_Deconvolution_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"
#include "WireCellData/GeomWire.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"


namespace WireCell2dToy{
  class uBooNEData2DDeconvolutionFDS : public WireCell::FrameDataSource
  {
  public:
    uBooNEData2DDeconvolutionFDS(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, WireCell::ChirpMap& umap, WireCell::ChirpMap& vmap, WireCell::ChirpMap& wmap, int nframes_total = -1, float time_offset_uv = 0, float time_offset_uw = 0, float overall_time_offset = 0);
    ~uBooNEData2DDeconvolutionFDS();

    virtual int size() const;
    virtual int jump(int frame_number);
    
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

    TGraph **gu_2D_g;
    TGraph **gv_2D_g;
    TGraph **gw_2D_g;
    
    TH2I *hu_2D_g;
    TH2I *hv_2D_g;
    TH2I *hw_2D_g;

    float scale_u_2d, scale_v_2d;
    
   
    void Deconvolute_2D(int plane);

    void restore_baseline(TH1F *h1);
    
  };
}

#endif
