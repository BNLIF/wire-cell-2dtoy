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
    /* std::vector <float>& get_uplane_rms(){return uplane_rms_g;}; */
    /* std::vector <float>& get_vplane_rms(){return vplane_rms_g;}; */
    
    std::vector <float>& get_wplane_rms(){return wplane_rms;};

    TH2I *Get_hu_gaus(){return hu_2D_g_gaus;};
    TH2I *Get_hv_gaus(){return hv_2D_g_gaus;};
    TH2I *Get_hw_gaus(){return hw_1D_g;};
    
    void set_flag_2D(int flag){flag_2D = flag;};

    void Save();
    
  private:
    int flag_2D;
    WireCell::FrameDataSource& fds;
    const WireCell::GeomDataSource& gds;
    int  max_frames;
    int nbin;

    float scale_u_1d, scale_v_1d;
    float scale_u_2d, scale_v_2d;

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
    
    void restore_baseline(TH1F *h1);
    double cal_rms(TH1F *h1, int chid);

    void Deconvolute_U_1D_c();
    void Deconvolute_V_1D_c();
    void Deconvolute_W_1D_g();

    void Deconvolute_U_2D_g(); // both with filter and without filter 
    void Deconvolute_V_1D_g(); // both with filter and without filter 
    void Deconvolute_V_2D_g(); // both with filter and without filter 

    void ROI_cal(TH1F *h1_1, TH1F *h2_1, TH1F *h3_1, TH1F *h4_1, TH1F *h5_1, Double_t threshold0,Double_t threshold2, TH1F *hresult, TH1F *hresult1, int flag_u);
    Double_t local_ave(TH1F *h1, Int_t bin, Int_t width);
    Int_t find_ROI_begin(TH1F *h1, Int_t bin, Double_t th);
    Int_t find_ROI_end(TH1F *h1, Int_t bin, Double_t th);

    // U plane  1D-U-c
    TH2I *hu_1D_c;
    TH2I *hu_1D_c_gaus;
    TGraph *gu_1D_c;

    // 2D-U-g-f, 2D-U-g
    TH2I *hu_2D_g_f;
    TH2I *hu_2D_g;
    TH2I *hu_2D_g_gaus;
    TGraph **gu_2D_g;
    
          
    //V plane 1-D-V-c
    TH2I *hv_1D_c;
    TH2I *hv_1D_c_gaus;
    TGraph *gv_1D_c;

    //1-D-V-g-f, 1-D-V-g
    /* TH2I *hv_1D_g_f; */
    /* TH2I *hv_1D_g; */
    /* TH2I *hv_1D_g_gaus; */
    /* TGraph *gv_1D_g; */

    TH2I *hv_2D_g_f; 
    TH2I *hv_2D_g;
    TH2I *hv_2D_g_gaus;
    TGraph **gv_2D_g;


    //W-plane  1-D-W-g
    TH2I *hw_1D_g;
    TH2I *hw_1D_g_gaus;
    TGraph *gw_1D_g;
    TGraph *gw_1D_c;

    
 

  };
}

#endif
