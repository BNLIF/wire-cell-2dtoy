#ifndef WIRECELL2dToy_uBooNE_Data_2D_Deconvolution_H
#define WIRECELL2dToy_uBooNE_Data_2D_Deconvolution_H

#include "WCPNav/FrameDataSource.h"
#include "WCPNav/GeomDataSource.h"
#include "WCPData/GeomWire.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TTree.h"


namespace WCP2dToy{
  class uBooNEData2DDeconvolutionFDS : public WCP::FrameDataSource
  {
  public:
    uBooNEData2DDeconvolutionFDS(WCP::FrameDataSource& fds, const WCP::GeomDataSource& gds, WCP::ChirpMap& umap, WCP::ChirpMap& vmap, WCP::ChirpMap& wmap, int nframes_total = -1, float time_offset_uv = 0, float time_offset_uw = 0, float overall_time_offset = 0, int flag_sim = 0);
    
    uBooNEData2DDeconvolutionFDS(TH2I *hu_decon, TH2I *hv_decon, TH2I *hw_decon, TTree *T_bad, const WCP::GeomDataSource& gds);

    ~uBooNEData2DDeconvolutionFDS();

    /* int get_run_no(){return run_no;}; */
    /* int get_subrun_no(){return subrun_no;}; */
    /* int get_event_no(){return event_no;}; */

    WCP::ChirpMap& get_u_cmap(){return umap;};
    WCP::ChirpMap& get_v_cmap(){return vmap;};
    WCP::ChirpMap& get_w_cmap(){return wmap;};
    
    TH2I* get_u_gaus(){return hu_2D_gg;};
    TH2I* get_v_gaus(){return hv_2D_gg;};
    TH2I* get_w_gaus(){return hw_2D_gg;};

    TH2I* get_u_wiener(){return hu_2D_g;};
    TH2I* get_v_wiener(){return hv_2D_g;};
    TH2I* get_w_wiener(){return hw_2D_g;};

    virtual int size() const;
    virtual int jump(int frame_number);
    void Clear();
    
  private:
    WCP::FrameDataSource *fds;
    const WCP::GeomDataSource& gds;
    int  max_frames;
    int nbin;

    /* int run_no, subrun_no, event_no; */


    WCP::ChirpMap umap;
    WCP::ChirpMap vmap;
    WCP::ChirpMap wmap;

    float time_offset_uv;
    float time_offset_uw;
    float overall_time_offset;

    int flag_sim;

    int nwire_u, nwire_v, nwire_w;

    TGraph **gu_2D_g;
    TGraph **gv_2D_g;
    TGraph **gw_2D_g;
    
    TH2I *hu_2D_g;
    TH2I *hv_2D_g;
    TH2I *hw_2D_g;

    TH2I *hu_2D_gg;
    TH2I *hv_2D_gg;
    TH2I *hw_2D_gg;

    float scale_u_2d, scale_v_2d;
    
    bool load_results_from_file;
   
    void Deconvolute_2D(int plane);

    void restore_baseline(TH1F *h1);

    /* float ref_resp_func[100]; */
    /* float u_resp_func[2400][100]; */
    /* float v_resp_func[2400][100]; */
    /* float w_resp_func[3456][100]; */
    
  };
}

#endif
