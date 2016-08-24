#ifndef WIRECELL2dToy_uBooNE_Data_ROI_H
#define WIRECELL2dToy_uBooNE_Data_ROI_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"

#include "TH2F.h"


namespace WireCell2dToy{
  class uBooNEDataROI 
  {
  public:
    uBooNEDataROI(WireCell::FrameDataSource& raw_fds, WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, WireCell::ChirpMap& umap, WireCell::ChirpMap& vmap, WireCell::ChirpMap& wmap);
    ~uBooNEDataROI();

    std::vector<std::pair<int,int>>& get_self_rois(int chid) {
      if (chid < nwire_u){
	return self_rois_u.at(chid);
      }else if (chid < nwire_u + nwire_v){
	return self_rois_v.at(chid - nwire_u);
      }else{
      	return self_rois_w.at(chid - nwire_u - nwire_v);
      }
    }

    std::vector<std::pair<int,int>>& get_others_rois(int chid) {
      if (chid < nwire_u){
	return others_rois_u.at(chid);
      }else if (chid < nwire_u + nwire_v){
	return others_rois_v.at(chid - nwire_u);
      }else{
      	return others_rois_w.at(chid - nwire_u - nwire_v);
      }
    }

    std::vector<std::pair<int,int>>& get_combined_rois(int chid) {
      if (chid < nwire_u){
	return combined_rois_u.at(chid);
      }else if (chid < nwire_u + nwire_v){
	return combined_rois_v.at(chid - nwire_u);
      }else{
      	return combined_rois_w.at(chid - nwire_u - nwire_v);
      }
    }

    std::vector <float>& get_uplane_rms(){return uplane_rms;};
    std::vector <float>& get_vplane_rms(){return vplane_rms;};
    std::vector <float>& get_wplane_rms(){return wplane_rms;};
    
  private:
    WireCell::FrameDataSource& fds;
    WireCell::FrameDataSource& raw_fds;
    const WireCell::GeomDataSource& gds;
    WireCell::ChirpMap& umap;
    WireCell::ChirpMap& vmap;
    WireCell::ChirpMap& wmap;

    int nwire_u;
    int nwire_v;
    int nwire_w;

    std::vector<std::vector<std::pair<int,int>>> self_rois_u;
    std::vector<std::vector<std::pair<int,int>>> self_rois_v;
    std::vector<std::vector<std::pair<int,int>>> self_rois_w;

    std::vector<std::vector<std::pair<int,int>>> others_rois_u;
    std::vector<std::vector<std::pair<int,int>>> others_rois_v;
    std::vector<std::vector<std::pair<int,int>>> others_rois_w;

    std::vector<std::vector<std::pair<int,int>>> combined_rois_u;
    std::vector<std::vector<std::pair<int,int>>> combined_rois_v;
    std::vector<std::vector<std::pair<int,int>>> combined_rois_w;

    std::vector<float> uplane_rms; // calibrated field response
    std::vector<float> vplane_rms; // calibrated field response
    std::vector<float> wplane_rms; // calibrated field response

    void restore_baseline(TH1F *h1);
    double cal_rms(TH1F *h1, int chid);
    void find_ROI_by_decon_itself(int th_factor = 5, int pad = 0);
    void find_ROI_by_raw_itself(int th_factor = 5, int pad = 0);
    void find_ROI_by_others();
    void merge_ROIs();

    

  };
}


#endif
