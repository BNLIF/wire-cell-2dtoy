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

    std::vector<std::pair<int,int>>& get_loose_rois(int chid) {
      if (chid < nwire_u){
	return loose_rois_u.at(chid);
      }else if (chid < nwire_u + nwire_v){
	return loose_rois_v.at(chid - nwire_u);
      }else{
      	return loose_rois_w.at(chid - nwire_u - nwire_v);
      }
    }

    /* std::vector<std::pair<int,int>>& get_others_rois(int chid) { */
    /*   if (chid < nwire_u){ */
    /* 	return others_rois_u.at(chid); */
    /*   }else if (chid < nwire_u + nwire_v){ */
    /* 	return others_rois_v.at(chid - nwire_u); */
    /*   }else{ */
    /*   	return others_rois_w.at(chid - nwire_u - nwire_v); */
    /*   } */
    /* } */

    /* std::vector<std::pair<int,int>>& get_combined_rois(int chid) { */
    /*   if (chid < nwire_u){ */
    /* 	return combined_rois_u.at(chid); */
    /*   }else if (chid < nwire_u + nwire_v){ */
    /* 	return combined_rois_v.at(chid - nwire_u); */
    /*   }else{ */
    /*   	return combined_rois_w.at(chid - nwire_u - nwire_v); */
    /*   } */
    /* } */

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

    std::vector<std::vector<std::pair<int,int>>> self_rois_u; // tight ROIs
    std::vector<std::vector<std::pair<int,int>>> self_rois_v; // tight ROIs
    std::vector<std::vector<std::pair<int,int>>> self_rois_w; // tight ROIs
    
    std::vector<std::vector<std::pair<int,int>>> loose_rois_u; // tight ROIs
    std::vector<std::vector<std::pair<int,int>>> loose_rois_v; // tight ROIs
    std::vector<std::vector<std::pair<int,int>>> loose_rois_w; // tight ROIs

    /* std::vector<std::vector<std::pair<int,int>>> others_rois_u; */
    /* std::vector<std::vector<std::pair<int,int>>> others_rois_v; */
    /* std::vector<std::vector<std::pair<int,int>>> others_rois_w; */
    /* std::vector<std::vector<std::pair<int,int>>> combined_rois_u; */
    /* std::vector<std::vector<std::pair<int,int>>> combined_rois_v; */
    /* std::vector<std::vector<std::pair<int,int>>> combined_rois_w; */

    std::vector<float> uplane_rms; // calibrated field response
    std::vector<float> vplane_rms; // calibrated field response
    std::vector<float> wplane_rms; // calibrated field response

    void create_ROI_connect_info(float asy = 0.1);
    void restore_baseline(TH1F *h1);
    double cal_rms(TH1F *h1, int chid);
    void find_ROI_by_decon_itself(int th_factor_ind = 5, int th_factor_col = 5, int pad = 0);
    void extend_ROI_self(int pad = 5);
    
    void find_ROI_loose(int rebin=6);
    
    Double_t local_ave(TH1F *h1, Int_t bin, Int_t width);
    Int_t find_ROI_end(TH1F*h1, Int_t bin, Double_t th = 0);
    Int_t find_ROI_begin(TH1F*h1, Int_t bin, Double_t th = 0);
    

    //void find_ROI_by_raw_itself(int th_factor = 5, int pad = 5);
    // void extend_ROI_others(int pad = 5);

    //void find_ROI_by_others();
    //void merge_ROIs();

    

  };
}


#endif
