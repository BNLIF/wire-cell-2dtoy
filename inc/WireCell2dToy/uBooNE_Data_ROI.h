#ifndef WIRECELL2dToy_uBooNE_Data_ROI_H
#define WIRECELL2dToy_uBooNE_Data_ROI_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"

#include "TH2F.h"


namespace WireCell2dToy{
  class uBooNEDataROI 
  {
  public:
    uBooNEDataROI(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, WireCell::ChirpMap& umap, WireCell::ChirpMap& vmap, WireCell::ChirpMap& wmap);
    ~uBooNEDataROI();

    std::vector<std::pair<int,int>>& get_self_rois(int chid) {
      if (chid < nwire_u){
	return self_rois_u.at(chid);
      }else if (chid < nwire_u + nwire_v){
	return self_rois_v.at(chid - nwire_u);
      }else if (chid < nwire_u + nwire_v + nwire_w){
      	return self_rois_w.at(chid - nwire_u - nwire_v);
      }
    }

  private:
    WireCell::FrameDataSource& fds;
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

    void restore_baseline(TH1F *h1);
    double cal_rms(TH1F *h1, int chid);
    void find_ROI_by_itself();
  };
}


#endif
