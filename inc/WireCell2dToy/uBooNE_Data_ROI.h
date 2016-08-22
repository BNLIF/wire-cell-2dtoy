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

  private:
    WireCell::FrameDataSource& fds;
    const WireCell::GeomDataSource& gds;
    WireCell::ChirpMap& umap;
    WireCell::ChirpMap& vmap;
    WireCell::ChirpMap& wmap;

    int nwire_u;
    int nwire_v;
    int nwire_w;

    void restore_baseline(TH1F *h1);
    double cal_rms(TH1F *h1, int chid);
    void find_ROI_by_itself();
  };
}


#endif
