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


    void restore_baseline(TH1F *h1);
    void find_ROI_by_itself();
  };
}


#endif
