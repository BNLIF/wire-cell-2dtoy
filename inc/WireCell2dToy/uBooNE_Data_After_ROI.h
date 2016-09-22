#ifndef WIRECELL2dToy_uBOONE_Data_After_ROI_H
#define WIRECELL2dToy_uBOONE_Data_After_ROI_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"
#include "WireCell2dToy/uBooNE_Data_ROI.h"
#include "WireCellData/SignalROI.h"
#include "TH1F.h"

namespace WireCell2dToy{
  class uBooNEDataAfterROI : public WireCell::FrameDataSource
  {
  public:
    uBooNEDataAfterROI(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, WireCell2dToy::uBooNEDataROI& rois, int rebin=6);
    ~uBooNEDataAfterROI();

    virtual int jump(int frame_number);
    virtual int size() const;

    void Clear();
  private:
    WireCell::FrameDataSource& fds;
    const WireCell::GeomDataSource& gds;
    WireCell2dToy::uBooNEDataROI& rois;
    int nwire_u;
    int nwire_v;
    int nwire_w;
    
    void generate_merge_ROIs();
    
    void CleanUpROIs();
    void BreakROIs();
    void BreakROI(WireCell::SignalROI *roi, float rms);
    void BreakROI1(WireCell::SignalROI *roi);
    void ShrinkROIs();
    void ShrinkROI(WireCell::SignalROI *roi);

    void unlink(WireCell::SignalROI *prev_roi, WireCell::SignalROI *next_roi);
    void link(WireCell::SignalROI *prev_roi, WireCell::SignalROI *next_roi);

    WireCell::SignalROIChList rois_u_tight;
    WireCell::SignalROIChList rois_v_tight;
    WireCell::SignalROIChList rois_w_tight;
    
    WireCell::SignalROIChList rois_u_loose;
    WireCell::SignalROIChList rois_v_loose;
    WireCell::SignalROIChList rois_w_loose;
    
    WireCell::SignalROIMap front_rois;
    WireCell::SignalROIMap back_rois;
    WireCell::SignalROIMap contained_rois;

    int rebin; 
  };
}

#endif
