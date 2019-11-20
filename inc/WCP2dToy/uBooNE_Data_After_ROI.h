#ifndef WIRECELL2dToy_uBOONE_Data_After_ROI_H
#define WIRECELL2dToy_uBOONE_Data_After_ROI_H

#include "WCPNav/FrameDataSource.h"
#include "WCPNav/GeomDataSource.h"
#include "WCP2dToy/uBooNE_Data_ROI.h"
#include "WCP2dToy/uBooNE_Data_2D_Deconvolution.h"
#include "WCPData/SignalROI.h"
#include "TH1F.h"

namespace WCP2dToy{
  class uBooNEDataAfterROI : public WCP::FrameDataSource
  {
  public:
    uBooNEDataAfterROI(WCP::FrameDataSource& fds, const WCP::GeomDataSource& gds, WCP2dToy::uBooNEDataROI& rois, int rebin=6);
    ~uBooNEDataAfterROI();

    virtual int jump(int frame_number);
    virtual int size() const;

    WCP::SignalROIChList& get_u_rois(){return rois_u_loose;};
    WCP::SignalROIChList& get_v_rois(){return rois_v_loose;};
    WCP::SignalROIChList& get_w_rois(){return rois_w_tight;};


    void Clear();
  private:
    WCP::FrameDataSource& fds;
    const WCP::GeomDataSource& gds;
    WCP2dToy::uBooNEDataROI& rois;
    int nwire_u;
    int nwire_v;
    int nwire_w;
    
    void generate_merge_ROIs();
    
    void CleanUpCollectionROIs();
    void CleanUpInductionROIs();

    void CleanUpROIs();
    void CheckROIs();
    
    void BreakROIs();
    void BreakROI(WCP::SignalROI *roi, float rms);
    void BreakROI1(WCP::SignalROI *roi);
    void ShrinkROIs();
    void ShrinkROI(WCP::SignalROI *roi);

    void ExtendROIs();

    void TestROIs();
    
    void unlink(WCP::SignalROI *prev_roi, WCP::SignalROI *next_roi);
    void link(WCP::SignalROI *prev_roi, WCP::SignalROI *next_roi);

    WCP::SignalROIChList rois_u_tight;
    WCP::SignalROIChList rois_v_tight;
    WCP::SignalROIChList rois_w_tight;
    
    WCP::SignalROIChList rois_u_loose;
    WCP::SignalROIChList rois_v_loose;
    WCP::SignalROIChList rois_w_loose;
    
    WCP::SignalROIMap front_rois;
    WCP::SignalROIMap back_rois;
    WCP::SignalROIMap contained_rois;

    int rebin; 
  };
}

#endif
