#ifndef WIRECELL2dToy_uBOONE_Data_After_ROI_gaus_H
#define WIRECELL2dToy_uBOONE_Data_After_ROI_gaus_H


#include "WireCell2dToy/uBooNE_Data_After_ROI.h"
#include "WireCell2dToy/uBooNE_Data_2D_Deconvolution.h"

namespace WireCell2dToy{
  class uBooNEDataAfterROI_Gaus : public WireCell::FrameDataSource
  {
  public: 
    uBooNEDataAfterROI_Gaus(uBooNEData2DDeconvolutionFDS* frame_data, uBooNEDataAfterROI* frame_rois);
    ~uBooNEDataAfterROI_Gaus();

    virtual int jump(int frame_number);
    virtual int size() const;

    void Clear();

  private:
    uBooNEData2DDeconvolutionFDS* frame_data;
    uBooNEDataAfterROI* frame_rois;
    
  };
}

#endif
