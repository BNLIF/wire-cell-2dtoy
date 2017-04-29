#include "WireCell2dToy/uBooNE_Data_After_ROI_gaus.h"

using namespace WireCell;

WireCell2dToy::uBooNEDataAfterROI_Gaus::uBooNEDataAfterROI_Gaus(WireCell2dToy::uBooNEData2DDeconvolutionFDS* frame_data, WireCell2dToy::uBooNEDataAfterROI* frame_rois)
  : frame_data(frame_data)
  , frame_rois(frame_rois)
{

}

WireCell2dToy::uBooNEDataAfterROI_Gaus::~uBooNEDataAfterROI_Gaus(){
  
}

void WireCell2dToy::uBooNEDataAfterROI_Gaus::Clear(){

}

int WireCell2dToy::uBooNEDataAfterROI_Gaus::size() const{
  return 1;
}

int WireCell2dToy::uBooNEDataAfterROI_Gaus::jump(int frame_number){
  Clear();
  if (frame.index == frame_number) {
    return frame_number;
  }
  frame.clear();

  // prepare to fill content
  TH2I *hu_data = frame_data->get_u_gaus();
  TH2I *hv_data = frame_data->get_v_gaus();  
  TH2I *hw_data = frame_data->get_w_gaus();

  // get ROIs  
  WireCell::SignalROIChList& rois_u = frame_rois->get_u_rois();
  WireCell::SignalROIChList& rois_v = frame_rois->get_v_rois();
  WireCell::SignalROIChList& rois_w = frame_rois->get_w_rois();

  
}
