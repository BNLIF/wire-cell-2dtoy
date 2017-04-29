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

  const Frame& frame1 = frame_data->get();
  size_t ntraces = frame1.traces.size();
  int nticks = frame_data->Get_Bins_Per_Frame();

  int bins_per_frame = frame_rois->Get_Bins_Per_Frame();
  int rebin = round(nticks/bins_per_frame);

  TH1F *htemp_signal = new TH1F("htemp_signal","htemp_signal",nticks,0,nticks);
  
  


  // prepare to fill content
  TH2I *hu_data = frame_data->get_u_gaus();
  TH2I *hv_data = frame_data->get_v_gaus();  
  TH2I *hw_data = frame_data->get_w_gaus();

  // get ROIs  
  WireCell::SignalROIChList& rois_u = frame_rois->get_u_rois();
  WireCell::SignalROIChList& rois_v = frame_rois->get_v_rois();
  WireCell::SignalROIChList& rois_w = frame_rois->get_w_rois();

  
  // load results back into data
  for (size_t ind=0; ind<ntraces; ++ind) {
    htemp_signal->Reset();
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();

    // do the baseline correction and load the results
    
    


    // save into frames ... 
    WireCell::Trace trace1;
    trace1.chid = chid;
    trace1.tbin = 0;
    trace1.charge.resize(bins_per_frame, 0.0);
    
    for (int i=0;i!=bins_per_frame;i++){
      float sum = 0;
      for (int j=0;j!=rebin;j++){
	sum += htemp_signal->GetBinContent(i*rebin+j+1);
      }
      trace1.charge.at(i) = sum;
    }
    frame.traces.push_back(trace1);
  }
  

  delete htemp_signal;
  frame.index = frame_number;
  return frame.index;
}
