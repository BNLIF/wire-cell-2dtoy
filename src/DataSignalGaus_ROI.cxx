#include "WireCell2dToy/DataSignalGaus_ROI.h"

using namespace WireCell;

WireCell2dToy::DataSignalGausROIFDS::DataSignalGausROIFDS(WireCell2dToy::DataSignalWienROIFDS& fds, int nframes_total)
  : fds(fds)
  , max_frames(nframes_total)
{
}

WireCell2dToy::DataSignalGausROIFDS::~DataSignalGausROIFDS(){
}

int WireCell2dToy::DataSignalGausROIFDS::size() const{
  return max_frames;
}

int WireCell2dToy::DataSignalGausROIFDS::jump(int frame_number){
  // fill the frame data ... 
  if (frame.index == frame_number) {
    return frame_number;
  }
  frame.clear();
  
  TH2I *hu_gaus = fds.Get_hu_gaus();
  TH2I *hv_gaus = fds.Get_hv_gaus();
  TH2I *hw_gaus = fds.Get_hw_gaus();

  int nticks = hu_gaus->GetNbinsY();
  int nwire_u = hu_gaus->GetNbinsX();
  int nwire_v = hv_gaus->GetNbinsX();
  int nwire_w = hw_gaus->GetNbinsX();
  
  for (int i=0;i!=nwire_u;i++){
    Trace t;
    t.chid = i;
    t.tbin = 0;
    t.charge.resize(nticks, 0.0);
    for (int j=0;j!= nticks; j++){
      t.charge.at(j) = hu_gaus->GetBinContent(j+1);
    }
    frame.traces.push_back(t);
  }
  
  for (int i=0;i!=nwire_v;i++){
    Trace t;
    t.chid = i+nwire_u;
    t.tbin = 0;
    t.charge.resize(nticks, 0.0);
    for (int j=0;j!= nticks; j++){
      t.charge.at(j) = hv_gaus->GetBinContent(j+1);
    }
    frame.traces.push_back(t);
  }

  for (int i=0;i!=nwire_w;i++){
    Trace t;
    t.chid = i+nwire_u+nwire_v;
    t.tbin = 0;
    t.charge.resize(nticks, 0.0);
    for (int j=0;j!= nticks; j++){
      t.charge.at(j) = hw_gaus->GetBinContent(j+1);
    }
    frame.traces.push_back(t);
  }


  frame.index = frame_number;
  return frame.index;
}
