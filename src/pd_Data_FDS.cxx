#include "WireCell2dToy/pd_Data_FDS.h"

using namespace WireCell;

WireCell2dToy::pdDataFDS::pdDataFDS(const WireCell::GeomDataSource& gds, TH2I *hu_decon, TH2I *hv_decon, TH2I *hw_decon, int eve_num)
  : gds(gds)
{
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  
  
  frame.clear();		// win or lose, we start anew

  frame.index =eve_num;
  bins_per_frame = hu_decon->GetNbinsY();
  // U plane
  for (size_t ind=0; ind < hu_decon->GetNbinsX(); ++ind) {
    WireCell::Trace trace;
    trace.chid = ind;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hu_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }
  
  // V plane
  for (size_t ind=0; ind < hv_decon->GetNbinsX(); ++ind) {
    WireCell::Trace trace;
    trace.chid = ind + nwire_u;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }

  // W plane
  for (size_t ind=0; ind < hw_decon->GetNbinsX(); ++ind) {
    WireCell::Trace trace;
    trace.chid = ind + nwire_u + nwire_v;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hw_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }

  //std::cout << frame.traces.size() << " " << bins_per_frame << std::endl;
}



WireCell2dToy::pdDataFDS::pdDataFDS(const WireCell::GeomDataSource& gds, TH2F *hu_decon, TH2F *hv_decon, TH2F *hw_decon, int eve_num)
  : gds(gds)
{
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  
  
  frame.clear();		// win or lose, we start anew

  frame.index =eve_num;
  bins_per_frame = hu_decon->GetNbinsY();
  // U plane
  for (size_t ind=0; ind < hu_decon->GetNbinsX(); ++ind) {
    WireCell::Trace trace;
    trace.chid = ind;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hu_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }
  
  // V plane
  for (size_t ind=0; ind < hv_decon->GetNbinsX(); ++ind) {
    WireCell::Trace trace;
    trace.chid = ind + nwire_u;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }

  // W plane
  for (size_t ind=0; ind < hw_decon->GetNbinsX(); ++ind) {
    WireCell::Trace trace;
    trace.chid = ind + nwire_u + nwire_v;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hw_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }

  //std::cout << frame.traces.size() << " " << bins_per_frame << std::endl;
}

void WireCell2dToy::pdDataFDS::refresh(TH2F *hu_decon, TH2F *hv_decon, TH2F *hw_decon, int eve_num){
    
  frame.clear();		// win or lose, we start anew

  frame.index =eve_num;
  bins_per_frame = hu_decon->GetNbinsY();
  // U plane
  for (size_t ind=0; ind < hu_decon->GetNbinsX(); ++ind) {
    WireCell::Trace trace;
    trace.chid = ind;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hu_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }
  
  // V plane
  for (size_t ind=0; ind < hv_decon->GetNbinsX(); ++ind) {
    WireCell::Trace trace;
    trace.chid = ind + nwire_u;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }

  // W plane
  for (size_t ind=0; ind < hw_decon->GetNbinsX(); ++ind) {
    WireCell::Trace trace;
    trace.chid = ind + nwire_u + nwire_v;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hw_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }
}


WireCell2dToy::pdDataFDS::~pdDataFDS(){
}

int WireCell2dToy::pdDataFDS::jump(int frame_number){
  return frame.index;
}

int WireCell2dToy::pdDataFDS::size() const{
  return 1;
}
