#include "WireCell2dToy/uBooNE_Data_After_ROI.h"

using namespace WireCell;


WireCell2dToy::uBooNEDataAfterROI::uBooNEDataAfterROI(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, WireCell2dToy::uBooNEDataROI& rois, int rebin)
  : fds(fds)
  , gds(gds)
  , rois(rois)
  , rebin(rebin)
{
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  bins_per_frame = int(fds.Get_Bins_Per_Frame()/rebin);
}


WireCell2dToy::uBooNEDataAfterROI::~uBooNEDataAfterROI(){
  
}

int WireCell2dToy::uBooNEDataAfterROI::size() const{
  return 1;
}


int WireCell2dToy::uBooNEDataAfterROI::jump(int frame_number){
  if (frame.index == frame_number) {
    return frame_number;
  }
  frame.clear();

  fds.jump(frame_number);
  const Frame& frame1 = fds.get();
  size_t ntraces = frame1.traces.size();
  int nticks = fds.Get_Bins_Per_Frame();

  std::vector<std::pair<int,int>>& uboone_rois = rois.get_loose_rois(0);

  TH1F *htemp_signal = new TH1F("htemp_signal","htemp_signal",nticks,0,nticks);
  TH1F *htemp_flag = new TH1F("htemp_flag","htemp_flag",nticks,0,nticks);
  
  for (size_t ind=0; ind<ntraces; ++ind) {
    htemp_signal->Reset();
    htemp_flag->Reset();
    
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();

    //copy signal
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp_signal->SetBinContent(tt,trace.charge.at(i));
    }
    
    //set ROIs
    uboone_rois = rois.get_loose_rois(chid);
   

    // std::cout << htemp_signal->GetBinContent(100) << " " << uboone_rois.size() << std::endl;
    
    for (int i=0;i!=uboone_rois.size();i++){
      for (int j=uboone_rois.at(i).first; j<=uboone_rois.at(i).second; j++){
	htemp_flag->SetBinContent(j+1,1);
      }
    }
    //remove outside ROIs
    for (int i=0;i!=htemp_signal->GetNbinsX();i++){
      if (htemp_flag->GetBinContent(i+1)==0)
	htemp_signal->SetBinContent(i+1,0);
    }

    //apply adaptive baseline ... 
    for (int i=0;i!=uboone_rois.size();i++){
      int start_bin = uboone_rois.at(i).first;
      int end_bin = uboone_rois.at(i).second;
      float start_content = htemp_signal->GetBinContent(start_bin+1);
      float end_content = htemp_signal->GetBinContent(end_bin+1);

      for (int j=start_bin;j<=end_bin;j++){
	float content = start_content + (j-start_bin) * 1.0 /(end_bin - start_bin) * (end_content-start_content);
	htemp_signal->SetBinContent(j+1,htemp_signal->GetBinContent(j+1) - content);
      }
    }

    // save into frames ... 
    WireCell::Trace trace1;
    trace1.chid = chid;
    trace1.tbin = 0;
    trace1.charge.resize(bins_per_frame, 0.0);
    
    for (int i=0;i!=bins_per_frame;i++){
      float sum = 0;
      for (int j=0;j!=rebin;j++){
	sum += htemp_signal->GetBinContent(i*rebin+j+1);
	//sum += htemp_flag->GetBinContent(i*rebin+j+1);
      }
      trace1.charge.at(i) = sum;
      //      std::cout << rebin << " " << i << " " << sum << std::endl;
    }
    frame.traces.push_back(trace1);
  }
  
  delete htemp_signal;
  delete htemp_flag;

  frame.index = frame_number;
  return frame.index;

}
