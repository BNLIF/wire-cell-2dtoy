#include "WireCell2dToy/ToySignalPre.h"
#include "WireCellData/GeomWire.h"
#include "TFile.h"
#include "TVirtualFFT.h"

using namespace WireCell;

WireCell2dToy::ToySignalPreFDS::ToySignalPreFDS(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, int bins_per_frame, int nframes_total)
  : fds(fds)
  , gds(gds)
  , bins_per_frame(bins_per_frame)
  , max_frames(nframes_total)
{  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  hu = new TH1F*[nwire_u];
  hv = new TH1F*[nwire_v];
  hw = new TH1F*[nwire_w];
  
  for (int i=0;i!=nwire_u;i++){
    hu[i] = new TH1F(Form("U_%d",i),Form("U_%d",i),bins_per_frame,0,bins_per_frame);
  }
  for (int i=0;i!=nwire_v;i++){
    hv[i] = new TH1F(Form("V_%d",i),Form("V_%d",i),bins_per_frame,0,bins_per_frame);
  }
  for (int i=0;i!=nwire_w;i++){
    hw[i] = new TH1F(Form("W_%d",i),Form("W_%d",i),bins_per_frame,0,bins_per_frame);
  }
}


int WireCell2dToy::ToySignalPreFDS::size() const{
  return max_frames;
}

void WireCell2dToy::ToySignalPreFDS::Save(){
  TFile *file = new TFile("temp_pre.root","RECREATE");
  for (int i=0;i!=nwire_u;i++){
    TH1F *huu = (TH1F*)hu[i]->Clone(Form("U1_%d",i));
  }
  for (int i=0;i!=nwire_v;i++){
    TH1F *hvv = (TH1F*)hv[i]->Clone(Form("V1_%d",i));
  }
  for (int i=0;i!=nwire_w;i++){
    TH1F *hww = (TH1F*)hw[i]->Clone(Form("W1_%d",i));
  }
  file->Write();
  file->Close();
}


int WireCell2dToy::ToySignalPreFDS::jump(int frame_number){
  return 0;
}


WireCell2dToy::ToySignalPreFDS::~ToySignalPreFDS(){
  for (int i=0;i!=nwire_u;i++){
    delete hu[i] ;
  }
  delete hu;
  for (int i=0;i!=nwire_v;i++){
    delete hv[i] ;
  }
  delete hv;
  for (int i=0;i!=nwire_w;i++){
    delete hw[i] ;
  }
  delete hw;

 
}
