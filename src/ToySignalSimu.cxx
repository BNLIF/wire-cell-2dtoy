#include "WireCell2dToy/ToySignalSimu.h"

#include "WireCellData/GeomWire.h"
#include "TFile.h"

using namespace WireCell;

WireCell2dToy::ToySignalSimuFDS::ToySignalSimuFDS(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds,
						  int bins_per_frame, int nframes_total)
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
  
#include "data.txt"

}

int WireCell2dToy::ToySignalSimuFDS::size() const{
  return max_frames;
}

void WireCell2dToy::ToySignalSimuFDS::Save(){
  TFile *file = new TFile("temp.root","RECREATE");
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

int WireCell2dToy::ToySignalSimuFDS::jump(int frame_number){
  // do simulation
  for (int i=0;i!=nwire_u;i++){
    hu[i]->Reset();
  }
  for (int i=0;i!=nwire_v;i++){
    hv[i]->Reset();
  }
  for (int i=0;i!=nwire_w;i++){
    hw[i]->Reset();
  }
  
  fds.jump(frame_number);
  
  //fill in the data ... 
  const Frame& frame = fds.get();
  size_t ntraces = frame.traces.size();
  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();

    
    TH1F *htemp;
    if (chid < nwire_u){
      htemp = hu[chid];
    }else if (chid < nwire_u + nwire_v){
      htemp = hv[chid - nwire_u];
    }else{
      htemp = hw[chid - nwire_u - nwire_v];
    }

    for (int j = 0; j!= nbins; j++){
      float charge = htemp->GetBinContent(tbin + 1 + j);
      charge += trace.charge.at(j);
      htemp->SetBinContent(tbin+1+j,charge);
    }
    //std::cout << chid << std::endl;
    // std::cout << nwire_u << " " << nwire_v << " " << nwire_w << std::endl;
  }
  
  // do FFT to convolute with response function
  
  // add in random noise
  
  // do FFT again to remove response function and apply filter

  // correct the baseline 

  // fill the frame data ... 


  return 0;
}

WireCell2dToy::ToySignalSimuFDS::~ToySignalSimuFDS(){
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
