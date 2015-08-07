#include "WireCell2dToy/ToySignalSimu.h"

#include "WireCellData/GeomWire.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "TRandom.h"
#include "TF1.h"

using namespace WireCell;

WireCell2dToy::ToySignalSimuFDS::ToySignalSimuFDS(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds,int bins_per_frame1, int nframes_total, float time_offset_uv, float time_offset_uw, int flag_random, float overall_time_offset, int overall_time_shift)
  : fds(fds)
  , gds(gds)
  , max_frames(nframes_total)
  , time_offset_uv(time_offset_uv)
  , time_offset_uw(time_offset_uw)
  , flag_random(flag_random)
{  
  bins_per_frame = bins_per_frame1;

  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  //test save
  hu1 = new TH1F*[nwire_u];
  hv1 = new TH1F*[nwire_v];
  hw1 = new TH1F*[nwire_w];
  
  for (int i=0;i!=nwire_u;i++){
    hu1[i] = new TH1F(Form("U_%d",i),Form("U_%d",i),bins_per_frame,0,bins_per_frame);
  }
  for (int i=0;i!=nwire_v;i++){
    hv1[i] = new TH1F(Form("V_%d",i),Form("V_%d",i),bins_per_frame,0,bins_per_frame);
  }
  for (int i=0;i!=nwire_w;i++){
    hw1[i] = new TH1F(Form("W_%d",i),Form("W_%d",i),bins_per_frame,0,bins_per_frame);
  }
  //test save
  
  hu = new TH1F("U","U",bins_per_frame,0,bins_per_frame);
  hv = new TH1F("V","V",bins_per_frame,0,bins_per_frame);
  hw = new TH1F("W","W",bins_per_frame,0,bins_per_frame);

  #include "data.txt"

  gu = new TGraph(5000,xu,yu);
  gv = new TGraph(5000,xv,yv);
  gw = new TGraph(5000,xw,yw);

  hur = new TH1F("hur","hur",bins_per_frame,0,bins_per_frame); // half us tick
  hvr = new TH1F("hvr","hvr",bins_per_frame,0,bins_per_frame); // half us tick
  hwr = new TH1F("hwr","hwr",bins_per_frame,0,bins_per_frame); // half us tick
  
  for (int i=0; i!=bins_per_frame; i++){  
    double time = hur->GetBinCenter(i+1)/2.-50 ;
    hur->SetBinContent(i+1,gu->Eval(time-overall_time_offset));
    hvr->SetBinContent(i+1,gv->Eval(time-time_offset_uv-overall_time_offset));
    hwr->SetBinContent(i+1,gw->Eval(time-time_offset_uw-overall_time_offset));
  } 
  
  delete gu;
  delete gv;
  delete gw;
  
}

int WireCell2dToy::ToySignalSimuFDS::size() const{
  return max_frames;
}

void WireCell2dToy::ToySignalSimuFDS::Save(){
  TFile *file = new TFile("temp_simu.root","RECREATE");
  //test save
  for (int i=0;i!=nwire_u;i++){
    TH1F *huu = (TH1F*)hu1[i]->Clone(Form("U1_%d",i));
  }
  for (int i=0;i!=nwire_v;i++){
    TH1F *hvv = (TH1F*)hv1[i]->Clone(Form("V1_%d",i));
  }
  for (int i=0;i!=nwire_w;i++){
    TH1F *hww = (TH1F*)hw1[i]->Clone(Form("W1_%d",i));
  }
  //test save
  file->Write();
  file->Close();
}

int WireCell2dToy::ToySignalSimuFDS::jump(int frame_number){
  // do simulation
  
  //test save
  for (int i=0;i!=nwire_u;i++){
    hu1[i]->Reset();
  }
  for (int i=0;i!=nwire_v;i++){
    hv1[i]->Reset();
  }
  for (int i=0;i!=nwire_w;i++){
    hw1[i]->Reset();
  }
  //test save
  
  if (frame.index == frame_number) {
    return frame_number;
  }
  
  fds.jump(frame_number);
  
   // start FFT to convolute with response function
  TVirtualFFT::SetTransform(0);
  TH1 *hm = 0;
  TH1 *hp = 0;
  TH1 *hmr;
  TH1 *hpr;

  TH1 *hmr_u = 0;
  TH1 *hpr_u = 0;
  TH1 *hmr_v = 0;
  TH1 *hpr_v = 0;
  TH1 *hmr_w = 0;
  TH1 *hpr_w = 0;
  
  double value_re[9600]; // hack for now
  double value_im[9600];
  int  n  = bins_per_frame;
  TVirtualFFT *ifft;
  TH1 *fb = 0;

  // add in random noise
  double noise[2]={0.48,0.6};
  for (int i=0;i!=2;i++){
    noise[i] = 7.8/4.7;
  }

  hmr_u = hur->FFT(0,"MAG");  
  hpr_u = hur->FFT(0,"PH");
  hmr_v = hvr->FFT(0,"MAG");  
  hpr_v = hvr->FFT(0,"PH");
  hmr_w = hwr->FFT(0,"MAG");  
  hpr_w = hwr->FFT(0,"PH");
  
  frame.clear();
  
  //fill in the data ... 
  const Frame& frame1 = fds.get();
  size_t ntraces = frame1.traces.size();
  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();

    
    TH1F *htemp;
    if (chid < nwire_u){
      //      htemp = hu;
      htemp = hu1[chid];
      hmr = hmr_u;
      hpr = hpr_u;
    }else if (chid < nwire_u + nwire_v){
      //htemp = hv;
      htemp = hv1[chid - nwire_u];
      hmr = hmr_v;
      hpr = hpr_v;
    }else{
      //htemp = hw;
      htemp = hw1[chid - nwire_u-nwire_v];
      hmr = hmr_w;
      hpr = hpr_w;
    }

    for (int j = 0; j!= nbins; j++){
      float charge = htemp->GetBinContent(tbin + 1 + j);
      charge += trace.charge.at(j);
      htemp->SetBinContent(tbin+1+j,charge);  
    }
    
    std::vector<double> vcharge;
    for (int j=0;j!=htemp->GetNbinsX();j++){
      vcharge.push_back(htemp->GetBinContent(j+1));
    }
    htemp->Reset();
    for (int j=0;j!=htemp->GetNbinsX();j++){
      int tt = j+1;//+3200;
      if (tt <= bins_per_frame)
	htemp->SetBinContent(tt,vcharge.at(j));
    }

    
    hm = htemp->FFT(0,"MAG");
    hp = htemp->FFT(0,"PH");
    for (int j=0;j!=bins_per_frame;j++){
      double rho = hm->GetBinContent(j+1)*hmr->GetBinContent(j+1);
      double phi = hp->GetBinContent(j+1) + hpr->GetBinContent(j+1);
      value_re[j] = rho*cos(phi)/bins_per_frame;
      value_im[j] = rho*sin(phi)/bins_per_frame;
    }
    ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    fb = TH1::TransformHisto(ifft,0,"Re");
    for (int j=0;j!=bins_per_frame;j++){
      int content;
      if (flag_random ==1){
	content = round(fb->GetBinContent(j+1) * 7.8*4096./2000. + gRandom->Gaus(0,noise[1]));
      }else{
	content = round(fb->GetBinContent(j+1) * 7.8*4096./2000.);
      }
      htemp->SetBinContent(j+1,content);
    }
    delete hm;
    delete hp;
    delete fb;
    delete ifft;

    Trace t;
    t.chid = chid;
    t.tbin = 0;
    t.charge.resize(bins_per_frame, 0.0);
    for (int j=0;j!=bins_per_frame;j++){
      t.charge.at(j) = htemp->GetBinContent(j+1);
    }
    frame.traces.push_back(t);
    //std::cout << chid << std::endl;
    // std::cout << nwire_u << " " << nwire_v << " " << nwire_w << std::endl;
  }
  

  
 

  // //U-plane first
 
  // for (int i=0;i!=nwire_u;i++){
    
  // }
  // //V-plane
  
  // for (int i=0;i!=nwire_v;i++){
  //   hm = hv[i]->FFT(hm,"MAG");
  //   hp = hv[i]->FFT(hp,"PH");
  //   for (int j=0;j!=bins_per_frame;j++){
  //     double rho = hm->GetBinContent(j+1)*hmr->GetBinContent(j+1);
  //     double phi = hp->GetBinContent(j+1) + hpr->GetBinContent(j+1);
  //     value_re[j] = rho*cos(phi)/bins_per_frame;
  //     value_im[j] = rho*sin(phi)/bins_per_frame;
  //   }
  //   ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
  //   ifft->SetPointsComplex(value_re,value_im);
  //   ifft->Transform();
  //   fb = TH1::TransformHisto(ifft,fb,"Re");
  //   for (int j=0;j!=bins_per_frame;j++){
  //     int content;
  //     if (flag_random == 1){
  // 	content = round(fb->GetBinContent(j+1)*7.8*4096./2000. + gRandom->Gaus(0,noise[1]));
  //     }else{
  // 	content = round(fb->GetBinContent(j+1)*7.8*4096./2000.);
  //     }
  //     hv[i]->SetBinContent(j+1,content);
  //   }
  // }
  // //W-plane
  
  // for (int i=0;i!=nwire_w;i++){
  //   hm = hw[i]->FFT(hm,"MAG");
  //   hp = hw[i]->FFT(hp,"PH");
  //   for (int j=0;j!=bins_per_frame;j++){
  //     double rho = hm->GetBinContent(j+1)*hmr->GetBinContent(j+1);
  //     double phi = hp->GetBinContent(j+1) + hpr->GetBinContent(j+1);
  //     value_re[j] = rho*cos(phi)/bins_per_frame;
  //     value_im[j] = rho*sin(phi)/bins_per_frame;
  //   }
  //   ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
  //   ifft->SetPointsComplex(value_re,value_im);
  //   ifft->Transform();
  //   fb = TH1::TransformHisto(ifft,fb,"Re");
  //   for (int j=0;j!=bins_per_frame;j++){
  //     int content;
  //     if (flag_random ==1){
  // 	content = round(fb->GetBinContent(j+1)*7.8*4096./2000. + gRandom->Gaus(0,noise[0]));
  //     }else{
  // 	content = round(fb->GetBinContent(j+1)*7.8*4096./2000. );
  //     }
  //     hw[i]->SetBinContent(j+1,content);
  //   }
  // }
  // // done with FFT


  // fill the frame data ... 
  
  
  // //U-plane
  // for (int i=0;i!=nwire_u;i++){
    
  // }
  // //V-plane
  // for (int i=0;i!=nwire_v;i++){
  //   Trace t;
  //   t.chid = i+nwire_u;
  //   t.tbin = 0;
  //   t.charge.resize(bins_per_frame, 0.0);
  //   for (int j=0;j!=bins_per_frame;j++){
  //     t.charge.at(j) = hv[i]->GetBinContent(j+1);
  //   }
  //   frame.traces.push_back(t);
  // }
  // //W-plane
  // for (int i=0;i!=nwire_w;i++){
  //   Trace t;
  //   t.chid = i + nwire_u + nwire_v;
  //   t.tbin = 0;
  //   t.charge.resize(bins_per_frame, 0.0);
  //   for (int j=0;j!=bins_per_frame;j++){
  //     t.charge.at(j) = hw[i]->GetBinContent(j+1);
  //   }
  //   frame.traces.push_back(t);
  // }
  
  
  frame.index = frame_number;
  return frame.index;
}

WireCell2dToy::ToySignalSimuFDS::~ToySignalSimuFDS(){
  //test save
  for (int i=0;i!=nwire_u;i++){
    delete hu1[i] ;
  }
   for (int i=0;i!=nwire_w;i++){
    delete hw1[i] ;
  }
   for (int i=0;i!=nwire_v;i++){
    delete hv1[i] ;
  }
  //test save

  delete hu;
  delete hv;
  delete hw;

  // delete gu;
  // delete gv;
  // delete gw;

  delete hur;
  delete hvr;
  delete hwr;
}
