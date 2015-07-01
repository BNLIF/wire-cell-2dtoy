#include "WireCell2dToy/ToySignalSimu.h"

#include "WireCellData/GeomWire.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "TRandom.h"
#include "TF1.h"

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

  gu = new TGraph(5000,xu,yu);
  gv = new TGraph(5000,xv,yv);
  gw = new TGraph(5000,xw,yw);

  hur = new TH1F("hur","hur",bins_per_frame,0,bins_per_frame); // half us tick
  hvr = new TH1F("hvr","hvr",bins_per_frame,0,bins_per_frame); // half us tick
  hwr = new TH1F("hwr","hwr",bins_per_frame,0,bins_per_frame); // half us tick
  
  for (int i=0; i!=bins_per_frame; i++){  
    double time = hur->GetBinCenter(i+1)/2.-50;
    hur->SetBinContent(i+1,gu->Eval(time));
    hvr->SetBinContent(i+1,gv->Eval(time));
    hwr->SetBinContent(i+1,gw->Eval(time));
  } 
  
  
  
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
  TVirtualFFT::SetTransform(0);
  TH1 *hm = 0;
  TH1 *hp = 0;
  TH1 *hmr = 0;
  TH1 *hpr = 0;

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

  //U-plane first
  hmr = hur->FFT(hmr,"MAG");  
  hpr = hur->FFT(hpr,"PH");
  for (int i=0;i!=nwire_u;i++){
    hm = hu[i]->FFT(hm,"MAG");
    hp = hu[i]->FFT(hp,"PH");
    for (int j=0;j!=bins_per_frame;j++){
      double rho = hm->GetBinContent(j+1)*hmr->GetBinContent(j+1);
      double phi = hp->GetBinContent(j+1) + hpr->GetBinContent(j+1);
      value_re[j] = rho*cos(phi)/bins_per_frame;
      value_im[j] = rho*sin(phi)/bins_per_frame;
    }
    ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    fb = TH1::TransformHisto(ifft,fb,"Re");
    for (int j=0;j!=bins_per_frame;j++){
      int content = fb->GetBinContent(j+1) * 7.8*4096./2000. + gRandom->Gaus(0,noise[1]);
      hu[i]->SetBinContent(j+1,content);
    }
  }
  //V-plane
  hmr = hvr->FFT(hmr,"MAG");  
  hpr = hvr->FFT(hpr,"PH");
  for (int i=0;i!=nwire_v;i++){
    hm = hv[i]->FFT(hm,"MAG");
    hp = hv[i]->FFT(hp,"PH");
    for (int j=0;j!=bins_per_frame;j++){
      double rho = hm->GetBinContent(j+1)*hmr->GetBinContent(j+1);
      double phi = hp->GetBinContent(j+1) + hpr->GetBinContent(j+1);
      value_re[j] = rho*cos(phi)/bins_per_frame;
      value_im[j] = rho*sin(phi)/bins_per_frame;
    }
    ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    fb = TH1::TransformHisto(ifft,fb,"Re");
    for (int j=0;j!=bins_per_frame;j++){
      int content = fb->GetBinContent(j+1)*7.8*4096./2000. + gRandom->Gaus(0,noise[1]);
      hv[i]->SetBinContent(j+1,content);
    }
  }
  //W-plane
  hmr = hwr->FFT(hmr,"MAG");  
  hpr = hwr->FFT(hpr,"PH");
  for (int i=0;i!=nwire_w;i++){
    hm = hw[i]->FFT(hm,"MAG");
    hp = hw[i]->FFT(hp,"PH");
    for (int j=0;j!=bins_per_frame;j++){
      double rho = hm->GetBinContent(j+1)*hmr->GetBinContent(j+1);
      double phi = hp->GetBinContent(j+1) + hpr->GetBinContent(j+1);
      value_re[j] = rho*cos(phi)/bins_per_frame;
      value_im[j] = rho*sin(phi)/bins_per_frame;
    }
    ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    fb = TH1::TransformHisto(ifft,fb,"Re");
    for (int j=0;j!=bins_per_frame;j++){
      int content = fb->GetBinContent(j+1)*7.8*4096./2000. + gRandom->Gaus(0,noise[0]);
      hw[i]->SetBinContent(j+1,content);
    }
  }
    
  
  // do FFT again to remove response function and apply filter

  TF1 *filter_u = new TF1("filter_u","(x>0.0)*gaus*exp(-0.5*pow(x/[3],[4]))");
  double par[5]={1.73/0.959301, 1.69, 1.55, 0.19, 3.75};
  filter_u->SetParameters(par);

  TF1 *filter_v = new TF1("filter_v","(x>0.0)*gaus*exp(-0.5*pow(x/[3],[4]))");
  double par1[5]={1.74/0.941034, 1.46, 1.33, 0.23, 4.89};
  filter_v->SetParameters(par1);

  TF1 *filter_w = new TF1("filter_y","(x>0.0)*[0]*exp(-0.5*(((x-[1])/[2])^2)^[3])");
  double par2[4]={1.03/0.995635, 0.08, 0.15, 2.17};
  filter_w->SetParameters(par2);

  //U-plane first
  hmr = hur->FFT(hmr,"MAG");  
  hpr = hur->FFT(hpr,"PH");
  for (int i=0;i!=nwire_u;i++){
    hm = hu[i]->FFT(hm,"MAG");
    hp = hu[i]->FFT(hp,"PH");
    for (int j=0;j!=bins_per_frame;j++){
      double freq;
      if (j < bins_per_frame/2.){
	freq = j/(1.*bins_per_frame)*2.;
      }else{
	freq = (bins_per_frame - j)/(1.*bins_per_frame)*2.;
      }

      double rho = hm->GetBinContent(j+1)/hmr->GetBinContent(j+1)*filter_u->Eval(freq);
      double phi = hp->GetBinContent(j+1) - hpr->GetBinContent(j+1);
      value_re[j] = rho*cos(phi)/bins_per_frame;
      value_im[j] = rho*sin(phi)/bins_per_frame;

      //std::cout << i << rho << " " << phi << " " << freq << std::endl;
    }
    ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    fb = TH1::TransformHisto(ifft,fb,"Re");
    for (int j=0;j!=bins_per_frame;j++){
      double content = fb->GetBinContent(j+1) /( 7.8*4096./2000.);
      hu[i]->SetBinContent(j+1,content);
    }
  }
  //V-plane
  hmr = hvr->FFT(hmr,"MAG");  
  hpr = hvr->FFT(hpr,"PH");
  for (int i=0;i!=nwire_v;i++){
    hm = hv[i]->FFT(hm,"MAG");
    hp = hv[i]->FFT(hp,"PH");
    for (int j=0;j!=bins_per_frame;j++){
      double freq;
      if (j < bins_per_frame/2.){
	freq = j/(1.*bins_per_frame)*2.;
      }else{
	freq = (bins_per_frame - j)/(1.*bins_per_frame)*2.;
      }

      double rho = hm->GetBinContent(j+1)/hmr->GetBinContent(j+1)*filter_v->Eval(freq);
      double phi = hp->GetBinContent(j+1) - hpr->GetBinContent(j+1);
      value_re[j] = rho*cos(phi)/bins_per_frame;
      value_im[j] = rho*sin(phi)/bins_per_frame;
    }
    ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    fb = TH1::TransformHisto(ifft,fb,"Re");
    for (int j=0;j!=bins_per_frame;j++){
      double content = fb->GetBinContent(j+1) /( 7.8*4096./2000.);
      hv[i]->SetBinContent(j+1,content);
    }
  }
  //W-plane
  hmr = hwr->FFT(hmr,"MAG");  
  hpr = hwr->FFT(hpr,"PH");
  for (int i=0;i!=nwire_w;i++){
    hm = hw[i]->FFT(hm,"MAG");
    hp = hw[i]->FFT(hp,"PH");
    for (int j=0;j!=bins_per_frame;j++){
      double freq;
      if (j < bins_per_frame/2.){
	freq = j/(1.*bins_per_frame)*2.;
      }else{
	freq = (bins_per_frame - j)/(1.*bins_per_frame)*2.;
      }

      double rho = hm->GetBinContent(j+1)/hmr->GetBinContent(j+1)*filter_w->Eval(freq);
      double phi = hp->GetBinContent(j+1) - hpr->GetBinContent(j+1);
      value_re[j] = rho*cos(phi)/bins_per_frame;
      value_im[j] = rho*sin(phi)/bins_per_frame;
    }
    ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    fb = TH1::TransformHisto(ifft,fb,"Re");
    for (int j=0;j!=bins_per_frame;j++){
      double content = fb->GetBinContent(j+1) /( 7.8*4096./2000.);
      hw[i]->SetBinContent(j+1,content);
    }
  }
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

  delete gu;
  delete gv;
  delete gw;

  delete hur;
  delete hvr;
  delete hwr;
}
