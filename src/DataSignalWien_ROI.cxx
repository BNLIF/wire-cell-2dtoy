#include "WireCell2dToy/DataSignalWien_ROI.h"

#include "WireCellData/GeomWire.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "TRandom.h"
#include "TF1.h"

using namespace WireCell;

WireCell2dToy::DataSignalWienROIFDS::DataSignalWienROIFDS(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds,WireCell::ChirpMap& umap, WireCell::ChirpMap& vmap, WireCell::ChirpMap& wmap, int bins_per_frame1, int nframes_total, float time_offset_uv, float time_offset_uw, float overall_time_offset)
  : fds(fds)
  , gds(gds)
  , max_frames(nframes_total)
  , time_offset_uv(time_offset_uv)
  , time_offset_uw(time_offset_uw)
  , overall_time_offset(overall_time_offset)
  , umap(umap)
  , vmap(vmap)
  , wmap(wmap)
{  
  bins_per_frame = bins_per_frame1;

  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  nbin = fds.Get_Bins_Per_Frame();

  uplane_rms.resize(nwire_u,0);
  vplane_rms.resize(nwire_v,0);
  wplane_rms.resize(nwire_w,0);

  uplane_rms_g.resize(nwire_u,0);
  vplane_rms_g.resize(nwire_v,0);
  wplane_rms_g.resize(nwire_w,0);

  hu_1D_c = new TH2F("hu_1D_c","hu_1D_c",nwire_u,0,nwire_u,bins_per_frame,0,bins_per_frame);
  hu_2D_g_f = new TH2F("hu_2D_g_f","hu_2D_g_f",nwire_u,0,nwire_u,bins_per_frame,0,bins_per_frame);
  hu_2D_g = new TH2F("hu_2D_g","hu_2D_g",nwire_u,0,nwire_u,bins_per_frame,0,bins_per_frame);
  
  hv_1D_c = new TH2F("hv_1D_c","hv_1D_c",nwire_v,0,nwire_v,bins_per_frame,0,bins_per_frame);
  hv_1D_g_f = new TH2F("hv_1D_g_f","hv_1D_g_f",nwire_v,0,nwire_v,bins_per_frame,0,bins_per_frame);
  hv_1D_g = new TH2F("hv_1D_g","hv_1D_g",nwire_v,0,nwire_v,bins_per_frame,0,bins_per_frame);

  hw_1D_g = new TH2F("hw_1D_g","hw_1D_g",nwire_w,0,nwire_w,bins_per_frame,0,bins_per_frame);

  #include "data_70_ROI.txt"  //70kV 2D deconvolution for U

  gu_1D_c = new TGraph(5000,u_1D_c_x, u_1D_c_y);
  gv_1D_c = new TGraph(5000,v_1D_c_x, v_1D_c_y);
  gw_1D_c = new TGraph(5000,w_1D_c_x, w_1D_c_y);
  

  gu_2D_g = new TGraph*[4];
  gu_2D_g[0] = new TGraph(5000,u_2D_g_0_x,u_2D_g_0_y);
  gu_2D_g[1] = new TGraph(5000,u_2D_g_1_x,u_2D_g_1_y);
  gu_2D_g[2] = new TGraph(5000,u_2D_g_2_x,u_2D_g_2_y);
  gu_2D_g[3] = new TGraph(5000,u_2D_g_3_x,u_2D_g_3_y);

  gv_1D_g = new TGraph(5000,v_1D_g_x,v_1D_g_y);
  gw_1D_g = new TGraph(5000,w_1D_g_x,w_1D_g_y);

  
}

WireCell2dToy::DataSignalWienROIFDS::~DataSignalWienROIFDS(){
  delete hu_1D_c;
  delete hu_2D_g_f;
  delete hu_2D_g;
  
  delete hv_1D_c;
  delete hv_1D_g_f;
  delete hv_1D_g;
  
  delete hw_1D_g;

  delete gu_1D_c;
  delete gv_1D_c;
  delete gw_1D_c;
 
  delete gv_1D_g;
  delete gw_1D_g;
  
  for (int i=0;i!=4;i++){
    delete gu_2D_g[i];
  }
  
  delete [] gu_2D_g;
}



int WireCell2dToy::DataSignalWienROIFDS::size() const{
  return max_frames;
}

void WireCell2dToy::DataSignalWienROIFDS::Save(){
}


void WireCell2dToy::DataSignalWienROIFDS::Deconvolute_U_2D_g(){
  const Frame& frame1 = fds.get();
  size_t ntraces = frame1.traces.size();

  const int nticks = nbin;
  double value_re[9600];
  double value_im[9600];

  int scale = nbin/bins_per_frame;

  // filter 
  TF1 *filter_u = new TF1("filter_u","(x>0.0)*exp(-0.5*pow(x/[0],[1]))");
  double par[2]={1.43555e+01/200.*2.,4.95096e+00};
  filter_u->SetParameters(par);

  TF1 *filter_low = new TF1("filter_low","1-exp(-pow(x/0.0045,2))");

  //response function ... 

  // deconvolution ... 

  delete filter_u;
  delete filter_low;
}

void WireCell2dToy::DataSignalWienROIFDS::Deconvolute_V_1D_g(){
  const Frame& frame1 = fds.get();
  size_t ntraces = frame1.traces.size();

  const int nticks = nbin;
  double value_re[9600];
  double value_im[9600];

  double value_re1[9600];
  double value_im1[9600];
  
  int scale = nbin/bins_per_frame;

  //filter ... 
  TF1 *filter_v = new TF1("filter_v","(x>0.0)*exp(-0.5*pow(x/[0],[1]))");
  double par1[2]={1.47404e+01/200.*2.,4.97667e+00};
  filter_v->SetParameters(par1);
  
  TF1 *filter_low = new TF1("filter_low","1-exp(-pow(x/0.0045,2))");

  //response function 
  TH1F *hvr = new TH1F("hvr1","hvr1",nbin,0,nbin); // half us tick
  float scale_v = 1.251/1.074*0.91*0.85;
  for (int i=0;i!=nbin;i++){
    double time  = hvr->GetBinCenter(i+1)/2.-50;
    double x = time-time_offset_uv;
    if (x > -35 && x  < 15){
      hvr->SetBinContent(i+1,gv_1D_g->Eval(x)/scale_v);
    }
  }
  
  TH1 *hmr_v = hvr->FFT(0,"MAG");
  TH1 *hpr_v = hvr->FFT(0,"PH");

  
  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    if (chid <nwire_u || chid >=nwire_u+nwire_v) continue;
    
    TH1F *htemp = new TH1F("htemp","htemp",nbin,0,nbin);
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp->SetBinContent(tt,trace.charge.at(i));
    }
    
    TH1 *hm = htemp->FFT(0,"MAG");
    TH1 *hp = htemp->FFT(0,"PH");
    
    for (int i=0;i!=nbin;i++){
      double freq;
      if (i< nbin/2.){
	freq = i/(1.*nbin)*2.;
      }else{
	freq = (nbin - i)/(1.*nbin)*2.;
      }
      double rho = hm->GetBinContent(i+1)/hmr_v->GetBinContent(i+1)*filter_v->Eval(freq);
      double phi = hp->GetBinContent(i+1) - hpr_v->GetBinContent(i+1);
      if (i==0) rho = 0;
      value_re[i] = rho*cos(phi)/nbin;
      value_im[i] = rho*sin(phi)/nbin;
      
      rho *= filter_low->Eval(freq);
      value_re1[i] = rho*cos(phi)/nbin;
      value_im1[i] = rho*sin(phi)/nbin;
    }

    TVirtualFFT *ifft = TVirtualFFT::FFT(1,&nbin,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    TH1 *fb = TH1::TransformHisto(ifft,0,"Re");
    
    // put results back
    for (int i=0;i!=nbin;i++){
      htemp->SetBinContent(i+1,fb->GetBinContent(i+1)/( 14.*1.1*4096./2000.));
    }
    delete ifft;
    delete fb;

    // correct baseline 
    restore_baseline(htemp);

    // put results back into the 2-D histogram
    for (int i=0;i!=bins_per_frame;i++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*i+k+1);
      }
      hv_1D_g->SetBinContent(chid+1-nwire_u,i+1,sum);
    }


    
    ifft = TVirtualFFT::FFT(1,&nbin,"C2R M K");
    ifft->SetPointsComplex(value_re1,value_im1);
    ifft->Transform();
    fb = TH1::TransformHisto(ifft,0,"Re");
    
    // put results back
    for (int i=0;i!=nbin;i++){
      htemp->SetBinContent(i+1,fb->GetBinContent(i+1)/( 14.*1.1*4096./2000.));
    }
    delete ifft;
    delete fb;

    // correct baseline 
    restore_baseline(htemp);

    // calculate RMS 
    double rms = cal_rms(htemp,chid);
    vplane_rms_g[chid-nwire_u] = rms*scale;
    
    // put results back into the 2-D histogram
    for (int i=0;i!=bins_per_frame;i++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*i+k+1);
      }
      hv_1D_g_f->SetBinContent(chid+1-nwire_u,i+1,sum);
    }

    

    delete htemp;
    delete hm;
    delete hp;
    
  }
  

  delete hvr;
  delete hmr_v;
  delete hpr_v;
  delete filter_v;
  delete filter_low;

}

void WireCell2dToy::DataSignalWienROIFDS::Deconvolute_W_1D_g(){
  const Frame& frame1 = fds.get();
  size_t ntraces = frame1.traces.size();

  const int nticks = nbin;
  double value_re[9600];
  double value_im[9600];
  int scale = nbin/bins_per_frame;

  // wiener filter for the W-plane Gaussian ... 
  TF1 *filter_w = new TF1("filter_y","(x>0.0)*exp(-0.5*pow(x/[0],[1]))");
  double par2[2]={1.45874e+01/200.*2.,5.02219e+00};
  filter_w->SetParameters(par2);
  
  //response function
  TH1F *hwr = new TH1F("hwr1","hwr1",nbin,0,nbin); // half us tick
  for (int i=0;i!=nbin;i++){
    double time  = hwr->GetBinCenter(i+1)/2.-50;
    double x = time-time_offset_uw;
    if (x > -35 && x  < 15){
      hwr->SetBinContent(i+1,gw_1D_g->Eval(x));
    }
  }
  
  TH1 *hmr_w = hwr->FFT(0,"MAG");
  TH1 *hpr_w = hwr->FFT(0,"PH");

   for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    if (chid < nwire_u+nwire_v) continue;
    
    TH1F *htemp = new TH1F("htemp","htemp",nbin,0,nbin);
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp->SetBinContent(tt,trace.charge.at(i));
    }
    
    TH1 *hm = htemp->FFT(0,"MAG");
    TH1 *hp = htemp->FFT(0,"PH");

    for (int i=0;i!=nbin;i++){
      double freq;
      if (i< nbin/2.){
	freq = i/(1.*nbin)*2.;
      }else{
	freq = (nbin - i)/(1.*nbin)*2.;
      }
      double rho = hm->GetBinContent(i+1)/hmr_w->GetBinContent(i+1)*filter_w->Eval(freq);
      double phi = hp->GetBinContent(i+1) - hpr_w->GetBinContent(i+1);
      if (i==0) rho = 0;
      
      value_re[i] = rho*cos(phi)/nbin;
      value_im[i] = rho*sin(phi)/nbin;
    }

    TVirtualFFT *ifft = TVirtualFFT::FFT(1,&nbin,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    TH1 *fb = TH1::TransformHisto(ifft,0,"Re");
    
    // put results back
    for (int i=0;i!=nbin;i++){
      htemp->SetBinContent(i+1,fb->GetBinContent(i+1)/( 14.*1.1*4096./2000.));
    }
    delete ifft;
    delete fb;

    // correct baseline 
    restore_baseline(htemp);

    // calculate RMS 
    double rms = cal_rms(htemp,chid);
    wplane_rms[chid-nwire_u-nwire_v] = rms*scale;
    
    // put results back into the 2-D histogram
    for (int i=0;i!=bins_per_frame;i++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*i+k+1);
      }
      hw_1D_g->SetBinContent(chid+1-nwire_u-nwire_v,i+1,sum);
    }
    delete htemp;
    delete hm;
    delete hp;
    
  }


  delete filter_w;
  delete hmr_w;
  delete hpr_w;
  delete hwr;
}


void WireCell2dToy::DataSignalWienROIFDS::Deconvolute_V_1D_c(){
  const Frame& frame1 = fds.get();
  size_t ntraces = frame1.traces.size();

  const int nticks = nbin;
  double value_re[9600];
  double value_im[9600];
  int scale = nbin/bins_per_frame;

  // wiener filter for V-plane
  TF1 *filter_v = new TF1("filter_v","(x>0.0)*exp(-0.5*pow(x/[0],[1]))");
  double par1[2]={1.31649e+01/200.*2.,4.67294e+00};
  filter_v->SetParameters(par1);

  //response function
  TH1F *hvr = new TH1F("hvr1","hvr1",nbin,0,nbin); // half us tick
  float scale_v = 1.251/1.074*0.91*0.85;
  for (int i=0;i!=nbin;i++){
    double time  = hvr->GetBinCenter(i+1)/2.-50;
    double x = time-time_offset_uv;
    if (x > -35 && x  < 15){
      hvr->SetBinContent(i+1,gv_1D_c->Eval(x)/scale_v);
    }
  }
  
  TH1 *hmr_v = hvr->FFT(0,"MAG");
  TH1 *hpr_v = hvr->FFT(0,"PH");

  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    if (chid <nwire_u || chid >=nwire_u+nwire_v) continue;
    
    TH1F *htemp = new TH1F("htemp","htemp",nbin,0,nbin);
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp->SetBinContent(tt,trace.charge.at(i));
    }
    
    TH1 *hm = htemp->FFT(0,"MAG");
    TH1 *hp = htemp->FFT(0,"PH");

    for (int i=0;i!=nbin;i++){
      double freq;
      if (i< nbin/2.){
	freq = i/(1.*nbin)*2.;
      }else{
	freq = (nbin - i)/(1.*nbin)*2.;
      }
      double rho = hm->GetBinContent(i+1)/hmr_v->GetBinContent(i+1)*filter_v->Eval(freq);
      double phi = hp->GetBinContent(i+1) - hpr_v->GetBinContent(i+1);
      if (i==0) rho = 0;
      
      value_re[i] = rho*cos(phi)/nbin;
      value_im[i] = rho*sin(phi)/nbin;
    }

    TVirtualFFT *ifft = TVirtualFFT::FFT(1,&nbin,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    TH1 *fb = TH1::TransformHisto(ifft,0,"Re");
    
    // put results back
    for (int i=0;i!=nbin;i++){
      htemp->SetBinContent(i+1,fb->GetBinContent(i+1)/( 14.*1.1*4096./2000.));
    }
    delete ifft;
    delete fb;

    // correct baseline 
    restore_baseline(htemp);

    // calculate RMS 
    double rms = cal_rms(htemp,chid);
    vplane_rms[chid-nwire_u] = rms*scale;
    
    // put results back into the 2-D histogram
    for (int i=0;i!=bins_per_frame;i++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*i+k+1);
      }
      hv_1D_c->SetBinContent(chid+1-nwire_u,i+1,sum);
    }
    delete htemp;
    delete hm;
    delete hp;
    
  }

  delete hvr;
  delete hmr_v;
  delete hpr_v;
  delete filter_v;
}


void WireCell2dToy::DataSignalWienROIFDS::Deconvolute_U_1D_c(){
  
  const Frame& frame1 = fds.get();
  size_t ntraces = frame1.traces.size();

  const int nticks = nbin;
  double value_re[9600];
  double value_im[9600];
  int scale = nbin/bins_per_frame;

  // wiener filter for U-plane
  TF1 *filter_u = new TF1("filter_u","(x>0.0)*exp(-0.5*pow(x/[0],[1]))");
  double par[2]={1.22781e+01/200.*2.,4.96159e+00};
  filter_u->SetParameters(par);
  
  //response function ...
  TH1F *hur = new TH1F("hur1","hur1",nbin,0,nbin); // half us tick
  float scale_u = 1.51/1.16*0.91*0.85;
  for (int i=0;i!=nbin;i++){
    double time  = hur->GetBinCenter(i+1)/2.-50;
    double x = time;
    if (x > -35 && x  < 15){
      hur->SetBinContent(i+1,gu_1D_c->Eval(x)/scale_u);
    }
  }
  TH1 *hmr_u = hur->FFT(0,"MAG");
  TH1 *hpr_u = hur->FFT(0,"PH");
  
  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    if (chid >=nwire_u) continue;
    
    TH1F *htemp = new TH1F("htemp","htemp",nbin,0,nbin);
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp->SetBinContent(tt,trace.charge.at(i));
    }
    
    TH1 *hm = htemp->FFT(0,"MAG");
    TH1 *hp = htemp->FFT(0,"PH");

    for (int i=0;i!=nbin;i++){
      double freq;
      if (i< nbin/2.){
	freq = i/(1.*nbin)*2.;
      }else{
	freq = (nbin - i)/(1.*nbin)*2.;
      }
      double rho = hm->GetBinContent(i+1)/hmr_u->GetBinContent(i+1)*filter_u->Eval(freq);
      double phi = hp->GetBinContent(i+1) - hpr_u->GetBinContent(i+1);
      if (i==0) rho = 0;
      
      value_re[i] = rho*cos(phi)/nbin;
      value_im[i] = rho*sin(phi)/nbin;
    }

    TVirtualFFT *ifft = TVirtualFFT::FFT(1,&nbin,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    TH1 *fb = TH1::TransformHisto(ifft,0,"Re");
    
    // put results back
    for (int i=0;i!=nbin;i++){
      htemp->SetBinContent(i+1,fb->GetBinContent(i+1)/( 14.*1.1*4096./2000.));
    }
    delete ifft;
    delete fb;

    // correct baseline 
    restore_baseline(htemp);

    // calculate RMS 
    double rms = cal_rms(htemp,chid);
    uplane_rms[chid] = rms*scale;
    
    // put results back into the 2-D histogram
    for (int i=0;i!=bins_per_frame;i++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*i+k+1);
      }
      hu_1D_c->SetBinContent(chid+1,i+1,sum);
    }


    delete htemp;
    delete hm;
    delete hp;
    
  }
  

  delete hur;
  delete hmr_u;
  delete hpr_u;
  delete filter_u;
}

double WireCell2dToy::DataSignalWienROIFDS::cal_rms(TH1F *htemp, int chid){
  //calculate rms, this is to be used for threshold purpose
  float rms = 0, rms1 = 0,rms2 = 0;
  int start=-1, end=-1;
  if (chid < nwire_u){
    if (umap.find(chid)!=umap.end()){
      start = umap[chid].first;
      end = umap[chid].second;
    }
  }else if (chid < nwire_u + nwire_v){
    if (vmap.find(chid-nwire_u)!=vmap.end()){
      start = vmap[chid-nwire_u].first;
      end = vmap[chid-nwire_u].second;
    }
  }else{
    if (wmap.find(chid-nwire_u-nwire_v)!=wmap.end()){
      start = wmap[chid-nwire_u-nwire_v].first;
      end = wmap[chid-nwire_u-nwire_v].second;
    }
  }
  
   // new method to calculate RMS
    int min1 =0,max1=0;
    for (int i=0;i!=htemp->GetNbinsX();i++){
      if (i < start || i > end){
    	if (htemp->GetBinContent(i+1)>max1)
    	  max1 = int(htemp->GetBinContent(i+1));
    	if (htemp->GetBinContent(i+1)<min1)
    	  min1 = int(htemp->GetBinContent(i+1));
      }
    }
    TH1F *h6 = new TH1F("h6","h6",int(max1-min1+1),min1,max1+1);
    for (int i=0;i!=htemp->GetNbinsX();i++){
      if (i < start || i > end){
    	h6->Fill(int(htemp->GetBinContent(i+1)));
      }
    }
    if (h6->GetSum()>0){
      //calculate 0.16, 0.84 percentile ...  
      double xq;
      xq = 0.16;
      double par[2];
      h6->GetQuantiles(1,&par[0],&xq);
      xq = 0.84;
      h6->GetQuantiles(1,&par[1],&xq);
      rms = (par[1]-par[0])/2.;
      
      //try to exclude signal
      rms2 = 0;
     for (int i=0;i!=htemp->GetNbinsX();i++){
       if (i < start || i > end){
    	 if (fabs(htemp->GetBinContent(i+1)) < 5.0*rms){
    	  rms1 += pow(htemp->GetBinContent(i+1),2);
    	  rms2 ++;
    	}
       }
     }
     if (rms2!=0){
       rms1 = sqrt(rms1/rms2);
     }else{
       rms1 = 0;   
     }
    }else{
      rms1 =0;
    }
    delete h6;

    return rms1;
}


void WireCell2dToy::DataSignalWienROIFDS::restore_baseline(TH1F *htemp){
  //correct baseline 
  double max = htemp->GetMaximum();
  double min = htemp->GetMinimum();
  int nbin_b = max - min;
  if (nbin_b ==0) nbin_b = 1;
  TH1F *h1 = new TH1F("h1","h1",nbin_b,min,max);
  for (int j=0;j!=nbin;j++){
    h1->Fill(htemp->GetBinContent(j+1));
  }
  float ped = h1->GetMaximumBin()*(max-min)/(nbin_b*1.) + min;
  float ave=0,ncount = 0;
  for (int j=0;j!=nbin;j++){
    if (fabs(htemp->GetBinContent(j+1)-ped)<400){
      ave +=htemp->GetBinContent(j+1);
      ncount ++;
      }
    }
    if (ncount==0) ncount=1;
    ave = ave/ncount;
    
    for (int j=0;j!=nbin;j++){
      double content = htemp->GetBinContent(j+1);
      content -= ave;
      htemp->SetBinContent(j+1,content);
    }
    delete h1;
}


int WireCell2dToy::DataSignalWienROIFDS::jump(int frame_number){
  // fill the frame data ... 
  if (frame.index == frame_number) {
    return frame_number;
  }
  frame.clear();
    
  fds.jump(frame_number);

  //std::cout << nbin << " " << bins_per_frame << std::endl; 

  TVirtualFFT::SetTransform(0);

  std::cout << "Deconvolution with calibrated field response for 1-D U Plane" << std::endl;
  Deconvolute_U_1D_c();
  std::cout << "Deconvolution with calibrated field response for 1-D V Plane" << std::endl;
  Deconvolute_V_1D_c();
  
  std::cout << "Deconvolution with garfield field response for 1-D W Plane" << std::endl;
  Deconvolute_W_1D_g();
 
  std::cout << "Deconvolution with garfield field response for 1-D V Plane" << std::endl;
  Deconvolute_V_1D_g();
 
  std::cout << "Load results back into frame" << std::endl;
  // load the results back into the frame ... 
  for (int i=0;i!=nwire_u;i++){
    Trace t;
    t.chid = i;
    t.tbin = 0;
    t.charge.resize(bins_per_frame, 0.0);
    for (int j=0;j!= bins_per_frame; j++){
      t.charge.at(j) = hv_1D_g->GetBinContent(i+1,j+1);
    }
    frame.traces.push_back(t);
  }

  for (int i=0;i!=nwire_v;i++){
    Trace t;
    t.chid = i+nwire_u;
    t.tbin = 0;
    t.charge.resize(bins_per_frame, 0.0);
    for (int j=0;j!= bins_per_frame; j++){
      t.charge.at(j) = hv_1D_g_f->GetBinContent(i+1,j+1);
    }
    frame.traces.push_back(t);
  }

  for (int i=0;i!=nwire_w;i++){
    Trace t;
    t.chid = i+nwire_u+nwire_v;
    t.tbin = 0;
    t.charge.resize(bins_per_frame, 0.0);
    for (int j=0;j!= bins_per_frame; j++){
      t.charge.at(j) = hw_1D_g->GetBinContent(i+1,j+1);
    }
    frame.traces.push_back(t);
  }

  
  frame.index = frame_number;
  return frame.index;
}


