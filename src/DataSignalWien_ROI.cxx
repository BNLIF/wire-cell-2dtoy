#include "WireCell2dToy/DataSignalWien_ROI.h"

#include "WireCellData/GeomWire.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "TRandom.h"
#include "TF1.h"
#include "TSpectrum.h"

#include <set>
#include <algorithm>

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
  
  hu_1D_c = new TH2I("hu_1D_c","hu_1D_c",nwire_u,0,nwire_u,bins_per_frame,0,bins_per_frame);
  hu_1D_c_gaus = new TH2I("hu_1D_c_gaus","hu_1D_c_gaus",nwire_u,0,nwire_u,bins_per_frame,0,bins_per_frame);

  hu_2D_g_f = new TH2I("hu_2D_g_f","hu_2D_g_f",nwire_u,0,nwire_u,bins_per_frame,0,bins_per_frame);
  hu_2D_g = new TH2I("hu_2D_g","hu_2D_g",nwire_u,0,nwire_u,bins_per_frame,0,bins_per_frame);
  hu_2D_g_gaus = new TH2I("hu_2D_g_gaus","hu_2D_g_gaus",nwire_u,0,nwire_u,bins_per_frame,0,bins_per_frame);
  
  hv_1D_c = new TH2I("hv_1D_c","hv_1D_c",nwire_v,0,nwire_v,bins_per_frame,0,bins_per_frame);
  hv_1D_c_gaus = new TH2I("hv_1D_c_gaus","hv_1D_c_gaus",nwire_v,0,nwire_v,bins_per_frame,0,bins_per_frame);
  
  hv_2D_g_f = new TH2I("hv_2D_g_f","hv_2D_g_f",nwire_v,0,nwire_v,bins_per_frame,0,bins_per_frame);
  hv_2D_g = new TH2I("hv_2D_g","hv_2D_g",nwire_v,0,nwire_v,bins_per_frame,0,bins_per_frame);
  hv_2D_g_gaus = new TH2I("hv_2D_g_gaus","hv_2D_g_gaus",nwire_v,0,nwire_v,bins_per_frame,0,bins_per_frame);

  hw_1D_g = new TH2I("hw_1D_g","hw_1D_g",nwire_w,0,nwire_w,bins_per_frame,0,bins_per_frame);
  hw_1D_g_gaus = new TH2I("hw_1D_g_gaus","hw_1D_g_gaus",nwire_w,0,nwire_w,bins_per_frame,0,bins_per_frame);

  #include "data_70_ROI.txt"  //70kV 2D deconvolution for U

  gu_1D_c = new TGraph(5000,u_1D_c_x, u_1D_c_y);
  gv_1D_c = new TGraph(5000,v_1D_c_x, v_1D_c_y);
  gw_1D_c = new TGraph(5000,w_1D_c_x, w_1D_c_y);
  

  gu_2D_g = new TGraph*[4];
  gu_2D_g[0] = new TGraph(5000,u_2D_g_0_x,u_2D_g_0_y);
  gu_2D_g[1] = new TGraph(5000,u_2D_g_1_x,u_2D_g_1_y);
  gu_2D_g[2] = new TGraph(5000,u_2D_g_2_x,u_2D_g_2_y);
  gu_2D_g[3] = new TGraph(5000,u_2D_g_3_x,u_2D_g_3_y);

  gv_2D_g = new TGraph*[4];
  gv_2D_g[0] = new TGraph(5000,v_2D_g_0_x,v_2D_g_0_y);
  gv_2D_g[1] = new TGraph(5000,v_2D_g_1_x,v_2D_g_1_y);
  gv_2D_g[2] = new TGraph(5000,v_2D_g_2_x,v_2D_g_2_y);
  gv_2D_g[3] = new TGraph(5000,v_2D_g_3_x,v_2D_g_3_y);

  // gv_1D_g = new TGraph(5000,v_1D_g_x,v_1D_g_y);
  

  gw_1D_g = new TGraph(5000,w_1D_g_x,w_1D_g_y);

  
}

WireCell2dToy::DataSignalWienROIFDS::~DataSignalWienROIFDS(){
  delete hu_1D_c;
  delete hu_1D_c_gaus;

  delete hu_2D_g_f;
  delete hu_2D_g;
  delete hu_2D_g_gaus;
  
  delete hv_1D_c;
  delete hv_1D_c_gaus;

  delete hv_2D_g_f;
  delete hv_2D_g;
  delete hv_2D_g_gaus;
  
  delete hw_1D_g;  
  delete hw_1D_g_gaus;

  delete gu_1D_c;
  delete gv_1D_c;
  delete gw_1D_c;
 
  //delete gv_1D_g;
  delete gw_1D_g;
  
  for (int i=0;i!=4;i++){
    delete gu_2D_g[i];
    delete gv_2D_g[i];
  }
  
  delete [] gu_2D_g;
  delete [] gv_2D_g;
}



int WireCell2dToy::DataSignalWienROIFDS::size() const{
  return max_frames;
}

void WireCell2dToy::DataSignalWienROIFDS::Save(){
}


void  WireCell2dToy::DataSignalWienROIFDS::Deconvolute_V_2D_g(){
  const Frame& frame1 = fds.get();
  size_t ntraces = frame1.traces.size();

  const int nticks = nbin;
  double value_re[9600];
  double value_im[9600];

  int scale = nbin/bins_per_frame;

  // filter 
  TF1 *filter_v = new TF1("filter_v","(x>0.0)*exp(-0.5*pow(x/[0],[1]))");
  double par[2]={1.47404e+01/200.*2.,4.97667e+00};
  filter_v->SetParameters(par);

  TF1 *filter_low = new TF1("filter_low","1-exp(-pow(x/0.0045,2))");

  TF1 *filter_wire = new TF1("filter_wire","exp(-0.5*pow(x/[0],2))");
  double par4[1] = {1.0/sqrt(3.1415926)*1.4};
  filter_wire->SetParameters(par4);
  //filter_wire->Draw();
  
  TF1 *filter_g = new TF1("filger_g","1./sqrt(2.*3.1415926)/[0]*exp(-x*x/2./[0]/[0])");
  double par3[1] = {3./2.1};
  filter_g->SetParameters(par3);
  
  TH1F *hfilter_time_gaus =new TH1F("hfilter_time_gaus","hfilter_time_gaus",nbin,0,nbin);
  for (int i=0;i!=nbin;i++){
    double xx = hfilter_time_gaus->GetBinCenter(i+1)/2.-nbin/4.;
    hfilter_time_gaus->SetBinContent(i+1,filter_g->Eval(xx));
  }
  hfilter_time_gaus->Scale(1./hfilter_time_gaus->GetSum());
  TH1 *hfilter_gaus = 0;
  hfilter_gaus = hfilter_time_gaus->FFT(hfilter_gaus,"MAG");


  // response function ... 
  TH1F *hvr = new TH1F("hvr1","hvr1",nbin,0,nbin); // half us tick
  const int nchannels = nwire_v;
  
  float scale_v = 1.251/1.074*0.91*0.85;
  double rho_res[7][nticks], phi_res[7][nticks];

  for (int j=0;j!=7;j++){
    TGraph *gtemp;
    if (j==0 || j==6){
      gtemp = gv_2D_g[3];
    }else if (j==1 || j==5){
      gtemp = gv_2D_g[2];
    }else if (j==2 || j==4){
      gtemp = gv_2D_g[1];
    }else if (j==3){
      gtemp = gv_2D_g[0];
    }
    hvr->Reset();
    for (int i=0; i!=nbin; i++){  
      double time = hvr->GetBinCenter(i+1)/2.-50 ;
      //*** scale factors for 70kV ***//
      double x = time - time_offset_uv -0.113; // 0.113 is the calibration number ... 
      if (x > -35 && x  < 15){
	if (j!=6 && j!=0 ){
	  if (gtemp->Eval(x)>0 ){
	    hvr->SetBinContent(i+1,gtemp->Eval(x)/scale_v); //70kV
	  }else{
	    hvr->SetBinContent(i+1,gtemp->Eval(x)/scale_v); //70kV
	  }
	}
      }else{
	hvr->SetBinContent(i+1,0);
      }
    }
    TH1 *hmr_v = hvr->FFT(0,"MAG");
    TH1 *hpr_v = hvr->FFT(0,"PH");
    
    for (Int_t i=0;i!=nticks;i++){
      rho_res[j][i] = hmr_v->GetBinContent(i+1);
      phi_res[j][i] = hpr_v->GetBinContent(i+1);
    }

    // std::cout << "abc " << j << std::endl;
    delete hmr_v;
    delete hpr_v;
  }


  // deconvolution ... 
  std::vector<std::vector<double>> rho_v, phi_v,result_re,result_im;
  for (int i=0;i!=nchannels;i++){
    std::vector<double> temp,temp1,temp2,temp3;
    temp.resize(nticks,0);
    rho_v.push_back(temp);
    temp1.resize(nticks,0);
    phi_v.push_back(temp1);
    temp2.resize(nticks,0);
    result_re.push_back(temp2);
    temp3.resize(nticks,0);
    result_im.push_back(temp3);
  }
  int tbin_save[nchannels];

  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    if (chid >=nwire_u+nwire_v || chid < nwire_u) continue;
    
    TH1F *htemp = new TH1F("htemp","htemp",nbin,0,nbin);
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp->SetBinContent(tt,trace.charge.at(i));
    }
    
    TH1 *hm = htemp->FFT(0,"MAG");
    TH1 *hp = htemp->FFT(0,"PH");

    for (Int_t j=0;j!=nticks;j++){
      rho_v[chid-nwire_u][j] = hm->GetBinContent(j+1);
      phi_v[chid-nwire_u][j] = hp->GetBinContent(j+1);
    }
    tbin_save[chid-nwire_u] = tbin;
    
    delete hm;
    delete hp;
    delete htemp;
  }

  double resp_re[nchannels], resp_im[nchannels];
  for (Int_t i=0;i!=nticks;i++){
    Double_t freq;
    if (i < nticks/2.){
      freq = i/(1.*nticks)*2.;
    }else{
      freq = (nticks - i)/(1.*nticks)*2.;
    }
    
    for (Int_t j=0;j!=nchannels;j++){
      value_re[j] = rho_v[j][i]*cos(phi_v[j][i]);
      value_im[j] = rho_v[j][i]*sin(phi_v[j][i]);
      if (j<7){
  	resp_re[j] = rho_res[j][i]*cos(phi_res[j][i]);
  	resp_im[j] = rho_res[j][i]*sin(phi_res[j][i]);
      }else{
  	resp_re[j] = 0.;
  	resp_im[j] = 0.;
      }
    }

    Int_t m=nchannels;
    
    TVirtualFFT *ifft = TVirtualFFT::FFT(1,&m,"C2CFORWARD M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    Double_t temp_re[nchannels],temp_im[nchannels];
    ifft->GetPointsComplex(temp_re,temp_im);
    
    ifft->SetPointsComplex(resp_re,resp_im);
    ifft->Transform();
    Double_t temp1_re[nchannels],temp1_im[nchannels];
    ifft->GetPointsComplex(temp1_re,temp1_im);
    
    Double_t temp2_re[nchannels],temp2_im[nchannels];
    for (Int_t j=0;j!=nchannels;j++){
      Double_t freq_wire;
      if (j < nchannels/2.){
	freq_wire = j/(1.*nchannels)*2.;
      }else{
	freq_wire = (nchannels -j)/(1.*nchannels)*2.;
      }
      if (temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j]>0){
  	temp2_re[j] = (temp_re[j]*temp1_re[j]+temp_im[j]*temp1_im[j])/m/
  	  (temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j]) * filter_wire->Eval(freq_wire);
  	temp2_im[j] = (temp_im[j]*temp1_re[j]-temp_re[j]*temp1_im[j])
  	  /m/(temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j]) * filter_wire->Eval(freq_wire);
      }else{
  	temp2_re[j] = 0;
  	temp2_im[j] = 0;
      }
    }
    
    TVirtualFFT *ifft3 = TVirtualFFT::FFT(1,&m,"C2CBACKWARD M K");
    ifft3->SetPointsComplex(temp2_re,temp2_im);
    ifft3->Transform();
    Double_t temp3_re[nchannels],temp3_im[nchannels];
    ifft3->GetPointsComplex(temp3_re,temp3_im);


    for (Int_t j=0;j!=nchannels;j++){
      Int_t shift = j - 3;
      if (shift <0) shift += nchannels;
      result_re[j][i] = temp3_re[shift]/nticks;//*filter_v->Eval(freq);
      result_im[j][i] = temp3_im[shift]/nticks;//*filter_v->Eval(freq);
    }

    delete ifft;
    delete ifft3;
  }
  
  for (Int_t chid=0;chid!=nchannels;chid++){
    int n = nticks;
    TVirtualFFT *ifft2 = TVirtualFFT::FFT(1,&n,"C2R M K");
    double temp_re[nticks],temp_im[nticks];
    for (int j=0;j!=nticks;j++){
      Double_t freq;
      if (j < nticks/2.){
	freq = j/(1.*nticks)*2.;
      }else{
	freq = (nticks - j)/(1.*nticks)*2.;
      }
      temp_re[j] = result_re[chid][j]*filter_v->Eval(freq);
      temp_im[j] = result_im[chid][j]*filter_v->Eval(freq);
    }
    ifft2->SetPointsComplex(temp_re,temp_im);
    ifft2->Transform();
    TH1 *fb = 0;
    fb = TH1::TransformHisto(ifft2,fb,"Re");
    
    TH1F *htemp = new TH1F("htemp","htemp",nbin,0,nbin);
    for (int i=0;i!=nbin;i++){
      htemp->SetBinContent(i+1,fb->GetBinContent(i+1)/( 14.*1.1*4096./2000.));
    }
    delete ifft2;
    delete fb;
    
    // correct baseline 
    restore_baseline(htemp);
    
    int start = -1, end = -1;
    if (vmap.find(chid)!=vmap.end()){
      start = vmap[chid].first;
      end = vmap[chid].second;
    }

    // put results back into the 2-D histogram
    for (int j=0;j!=bins_per_frame;j++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*j+k+1);
      }
      int bin = j+100;
      if (bin >=start && bin <=end){
      }else{
	hv_2D_g->SetBinContent(chid+1,j+1,int(sum));
      }
    }
    
    
    TVirtualFFT *ifft4 = TVirtualFFT::FFT(1,&n,"C2R M K");
    for (int j=0;j!=nticks;j++){
      Double_t freq;
      if (j < nticks/2.){
	freq = j/(1.*nticks)*2.;
      }else{
	freq = (nticks - j)/(1.*nticks)*2.;
      }

      temp_re[j] = result_re[chid][j]*filter_v->Eval(freq)* filter_low->Eval(freq); 
      temp_im[j] = result_im[chid][j]*filter_v->Eval(freq)* filter_low->Eval(freq); 
    }
    ifft4->SetPointsComplex(temp_re,temp_im);
    ifft4->Transform();
    TH1 *fb1 = 0;
    fb1 = TH1::TransformHisto(ifft4,fb1,"Re");
    
    for (int i=0;i!=nbin;i++){
      htemp->SetBinContent(i+1,fb1->GetBinContent(i+1)/( 14.*1.1*4096./2000.));
    }
    delete ifft4;
    delete fb1;
    
    // correct baseline 
    restore_baseline(htemp);
    
    // calculate RMS 
    double rms = cal_rms(htemp,chid+nwire_u);
    vplane_rms_g[chid] = rms*scale;

    // put results back into the 2-D histogram
    for (int j=0;j!=bins_per_frame;j++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*j+k+1);
      }
      int bin = j+100;
      if (bin >=start && bin <=end){
      }else{
	hv_2D_g_f->SetBinContent(chid+1,j+1,int(sum));
      }
    }


    ifft4 = TVirtualFFT::FFT(1,&n,"C2R M K");
    for (int j=0;j!=nticks;j++){
      temp_re[j] = result_re[chid][j]*hfilter_gaus->GetBinContent(j+1);
      temp_im[j] = result_im[chid][j]*hfilter_gaus->GetBinContent(j+1);
    }
    temp_re[0] = 0.;
    temp_im[0] = 0.;
    
    ifft4->SetPointsComplex(temp_re,temp_im);
    ifft4->Transform();
    fb1 = 0;
    fb1 = TH1::TransformHisto(ifft4,fb1,"Re");
    
    for (int i=0;i!=nbin;i++){
      htemp->SetBinContent(i+1,fb1->GetBinContent(i+1)/( 14.*1.1*4096./2000.));
    }
    delete ifft4;
    delete fb1;
    
    // correct baseline 
    restore_baseline(htemp);

    // put results back into the 2-D histogram
    for (int j=0;j!=bins_per_frame;j++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*j+k+1);
      }
      int bin = j+100;
      if (bin >=start && bin <=end){
      }else{
	hv_2D_g_gaus->SetBinContent(chid+1,j+1,int(sum));
      }
    }



    delete htemp;    
  }


  delete hfilter_time_gaus;
  delete filter_g;
  delete hfilter_gaus; 

  delete hvr;
  delete filter_v;
  delete filter_low;
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

  TF1 *filter_wire = new TF1("filter_wire","exp(-0.5*pow(x/[0],2))");
  double par4[1] = {1.0/sqrt(3.1415926)*1.4};
  filter_wire->SetParameters(par4);
  // filter_wire->Draw();
  
  

  TF1 *filter_g = new TF1("filger_g","1./sqrt(2.*3.1415926)/[0]*exp(-x*x/2./[0]/[0])");
  double par3[1] = {3./2.1};
  filter_g->SetParameters(par3);
  
  TH1F *hfilter_time_gaus =new TH1F("hfilter_time_gaus","hfilter_time_gaus",nbin,0,nbin);
  for (int i=0;i!=nbin;i++){
    double xx = hfilter_time_gaus->GetBinCenter(i+1)/2.-nbin/4.;
    hfilter_time_gaus->SetBinContent(i+1,filter_g->Eval(xx));
  }
  hfilter_time_gaus->Scale(1./hfilter_time_gaus->GetSum());
  TH1 *hfilter_gaus = 0;
  hfilter_gaus = hfilter_time_gaus->FFT(hfilter_gaus,"MAG");


  // response function ... 
  TH1F *hur = new TH1F("hur1","hur1",nbin,0,nbin); // half us tick
  const int nchannels = nwire_u;
  float scale_u = 1.51/1.16*0.91*0.85;
  double rho_res[7][nticks], phi_res[7][nticks];

  for (int j=0;j!=7;j++){
    TGraph *gtemp;
    if (j==0 || j==6){
      gtemp = gu_2D_g[3];
    }else if (j==1 || j==5){
      gtemp = gu_2D_g[2];
    }else if (j==2 || j==4){
      gtemp = gu_2D_g[1];
    }else if (j==3){
      gtemp = gu_2D_g[0];
    }
    for (int i=0; i!=nbin; i++){  
      double time = hur->GetBinCenter(i+1)/2.-50 ;
      //*** scale factors for 70kV ***//
      double x = time ;
      if (x > -35 && x  < 15){
	if (gtemp->Eval(x)>0 ){
	  hur->SetBinContent(i+1,gtemp->Eval(x)/scale_u); //70kV
	}else{
	  hur->SetBinContent(i+1,gtemp->Eval(x)/scale_u); //70kV
	}
      }else{
	hur->SetBinContent(i+1,0);
      }
    }
    TH1 *hmr_u = hur->FFT(0,"MAG");
    TH1 *hpr_u = hur->FFT(0,"PH");
    
    for (Int_t i=0;i!=nticks;i++){
      rho_res[j][i] = hmr_u->GetBinContent(i+1);
      phi_res[j][i] = hpr_u->GetBinContent(i+1);
    }

    // std::cout << "abc " << j << std::endl;
    delete hmr_u;
    delete hpr_u;
  }


  // deconvolution ... 
  std::vector<std::vector<double>> rho_u, phi_u,result_re,result_im;
  for (int i=0;i!=nchannels;i++){
    std::vector<double> temp,temp1,temp2,temp3;
    temp.resize(nticks,0);
    rho_u.push_back(temp);
    temp1.resize(nticks,0);
    phi_u.push_back(temp1);
    temp2.resize(nticks,0);
    result_re.push_back(temp2);
    temp3.resize(nticks,0);
    result_im.push_back(temp3);
  }
  int tbin_save[nchannels];

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

    for (Int_t j=0;j!=nticks;j++){
      rho_u[chid][j] = hm->GetBinContent(j+1);
      phi_u[chid][j] = hp->GetBinContent(j+1);
    }
    tbin_save[chid] = tbin;
    
    delete hm;
    delete hp;
    delete htemp;
  }

  double resp_re[nchannels], resp_im[nchannels];
  for (Int_t i=0;i!=nticks;i++){
    Double_t freq;
    if (i < nticks/2.){
      freq = i/(1.*nticks)*2.;
    }else{
      freq = (nticks - i)/(1.*nticks)*2.;
    }
    
    for (Int_t j=0;j!=nchannels;j++){
      value_re[j] = rho_u[j][i]*cos(phi_u[j][i]);
      value_im[j] = rho_u[j][i]*sin(phi_u[j][i]);
      if (j<7){
  	resp_re[j] = rho_res[j][i]*cos(phi_res[j][i]);
  	resp_im[j] = rho_res[j][i]*sin(phi_res[j][i]);
      }else{
  	resp_re[j] = 0.;
  	resp_im[j] = 0.;
      }
    }

    Int_t m=nchannels;
    
    TVirtualFFT *ifft = TVirtualFFT::FFT(1,&m,"C2CFORWARD M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    Double_t temp_re[nchannels],temp_im[nchannels];
    ifft->GetPointsComplex(temp_re,temp_im);
    
    ifft->SetPointsComplex(resp_re,resp_im);
    ifft->Transform();
    Double_t temp1_re[nchannels],temp1_im[nchannels];
    ifft->GetPointsComplex(temp1_re,temp1_im);
    
    Double_t temp2_re[nchannels],temp2_im[nchannels];
    for (Int_t j=0;j!=nchannels;j++){
      Double_t freq_wire;
      if (j < nchannels/2.){
	freq_wire = j/(1.*nchannels)*2.;
      }else{
	freq_wire = (nchannels -j)/(1.*nchannels)*2.;
      }


      if (temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j]>0){
  	temp2_re[j] = (temp_re[j]*temp1_re[j]+temp_im[j]*temp1_im[j])/m/
  	  (temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j])*filter_wire->Eval(freq_wire);
  	temp2_im[j] = (temp_im[j]*temp1_re[j]-temp_re[j]*temp1_im[j])
  	  /m/(temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j])*filter_wire->Eval(freq_wire);
      }else{
  	temp2_re[j] = 0;
  	temp2_im[j] = 0;
      }
    }
    
    TVirtualFFT *ifft3 = TVirtualFFT::FFT(1,&m,"C2CBACKWARD M K");
    ifft3->SetPointsComplex(temp2_re,temp2_im);
    ifft3->Transform();
    Double_t temp3_re[nchannels],temp3_im[nchannels];
    ifft3->GetPointsComplex(temp3_re,temp3_im);


    for (Int_t j=0;j!=nchannels;j++){
      Int_t shift = j - 3;
      if (shift <0) shift += nchannels;
      result_re[j][i] = temp3_re[shift]/nticks;//*filter_u->Eval(freq);
      result_im[j][i] = temp3_im[shift]/nticks;//*filter_u->Eval(freq);
    }

    delete ifft;
    delete ifft3;
  }
  
  for (Int_t chid=0;chid!=nchannels;chid++){
    int n = nticks;
    TVirtualFFT *ifft2 = TVirtualFFT::FFT(1,&n,"C2R M K");
    double temp_re[nticks],temp_im[nticks];
    for (int j=0;j!=nticks;j++){
      Double_t freq;
      if (j < nticks/2.){
	freq = j/(1.*nticks)*2.;
      }else{
	freq = (nticks - j)/(1.*nticks)*2.;
      }
      temp_re[j] = result_re[chid][j]*filter_u->Eval(freq);
      temp_im[j] = result_im[chid][j]*filter_u->Eval(freq);
    }
    ifft2->SetPointsComplex(temp_re,temp_im);
    ifft2->Transform();
    TH1 *fb = 0;
    fb = TH1::TransformHisto(ifft2,fb,"Re");
    
    TH1F *htemp = new TH1F("htemp","htemp",nbin,0,nbin);
    for (int i=0;i!=nbin;i++){
      htemp->SetBinContent(i+1,fb->GetBinContent(i+1)/( 14.*1.1*4096./2000.));
    }
    delete ifft2;
    delete fb;
    
    // correct baseline 
    restore_baseline(htemp);
    
    int start = -1, end = -1;
    if (umap.find(chid)!=umap.end()){
      start = umap[chid].first;
      end = umap[chid].second;
    }

    // put results back into the 2-D histogram
    for (int j=0;j!=bins_per_frame;j++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*j+k+1);
      }
      int bin = j+100;
      if (bin >=start && bin <=end){
      }else{
	hu_2D_g->SetBinContent(chid+1,j+1,int(sum));
      }
    }
    
    
    TVirtualFFT *ifft4 = TVirtualFFT::FFT(1,&n,"C2R M K");
    for (int j=0;j!=nticks;j++){
      Double_t freq;
      if (j < nticks/2.){
	freq = j/(1.*nticks)*2.;
      }else{
	freq = (nticks - j)/(1.*nticks)*2.;
      }

      temp_re[j] = result_re[chid][j]*filter_u->Eval(freq)* filter_low->Eval(freq); 
      temp_im[j] = result_im[chid][j]*filter_u->Eval(freq)* filter_low->Eval(freq); 
    }
    ifft4->SetPointsComplex(temp_re,temp_im);
    ifft4->Transform();
    TH1 *fb1 = 0;
    fb1 = TH1::TransformHisto(ifft4,fb1,"Re");
    
    for (int i=0;i!=nbin;i++){
      htemp->SetBinContent(i+1,fb1->GetBinContent(i+1)/( 14.*1.1*4096./2000.));
    }
    delete ifft4;
    delete fb1;
    
    // correct baseline 
    restore_baseline(htemp);
    
    // calculate RMS 
    double rms = cal_rms(htemp,chid);
    uplane_rms_g[chid] = rms*scale;

    // put results back into the 2-D histogram
    for (int j=0;j!=bins_per_frame;j++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*j+k+1);
      }
      int bin = j+100;
      if (bin >=start && bin <=end){
      }else{
	hu_2D_g_f->SetBinContent(chid+1,j+1,int(sum));
      }
    }


    ifft4 = TVirtualFFT::FFT(1,&n,"C2R M K");
    for (int j=0;j!=nticks;j++){
      temp_re[j] = result_re[chid][j]*hfilter_gaus->GetBinContent(j+1);
      temp_im[j] = result_im[chid][j]*hfilter_gaus->GetBinContent(j+1);
    }
    temp_re[0] = 0.;
    temp_im[0] = 0.;
    
    ifft4->SetPointsComplex(temp_re,temp_im);
    ifft4->Transform();
    fb1 = 0;
    fb1 = TH1::TransformHisto(ifft4,fb1,"Re");
    
    for (int i=0;i!=nbin;i++){
      htemp->SetBinContent(i+1,fb1->GetBinContent(i+1)/( 14.*1.1*4096./2000.));
    }
    delete ifft4;
    delete fb1;
    
    // correct baseline 
    restore_baseline(htemp);

    // put results back into the 2-D histogram
    for (int j=0;j!=bins_per_frame;j++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*j+k+1);
      }
      int bin = j+100;
      if (bin >=start && bin <=end){
      }else{
	hu_2D_g_gaus->SetBinContent(chid+1,j+1,int(sum));
      }
    }



    delete htemp;    
  }


  delete hfilter_time_gaus;
  delete filter_g;
  delete hfilter_gaus; 

  delete hur;
  delete filter_u;
  delete filter_low;
}


void WireCell2dToy::DataSignalWienROIFDS::ROI_cal(TH1F *h1_1, TH1F *h2_1, TH1F *h3_1, TH1F *h4_1, TH1F *h5_1, Double_t threshold0,Double_t threshold2, TH1F *hresult, TH1F *hresult1, int flag_u = 0){
  
    Double_t th = threshold2*2.0; // do 2 sigma
    Double_t th1 = threshold0*3.6; // do 3.6 sigma
    Double_t th2 = threshold2*4.0; // do three sigma ... 

    //std::cout << th << " " << th1 << " " << th2 << std::endl;
    //h1_1 1D_c
    //h2_1 2D_g
    //h3_1 2D_g_f
    //h4_1 1D gaussian filter
    //h5_1 2D gaussian filter 

    //calibrate the 1D histogram (temp fix)
    h2_1->Reset();
    h3_1->Reset();
    h5_1->Reset();
    
    //calibrate the 2D histgraom (temp fix)
    // h1_1->Reset();
    // h4_1->Reset();

    std::vector<std::pair <Int_t,Int_t> > ROIs;

    for (Int_t i=0;i<h3_1->GetNbinsX();i++){
      Double_t content = h3_1->GetBinContent(i+1);
      Double_t next_content;
      
      if (i!=h3_1->GetNbinsX()-1){
  	next_content = h3_1->GetBinContent(i+2);
      }else{
  	next_content = h3_1->GetBinContent(1);
      }
      
      if ((content > th )){
  	// go to find the beginning of ROI
  	Int_t begin;
  	begin = find_ROI_begin(h3_1,i, th*0.8) ;
	
  	Int_t end;
  	// go to find the end of ROI 
  	end = find_ROI_end(h3_1,i, th*0.8) ;
	
  	// save ROI 
  	ROIs.push_back(std::make_pair(begin,end));
  	//    std::cout << begin << " " << end << endl;
	
  	// reset the beginning of the search for next ROI
  	if (end <h3_1->GetNbinsX()){
  	  i = end;
  	}else{
  	  i=h3_1->GetNbinsX();
  	}
      }
    }
    
    // std::cout << ROIs.size() << std::endl;

    std::vector<std::pair <Int_t,Int_t> > ROIs_1;
    for (Int_t i=0;i<h1_1->GetNbinsX();i++){
      Double_t content = h1_1->GetBinContent(i+1);
      
      if ((content > th1 )){
  	// go to find the beginning of ROI
  	Int_t begin;
  	begin = find_ROI_begin(h1_1,i, th*0.25) ;
	
  	Int_t end;
  	// go to find the end of ROI 
  	end = find_ROI_end(h1_1,i, th*0.25) ;
	
  	// save ROI 
  	ROIs_1.push_back(std::make_pair(begin,end));
  	//    std::cout << begin << " " << end << endl;
	
  	// reset the beginning of the search for next ROI
  	if (end <h1_1->GetNbinsX()){
  	  i = end;
  	}else{
  	  i=h1_1->GetNbinsX();
  	}
      }
    }
    
    //std::cout << ROIs_1.size() << std::endl;

    TH1F *h2_3 = (TH1F*)h2_1->Clone("h2_3");
    h2_3->Reset();
        
    for (Int_t i=0;i!=ROIs_1.size();i++){
      Int_t begin = ROIs_1.at(i).first;
      Int_t end = ROIs_1.at(i).second;
      
      Double_t content_begin;
      Double_t content_end;
      
      if (begin <0){
  	content_begin = h2_1->GetBinContent(begin+1+h2_1->GetNbinsX());
      }else{
  	content_begin = h2_1->GetBinContent(begin+1);
      }
      
      if (end >= h2_1->GetNbinsX()){
  	content_end = h2_1->GetBinContent(end+1-h2_1->GetNbinsX());
      }else{
  	content_end = h2_1->GetBinContent(end+1);
      }
      
      TH1F *htemp = new TH1F("htemp","htemp",end-begin+1,begin,end+1);

      for (Int_t j=begin;j<=end;j++){
  	Double_t content_current;
  	if (j <0){
  	  content_current = h2_1->GetBinContent(j+1+h2_1->GetNbinsX());
  	}else if (j >= h2_1->GetNbinsX()){
  	  content_current = h2_1->GetBinContent(j+1 - h2_1->GetNbinsX());
  	}else{
  	  content_current = h2_1->GetBinContent(j+1);
  	}
	
  	content_current -=  content_begin + (content_end-content_begin) * (j-begin)/(end-begin);
	
	// if (content_current <= th2) content_current = 0;

	// if (content_current < h1_1->GetBinContent(j+1))
	//   content_current = h1_1->GetBinContent(j+1);

  	// h2_3->SetBinContent(j+1,content_current);

	htemp->SetBinContent(j-begin+1,content_current);
      }

       TSpectrum *s = new TSpectrum(100);
       Int_t nfound = s->Search(htemp,2,"nobackground new",0.1);
       
       if (nfound >1){
	 //cout << htemp->GetNbinsX() << " " << nfound << " " << begin << " " << end << endl;
	 Int_t npeaks = s->GetNPeaks();
	 Double_t *peak_pos = s->GetPositionX();
	 
	 // const int size = npeaks;
	 // const int size_1 = npeaks+1;
	 int order_peak_pos[105];
	 
	 for (Int_t j=0;j!=npeaks;j++){
	   order_peak_pos[j] = *(peak_pos+j);
	 }
	 std::sort(order_peak_pos,order_peak_pos + npeaks);
	 
	 // for (Int_t j=0;j!=npeaks;j++){
	 // 	std::cout << order_peak_pos[j] << std::endl;
	 // }
	 
	 Float_t valley_pos[25];
	 valley_pos[0] = begin;
	 
	 for (Int_t j=0;j!=npeaks-1;j++){
	   Float_t min = 1e9;
	   for (Int_t k = order_peak_pos[j]-begin; k< order_peak_pos[j+1]-begin;k++){
	     
	     if (htemp->GetBinContent(k+1) < min){
	       min = htemp->GetBinContent(k+1);
	       valley_pos[j+1] = k+begin;
	     }
	   }
	   //	std::cout << valley_pos[j+1] << std::endl;
	   //std::cout << *(peak_pos+j) << std::endl;
	 }
	 valley_pos[npeaks] = end;
	 
	 std::set<int> saved_boundaries;
	 
	 for (Int_t j=0;j!=npeaks;j++){
	   int flag = 0;
	   Int_t start_pos = valley_pos[j];
	   Double_t start_content = htemp->GetBinContent(valley_pos[j]-begin+1);
	   Int_t end_pos = valley_pos[j+1];
	   Double_t end_content = htemp->GetBinContent(valley_pos[j+1]-begin+1);
	   Int_t Peak_pos = order_peak_pos[j];
	   Double_t peak_content = htemp->GetBinContent(order_peak_pos[j]-begin+1);
	   
	   
	   // std::cout << start_pos << " " << start_content << " " << Peak_pos << " " << peak_content
	   // 	  << " " << end_pos << " " << peak_content << std::endl;
	   
	   if ((start_content > th || end_content > th) &&
	       start_content < peak_content /2. && end_content < peak_content/2.) flag = 1;
	   // deal with the small peaks ... 
	   if (peak_content > th2 && peak_content < th2 * 1.5) flag =1;
	   
	   if (flag==1){
	     saved_boundaries.insert(start_pos);
	     saved_boundaries.insert(end_pos);
	     //	  std::cout << start_pos << " " << start_content << " " << end_pos << 
	     //  " " << end_content << std::endl;
	   }
	   //	std::cout << valley_pos[j] << " " << *(peak_pos+j) << " " << valley_pos[j+1] <<std::endl;
	 }
	 
	 //   std::cout << saved_boundaries.size() << std::endl;
	 for (Int_t j=0;j!=npeaks;j++){
	   int flag = 0;
	   Int_t start_pos = valley_pos[j];
	   Double_t start_content = htemp->GetBinContent(valley_pos[j]-begin+1);
	   Int_t end_pos = valley_pos[j+1];
	   Double_t end_content = htemp->GetBinContent(valley_pos[j+1]-begin+1);
	   Int_t Peak_pos = order_peak_pos[j];
	   Double_t peak_content = htemp->GetBinContent(order_peak_pos[j]-begin+1);
	   
	   if (saved_boundaries.find(start_pos) != saved_boundaries.end() ||
	       saved_boundaries.find(end_pos) != saved_boundaries.end()){
	     
	     for (Int_t k = start_pos; k!=end_pos+1;k++){
	       Double_t temp_content = htemp->GetBinContent(k-begin+1) - (start_content + (end_content-start_content) * (k-start_pos) / (end_pos - start_pos));
	       htemp->SetBinContent(k-begin+1,temp_content);
	     }
	     //std::cout << "Adaptive Baseline " << start_pos << " " << end_pos << std::endl;
	   }
	 }
	 //htemp->Draw();
	 //TH1 *hb = s->Background(htemp,20,"same");
	 //for (Int_t j=0;j!=hb->GetNbinsX();j++){
	 //h2_5->SetBinContent(j+begin+1,hb->GetBinContent(j+1));
	 //}
	 //l3->Draw("same");
	 //hb->Draw("same");
      //hb->SetLineColor(4);
      // c1->Update();
      // int abc;
      // cin >> abc;
      //delete hb;
       }
       
       //int flag = 0;
       for (Int_t j=0;j!=htemp->GetNbinsX();j++){
	 if (htemp->GetBinContent(j+1) >= th2){
	   //flag  = 1;
	   h2_3->SetBinContent(j+begin+1,htemp->GetBinContent(j+1));
	 }else{
	   if (h1_1->GetBinContent(j+begin+1) > th1)
	     h2_3->SetBinContent(j+begin+1,h1_1->GetBinContent(j+begin+1));
	 }
       }
       // if (flag == 0){
       // 	 for (Int_t j=0;j!=htemp->GetNbinsX();j++){
       // 	    h2_3->SetBinContent(j+begin+1,h1_1->GetBinContent(j+begin+1));
       // 	 }
       // }
       
       
       delete htemp;
       delete s;

    }
    
    
    TH1F *h2_2 = (TH1F*)h2_1->Clone("h2_2");
    h2_2->Reset();
    
    
    for (Int_t i=0;i!=ROIs.size();i++){
      Int_t begin = ROIs.at(i).first;
      Int_t end = ROIs.at(i).second;
      
      Double_t content_begin;
      Double_t content_end;
      
      if (begin <0){
  	content_begin = h2_1->GetBinContent(begin+1+h2_1->GetNbinsX());
      }else{
  	content_begin = h2_1->GetBinContent(begin+1);
      }
      
      if (end >= h2_1->GetNbinsX()){
  	content_end = h2_1->GetBinContent(end+1-h2_1->GetNbinsX());
      }else{
  	content_end = h2_1->GetBinContent(end+1);
      }
      
      TH1F *htemp = new TH1F("htemp","htemp",end-begin+1,begin,end+1);

      for (Int_t j=begin;j<=end;j++){
  	Double_t content_current;
  	if (j <0){
  	  content_current = h2_1->GetBinContent(j+1+h2_1->GetNbinsX());
  	}else if (j >= h2_1->GetNbinsX()){
  	  content_current = h2_1->GetBinContent(j+1 - h2_1->GetNbinsX());
  	}else{
  	  content_current = h2_1->GetBinContent(j+1);
  	}
	
  	content_current -=  content_begin + (content_end-content_begin) * (j-begin)/(end-begin);
  	// if(content_current > th2)
  	//   h2_2->SetBinContent(j+1,content_current);

	htemp->SetBinContent(j-begin+1,content_current);
      }

      TSpectrum *s = new TSpectrum(100);
      Int_t nfound = s->Search(htemp,2,"nobackground new",0.1);
      
      if (nfound >1 ){
	//cout << htemp->GetNbinsX() << " " << nfound << " " << begin << " " << end << endl;
	Int_t npeaks = s->GetNPeaks();
	Double_t *peak_pos = s->GetPositionX();
	
	// const int size = npeaks;
	// const int size_1 = npeaks+1;
	int order_peak_pos[105];
	
	for (Int_t j=0;j!=npeaks;j++){
	  order_peak_pos[j] = *(peak_pos+j);
	}
	std::sort(order_peak_pos,order_peak_pos + npeaks);
	
	// for (Int_t j=0;j!=npeaks;j++){
	// 	std::cout << order_peak_pos[j] << std::endl;
	// }
	
	Float_t valley_pos[25];
	valley_pos[0] = begin;
	
	for (Int_t j=0;j!=npeaks-1;j++){
	  Float_t min = 1e9;
	  for (Int_t k = order_peak_pos[j]-begin; k< order_peak_pos[j+1]-begin;k++){
	    
	    if (htemp->GetBinContent(k+1) < min){
	      min = htemp->GetBinContent(k+1);
	      valley_pos[j+1] = k+begin;
	    }
	  }
	  //	std::cout << valley_pos[j+1] << std::endl;
	  //std::cout << *(peak_pos+j) << std::endl;
	}
	valley_pos[npeaks] = end;
	
	std::set<int> saved_boundaries;
	
	for (Int_t j=0;j!=npeaks;j++){
	  int flag = 0;
	  Int_t start_pos = valley_pos[j];
	  Double_t start_content = htemp->GetBinContent(valley_pos[j]-begin+1);
	  Int_t end_pos = valley_pos[j+1];
	  Double_t end_content = htemp->GetBinContent(valley_pos[j+1]-begin+1);
	  Int_t Peak_pos = order_peak_pos[j];
	  Double_t peak_content = htemp->GetBinContent(order_peak_pos[j]-begin+1);
	  
	  if ((start_content >= th2 || end_content >= th2) &&
	      start_content < peak_content /2. && end_content < peak_content/2.) flag = 1;
	  
	  // deal with the small peaks ... 
	  if (peak_content > th2 && peak_content < th2 * 1.5) flag =1;
	  
	  if (flag==1){
	    saved_boundaries.insert(start_pos);
	    saved_boundaries.insert(end_pos);
	    //	  std::cout << start_pos << " " << start_content << " " << end_pos << 
	    //  " " << end_content << std::endl;
	  }
	  //	std::cout << valley_pos[j] << " " << *(peak_pos+j) << " " << valley_pos[j+1] <<std::endl;
	}
	
	//  std::cout << saved_boundaries.size() << std::endl;
	for (Int_t j=0;j!=npeaks;j++){
	  int flag = 0;
	  Int_t start_pos = valley_pos[j];
	  Double_t start_content = htemp->GetBinContent(valley_pos[j]-begin+1);
	  Int_t end_pos = valley_pos[j+1];
	  Double_t end_content = htemp->GetBinContent(valley_pos[j+1]-begin+1);
	  Int_t Peak_pos = order_peak_pos[j];
	  Double_t peak_content = htemp->GetBinContent(order_peak_pos[j]-begin+1);
	  
	  if (saved_boundaries.find(start_pos) != saved_boundaries.end() ||
	      saved_boundaries.find(end_pos) != saved_boundaries.end()){
	    
	    for (Int_t k = start_pos; k!=end_pos+1;k++){
	      Double_t temp_content = htemp->GetBinContent(k-begin+1) - (start_content + (end_content-start_content) * (k-start_pos) / (end_pos - start_pos));
	      htemp->SetBinContent(k-begin+1,temp_content);
	    }
	    //std::cout << "Adaptive Baseline " << start_pos << " " << end_pos << std::endl;
	  }
	}
	//htemp->Draw();
	//TH1 *hb = s->Background(htemp,20,"same");
	//for (Int_t j=0;j!=hb->GetNbinsX();j++){
	//h2_5->SetBinContent(j+begin+1,hb->GetBinContent(j+1));
	//}
	//l3->Draw("same");
	//hb->Draw("same");
	//hb->SetLineColor(4);
	// c1->Update();
	// int abc;
	// cin >> abc;
	//delete hb;
      }
      
      
      for (Int_t j=0;j!=htemp->GetNbinsX();j++){
	if (htemp->GetBinContent(j+1) > th2)
	  h2_2->SetBinContent(j+begin+1,htemp->GetBinContent(j+1));
      }
      
      
      delete htemp;
      delete s;
      
    }
    
    for (Int_t i=0;i!=h2_2->GetNbinsX();i++){
      Double_t content1 = h2_2->GetBinContent(i+1);
      Double_t content2 = h2_3->GetBinContent(i+1);
    
      if (content1 > content2){
   	hresult->SetBinContent(i+1,content1);
      }else{
   	hresult->SetBinContent(i+1,content2);
      }
    }

    delete h2_3;
    delete h2_2;


    // save the Gaussian results ... 

    TH1F *h2_4 = (TH1F*)h2_1->Clone("h2_4");
    h2_4->Reset();
        
    for (Int_t i=0;i!=ROIs_1.size();i++){
      Int_t begin = ROIs_1.at(i).first;
      Int_t end = ROIs_1.at(i).second;
      
      Double_t content_begin;
      Double_t content_end;
      
      if (begin <0){
  	content_begin = h5_1->GetBinContent(begin+1+h5_1->GetNbinsX());
      }else{
  	content_begin = h5_1->GetBinContent(begin+1);
      }
      
      if (end >= h5_1->GetNbinsX()){
  	content_end = h5_1->GetBinContent(end+1-h5_1->GetNbinsX());
      }else{
  	content_end = h5_1->GetBinContent(end+1);
      }
      
      for (Int_t j=begin;j<=end;j++){
  	Double_t content_current;
  	if (j <0){
  	  content_current = h5_1->GetBinContent(j+1+h5_1->GetNbinsX());
  	}else if (j >= h5_1->GetNbinsX()){
  	  content_current = h5_1->GetBinContent(j+1 - h5_1->GetNbinsX());
  	}else{
  	  content_current = h5_1->GetBinContent(j+1);
  	}
	
  	content_current -=  content_begin + (content_end-content_begin) * (j-begin)/(end-begin);
	
	if (content_current < h4_1->GetBinContent(j+1) && flag_u ==1)
	  content_current = h4_1->GetBinContent(j+1);

  	h2_4->SetBinContent(j+1,content_current);
      }
    }
    
    
    TH1F *h2_5 = (TH1F*)h2_1->Clone("h2_5");
    h2_5->Reset();
    
    
    for (Int_t i=0;i!=ROIs.size();i++){
      Int_t begin = ROIs.at(i).first;
      Int_t end = ROIs.at(i).second;
      
      Double_t content_begin;
      Double_t content_end;
      
      if (begin <0){
  	content_begin = h5_1->GetBinContent(begin+1+h5_1->GetNbinsX());
      }else{
  	content_begin = h5_1->GetBinContent(begin+1);
      }
      
      if (end >= h5_1->GetNbinsX()){
  	content_end = h5_1->GetBinContent(end+1-h5_1->GetNbinsX());
      }else{
  	content_end = h5_1->GetBinContent(end+1);
      }
      
      for (Int_t j=begin;j<=end;j++){
  	Double_t content_current;
  	if (j <0){
  	  content_current = h5_1->GetBinContent(j+1+h5_1->GetNbinsX());
  	}else if (j >= h5_1->GetNbinsX()){
  	  content_current = h5_1->GetBinContent(j+1 - h5_1->GetNbinsX());
  	}else{
  	  content_current = h5_1->GetBinContent(j+1);
  	}
	
  	content_current -=  content_begin + (content_end-content_begin) * (j-begin)/(end-begin);
	h2_5->SetBinContent(j+1,content_current);
      }
    }
    
    for (Int_t i=0;i!=h2_4->GetNbinsX();i++){
      Double_t content1 = h2_4->GetBinContent(i+1);
      Double_t content2 = h2_5->GetBinContent(i+1);
      if (content1 > content2){
   	hresult1->SetBinContent(i+1,content1);
      }else{
   	hresult1->SetBinContent(i+1,content2);
      }
    }

    delete h2_4;
    delete h2_5;


}




Int_t WireCell2dToy::DataSignalWienROIFDS::find_ROI_end(TH1F *h1, Int_t bin, Double_t th=0){
  Int_t end = bin;
  Double_t content = h1->GetBinContent(end+1);
  while(content>th){
    end ++;
    
    if (end >=h1->GetNbinsX()){
      content = h1->GetBinContent(end - h1->GetNbinsX()+1);
    }else{
      content = h1->GetBinContent(end+1);
    }
  }

  while(local_ave(h1,end+1,1) < local_ave(h1,end,1)){
    end++;
  } 
  return end;
}

Int_t WireCell2dToy::DataSignalWienROIFDS::find_ROI_begin(TH1F *h1, Int_t bin, Double_t th=0){
  // find the first one before bin and is below threshold ... 
  Int_t begin = bin;
  Double_t content = h1->GetBinContent(begin+1);
  while(content > th){
    begin --;

    if (begin <0){
      content = h1->GetBinContent(begin+h1->GetNbinsX()+1);
    }else{
      content = h1->GetBinContent(begin+1);
    }

  }
  
  // calculate the local average
  // keep going and find the minimum
  while( local_ave(h1,begin-1,1) < local_ave(h1,begin,1)){
    //cout << "X " << begin << " " <<local_ave(h1,begin,1)  << " " << local_ave(h1,begin-1,1)  << endl;
    begin --;
  }
  //  Double_t current_ave = local_ave(h1, begin, 1);
  //  cout << bin << " " << begin << " " << content << " " << current_ave << endl;
  
  return begin;
}

Double_t WireCell2dToy::DataSignalWienROIFDS::local_ave(TH1F *h1, Int_t bin, Int_t width){
  Double_t sum1 = 0;
  Double_t sum2 = 0;
  
  for (Int_t i=-width;i<width+1;i++){
    Int_t current_bin = bin + i;

    while (current_bin <0)
      current_bin += h1->GetNbinsX();
    while (current_bin >= h1->GetNbinsX())
      current_bin -= h1->GetNbinsX();
    
    sum1 += h1->GetBinContent(current_bin+1);
    sum2 ++;
  }

  if (sum2>0){
    return sum1/sum2;
  }else{
    return 0;
  }
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

  TF1 *filter_g = new TF1("filger_g","1./sqrt(2.*3.1415926)/[0]*exp(-x*x/2./[0]/[0])");
  double par3[1] = {3./2.1};
  filter_g->SetParameters(par3);
  
  TH1F *hfilter_time_gaus =new TH1F("hfilter_time_gaus","hfilter_time_gaus",nbin,0,nbin);
  for (int i=0;i!=nbin;i++){
    double xx = hfilter_time_gaus->GetBinCenter(i+1)/2.-nbin/4.;
    hfilter_time_gaus->SetBinContent(i+1,filter_g->Eval(xx));
  }
  hfilter_time_gaus->Scale(1./hfilter_time_gaus->GetSum());
  TH1 *hfilter_gaus = 0;
  hfilter_gaus = hfilter_time_gaus->FFT(hfilter_gaus,"MAG");


  //response function 
  TH1F *hvr = new TH1F("hvr1","hvr1",nbin,0,nbin); // half us tick
  float scale_v = 1.251/1.074*0.91*0.85;
  for (int i=0;i!=nbin;i++){
    double time  = hvr->GetBinCenter(i+1)/2.-50;
    double x = time-time_offset_uv-0.113;
    if (x > -35 && x  < 15){
      if (gv_2D_g[0]->Eval(x)>0){
	hvr->SetBinContent(i+1,gv_2D_g[0]->Eval(x)/scale_v);
      }else{
	hvr->SetBinContent(i+1,gv_2D_g[0]->Eval(x)/scale_v);
      }
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

    int start=-1, end=-1;
    if (vmap.find(chid-nwire_u)!=vmap.end()){
      start = vmap[chid-nwire_u].first;
      end = vmap[chid-nwire_u].second;
    }
    // put results back into the 2-D histogram
    for (int i=0;i!=bins_per_frame;i++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*i+k+1);
      }
      int bin = i+100;
      if (bin >=start && bin <=end){
      }else{
	hv_2D_g->SetBinContent(chid+1-nwire_u,i+1,int(sum));
      }
    }

    
    for (int i=0;i!=nbin;i++){
      
      double rho = hm->GetBinContent(i+1)/hmr_v->GetBinContent(i+1)*hfilter_gaus->GetBinContent(i+1);
      double phi = hp->GetBinContent(i+1) - hpr_v->GetBinContent(i+1);
      if (i==0) rho = 0;
      value_re[i] = rho*cos(phi)/nbin;
      value_im[i] = rho*sin(phi)/nbin;

    }

    ifft = TVirtualFFT::FFT(1,&nbin,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
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

    // put results back into the 2-D histogram
    for (int i=0;i!=bins_per_frame;i++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*i+k+1);
      }
      int bin = i+100;
      if (bin >=start && bin <=end){
      }else{
	hv_2D_g_gaus->SetBinContent(chid+1-nwire_u,i+1,int(sum));
      }
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
      int bin = i+100;
      if (bin >=start && bin <=end){
      }else{
	hv_2D_g_f->SetBinContent(chid+1-nwire_u,i+1,int(sum));
      }
    }

    

    delete htemp;
    delete hm;
    delete hp;
    
  }
  
  delete hfilter_time_gaus;
  delete filter_g;
  delete hfilter_gaus; 

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
  
  TF1 *filter_g = new TF1("filger_g","1./sqrt(2.*3.1415926)/[0]*exp(-x*x/2./[0]/[0])");
  double par3[1] = {3./2.1};
  filter_g->SetParameters(par3);
  
  TH1F *hfilter_time_gaus =new TH1F("hfilter_time_gaus","hfilter_time_gaus",nbin,0,nbin);
  for (int i=0;i!=nbin;i++){
    double xx = hfilter_time_gaus->GetBinCenter(i+1)/2.-nbin/4.;
    hfilter_time_gaus->SetBinContent(i+1,filter_g->Eval(xx));
  }
  hfilter_time_gaus->Scale(1./hfilter_time_gaus->GetSum());
  TH1 *hfilter_gaus = 0;
  hfilter_gaus = hfilter_time_gaus->FFT(hfilter_gaus,"MAG");
  


  //response function
  TH1F *hwr = new TH1F("hwr1","hwr1",nbin,0,nbin); // half us tick
  for (int i=0;i!=nbin;i++){
    double time  = hwr->GetBinCenter(i+1)/2.-50;
    double x = time-time_offset_uw + 0.803;
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
    
    int start=-1, end=-1;
    if (wmap.find(chid-nwire_u-nwire_v)!=wmap.end()){
      start = wmap[chid-nwire_u-nwire_v].first;
      end = wmap[chid-nwire_u-nwire_v].second;
    }
    
    // put results back into the 2-D histogram
    for (int i=0;i!=bins_per_frame;i++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*i+k+1);
      }
      int bin = i+100;
      if (bin >=start && bin <=end){
      }else{
	hw_1D_g->SetBinContent(chid+1-nwire_u-nwire_v,i+1,int(sum));
      }
    }

    for (int i=0;i!=nbin;i++){
      
      double rho = hm->GetBinContent(i+1)/hmr_w->GetBinContent(i+1)*hfilter_gaus->GetBinContent(i+1);
      double phi = hp->GetBinContent(i+1) - hpr_w->GetBinContent(i+1);
      if (i==0) rho = 0;
      
      value_re[i] = rho*cos(phi)/nbin;
      value_im[i] = rho*sin(phi)/nbin;
    }

    ifft = TVirtualFFT::FFT(1,&nbin,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
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
    
    // put results back into the 2-D histogram
    for (int i=0;i!=bins_per_frame;i++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*i+k+1);
      }
      int bin = i+100;
      if (bin >=start && bin <=end){
      }else{
	hw_1D_g_gaus->SetBinContent(chid+1-nwire_u-nwire_v,i+1,int(sum));
      }
    }

    
    delete htemp;
    delete hm;
    delete hp;
    
  }

   delete hfilter_time_gaus;
  delete filter_g;
  delete hfilter_gaus; 
  


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

  TF1 *filter_g = new TF1("filger_g","1./sqrt(2.*3.1415926)/[0]*exp(-x*x/2./[0]/[0])");
  double par3[1] = {3.5/2.1};
  filter_g->SetParameters(par3);
  
  TH1F *hfilter_time_gaus =new TH1F("hfilter_time_gaus","hfilter_time_gaus",nbin,0,nbin);
  for (int i=0;i!=nbin;i++){
    double xx = hfilter_time_gaus->GetBinCenter(i+1)/2.-nbin/4.;
    hfilter_time_gaus->SetBinContent(i+1,filter_g->Eval(xx));
  }
  hfilter_time_gaus->Scale(1./hfilter_time_gaus->GetSum());
  TH1 *hfilter_gaus = 0;
  hfilter_gaus = hfilter_time_gaus->FFT(hfilter_gaus,"MAG");


  //response function
  TH1F *hvr = new TH1F("hvr1","hvr1",nbin,0,nbin); // half us tick
  float scale_v = 1.251/1.074*0.91*0.85;
  for (int i=0;i!=nbin;i++){
    double time  = hvr->GetBinCenter(i+1)/2.-50;
    double x = time-time_offset_uv + 0.887;
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
    
    int start=-1, end=-1;
    if (vmap.find(chid-nwire_u)!=vmap.end()){
      start = vmap[chid-nwire_u].first;
      end = vmap[chid-nwire_u].second;
    }
    
    // put results back into the 2-D histogram
    for (int i=0;i!=bins_per_frame;i++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*i+k+1);
      }
      int bin = i+100;
      if (bin >=start && bin <=end){
      }else{
	hv_1D_c->SetBinContent(chid+1-nwire_u,i+1,int(sum));
      }
    }

    // put in Gaussian results
     for (int i=0;i!=nbin;i++){
      
      double rho = hm->GetBinContent(i+1)/hmr_v->GetBinContent(i+1)*hfilter_gaus->GetBinContent(i+1);
      double phi = hp->GetBinContent(i+1) - hpr_v->GetBinContent(i+1);
      if (i==0) rho = 0;
      
      value_re[i] = rho*cos(phi)/nbin;
      value_im[i] = rho*sin(phi)/nbin;
    }

    ifft = TVirtualFFT::FFT(1,&nbin,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
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

   
    
    // put results back into the 2-D histogram
    for (int i=0;i!=bins_per_frame;i++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*i+k+1);
      }
      int bin = i+100;
      if (bin >=start && bin <=end){
      }else{
	hv_1D_c_gaus->SetBinContent(chid+1-nwire_u,i+1,int(sum));
      }
    }



    delete htemp;
    delete hm;
    delete hp;
    
  }

  delete hfilter_time_gaus;
  delete filter_g;
  delete hfilter_gaus; 
  
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
  
  TF1 *filter_g = new TF1("filger_g","1./sqrt(2.*3.1415926)/[0]*exp(-x*x/2./[0]/[0])");
  double par3[1] = {3.5/2.1};
  filter_g->SetParameters(par3);
  
  TH1F *hfilter_time_gaus =new TH1F("hfilter_time_gaus","hfilter_time_gaus",nbin,0,nbin);
  for (int i=0;i!=nbin;i++){
    double xx = hfilter_time_gaus->GetBinCenter(i+1)/2.-nbin/4.;
    hfilter_time_gaus->SetBinContent(i+1,filter_g->Eval(xx));
  }
  hfilter_time_gaus->Scale(1./hfilter_time_gaus->GetSum());
  TH1 *hfilter_gaus = 0;
  hfilter_gaus = hfilter_time_gaus->FFT(hfilter_gaus,"MAG");


  //response function ...
  TH1F *hur = new TH1F("hur1","hur1",nbin,0,nbin); // half us tick
  float scale_u = 1.51/1.16*0.91*0.85;
  for (int i=0;i!=nbin;i++){
    double time  = hur->GetBinCenter(i+1)/2.-50;
    double x = time + 0.1;
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
    
    int start=-1, end=-1;
    if (umap.find(chid)!=umap.end()){
      start = umap[chid].first;
      end = umap[chid].second;
    }

    // put results back into the 2-D histogram
    for (int i=0;i!=bins_per_frame;i++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*i+k+1);
      }
      int bin = i+100;
      if (bin >=start && bin <=end){
      }else{
	hu_1D_c->SetBinContent(chid+1,i+1,int(sum));
      }
    }

    // put in the Gaussian results ... 
    for (int i=0;i!=nbin;i++){
      double rho = hm->GetBinContent(i+1)/hmr_u->GetBinContent(i+1)*hfilter_gaus->GetBinContent(i+1);
      double phi = hp->GetBinContent(i+1) - hpr_u->GetBinContent(i+1);
      if (i==0) rho = 0;
      
      value_re[i] = rho*cos(phi)/nbin;
      value_im[i] = rho*sin(phi)/nbin;
    }

    ifft = TVirtualFFT::FFT(1,&nbin,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
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
    
    // put results back into the 2-D histogram
    for (int i=0;i!=bins_per_frame;i++){
      double sum = 0;
      for (int k=0;k!=scale;k++){
	sum += htemp->GetBinContent(scale*i+k+1);
      }
      int bin = i+100;
      if (bin >=start && bin <=end){
      }else{
	hu_1D_c_gaus->SetBinContent(chid+1,i+1,int(sum));
      }
    }
    

    delete htemp;
    delete hm;
    delete hp;
    
  }
  

  delete hfilter_time_gaus;
  delete filter_g;
  delete hfilter_gaus; 

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
 
  std::cout << "Deconvolution with garfield field response for 2-D V Plane" << std::endl;
  //Deconvolute_V_1D_g();
  Deconvolute_V_2D_g();

  std::cout << "Deconvolution with garfield field response for 2-D U Plane" << std::endl;
  Deconvolute_U_2D_g();
 
  // 1: 1D_c, 2: 2D_g, 3: 2D_g_filter, 0: ROI
  int flag_save = 0;
  
  
  std::cout << "Load results back into frame" << std::endl;
  // load the results back into the frame ... 
  for (int i=0;i!=nwire_u;i++){
    Trace t;
    t.chid = i;
    t.tbin = 0;
    t.charge.resize(bins_per_frame, 0.0);
    int nticks = bins_per_frame;
    
    TH1F *h1_1 = new TH1F("h1_1","h1_1",nticks,0,nticks);
    TH1F *h2_1 = new TH1F("h2_1","h2_1",nticks,0,nticks);
    TH1F *h3_1 = new TH1F("h3_1","h3_1",nticks,0,nticks);
    TH1F *h4_1 = new TH1F("h4_1","h4_1",nticks,0,nticks);
    TH1F *h5_1 = new TH1F("h5_1","h5_1",nticks,0,nticks);
    TH1F *hresult = new TH1F("hresult","hresult",nticks,0,nticks);
    TH1F *hresult1 = new TH1F("hresult1","hresult1",nticks,0,nticks);

    Double_t th1 = uplane_rms.at(i);
    Double_t th2 = uplane_rms_g.at(i);

    

    for (Int_t j=0;j!=nticks;j++){
      Double_t content = hu_1D_c->GetBinContent(i+1,j+1);
      h1_1->SetBinContent(j+1,content);
    
      content = hu_2D_g->GetBinContent(i+1,j+1);
      h2_1->SetBinContent(j+1,content);
      
      content = hu_2D_g_f->GetBinContent(i+1,j+1);
      h3_1->SetBinContent(j+1,content);

      content = hu_1D_c_gaus->GetBinContent(i+1,j+1);
      h4_1->SetBinContent(j+1,content);
      
      content = hu_2D_g_gaus->GetBinContent(i+1,j+1);
      h5_1->SetBinContent(j+1,content);
    }
    
    ROI_cal(h1_1,h2_1,h3_1,h4_1,h5_1,th1,th2,hresult,hresult1,1);

    for (Int_t j=0;j!=nticks;j++){
      hu_2D_g_gaus->SetBinContent(i+1,j+1,hresult1->GetBinContent(j+1));
    }


    for (int j=0;j!= bins_per_frame; j++){
      if (flag_save == 1){
	// save 1D_c result
	t.charge.at(j) = h1_1->GetBinContent(j+1);
      }else if (flag_save == 2){
	// save 2D_g result
	t.charge.at(j) = h2_1->GetBinContent(j+1);
      }else if (flag_save == 3){
	// save 2D_g_filter result
	t.charge.at(j) = h3_1->GetBinContent(j+1);
      }else if (flag_save == 0){
	//normal result
	t.charge.at(j) = hresult->GetBinContent(j+1);//hu_2D_g_f->GetBinContent(i+1,j+1);
      }else if (flag_save ==4){
	t.charge.at(j) = h4_1->GetBinContent(j+1);
      }else if (flag_save ==5){
	t.charge.at(j) = h5_1->GetBinContent(j+1);
      }else if (flag_save ==6){
	t.charge.at(j) = hresult1->GetBinContent(j+1);
      }
    }
    frame.traces.push_back(t);

    delete h1_1;
    delete h2_1;
    delete h3_1;
    delete h4_1;
    delete h5_1;
    delete hresult;
    delete hresult1;
  }

  for (int i=0;i!=nwire_v;i++){
    Trace t;
    t.chid = i+nwire_u;
    t.tbin = 0;
    t.charge.resize(bins_per_frame, 0.0);

    int nticks = bins_per_frame;
    
    TH1F *h1_1 = new TH1F("h1_1","h1_1",nticks,0,nticks);
    TH1F *h2_1 = new TH1F("h2_1","h2_1",nticks,0,nticks);
    TH1F *h3_1 = new TH1F("h3_1","h3_1",nticks,0,nticks);
    TH1F *h4_1 = new TH1F("h4_1","h4_1",nticks,0,nticks);
    TH1F *h5_1 = new TH1F("h5_1","h5_1",nticks,0,nticks);
    TH1F *hresult = new TH1F("hresult","hresult",nticks,0,nticks);
    TH1F *hresult1 = new TH1F("hresult1","hresult1",nticks,0,nticks);
    
    Double_t th1 = vplane_rms.at(i);
    Double_t th2 = vplane_rms_g.at(i);

    for (Int_t j=0;j!=nticks;j++){
      Double_t content = hv_1D_c->GetBinContent(i+1,j+1);
      h1_1->SetBinContent(j+1,content);
    
      content = hv_2D_g->GetBinContent(i+1,j+1);
      h2_1->SetBinContent(j+1,content);
      
      content = hv_2D_g_f->GetBinContent(i+1,j+1);
      h3_1->SetBinContent(j+1,content);

      content = hv_1D_c_gaus->GetBinContent(i+1,j+1);
      h4_1->SetBinContent(j+1,content);
    
      content = hv_2D_g_gaus->GetBinContent(i+1,j+1);
      h5_1->SetBinContent(j+1,content);
    }
    
    ROI_cal(h1_1,h2_1,h3_1,h4_1,h5_1,th1,th2,hresult,hresult1,1);

    for (Int_t j=0;j!=nticks;j++){
      hv_2D_g_gaus->SetBinContent(i+1,j+1,hresult1->GetBinContent(j+1));
    }

    for (int j=0;j!= bins_per_frame; j++){
      if (flag_save == 1){
	t.charge.at(j) = h1_1->GetBinContent(j+1);
      }else if (flag_save == 2){
	t.charge.at(j) = h2_1->GetBinContent(j+1);
      }else if (flag_save == 3){
	t.charge.at(j) = h3_1->GetBinContent(j+1);
      }else if (flag_save == 0){
	t.charge.at(j) = hresult->GetBinContent(j+1);//hu_2D_g->GetBinContent(i+1,j+1);
      }else if (flag_save ==4){
	t.charge.at(j) = h4_1->GetBinContent(j+1);
      }else if (flag_save ==5){
	t.charge.at(j) = h5_1->GetBinContent(j+1);
      }else if (flag_save ==6){
	t.charge.at(j) = hresult1->GetBinContent(j+1);
      }
    }
    frame.traces.push_back(t);

    delete h1_1;
    delete h2_1;
    delete h3_1;
    delete h4_1;
    delete h5_1;
    delete hresult;
    delete hresult1;
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


