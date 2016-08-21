#include "WireCell2dToy/uBooNE_Data_2D_Deconvolution.h"

#include "WireCellData/GeomWire.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "TRandom.h"
#include "TF1.h"
#include "TSpectrum.h"

#include <set>
#include <algorithm>

using namespace WireCell;


WireCell2dToy::uBooNEData2DDeconvolutionFDS::uBooNEData2DDeconvolutionFDS(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds,WireCell::ChirpMap& umap, WireCell::ChirpMap& vmap, WireCell::ChirpMap& wmap, int nframes_total, float time_offset_uv, float time_offset_uw, float overall_time_offset)
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
  bins_per_frame = fds.Get_Bins_Per_Frame();
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  hu_2D_g = new TH2I("hu_2D_g","hu_2D_g",nwire_u,0,nwire_u,bins_per_frame,0,bins_per_frame);
  hv_2D_g = new TH2I("hv_2D_g","hv_2D_g",nwire_v,0,nwire_v,bins_per_frame,0,bins_per_frame);
  hw_2D_g = new TH2I("hw_2D_g","hw_2D_g",nwire_w,0,nwire_w,bins_per_frame,0,bins_per_frame);

  #include "data_70_2D_11.txt"
  
  gu_2D_g = new TGraph*[11];
  gu_2D_g[0] = new TGraph(5000,u_2D_g_0_x,u_2D_g_0_y);
  gu_2D_g[1] = new TGraph(5000,u_2D_g_1_x,u_2D_g_1_y);
  gu_2D_g[2] = new TGraph(5000,u_2D_g_2_x,u_2D_g_2_y);
  gu_2D_g[3] = new TGraph(5000,u_2D_g_3_x,u_2D_g_3_y);
  gu_2D_g[4] = new TGraph(5000,u_2D_g_4_x,u_2D_g_4_y);
  gu_2D_g[5] = new TGraph(5000,u_2D_g_5_x,u_2D_g_5_y);
  gu_2D_g[6] = new TGraph(5000,u_2D_g_6_x,u_2D_g_6_y);
  gu_2D_g[7] = new TGraph(5000,u_2D_g_7_x,u_2D_g_7_y);
  gu_2D_g[8] = new TGraph(5000,u_2D_g_8_x,u_2D_g_8_y);
  gu_2D_g[9] = new TGraph(5000,u_2D_g_9_x,u_2D_g_9_y);
  gu_2D_g[10] = new TGraph(5000,u_2D_g_10_x,u_2D_g_10_y);

  gv_2D_g = new TGraph*[11];
  gv_2D_g[0] = new TGraph(5000,v_2D_g_0_x,v_2D_g_0_y);
  gv_2D_g[1] = new TGraph(5000,v_2D_g_1_x,v_2D_g_1_y);
  gv_2D_g[2] = new TGraph(5000,v_2D_g_2_x,v_2D_g_2_y);
  gv_2D_g[3] = new TGraph(5000,v_2D_g_3_x,v_2D_g_3_y);
  gv_2D_g[4] = new TGraph(5000,v_2D_g_4_x,v_2D_g_4_y);
  gv_2D_g[5] = new TGraph(5000,v_2D_g_5_x,v_2D_g_5_y);
  gv_2D_g[6] = new TGraph(5000,v_2D_g_6_x,v_2D_g_6_y);
  gv_2D_g[7] = new TGraph(5000,v_2D_g_7_x,v_2D_g_7_y);
  gv_2D_g[8] = new TGraph(5000,v_2D_g_8_x,v_2D_g_8_y);
  gv_2D_g[9] = new TGraph(5000,v_2D_g_9_x,v_2D_g_9_y);
  gv_2D_g[10] = new TGraph(5000,v_2D_g_10_x,v_2D_g_10_y);

  
  gw_2D_g = new TGraph*[11];
  gw_2D_g[0] = new TGraph(5000,w_2D_g_0_x,w_2D_g_0_y);
  gw_2D_g[1] = new TGraph(5000,w_2D_g_1_x,w_2D_g_1_y);
  gw_2D_g[2] = new TGraph(5000,w_2D_g_2_x,w_2D_g_2_y);
  gw_2D_g[3] = new TGraph(5000,w_2D_g_3_x,w_2D_g_3_y);
  gw_2D_g[4] = new TGraph(5000,w_2D_g_4_x,w_2D_g_4_y);
  gw_2D_g[5] = new TGraph(5000,w_2D_g_5_x,w_2D_g_5_y);
  gw_2D_g[6] = new TGraph(5000,w_2D_g_6_x,w_2D_g_6_y);
  gw_2D_g[7] = new TGraph(5000,w_2D_g_7_x,w_2D_g_7_y);
  gw_2D_g[8] = new TGraph(5000,w_2D_g_8_x,w_2D_g_8_y);
  gw_2D_g[9] = new TGraph(5000,w_2D_g_9_x,w_2D_g_9_y);
  gw_2D_g[10] = new TGraph(5000,w_2D_g_10_x,w_2D_g_10_y);

  scale_u_2d = 1.0;
  scale_v_2d = 1.0;
}

WireCell2dToy::uBooNEData2DDeconvolutionFDS::~uBooNEData2DDeconvolutionFDS()
{
  delete hu_2D_g;
  delete hv_2D_g;
  delete hw_2D_g;

  for (int i=0;i!=11;i++){
    delete gu_2D_g[i];
    delete gv_2D_g[i];
    delete gw_2D_g[i];
  }
  
  delete [] gu_2D_g;
  delete [] gv_2D_g;
  delete [] gw_2D_g;
}

int WireCell2dToy::uBooNEData2DDeconvolutionFDS::jump(int frame_number){
  if (frame.index == frame_number) {
    return frame_number;
  }
  frame.clear();
  fds.jump(frame_number);

  TVirtualFFT::SetTransform(0);


  std::cout << "Deconvolution with garfield field response for 2-D W Plane" << std::endl;
  Deconvolute_W();
 
  std::cout << "Deconvolution with garfield field response for 2-D V Plane" << std::endl;
  Deconvolute_V();

  std::cout << "Deconvolution with garfield field response for 2-D U Plane" << std::endl;
  Deconvolute_U();


  // load data to be added ... 
  for (int i=0;i!=nwire_u;i++){
    Trace t;
    t.chid = i;
    t.tbin = 0;
    t.charge.resize(bins_per_frame, 0.0);
    for (int j=0;j!= bins_per_frame; j++){
      t.charge.at(j) = hu_2D_g->GetBinContent(i+1,j+1);
    }
    frame.traces.push_back(t);
  }

  for (int i=0;i!=nwire_v;i++){
    Trace t;
    t.chid = i+nwire_u;
    t.tbin = 0;
    t.charge.resize(bins_per_frame, 0.0);
    for (int j=0;j!= bins_per_frame; j++){
      t.charge.at(j) = hv_2D_g->GetBinContent(i+1,j+1);
    }
    frame.traces.push_back(t);
  }

  for (int i=0;i!=nwire_w;i++){
    Trace t;
    t.chid = i + nwire_u + nwire_v;
    t.tbin = 0;
    t.charge.resize(bins_per_frame, 0.0);
    for (int j=0;j!= bins_per_frame; j++){
      t.charge.at(j) = hw_2D_g->GetBinContent(i+1,j+1);
    }
    frame.traces.push_back(t);
  }
  

  frame.index = frame_number;
  return frame.index;
}

int WireCell2dToy::uBooNEData2DDeconvolutionFDS::size() const{
  return max_frames;
}

void  WireCell2dToy::uBooNEData2DDeconvolutionFDS::Deconvolute_U(){
  const Frame& frame1 = fds.get();
  size_t ntraces = frame1.traces.size();

  const int nticks = bins_per_frame;
  int nbin = nticks;
  double value_re[nticks];
  double value_im[nticks];

  TF1 *filter_u = new TF1("filter_u","(x>0.0)*exp(-0.5*pow(x/[0],[1]))");
  double par[2]={1.43555e+01/200.*2.,4.95096e+00};
  filter_u->SetParameters(par);

  TF1 *filter_wire = new TF1("filter_wire","exp(-0.5*pow(x/[0],2))");
  double par4[1] = {1.0/sqrt(3.1415926)*1.4};
  filter_wire->SetParameters(par4);

  TF1 *filter_low = new TF1("filter_low","1-exp(-pow(x/0.02,2))");

  TH1F *hur = new TH1F("hur1","hur1",nbin,0,nbin); // half us tick
  const int nchannels = nwire_u;
  float scale_u = scale_u_2d;
  double rho_res[21][nticks], phi_res[21][nticks];

   for (int j=0;j!=21;j++){ 
    TGraph *gtemp;
    if (j==0 || j==20){
      gtemp = gu_2D_g[10];
    }else if (j==1 || j==19){
      gtemp = gu_2D_g[9];
    }else if (j==2 || j==18){
      gtemp = gu_2D_g[8];
    }else if (j==3 || j==17){
      gtemp = gu_2D_g[7];
    }else if (j==4 || j==16){
      gtemp = gu_2D_g[6];
    }else if (j==5 || j==15){
      gtemp = gu_2D_g[5];
    }else if (j==6 || j==14){
      gtemp = gu_2D_g[4];
    }else if (j==7 || j==13){
      gtemp = gu_2D_g[3];
    }else if (j==8 || j==12){
      gtemp = gu_2D_g[2];
    }else if (j==9 || j==11){
      gtemp = gu_2D_g[1];
    }else if (j==10){
      gtemp = gu_2D_g[0];
    }

    for (int i=0; i!=nbin; i++){  
      double time = hur->GetBinCenter(i+1)/2.-90 ;
      //*** scale factors for 70kV ***//
      double x = time ;
      if (x > -84 && x  < 15.8){
	hur->SetBinContent(i+1,gtemp->Eval(x)/scale_u); //70kV
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
       if (j<21){
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
       Int_t shift = j - 10;
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
       temp_re[j] = result_re[chid][j]*filter_u->Eval(freq) * filter_low->Eval(freq);
       temp_im[j] = result_im[chid][j]*filter_u->Eval(freq) * filter_low->Eval(freq);
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
       int bin = j+180;
       if (bin >= bins_per_frame) bin -= bins_per_frame;
       if (bin >=start && bin <=end){
       }else{
	 hu_2D_g->SetBinContent(chid+1,bin+1,int(htemp->GetBinContent(j+1)));
       }
     }
     
     delete htemp;
   }
     
   
   delete hur;
   delete filter_u;
   delete filter_low;
   delete filter_wire;
   
}

void  WireCell2dToy::uBooNEData2DDeconvolutionFDS::Deconvolute_V(){
  
}

void  WireCell2dToy::uBooNEData2DDeconvolutionFDS::Deconvolute_W(){
  
}



void WireCell2dToy::uBooNEData2DDeconvolutionFDS::restore_baseline(TH1F *htemp){
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
