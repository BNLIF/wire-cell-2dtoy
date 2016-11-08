#include "WireCell2dToy/ToySignalGaus.h"
#include "WireCell2dToy/ToySignalSimu.h"
#include "WireCellData/GeomWire.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TMatrixDSparse.h"
#include "TMatrixD.h"
#include "TVectorD.h"

using namespace WireCell;

WireCell2dToy::ToySignalGausFDS::ToySignalGausFDS(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, int bins_per_frame1, int nframes_total, float time_offset_uv, float time_offset_uw, float overall_time_offset)
  : fds(fds)
  , gds_flag(0)
  , gds(&gds)
  , dgds(0)
  , max_frames(nframes_total)
  , time_offset_uv(time_offset_uv)
  , time_offset_uw(time_offset_uw)
  , overall_time_offset(overall_time_offset)
{  
  bins_per_frame = bins_per_frame1;

  // GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  // GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  // GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  // nwire_u = wires_u.size();
  // nwire_v = wires_v.size();
  // nwire_w = wires_w.size();

  nbin = fds.Get_Bins_Per_Frame();

  // hu = new TH1F*[nwire_u];
  // hv = new TH1F*[nwire_v];
  // hw = new TH1F*[nwire_w];

  // for (int i=0;i!=nwire_u;i++){
  //   hu[i] = new TH1F(Form("U3_%d",i),Form("U3_%d",i),nbin,0,nbin);
  // }
  // for (int i=0;i!=nwire_v;i++){
  //   hv[i] = new TH1F(Form("V3_%d",i),Form("V3_%d",i),nbin,0,nbin);
  // }
  // for (int i=0;i!=nwire_w;i++){
  //   hw[i] = new TH1F(Form("W3_%d",i),Form("W3_%d",i),nbin,0,nbin);
  // }
  
  hu = new TH1F("U3","U3",nbin,0,nbin);
  hv = new TH1F("V3","V3",nbin,0,nbin);
  hw = new TH1F("W3","W3",nbin,0,nbin);

  
  //define filters
  filter_g = new TF1("filger_g","1./sqrt(2.*3.1415926)/[0]*exp(-x*x/2./[0]/[0])");
  double par3[1] = {2./2.2};
  filter_g->SetParameters(par3);
  
  hfilter_time_gaus =new TH1F("hfilter_time_gaus","hfilter_time_gaus",nbin,0,nbin);
  for (int i=0;i!=nbin;i++){
    double xx = hfilter_time_gaus->GetBinCenter(i+1)/2.-nbin/4.;
    hfilter_time_gaus->SetBinContent(i+1,filter_g->Eval(xx));
  }
  hfilter_time_gaus->Scale(1./hfilter_time_gaus->GetSum());
  hfilter_gaus = 0;
  hfilter_gaus = hfilter_time_gaus->FFT(hfilter_gaus,"MAG");
  

  //get in the response function ... 
  #include "data.txt"

  gu = new TGraph(5000,xu,yu);
  gv = new TGraph(5000,xv,yv);
  gw = new TGraph(5000,xw,yw);

  hur = new TH1F("hur2","hur2",nbin,0,nbin); // half us tick
  hvr = new TH1F("hvr2","hvr2",nbin,0,nbin); // half us tick
  hwr = new TH1F("hwr2","hwr2",nbin,0,nbin); // half us tick
  
  for (int i=0; i!=nbin; i++){  
    double time = hur->GetBinCenter(i+1)/2.-50 ;
    hur->SetBinContent(i+1,gu->Eval(time-overall_time_offset));
    hvr->SetBinContent(i+1,gv->Eval(time-time_offset_uv-overall_time_offset));
    hwr->SetBinContent(i+1,gw->Eval(time-time_offset_uw-overall_time_offset));
  } 
  
  hmr_u = 0;
  hmr_v = 0;
  hmr_w = 0;
  hpr_u = 0;
  hpr_v = 0;
  hpr_w = 0;
}

WireCell2dToy::ToySignalGausFDS::ToySignalGausFDS(WireCell::FrameDataSource& fds, const WireCell::DetectorGDS& gds, int bins_per_frame1, int nframes_total, float time_offset_uv, float time_offset_uw, float overall_time_offset)
  : fds(fds)
  , gds_flag(1)
  , gds(0)
  , dgds(&gds)
  , max_frames(nframes_total)
  , time_offset_uv(time_offset_uv)
  , time_offset_uw(time_offset_uw)
  , overall_time_offset(overall_time_offset)
{  
  bins_per_frame = bins_per_frame1;

  // GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  // GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  // GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  // nwire_u = wires_u.size();
  // nwire_v = wires_v.size();
  // nwire_w = wires_w.size();

  nbin = fds.Get_Bins_Per_Frame();

  // hu = new TH1F*[nwire_u];
  // hv = new TH1F*[nwire_v];
  // hw = new TH1F*[nwire_w];

  // for (int i=0;i!=nwire_u;i++){
  //   hu[i] = new TH1F(Form("U3_%d",i),Form("U3_%d",i),nbin,0,nbin);
  // }
  // for (int i=0;i!=nwire_v;i++){
  //   hv[i] = new TH1F(Form("V3_%d",i),Form("V3_%d",i),nbin,0,nbin);
  // }
  // for (int i=0;i!=nwire_w;i++){
  //   hw[i] = new TH1F(Form("W3_%d",i),Form("W3_%d",i),nbin,0,nbin);
  // }
  
  hu = new TH1F("U3","U3",nbin,0,nbin);
  hv = new TH1F("V3","V3",nbin,0,nbin);
  hw = new TH1F("W3","W3",nbin,0,nbin);

  
  //define filters
  filter_g = new TF1("filger_g","1./sqrt(2.*3.1415926)/[0]*exp(-x*x/2./[0]/[0])");
  double par3[1] = {2./2.2};
  filter_g->SetParameters(par3);
  
  hfilter_time_gaus =new TH1F("hfilter_time_gaus","hfilter_time_gaus",nbin,0,nbin);
  for (int i=0;i!=nbin;i++){
    double xx = hfilter_time_gaus->GetBinCenter(i+1)/2.-nbin/4.;
    hfilter_time_gaus->SetBinContent(i+1,filter_g->Eval(xx));
  }
  hfilter_time_gaus->Scale(1./hfilter_time_gaus->GetSum());
  hfilter_gaus = 0;
  hfilter_gaus = hfilter_time_gaus->FFT(hfilter_gaus,"MAG");
  

  //get in the response function ... 
  #include "data.txt"

  gu = new TGraph(5000,xu,yu);
  gv = new TGraph(5000,xv,yv);
  gw = new TGraph(5000,xw,yw);

  hur = new TH1F("hur2","hur2",nbin,0,nbin); // half us tick
  hvr = new TH1F("hvr2","hvr2",nbin,0,nbin); // half us tick
  hwr = new TH1F("hwr2","hwr2",nbin,0,nbin); // half us tick
  
  for (int i=0; i!=nbin; i++){  
    double time = hur->GetBinCenter(i+1)/2.-50 ;
    hur->SetBinContent(i+1,gu->Eval(time-overall_time_offset));
    hvr->SetBinContent(i+1,gv->Eval(time-time_offset_uv-overall_time_offset));
    hwr->SetBinContent(i+1,gw->Eval(time-time_offset_uw-overall_time_offset));
  } 
  
  hmr_u = 0;
  hmr_v = 0;
  hmr_w = 0;
  hpr_u = 0;
  hpr_v = 0;
  hpr_w = 0;
}




int WireCell2dToy::ToySignalGausFDS::size() const{
  return max_frames;
}

void WireCell2dToy::ToySignalGausFDS::Save(){
  TFile *file = new TFile("temp_gaus.root","RECREATE");
  // for (int i=0;i!=nwire_u;i++){
  //   TH1F *huu = (TH1F*)hu[i]->Clone(Form("U1_%d",i));
  // }
  // for (int i=0;i!=nwire_v;i++){
  //   TH1F *hvv = (TH1F*)hv[i]->Clone(Form("V1_%d",i));
  // }
  // for (int i=0;i!=nwire_w;i++){
  //   TH1F *hww = (TH1F*)hw[i]->Clone(Form("W1_%d",i));
  // }

  // TH1F *hftg = (TH1F*)hfilter_time_gaus->Clone("hftg");
  // TH1 *hfg = (TH1*)hfilter_gaus->Clone("hfg");
  
 
  
  file->Write();
  file->Close();
}


int WireCell2dToy::ToySignalGausFDS::jump(int frame_number){
  // fill the frame data ... 
  if (frame.index == frame_number) {
    return frame_number;
  }
  frame.clear();
  int scale = nbin/bins_per_frame;

  // std::cout << "test1 " << std::endl;
  fds.jump(frame_number);
  //std::cout << "test2 " << std::endl;
  const Frame& frame1 = fds.get();
  //std::cout << "test3 " << std::endl;
  size_t ntraces = frame1.traces.size();
  //std::cout << "test4 " << std::endl;

  TVirtualFFT::SetTransform(0);
  
  
  hmr_u = hur->FFT(0,"MAG");
  hmr_v = hvr->FFT(0,"MAG");
  hmr_w = hwr->FFT(0,"MAG");

  hpr_u = hur->FFT(0,"PH");
  hpr_v = hvr->FFT(0,"PH");
  hpr_w = hwr->FFT(0,"PH");
  
  TH1 *hm = 0;
  TH1 *hp = 0;
  double value_re[9600];
  double value_im[9600];
  TVirtualFFT *ifft;
  TH1 *fb = 0;
  int n = nbin;
  TH1 *hmr, *hpr;

  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    
    //std::cout << tbin << " " << chid << " " << nbins << std::endl;

    TH1F *htemp;
    TF1 *filter;
    
    if (gds_flag == 0){
      WirePlaneType_t plane = gds->by_channel(chid).at(0)->plane();
      if (plane == WirePlaneType_t(0)){
	htemp = hu;
	hmr = hmr_u;
	hpr = hpr_u;
      }else if (plane == WirePlaneType_t(1)){
	htemp = hv;
	hmr = hmr_v;
	hpr = hpr_v;
      }else if (plane == WirePlaneType_t(2)){
	htemp= hw;
	hmr = hmr_w;
	hpr = hpr_w;
      }
    }else{
      WirePlaneType_t plane = dgds->by_channel(chid).at(0)->plane();
      if (plane == WirePlaneType_t(0)){
	htemp = hu;
	hmr = hmr_u;
	hpr = hpr_u;
      }else if (plane == WirePlaneType_t(1)){
	htemp = hv;
	hmr = hmr_v;
	hpr = hpr_v;
      }else if (plane == WirePlaneType_t(2)){
	htemp= hw;
	hmr = hmr_w;
	hpr = hpr_w;
      }
    }
    // if (chid < nwire_u){
    //   htemp = hu;
    //   hmr = hmr_u;
    //   hpr = hpr_u;
    // }else if (chid < nwire_u + nwire_v){
    //   htemp = hv;
    //   hmr = hmr_v;
    //   hpr = hpr_v;
    // }else{
    //   htemp = hw;
    //   hmr = hmr_w;
    //   hpr = hpr_w;
    // }

    

    htemp->Reset();

    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;//  + overall_time_shift * 2;
      // if (tt > nbin) tt -= nbin;
      // if (tt < 1 ) tt += nbin;
      htemp->SetBinContent(tt,trace.charge.at(i));
    }
        
    hm = htemp->FFT(0,"MAG");
    hp = htemp->FFT(0,"PH");
    
    for (int i=0;i!=nbin;i++){
      double rho = hm->GetBinContent(i+1)/hmr->GetBinContent(i+1)*hfilter_gaus->GetBinContent(i+1);
      double phi = hp->GetBinContent(i+1) - hpr->GetBinContent(i+1);

      if (i==0) rho = 0;
      value_re[i] = rho*cos(phi)/nbin;
      value_im[i] = rho*sin(phi)/nbin;
    }
    
    ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    fb = TH1::TransformHisto(ifft,0,"Re");

    for (int i=0;i!=nbin;i++){
      htemp->SetBinContent(i+1,fb->GetBinContent(i+1)/( 7.8*4096./2000.));
    }
    delete hm;
    delete hp;
    delete ifft;
    delete fb;
    

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
    
    //save into frame

    Trace t;
    t.chid = chid;
    t.tbin = tbin;
    t.charge.resize(bins_per_frame, 0.0);
    for (int j=0;j!=bins_per_frame;j++){
      t.charge.at(j) = 0;
      for (int k=0;k!=scale;k++){
	t.charge.at(j) += htemp->GetBinContent(scale*j+k+1);
      }
    }
    frame.traces.push_back(t);
   
  }

  delete hmr_u;
  delete hmr_v;
  delete hmr_w;
  delete hpr_u;
  delete hpr_v;
  delete hpr_w;
  
  

  frame.index = frame_number;
  return frame.index;
}


WireCell2dToy::ToySignalGausFDS::~ToySignalGausFDS(){
  // for (int i=0;i!=nwire_u;i++){
  //   delete hu[i] ;
  // }
  delete hu;
  // for (int i=0;i!=nwire_v;i++){
  //   delete hv[i] ;
  // }
  delete hv;
  // for (int i=0;i!=nwire_w;i++){
  //   delete hw[i] ;
  // }
  delete hw;
  
  delete filter_g;
  delete hfilter_time_gaus;
  delete hfilter_gaus;
  
  delete hur;
  delete hvr;
  delete hwr;
 
  delete gu;
  delete gv;
  delete gw;
 
}
