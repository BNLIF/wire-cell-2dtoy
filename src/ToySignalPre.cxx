#include "WireCell2dToy/ToySignalPre.h"
#include "WireCell2dToy/ToySignalSimu.h"
#include "WireCellData/GeomWire.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TMatrixDSparse.h"
#include "TMatrixD.h"
#include "TVectorD.h"

using namespace WireCell;

WireCell2dToy::ToySignalPreFDS::ToySignalPreFDS(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, int bins_per_frame1, int nframes_total)
  : fds(fds)
  , gds(gds)
  , max_frames(nframes_total)
{  
  bins_per_frame = bins_per_frame1;

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
    hu[i] = new TH1F(Form("U3_%d",i),Form("U3_%d",i),bins_per_frame,0,bins_per_frame);
  }
  for (int i=0;i!=nwire_v;i++){
    hv[i] = new TH1F(Form("V3_%d",i),Form("V3_%d",i),bins_per_frame,0,bins_per_frame);
  }
  for (int i=0;i!=nwire_w;i++){
    hw[i] = new TH1F(Form("W3_%d",i),Form("W3_%d",i),bins_per_frame,0,bins_per_frame);
  }
  

  //define filter
  TF1 *filter_u = new TF1("filter_u","(x>0.0)*gaus*exp(-0.5*pow(x/[3],[4]))");
  double par[5]={1.73/0.959301, 1.69, 1.55, 0.19, 3.75};
  filter_u->SetParameters(par);

  TF1 *filter_v = new TF1("filter_v","(x>0.0)*gaus*exp(-0.5*pow(x/[3],[4]))");
  double par1[5]={1.74/0.941034, 1.46, 1.33, 0.23, 4.89};
  filter_v->SetParameters(par1);

  TF1 *filter_w = new TF1("filter_y","(x>0.0)*[0]*exp(-0.5*(((x-[1])/[2])^2)^[3])");
  double par2[4]={1.03/0.995635, 0.08, 0.15, 2.17};
  filter_w->SetParameters(par2);


  //move filter to time domain
  int nbin = fds.Get_Bins_Per_Frame();
  hfilter_u = new TH1F("hfilter_u","hfilter_u",nbin,0,nbin);
  hfilter_v = new TH1F("hfilter_v","hfilter_v",nbin,0,nbin);
  hfilter_w = new TH1F("hfilter_w","hfilter_w",nbin,0,nbin);
  
  for (Int_t i=0;i!=nbin;i++){
    Double_t frequency = (nbin/2.-fabs(i-nbin/2.))*2./nbin;
    hfilter_u->SetBinContent(i+1,filter_u->Eval(frequency));
    hfilter_v->SetBinContent(i+1,filter_v->Eval(frequency));
    hfilter_w->SetBinContent(i+1,filter_w->Eval(frequency));
  }

  hfilter_time_u = new TH1F("hfilter_time_u","hfilter_time_u",nbin*2,-nbin/4.,nbin/4.);
  hfilter_time_v = new TH1F("hfilter_time_v","hfilter_time_v",nbin*2,-nbin/4.,nbin/4.);
  hfilter_time_w = new TH1F("hfilter_time_w","hfilter_time_w",nbin*2,-nbin/4.,nbin/4.);
  
  double value_re[9600],value_im[9600];
  Int_t n = nbin;
  TVirtualFFT *ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
  TH1 *fb;
  double baseline ;

  //U-plane
  for (Int_t i=0;i!=nbin;i++){
    Double_t filter = hfilter_u->GetBinContent(i+1);
    value_re[i] = filter/nbin;
    value_im[i] = 0;
  }
  ifft->SetPointsComplex(value_re,value_im);
  ifft->Transform();  
  fb = 0;
  fb = TH1::TransformHisto(ifft,fb,"Re");
  
  baseline = -1./nbin;
  for (int i=0;i!=nbin;i++){
    Double_t content = fb->GetBinContent(i+1);
    content -=baseline;
    Int_t newbin = i + nbin/2;
    if (newbin > nbin) newbin -=nbin;
    hfilter_time_u->SetBinContent(2*newbin,content/2.);
    hfilter_time_u->SetBinContent(2*newbin+1,content/2.);
  }
  // V-plane
  for (Int_t i=0;i!=nbin;i++){
    Double_t filter = hfilter_v->GetBinContent(i+1);
    value_re[i] = filter/nbin;
    value_im[i] = 0;
  }
  ifft->SetPointsComplex(value_re,value_im);
  ifft->Transform();  
  fb = 0;
  fb = TH1::TransformHisto(ifft,fb,"Re");
  
  baseline = -1./nbin;
  for (int i=0;i!=nbin;i++){
    Double_t content = fb->GetBinContent(i+1);
    content -=baseline;
    Int_t newbin = i + nbin/2;
    if (newbin > nbin) newbin -=nbin;
    hfilter_time_v->SetBinContent(2*newbin,content/2.);
    hfilter_time_v->SetBinContent(2*newbin+1,content/2.);
  }
  // W-plane
  for (Int_t i=0;i!=nbin;i++){
    Double_t filter = hfilter_w->GetBinContent(i+1);
    value_re[i] = filter/nbin;
    value_im[i] = 0;
  }
  ifft->SetPointsComplex(value_re,value_im);
  ifft->Transform();  
  fb = 0;
  fb = TH1::TransformHisto(ifft,fb,"Re");
  
  baseline = -1./nbin;
  for (int i=0;i!=nbin;i++){
    Double_t content = fb->GetBinContent(i+1);
    content -=baseline;
    Int_t newbin = i + nbin/2;
    if (newbin > nbin) newbin -=nbin;
    hfilter_time_w->SetBinContent(2*newbin,content/2.);
    hfilter_time_w->SetBinContent(2*newbin+1,content/2.);
  }
  
  double xx[1000],yy[1000];
  for (int i=0;i!=1000;i++){
    xx[i] = hfilter_time_u->GetBinCenter(nbin+i-500);
    yy[i] = hfilter_time_u->GetBinContent(nbin+i-500);
  }
  gu = new TGraph(1000,xx,yy);
  
  for (int i=0;i!=1000;i++){
    xx[i] = hfilter_time_v->GetBinCenter(nbin+i-500);
    yy[i] = hfilter_time_v->GetBinContent(nbin+i-500);
  }
  gv = new TGraph(1000,xx,yy);
  
  for (int i=0;i!=1000;i++){
    xx[i] = hfilter_time_w->GetBinCenter(nbin+i-500);
    yy[i] = hfilter_time_w->GetBinContent(nbin+i-500);
  }
  gw = new TGraph(1000,xx,yy);
  
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

  TH1F *hfilter_uu = (TH1F*)hfilter_time_u->Clone("UU");
  TH1F *hfilter_vv = (TH1F*)hfilter_time_v->Clone("VV");
  TH1F *hfilter_ww = (TH1F*)hfilter_time_w->Clone("WW");

  gu->Write("gu");
  gv->Write("gv");
  gw->Write("gw");
  
  file->Write();
  file->Close();
}


int WireCell2dToy::ToySignalPreFDS::jump(int frame_number){
  // fill the frame data ... 
  frame.clear();

  int n = bins_per_frame;
  int m = fds.Get_Bins_Per_Frame();
  int scale = m/n;

  TMatrixDSparse A_u(m,n);
  TMatrixDSparse AT_u(n,m);
  TMatrixD MB_u(n,n);
  TMatrixD FM_u(n,m);

  double x[40],y[40];
  //treat filter ... 
  for (int i=0;i!=40;i++){
    x[i] = (i-20+0.5)*0.5;
    y[i] =  gu->Eval(x[i]+0.5/4.);
    y[i] += gu->Eval(x[i]-0.5/4.);
  }
  // fill matrix
  for (int i=0;i!=n;i++){
    for (int j=0;j!=40;j++){
      Int_t row = i;
      Int_t column = scale*i-20+j+1;
      if (column >=0 && column < m){
	A_u(column,row) = y[j];
      }
    }
  }
  //operate matrix and 

  //std::cout << " a" << std::endl;
  AT_u.Transpose(A_u);
  //std::cout << " b" << std::endl;
  MB_u = (AT_u) * (A_u);
  std::cout << " Start to Invert U Matrix" << std::endl;
  MB_u.Invert();
  std::cout << " Finish to Invert U Matrix" << std::endl;
  FM_u = (MB_u) * (AT_u) ;
  //std::cout << " e" << std::endl;
  

  //V-plane
  TMatrixDSparse A_v(m,n);
  TMatrixDSparse AT_v(n,m);
  TMatrixD MB_v(n,n);
  TMatrixD FM_v(n,m);
  
   //treat filter ... 
  for (int i=0;i!=40;i++){
    x[i] = (i-20+0.5)*0.5;
    y[i] =  gv->Eval(x[i]+0.5/4.);
    y[i] += gv->Eval(x[i]-0.5/4.);
  }
  // fill matrix
  for (int i=0;i!=n;i++){
    for (int j=0;j!=40;j++){
      Int_t row = i;
      Int_t column = scale*i-20+j+1;
      if (column >=0 && column < m){
	A_v(column,row) = y[j];
      }
    }
  }
  //operate matrix and 
  AT_v.Transpose(A_v);
  MB_v = (AT_v) * (A_v);
  std::cout << " Start to Invert V Matrix" << std::endl;
  MB_v.Invert();
  std::cout << " Finish to Invert V Matrix" << std::endl;
  FM_v = (MB_v) * (AT_v) ;

  

  TMatrixDSparse A_w(m,n);
  TMatrixDSparse AT_w(n,m);
  TMatrixD MB_w(n,n);
  TMatrixD FM_w(n,m);

   //treat filter ... 
  for (int i=0;i!=40;i++){
    x[i] = (i-20+0.5)*0.5;
    y[i] =  gw->Eval(x[i]+0.5/4.);
    y[i] += gw->Eval(x[i]-0.5/4.);
  }
  // fill matrix
  for (int i=0;i!=n;i++){
    for (int j=0;j!=40;j++){
      Int_t row = i;
      Int_t column = scale*i-20+j+1;
      if (column >=0 && column < m){
	A_w(column,row) = y[j];
      }
    }
  }
  //operate matrix and 
  AT_w.Transpose(A_w);
  MB_w = (AT_w) * (A_w);
  std::cout << " Start to Invert W Matrix" << std::endl;
  MB_w.Invert();
  std::cout << " End to Invert W Matrix" << std::endl;
  FM_w = (MB_w) * (AT_w) ;

  

  const Frame& frame1 = fds.get();
  size_t ntraces = frame1.traces.size();


  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();

    TVectorD rx(n);
    TVectorD ry(m);
    
    // fill raw data ... 
    for (int i=0;i!=m;i++){
      ry[i] = trace.charge.at(i);
    }

    //if (chid ==447) ry.Print();

    TH1F *htemp;
    TGraph *gfilter;
    if (chid < nwire_u){
      htemp = hu[chid];
      rx = FM_u * ry;
    }else if (chid < nwire_u + nwire_v){
      htemp = hv[chid - nwire_u];
      rx = FM_v * ry;
    }else{
      htemp = hw[chid - nwire_u - nwire_v];
      rx = FM_w * ry;
    }
    for (int i = 0;i!=n;i++){
      htemp->SetBinContent(i+1,rx[i]);
    }
    
    Trace t;
    t.chid = chid;
    t.tbin = tbin;
    t.charge.resize(bins_per_frame, 0.0);
    for (int j=0;j!=bins_per_frame;j++){
      t.charge.at(j) = htemp->GetBinContent(j+1);
    }
    frame.traces.push_back(t);

  }


  
  frame.index = frame_number;
  return frame.index;
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

  delete hfilter_u;
  delete hfilter_v;
  delete hfilter_w;

  delete hfilter_time_u;
  delete hfilter_time_v;
  delete hfilter_time_w;
  
  delete gu;
  delete gv;
  delete gw;
  
 
}
