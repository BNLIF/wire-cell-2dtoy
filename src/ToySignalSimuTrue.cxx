#include "WCP2dToy/ToySignalSimuTrue.h"

#include "WCPData/GeomWire.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "TRandom.h"
#include "TF1.h"

using namespace WCP;

WCP2dToy::ToySignalSimuTrueFDS::ToySignalSimuTrueFDS(WCP::FrameDataSource& fds1, const WCP::GeomDataSource& gds,
							  int bins_per_frame1, int nframes_total, int flag_smear)
  : fds(&fds1)
  , max_frames(nframes_total)
  , flag_smear(flag_smear)
  , gds(&gds)
  , dgds(0)
  , gds_flag(0)
{  
  bins_per_frame = bins_per_frame1;

  // GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  // GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  // GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  // nwire_u = wires_u.size();
  // nwire_v = wires_v.size();
  // nwire_w = wires_w.size();
  
  nbin = fds1.Get_Bins_Per_Frame();
 
  // hu = new TH1F*[nwire_u];
  // hv = new TH1F*[nwire_v];
  // hw = new TH1F*[nwire_w];
  
  // for (int i=0;i!=nwire_u;i++){
  //   hu[i] = new TH1F(Form("U5_%d",i),Form("U5_%d",i),nbin,0,nbin);
  // }
  // for (int i=0;i!=nwire_v;i++){
  //   hv[i] = new TH1F(Form("V5_%d",i),Form("V5_%d",i),nbin,0,nbin);
  // }
  // for (int i=0;i!=nwire_w;i++){
  //   hw[i] = new TH1F(Form("W5_%d",i),Form("W5_%d",i),nbin,0,nbin);
  // }

  hu = new TH1F("U5","U5",nbin,0,nbin);
  hv = new TH1F("V5","V5",nbin,0,nbin);
  hw = new TH1F("W5","W5",nbin,0,nbin);
  
  //define filters
  filter_g = new TF1("filger_g","1./sqrt(2.*3.1415926)/[0]*exp(-x*x/2./[0]/[0])");
  double par3[1] = {2./2.2};
  filter_g->SetParameters(par3);
  
  hfilter_time_gaus =new TH1F("hfilter_time_gaus1","hfilter_time_gaus1",nbin,0,nbin);
  for (int i=0;i!=nbin;i++){
    double xx = hfilter_time_gaus->GetBinCenter(i+1)/2.-nbin/4.;
    hfilter_time_gaus->SetBinContent(i+1,filter_g->Eval(xx));
  }
  hfilter_time_gaus->Scale(1./hfilter_time_gaus->GetSum());
  hfilter_gaus = 0;
  hfilter_gaus = hfilter_time_gaus->FFT(0,"MAG");
  
}


WCP2dToy::ToySignalSimuTrueFDS::ToySignalSimuTrueFDS(WCP::FrameDataSource& fds1, const WCP::DetectorGDS& gds,
							  int bins_per_frame1, int nframes_total, int flag_smear)
  : fds(&fds1)
  , max_frames(nframes_total)
  , flag_smear(flag_smear)
  , dgds(&gds)
  , gds(0)
  , gds_flag(1)
{  
  bins_per_frame = bins_per_frame1;

  // GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  // GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  // GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  // nwire_u = wires_u.size();
  // nwire_v = wires_v.size();
  // nwire_w = wires_w.size();
  
  nbin = fds1.Get_Bins_Per_Frame();
  
  // hu = new TH1F*[nwire_u];
  // hv = new TH1F*[nwire_v];
  // hw = new TH1F*[nwire_w];
  
  // for (int i=0;i!=nwire_u;i++){
  //   hu[i] = new TH1F(Form("U5_%d",i),Form("U5_%d",i),nbin,0,nbin);
  // }
  // for (int i=0;i!=nwire_v;i++){
  //   hv[i] = new TH1F(Form("V5_%d",i),Form("V5_%d",i),nbin,0,nbin);
  // }
  // for (int i=0;i!=nwire_w;i++){
  //   hw[i] = new TH1F(Form("W5_%d",i),Form("W5_%d",i),nbin,0,nbin);
  // }

  hu = new TH1F("U5","U5",nbin,0,nbin);
  hv = new TH1F("V5","V5",nbin,0,nbin);
  hw = new TH1F("W5","W5",nbin,0,nbin);
  
  //define filters
  filter_g = new TF1("filger_g","1./sqrt(2.*3.1415926)/[0]*exp(-x*x/2./[0]/[0])");
  double par3[1] = {2./2.2};
  filter_g->SetParameters(par3);
  
  hfilter_time_gaus =new TH1F("hfilter_time_gaus1","hfilter_time_gaus1",nbin,0,nbin);
  for (int i=0;i!=nbin;i++){
    double xx = hfilter_time_gaus->GetBinCenter(i+1)/2.-nbin/4.;
    hfilter_time_gaus->SetBinContent(i+1,filter_g->Eval(xx));
  }
  hfilter_time_gaus->Scale(1./hfilter_time_gaus->GetSum());
  hfilter_gaus = 0;
  hfilter_gaus = hfilter_time_gaus->FFT(0,"MAG");
  
}




int WCP2dToy::ToySignalSimuTrueFDS::size() const{
  return max_frames;
}

void WCP2dToy::ToySignalSimuTrueFDS::Save(){
  TFile *file = new TFile("temp_true.root","RECREATE");
  // for (int i=0;i!=nwire_u;i++){
  //   TH1F *huu = (TH1F*)hu[i]->Clone(Form("U1_%d",i));
  // }
  // for (int i=0;i!=nwire_v;i++){
  //   TH1F *hvv = (TH1F*)hv[i]->Clone(Form("V1_%d",i));
  // }
  // for (int i=0;i!=nwire_w;i++){
  //   TH1F *hww = (TH1F*)hw[i]->Clone(Form("W1_%d",i));
  // }

    
  file->Write();
  file->Close();
}

int WCP2dToy::ToySignalSimuTrueFDS::jump(int frame_number){
  // do simulation
  // for (int i=0;i!=nwire_u;i++){
  //   hu[i]->Reset();
  // }
  // for (int i=0;i!=nwire_v;i++){
  //   hv[i]->Reset();
  // }
  // for (int i=0;i!=nwire_w;i++){
  //   hw[i]->Reset();
  // }
  
  if (frame.index == frame_number) {
    return frame_number;
  }

  // start FFT to convolute with response function
  TVirtualFFT::SetTransform(0);
  TH1 *hm = 0;
  TH1 *hp = 0;
    
  double value_re[9600]; // hack for now
  double value_im[9600];
  int  n  = nbin;
  TVirtualFFT *ifft=0;
  TH1 *fb = 0;

  frame.clear();
  int scale = nbin/bins_per_frame;
  
  fds->jump(frame_number);
  
  // std::cout << "Xin1 " << " " << nbin << " " << bins_per_frame << " " << scale << std::endl; 

  //fill in the data ... 
  const Frame& frame1 = fds->get();
  size_t ntraces = frame1.traces.size();
  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    
    
    TH1F *htemp;
    if (gds_flag == 0){
      // regular gds
      
      WirePlaneType_t plane = gds->by_channel(chid).at(0)->plane();
      if (plane == WirePlaneType_t(0)){
	htemp = hu;
      }else if (plane == WirePlaneType_t(1)){
	htemp = hv;
      }else if (plane == WirePlaneType_t(2)){
	htemp= hw;
      }
    }else{
      // detector gds 
      //dgds->by_channel(chid);
      WirePlaneType_t plane = dgds->by_channel(chid).at(0)->plane();
      if (plane == WirePlaneType_t(0)){
	htemp = hu;
      }else if (plane == WirePlaneType_t(1)){
	htemp = hv;
      }else if (plane == WirePlaneType_t(2)){
	htemp= hw;
      }
    }
    // if (chid < nwire_u){
    //   //htemp = hu[chid];
    //   htemp = hu;
    // }else if (chid < nwire_u + nwire_v){
    //   //htemp = hv[chid - nwire_u];
    //   htemp = hv;
    // }else{
    //   //htemp = hw[chid - nwire_u - nwire_v];
    //   htemp = hw;
    // }
    
    htemp->Reset();
    
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
      if (tt <= nbin)
	htemp->SetBinContent(tt,vcharge.at(j));
    }
    
    if (flag_smear == 1){
      hm = htemp->FFT(0,"MAG");
      hp = htemp->FFT(0,"PH");
      for (int j=0;j!=nbin;j++){
	double rho = hm->GetBinContent(j+1)*hfilter_gaus->GetBinContent(j+1);
	double phi = hp->GetBinContent(j+1);
	value_re[j] = rho*cos(phi)/nbin;
	value_im[j] = rho*sin(phi)/nbin;
      }
      ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
      ifft->SetPointsComplex(value_re,value_im);
      ifft->Transform();
      fb = TH1::TransformHisto(ifft,0,"Re");
      for (int j=0;j!=nbin;j++){
	int content = fb->GetBinContent(j+1) ;
	htemp->SetBinContent(j+1,content);
      }
      
      
    }


    //std::cout << "Xin2 " << std::endl; 

    delete hm;
    delete hp;
    if (ifft!=0)
      delete ifft;
    delete fb;
    
    Trace t;
    t.chid = chid;
    t.tbin = 0;
    t.charge.resize(bins_per_frame, 0.0);
    for (int j=0;j!=bins_per_frame;j++){
      t.charge.at(j) = 0;
      for (int k=0;k!=scale;k++){
	t.charge.at(j) += htemp->GetBinContent(scale*j+k+1);
      }
    }
    frame.traces.push_back(t);
    
    //std::cout << "Xin3 " << std::endl; 

    //std::cout << chid << std::endl;
    // std::cout << nwire_u << " " << nwire_v << " " << nwire_w << std::endl;
  }
  
  
  
 
  
  
  // //U-plane first
  // for (int i=0;i!=nwire_v;i++){
    
  // }
  // //V-plane
  // for (int i=0;i!=nwire_w;i++){
  //   hm = hw[i]->FFT(hm,"MAG");
  //   hp = hw[i]->FFT(hp,"PH");
  //   for (int j=0;j!=nbin;j++){
  //     double rho = hm->GetBinContent(j+1)*hfilter_gaus->GetBinContent(j+1);
  //     double phi = hp->GetBinContent(j+1);
  //     value_re[j] = rho*cos(phi)/nbin;
  //     value_im[j] = rho*sin(phi)/nbin;
  //   }
  //   ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
  //   ifft->SetPointsComplex(value_re,value_im);
  //   ifft->Transform();
  //   fb = TH1::TransformHisto(ifft,fb,"Re");
  //   for (int j=0;j!=nbin;j++){
  //     int content = fb->GetBinContent(j+1) ;
  //     hw[i]->SetBinContent(j+1,content);
  //   }
  // }

  // //W-plane
  // for (int i=0;i!=nwire_u;i++){
  //   hm = hu[i]->FFT(hm,"MAG");
  //   hp = hu[i]->FFT(hp,"PH");
  //   for (int j=0;j!=nbin;j++){
  //     double rho = hm->GetBinContent(j+1)*hfilter_gaus->GetBinContent(j+1);
  //     double phi = hp->GetBinContent(j+1);
  //     value_re[j] = rho*cos(phi)/nbin;
  //     value_im[j] = rho*sin(phi)/nbin;
  //   }
  //   ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
  //   ifft->SetPointsComplex(value_re,value_im);
  //   ifft->Transform();
  //   fb = TH1::TransformHisto(ifft,fb,"Re");
  //   for (int j=0;j!=nbin;j++){
  //     int content = fb->GetBinContent(j+1) ;
  //     hu[i]->SetBinContent(j+1,content);
  //   }
  // }



  // done with FFT


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
  //     t.charge.at(j)=0;
  //     for (int k=0;k!=scale;k++){
  // 	t.charge.at(j) += hv[i]->GetBinContent(scale*j+k+1);
  //     }
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
  //     t.charge.at(j) = 0;
  //     for (int k=0;k!=scale;k++){
  // 	t.charge.at(j) += hw[i]->GetBinContent(scale*j+k+1);
  //     }
  //   }
  //   frame.traces.push_back(t);
  // }
  
  
  frame.index = frame_number;
  return frame.index;
}

WCP2dToy::ToySignalSimuTrueFDS::~ToySignalSimuTrueFDS(){
  
  fds = 0;
  

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

  //frame.clear();

  delete filter_g;
  delete hfilter_time_gaus;
  delete hfilter_gaus;
}
