#include "WireCell2dToy/DataSignalWien_2D.h"

#include "WireCellData/GeomWire.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "TRandom.h"
#include "TF1.h"

using namespace WireCell;

WireCell2dToy::DataSignalWien2DFDS::DataSignalWien2DFDS(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds,WireCell::ChirpMap& umap, WireCell::ChirpMap& vmap, WireCell::ChirpMap& wmap, int bins_per_frame1, int nframes_total, float time_offset_uv, float time_offset_uw, float overall_time_offset)
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

  hu = new TH1F("U4","U4",nbin,0,nbin);
  hv = new TH1F("V4","V4",nbin,0,nbin);
  hw = new TH1F("W4","W4",nbin,0,nbin);

  #include "data_70_2.txt"  //70kV 2D deconvolution for U
  
  gu = new TGraph*[4];
  gu[0] = new TGraph(5000,xu,yu0);
  gu[1] = new TGraph(5000,xu,yu1);
  gu[2] = new TGraph(5000,xu,yu2);
  gu[3] = new TGraph(5000,xu,yu3);
  

  gv = new TGraph(5000,xv,yv);
  gw = new TGraph(5000,xw,yw);
}

int WireCell2dToy::DataSignalWien2DFDS::size() const{
  return max_frames;
}

void WireCell2dToy::DataSignalWien2DFDS::Save(){
  // TFile *file = new TFile("temp_wien.root","RECREATE");
  // file->Write();
  // file->Close();
}

int WireCell2dToy::DataSignalWien2DFDS::jump(int frame_number){
  // fill the frame data ... 
  if (frame.index == frame_number) {
    return frame_number;
  }
  frame.clear();
  int scale = nbin/bins_per_frame;

  
  fds.jump(frame_number);
  const Frame& frame1 = fds.get();
  size_t ntraces = frame1.traces.size();

  TVirtualFFT::SetTransform(0);
  hur = new TH1F("hur1","hur1",nbin,0,nbin); // half us tick
  hvr = new TH1F("hvr1","hvr1",nbin,0,nbin); // half us tick
  hwr = new TH1F("hwr1","hwr1",nbin,0,nbin); // half us tick
  
  const int nchannels = nwire_u;
  const int nticks = nbin;
  
  // do U-plane
  float scale_u = 1.51/1.16*0.91*0.85;
  float scale_v = 1.251/1.074*0.91*0.85;
  double rho_res[7][nticks], phi_res[7][nticks];

  for (int j=0;j!=7;j++){
    TGraph *gtemp;
    if (j==0 || j==6){
      gtemp = gu[3];
    }else if (j==1 || j==5){
      gtemp = gu[2];
    }else if (j==2 || j==4){
      gtemp = gu[1];
    }else if (j==3){
      gtemp = gu[0];
    }
    for (int i=0; i!=nbin; i++){  
      double time = hur->GetBinCenter(i+1)/2.-50 ;
      //*** scale factors for 70kV ***//
      double x = time;
      if (x > -35 && x  < 15){
	if (gtemp->Eval(x)>0 && x<0){
	  hur->SetBinContent(i+1,gtemp->Eval(x)/scale_u*1.0); //70kV
	}else{
	  hur->SetBinContent(i+1,gtemp->Eval(x)/scale_u); //70kV
	}
      }else{
	hur->SetBinContent(i+1,0);
      }
    }
    hmr_u = hur->FFT(0,"MAG");
    hpr_u = hur->FFT(0,"PH");
    
    for (Int_t i=0;i!=nticks;i++){
      rho_res[j][i] = hmr_u->GetBinContent(i+1);
      phi_res[j][i] = hpr_u->GetBinContent(i+1);
    }

    // std::cout << "abc " << j << std::endl;
    
    delete hmr_u;
    delete hpr_u;
  }

  for (int i=0; i!=nbin; i++){  
    double time = hur->GetBinCenter(i+1)/2.-50 ;
    double x = time;
    x = time-time_offset_uv;  //70kV
    if (x > -35 && x  < 15){
      if (gv->Eval(x)>0 && x<0){
	hvr->SetBinContent(i+1,gv->Eval(x)/scale_v*1.0);  //70kV
      }else{
	hvr->SetBinContent(i+1,gv->Eval(x)/scale_v);  //70kV
      }
    }
    x = time-time_offset_uw;  //70kV
    if (x > -35 && x  < 15){
      hwr->SetBinContent(i+1,gw->Eval(x));
    } 
  } 
  hmr_v = hvr->FFT(0,"MAG");
  hmr_w = hwr->FFT(0,"MAG");
  hpr_v = hvr->FFT(0,"PH");
  hpr_w = hwr->FFT(0,"PH");
    

  //  std::cout << "abc1 " << std::endl;

  TH1 *hm = 0;
  TH1 *hp = 0;
  double value_re[9600];
  double value_im[9600];
  TVirtualFFT *ifft;
  TH1 *fb = 0;
  int n = nbin;
  TH1 *hmr, *hpr;
  
  // do FFT again to remove response function and apply filter

  TF1 *filter_u = new TF1("filter_u","(x>0.0)*gaus*exp(-0.5*pow(x/[3],[4]))");
  double par[5]={1.73/0.959301, 1.69, 1.55, 0.19, 3.75};
  filter_u->SetParameters(par);

  TF1 *filter_v = new TF1("filter_v","(x>0.0)*gaus*exp(-0.5*pow(x/[3],[4]))");
  double par1[5]={1.74/0.941034, 1.46, 1.33, 0.23, 4.89};
  filter_v->SetParameters(par1);

  TF1 *filter_w = new TF1("filter_y","(x>0.0)*[0]*exp(-0.5*pow(pow((x-[1])/[2],2),[3]))");
  double par2[4]={1.03/0.995635, 0.08, 0.15, 2.17};
  filter_w->SetParameters(par2);
  
  std::vector<std::vector<double>> rho_u, phi_u,result_re,result_im;
  //double result_re[nchannels][nticks],result_im[nchannels][nticks];
  //float rho_u[nchannels][nticks],phi_u[nchannels][nticks];
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
    
    TH1F *htemp;
    TF1 *filter;
    int flag_2d = 0; // do 2-D deconvolution?
    if (chid < nwire_u){
      htemp = hu;
      hmr = hmr_u;
      hpr = hpr_u;
      filter = filter_u;
      flag_2d = 1;
    }else if (chid < nwire_u + nwire_v){
      htemp = hv;
      hmr = hmr_v;
      hpr = hpr_v;
      filter = filter_v;
    }else{
      htemp = hw;
      hmr = hmr_w;
      hpr = hpr_w;
      filter = filter_w;
    }
    htemp->Reset();
    
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp->SetBinContent(tt,trace.charge.at(i));
    }    
    hm = htemp->FFT(0,"MAG");
    hp = htemp->FFT(0,"PH");
    if (flag_2d ==1){
      // save 
      for (Int_t j=0;j!=nticks;j++){
  	rho_u[chid][j] = hm->GetBinContent(j+1);
  	phi_u[chid][j] = hp->GetBinContent(j+1);
      }
      tbin_save[chid] = tbin;

      delete hm;
      delete hp;
    }else{
      for (int i=0;i!=nbin;i++){
  	double freq;
  	if (i < nbin/2.){
  	  freq = i/(1.*nbin)*2.;
  	}else{
  	  freq = (nbin - i)/(1.*nbin)*2.;
  	}
	
  	double rho = hm->GetBinContent(i+1)/hmr->GetBinContent(i+1)*filter->Eval(freq);
	//	if (freq < 0.0045) rho = 0;
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
  	htemp->SetBinContent(i+1,fb->GetBinContent(i+1)/( 14.*4096./2000.));
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
      
      
      if (chid < nwire_u){
  	uplane_rms[chid] = rms1*scale;
      }else if (chid < nwire_u + nwire_v){
  	vplane_rms[chid - nwire_u] = rms1*scale;
      }else{
  	wplane_rms[chid - nwire_u - nwire_v] = rms1*scale;
      }
    }
  }
   
  std::cout << "2D deconvolution for U-plane" << std::endl;
  // DO 2D deconvolution for U-plane 
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
      if (temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j]>0){
  	temp2_re[j] = (temp_re[j]*temp1_re[j]+temp_im[j]*temp1_im[j])/m/
  	  (temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j]);
  	temp2_im[j] = (temp_im[j]*temp1_re[j]-temp_re[j]*temp1_im[j])
  	  /m/(temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j]);
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

      //     if (freq > 0.0045){
      result_re[j][i] = temp3_re[shift]/nticks*filter_u->Eval(freq); //filter out the high frequency in time ... 
      result_im[j][i] = temp3_im[shift]/nticks*filter_u->Eval(freq);
      // }else{
    //   	result_re[j][i] = 0.;
    //   	result_im[j][i] = 0.;
    // }
    }
    
    delete ifft;
    delete ifft3;
  }

  // put results back ... 

  for (Int_t i=0;i!=nchannels;i++){
    int n = nticks;
    TVirtualFFT *ifft2 = TVirtualFFT::FFT(1,&n,"C2R M K");
    double temp_re[nticks],temp_im[nticks];
    for (int j=0;j!=nticks;j++){
      temp_re[j] = result_re[i][j];
      temp_im[j] = result_im[i][j];
    }
    ifft2->SetPointsComplex(temp_re,temp_im);
    ifft2->Transform();
    TH1 *fb = 0;
    fb = TH1::TransformHisto(ifft2,fb,"Re");

    // 
    TH1F *htemp = hu;
    for (int i=0;i!=nbin;i++){
      htemp->SetBinContent(i+1,fb->GetBinContent(i+1)/( 14.*4096./2000.));
    }
    
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
    int chid = i;
  
    //save into frame
    Trace t;
    t.chid = chid;
    t.tbin = tbin_save[chid];
    t.charge.resize(bins_per_frame, 0.0);
    for (int j=0;j!=bins_per_frame;j++){
      t.charge.at(j) = 0;
      for (int k=0;k!=scale;k++){
  	t.charge.at(j) += htemp->GetBinContent(scale*j+k+1);
      }
    }
    frame.traces.push_back(t);
    
    //calculate rms, this is to be used for threshold purpose
    float rms = 0, rms1 = 0,rms2 = 0;
    int start=-1, end=-1;
    if (chid < nwire_u){
      if (umap.find(chid)!=umap.end()){
  	start = umap[chid].first;
  	end = umap[chid].second;
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
      
      
    if (chid < nwire_u){
      uplane_rms[chid] = rms1*scale;
    }
    
    delete fb;
    delete ifft2;
  }

  //  std::cout << "abc: " << std::endl;


  delete filter_u;
  delete filter_v;
  delete filter_w;
  //  delete hmr_u;
  delete hmr_v;
  delete hmr_w;
  //delete hpr_u;
  delete hpr_v;
  delete hpr_w;
  
  frame.index = frame_number;
  return frame.index;
}

WireCell2dToy::DataSignalWien2DFDS::~DataSignalWien2DFDS(){
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
  
  for (int i=0;i!=4;i++){
  delete gu[i];
  } 
  delete [] gu;
  delete gv;
  delete gw;

  delete hur;
  delete hvr;
  delete hwr;
}
