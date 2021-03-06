#include "WCP2dToy/DataSignalWien.h"

#include "WCPData/GeomWire.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "TRandom.h"
#include "TF1.h"

using namespace WCP;

WCP2dToy::DataSignalWienFDS::DataSignalWienFDS(WCP::FrameDataSource& fds, const WCP::GeomDataSource& gds,WCP::ChirpMap& umap, WCP::ChirpMap& vmap, WCP::ChirpMap& wmap, int bins_per_frame1, int nframes_total, float time_offset_uv, float time_offset_uw, float overall_time_offset)
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

  // hu = new TH1F*[nwire_u];
  // hv = new TH1F*[nwire_v];
  // hw = new TH1F*[nwire_w];
  
  // for (int i=0;i!=nwire_u;i++){
  //   hu[i] = new TH1F(Form("U4_%d",i),Form("U4_%d",i),nbin,0,nbin);
  // }
  // for (int i=0;i!=nwire_v;i++){
  //   hv[i] = new TH1F(Form("V4_%d",i),Form("V4_%d",i),nbin,0,nbin);
  // }
  // for (int i=0;i!=nwire_w;i++){
  //   hw[i] = new TH1F(Form("W4_%d",i),Form("W4_%d",i),nbin,0,nbin);
  // }

  hu = new TH1F("U4","U4",nbin,0,nbin);
  hv = new TH1F("V4","V4",nbin,0,nbin);
  hw = new TH1F("W4","W4",nbin,0,nbin);

  //#include "data.txt"   //58kV
  #include "data_70.txt"  //70kV
  //#include "data_70_1.txt"  //70kV with Xin's tune on V and Y
  
  gu = new TGraph(5000,xu,yu);
  gv = new TGraph(5000,xv,yv);
  gw = new TGraph(5000,xw,yw);

  hur = new TH1F("hur1","hur1",nbin,0,nbin); // half us tick
  hvr = new TH1F("hvr1","hvr1",nbin,0,nbin); // half us tick
  hwr = new TH1F("hwr1","hwr1",nbin,0,nbin); // half us tick
  
  for (int i=0; i!=nbin; i++){  
    double time = hur->GetBinCenter(i+1)/2.-50 ;
    // hur->SetBinContent(i+1,gu->Eval(time-overall_time_offset));
    // hvr->SetBinContent(i+1,gv->Eval(time-time_offset_uv-overall_time_offset));
    // hwr->SetBinContent(i+1,gw->Eval(time-time_offset_uw-overall_time_offset));

    //*** scale factors for 58kV ***//
    //float scale = 0.86; 
    //float scale_u = 1.51/1.16*0.91;
    //float scale_v = 1.251/1.074*0.91;
    //double x = scale*(time-overall_time_offset);

    //*** scale factors for 70kV ***//
    float scale_u = 1.51/1.16*0.91*0.85;
    float scale_v = 1.251/1.074*0.91*0.85;
    double x = time;

    if (x > -35 && x  < 15){
      if (gu->Eval(x) > 0 && x < 0){
	//hur->SetBinContent(i+1,gu->Eval(x)*0.5/0.8/scale_u); //58kV
	hur->SetBinContent(i+1,gu->Eval(x)*0.6/scale_u); //70kV
      }else{
	//hur->SetBinContent(i+1,gu->Eval(x)/0.8/scale_u); //58kV
	hur->SetBinContent(i+1,gu->Eval(x)/scale_u); //70kV
      }
    }

    //x = scale*(time-time_offset_uv-overall_time_offset);  //58kV
    x = time-time_offset_uv;  //70kV

    if (x > -35 && x  < 15){
      if (gv->Eval(x) > 0 && x < 0){
	//hvr->SetBinContent(i+1,gv->Eval(x)*0.6/scale_v);  //58kV
	hvr->SetBinContent(i+1,gv->Eval(x)*0.7/scale_v); //70kV
      }else{
	//hvr->SetBinContent(i+1,gv->Eval(x)/scale_v);  //58kV
	hvr->SetBinContent(i+1,gv->Eval(x)/scale_v);  //70kV
      }
    }
    
    //x = scale*(time-time_offset_uw-overall_time_offset); //58kV
    x = time-time_offset_uw;  //70kV

    if (x > -35 && x  < 15){
      hwr->SetBinContent(i+1,gw->Eval(x));
    } 
  } 
  
  hmr_u = 0;
  hmr_v = 0;
  hmr_w = 0;
  hpr_u = 0;
  hpr_v = 0;
  hpr_w = 0;
  
}

int WCP2dToy::DataSignalWienFDS::size() const{
  return max_frames;
}

void WCP2dToy::DataSignalWienFDS::Save(){
  TFile *file = new TFile("temp_wien.root","RECREATE");
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

int WCP2dToy::DataSignalWienFDS::jump(int frame_number){
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
  
  // do FFT again to remove response function and apply filter

  TF1 *filter_u = new TF1("filter_u","(x>0.0)*exp(-0.5*pow(x/[0],[1]))");
  double par[2]={1.22781e+01/200.*2.,4.96159e+00};
  filter_u->SetParameters(par);

  TF1 *filter_v = new TF1("filter_v","(x>0.0)*exp(-0.5*pow(x/[0],[1]))");
  double par1[2]={1.31649e+01/200.*2.,4.67294e+00};
  filter_v->SetParameters(par1);

  TF1 *filter_w = new TF1("filter_y","(x>0.0)*exp(-0.5*pow(x/[0],[1]))");
  double par2[2]={1.26250e+01/200.*2.,4.48298e+00};
  filter_w->SetParameters(par2);


  // TF1 *filter_u = new TF1("filter_u","(x>0.0)*gaus*exp(-0.5*pow(x/[3],[4]))");
  // double par[5]={1.73/0.959301, 1.69, 1.55, 0.19, 3.75};
  // filter_u->SetParameters(par);
  
  // TF1 *filter_v = new TF1("filter_v","(x>0.0)*gaus*exp(-0.5*pow(x/[3],[4]))");
  // double par1[5]={1.74/0.941034, 1.46, 1.33, 0.23, 4.89};
  // filter_v->SetParameters(par1);

  // TF1 *filter_w = new TF1("filter_y","(x>0.0)*[0]*exp(-0.5*pow(pow((x-[1])/[2],2),[3]))");
  // double par2[4]={1.03/0.995635, 0.08, 0.15, 2.17};
  // filter_w->SetParameters(par2);



      // double rho = hm->GetBinContent(j+1)/hmr->GetBinContent(j+1)*filter_u->Eval(freq);
      
  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    
    //std::cout << tbin << " " << chid << " " << nbins << std::endl;

    TH1F *htemp;
    TF1 *filter;
    
    if (chid < nwire_u){
      htemp = hu;//[chid];
      hmr = hmr_u;
      hpr = hpr_u;
      filter = filter_u;
    }else if (chid < nwire_u + nwire_v){
      htemp = hv;//[chid - nwire_u];
      hmr = hmr_v;
      hpr = hpr_v;
      filter = filter_v;
    }else{
      htemp = hw;//[chid - nwire_u - nwire_v];
      hmr = hmr_w;
      hpr = hpr_w;
      filter = filter_w;
    }
    htemp->Reset();
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;// + overall_time_shift * 2;
      // if (tt > nbin) tt -= nbin;
      // if (tt < 1 ) tt += nbin;
      htemp->SetBinContent(tt,trace.charge.at(i));
    }
        
    hm = htemp->FFT(0,"MAG");
    hp = htemp->FFT(0,"PH");
    
    for (int i=0;i!=nbin;i++){
      
      double freq;
      if (i < nbin/2.){
       	freq = i/(1.*nbin)*2.;
      }else{
       	freq = (nbin - i)/(1.*nbin)*2.;
      }

      double rho = hm->GetBinContent(i+1)/hmr->GetBinContent(i+1)*filter->Eval(freq);
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
        
    // // old method to calculate RMS
    // for (int i=0;i!=htemp->GetNbinsX();i++){
    //   if (i < start || i > end){
    //  	rms += pow(htemp->GetBinContent(i+1),2);
    //  	rms2 ++;
    //   }
    // }
    // if (rms2!=0){
    //   rms = sqrt(rms/rms2);
    // }else{
    //   rms = 0;   
    // }

    // //try to exclude signal
    // rms2 = 0;
    //  for (int i=0;i!=htemp->GetNbinsX();i++){
    //   if (i < start || i > end){
    // 	if (fabs(htemp->GetBinContent(i+1)) < 5.0*rms){
    // 	  rms1 += pow(htemp->GetBinContent(i+1),2);
    // 	  rms2 ++;
    // 	}
    //   }
    // }
    // if (rms2!=0){
    //   rms1 = sqrt(rms1/rms2);
    // }else{
    //   rms1 = 0;   
    // }
    
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
 
  delete filter_u;
  delete filter_v;
  delete filter_w;

  delete hmr_u;
  delete hmr_v;
  delete hmr_w;
  delete hpr_u;
  delete hpr_v;
  delete hpr_w;
  
  frame.index = frame_number;
  return frame.index;
}

WCP2dToy::DataSignalWienFDS::~DataSignalWienFDS(){
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

  delete gu;
  delete gv;
  delete gw;

  delete hur;
  delete hvr;
  delete hwr;
}
