#include "WireCell2dToy/uBooNE_Data_ROI.h"

#include "WireCellData/GeomWire.h"
#include "TVirtualFFT.h"
#include "TF1.h"

using namespace WireCell;

WireCell2dToy::uBooNEDataROI::uBooNEDataROI(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, WireCell::ChirpMap& umap, WireCell::ChirpMap& vmap, WireCell::ChirpMap& wmap)
  : fds(fds)
  , gds(gds)
  , umap(umap)
  , vmap(vmap)
  , wmap(wmap)
{
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  find_ROI_by_itself();
}

WireCell2dToy::uBooNEDataROI::~uBooNEDataROI()
{
  
}

void WireCell2dToy::uBooNEDataROI::find_ROI_by_itself(){
  

  const int nbins = fds.Get_Bins_Per_Frame();

  // make a histogram for induction planes and fold it with a low-frequency stuff
  // if collection, just fill the histogram ... 
  TH1F *hresult = new TH1F("hresult","hresult",nbins,0,nbins);
  TF1 *filter_low = new TF1("filter_low","1-exp(-pow(x/0.02,2))");

  // load the data and do the convolution ... 
  const Frame& frame1 = fds.get();
  size_t ntraces = frame1.traces.size();
  double value_re[nbins];
  double value_im[nbins];
  
  for (int i=0;i!=ntraces;i++){
    const Trace& trace = frame1.traces[i];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nticks = trace.charge.size();
    hresult->Reset();
    for (int j=0;j!=nticks;j++){
      hresult->SetBinContent(j+1,trace.charge.at(j));
    }
    
    if (chid < nwire_u + nwire_v){
      // induction planes fold with low frequench stuff
      TH1 *hm = hresult->FFT(0,"MAG");
      TH1 *hp = hresult->FFT(0,"PH");
      TVirtualFFT *ifft2 = TVirtualFFT::FFT(1,&nticks,"C2R M K");

      for (int j=0;j!=nticks;j++){
	Double_t freq;
	if (j < nticks/2.){
	  freq = j/(1.*nticks)*2.;
	}else{
	  freq = (nticks - j)/(1.*nticks)*2.;
	}
	
	value_re[j] = hm->GetBinContent(j+1) * cos(hp->GetBinContent(j+1)) / nticks * filter_low->Eval(freq);
	value_im[j] = hm->GetBinContent(j+1) * sin(hp->GetBinContent(j+1)) / nticks * filter_low->Eval(freq);
      }
      ifft2->SetPointsComplex(value_re,value_im);
      ifft2->Transform();
      TH1 *fb = 0;
      fb = TH1::TransformHisto(ifft2,fb,"Re");
      
      for (int j=0;j!=nticks;j++){
	hresult->SetBinContent(j+1,fb->GetBinContent(j+1));
      }

      delete fb;
      delete hm;
      delete hp;
      delete ifft2;
    }
    restore_baseline(hresult);
    
    //std::cout << chid << " " << cal_rms(hresult,chid) << std::endl;

    // now find ROI ... 


  }

  
  delete hresult;
  delete filter_low;
}



void WireCell2dToy::uBooNEDataROI::restore_baseline(TH1F *htemp){
  //correct baseline 
  double max = htemp->GetMaximum();
  double min = htemp->GetMinimum();
  int nbin_b = max - min;
  if (nbin_b ==0) nbin_b = 1;
  TH1F *h1 = new TH1F("h1","h1",nbin_b,min,max);
  int nbin = htemp->GetNbinsX();
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



double WireCell2dToy::uBooNEDataROI::cal_rms(TH1F *htemp, int chid){
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
