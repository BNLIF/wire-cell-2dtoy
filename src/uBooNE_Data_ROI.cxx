#include "WireCell2dToy/uBooNE_Data_ROI.h"

#include "WireCellData/GeomWire.h"
#include "TVirtualFFT.h"
#include "TF1.h"

using namespace WireCell;

Double_t WireCell2dToy::uBooNEDataROI::local_ave(TH1F *h1, Int_t bin, Int_t width){
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


void WireCell2dToy::uBooNEDataROI::Clear(){
  self_rois_u.clear();
  self_rois_v.clear();
  self_rois_w.clear();
  
  loose_rois_u.clear();
  loose_rois_v.clear();
  loose_rois_w.clear();

  uplane_rms.clear();
  vplane_rms.clear();
  wplane_rms.clear();
}


Int_t WireCell2dToy::uBooNEDataROI::find_ROI_end(TH1F *h1, Int_t bin, Double_t th){
  Int_t end = bin;
  Double_t content = h1->GetBinContent(end+1);
  while(content>th){
    end ++;
    if (end >=h1->GetNbinsX()){
      content = h1->GetBinContent(end - h1->GetNbinsX()+1);
    }else{
      content = h1->GetBinContent(end+1);
    }
    if (end == h1->GetNbinsX()) break;
  }

  while(local_ave(h1,end+1,1) < local_ave(h1,end,1)){
    end++;
    if (end == h1->GetNbinsX()) break;
  } 
  return end;
}

Int_t WireCell2dToy::uBooNEDataROI::find_ROI_begin(TH1F *h1, Int_t bin, Double_t th){
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
    if (begin == 0) break;
  }
  
  // calculate the local average
  // keep going and find the minimum
  while( local_ave(h1,begin-1,1) < local_ave(h1,begin,1)){
    begin --;
    if (begin == 0) break;
  }
  
  return begin;
}



WireCell2dToy::uBooNEDataROI::uBooNEDataROI(WireCell::FrameDataSource& raw_fds,WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, WireCell::ChirpMap& umap, WireCell::ChirpMap& vmap, WireCell::ChirpMap& wmap, std::set<int>& lf_noisy_channels)
  : fds(fds)
  , raw_fds(raw_fds)
  , gds(gds)
  , umap(umap)
  , vmap(vmap)
  , wmap(wmap)
  , lf_noisy_channels(lf_noisy_channels)
  , nwire_u(0)
  , nwire_v(0)
  , nwire_w(0)
{
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  // std::cout << umap.size() << " " << vmap.size() << " " << wmap.size() <<  " " << umap[880].first << " " << umap[880].second << std::endl;

  self_rois_u.resize(nwire_u);
  self_rois_v.resize(nwire_v);
  self_rois_w.resize(nwire_w);

  loose_rois_u.resize(nwire_u);
  loose_rois_v.resize(nwire_v);
  loose_rois_w.resize(nwire_w);


  // others_rois_u.resize(nwire_u);
  // others_rois_v.resize(nwire_v);
  // others_rois_w.resize(nwire_w);
  
  // combined_rois_u.resize(nwire_u);
  // combined_rois_v.resize(nwire_v);
  // combined_rois_w.resize(nwire_w);

  uplane_rms.resize(nwire_u);
  vplane_rms.resize(nwire_v);
  wplane_rms.resize(nwire_w);
    

  std::cout << "Finding ROI based on decon itself " << std::endl;
  find_ROI_by_decon_itself(3.6,5.0,0); // 5 sigma, with zero padding
  //std::cout << self_rois_u.size() << " " << self_rois_v.size() << " " << self_rois_w.size() << std::endl;
  std::cout << "Extend Self ROIs" << std::endl;
  extend_ROI_self(5);
  std::cout << "Create more Self ROIs based on connectivity" << std::endl;
  create_ROI_connect_info();
  

  std::cout << "Finding Loose ROI" << std::endl;
  find_ROI_loose();
  std::cout << "Extend Loose ROI" << std::endl;
  extend_ROI_loose();



  //std::cout << "Fidning ROI based on raw itself " << std::endl;
  //find_ROI_by_raw_itself(3.6,5); // 3 sigma with 5 (half) padding
  //std::cout << "Finding ROI based on other two planes " << std::endl;
  //find_ROI_by_others();
  //std::cout << "Extend Two-plane ROIs" << std::endl;
  //extend_ROI_others(5);
  //two plane tiling is not working for  this configuration ... 

  //std::cout << "Merge ROIs " << std::endl;
  //merge_ROIs();
}

void WireCell2dToy::uBooNEDataROI::extend_ROI_loose(){
  // compare the loose one with tight one 
  for(int i=0;i!=nwire_u;i++){
    std::vector<std::pair<int,int>> temp_rois;
    for (int j=0;j!=loose_rois_u.at(i).size();j++){
      int start = loose_rois_u.at(i).at(j).first;
      int end = loose_rois_u.at(i).at(j).second;
      for (int k=0;k!=self_rois_u.at(i).size();k++){
	int temp_start = self_rois_u.at(i).at(k).first;
	int temp_end = self_rois_u.at(i).at(k).second;
	if (start > temp_start && start < temp_end)
	  start = temp_start;
	// loop through all the tight one to examine start
	if (end > temp_start && end < temp_end)
	  end = temp_end; 
	// loop through all the tight one to examine the end
      }
      if (temp_rois.size()==0){
	temp_rois.push_back(std::make_pair(start,end));
      }else{
	if (start < temp_rois.back().second){
	  temp_rois.back().second = end;
	}else{
	  temp_rois.push_back(std::make_pair(start,end));
	}
      }
    }
    loose_rois_u.at(i) = temp_rois;
  }


  for(int i=0;i!=nwire_v;i++){
    std::vector<std::pair<int,int>> temp_rois;
    for (int j=0;j!=loose_rois_v.at(i).size();j++){
      int start = loose_rois_v.at(i).at(j).first;
      int end = loose_rois_v.at(i).at(j).second;
      for (int k=0;k!=self_rois_v.at(i).size();k++){
	int temp_start = self_rois_v.at(i).at(k).first;
	int temp_end = self_rois_v.at(i).at(k).second;
	if (start > temp_start && start < temp_end)
	  start = temp_start;
	// loop through all the tight one to examine start
	if (end > temp_start && end < temp_end)
	  end = temp_end; 
	// loop through all the tight one to examine the end
      }
      if (temp_rois.size()==0){
	temp_rois.push_back(std::make_pair(start,end));
      }else{
	if (start < temp_rois.back().second){
	  temp_rois.back().second = end;
	}else{
	  temp_rois.push_back(std::make_pair(start,end));
	}
      }
    }
    loose_rois_v.at(i) = temp_rois;
  }
  

}

void WireCell2dToy::uBooNEDataROI::create_ROI_connect_info(float asy){
  
  // u 
  for (int i=0;i!=nwire_u-2;i++){
    for (int j=0; j!=self_rois_u.at(i).size();j++){
      int start1 = self_rois_u.at(i).at(j).first;
      int end1 = self_rois_u.at(i).at(j).second;
      int length1 = end1-start1+1;
      for (int k=0; k!=self_rois_u.at(i+2).size();k++){
	int start2 = self_rois_u.at(i+2).at(k).first;
	int end2 = self_rois_u.at(i+2).at(k).second;
	int length2 = end2 - start2 + 1;
	if ( fabs(length2 - length1) < (length2 + length1) * asy){
	  int start3 = (start1+start2)/2.;
	  int end3 = (end1+end2)/2.;
	  if (start3 < end3 && start3 <= end1 && start3 <=end2 && end3 >= start1 && end3 >=start2){
	    // go through existing ones to make sure there is no overlap
	    int flag = 0; 
	    for (int i1 = 0; i1!=self_rois_u.at(i+1).size();i1++){
	      int max_start = start3;
	      if (self_rois_u.at(i+1).at(i1).first > max_start)
		max_start = self_rois_u.at(i+1).at(i1).first;
	      int min_end = end3;
	      if (self_rois_u.at(i+1).at(i1).second < min_end)
		min_end = self_rois_u.at(i+1).at(i1).second ;
	      if (max_start < min_end){
		flag = 1;
		break;
	      }
	    }
	    if (flag == 0)
	      self_rois_u.at(i+1).push_back(std::make_pair(start3,end3));
	  }
	}
      } 
    }
  }
  // v
  for (int i=0;i!=nwire_v-2;i++){
    for (int j=0; j!=self_rois_v.at(i).size();j++){
      int start1 = self_rois_v.at(i).at(j).first;
      int end1 = self_rois_v.at(i).at(j).second;
      int length1 = end1-start1+1;
      for (int k=0; k!=self_rois_v.at(i+2).size();k++){
	int start2 = self_rois_v.at(i+2).at(k).first;
	int end2 = self_rois_v.at(i+2).at(k).second;
	int length2 = end2 - start2 + 1;
	if ( fabs(length2 - length1) < (length2 + length1) * asy){
	  int start3 = (start1+start2)/2.;
	  int end3 = (end1+end2)/2.;
	  if (start3 < end3 && start3 <= end1 && start3 <=end2 && end3 >= start1 && end3 >=start2){
	    // go through existing ones to make sure there is no overlap
	    int flag = 0; 
	    for (int i1 = 0; i1!=self_rois_v.at(i+1).size();i1++){
	      int max_start = start3;
	      if (self_rois_v.at(i+1).at(i1).first > max_start)
		max_start = self_rois_v.at(i+1).at(i1).first;
	      int min_end = end3;
	      if (self_rois_v.at(i+1).at(i1).second < min_end)
		min_end = self_rois_v.at(i+1).at(i1).second ;
	      if (max_start < min_end){
		flag = 1;
		break;
	      }
	    }
	    if (flag == 0)
	      self_rois_v.at(i+1).push_back(std::make_pair(start3,end3));
	  }
	}
      } 
    }
  }
  // w?
  for (int i=0;i!=nwire_w-2;i++){
    for (int j=0; j!=self_rois_w.at(i).size();j++){
      int start1 = self_rois_w.at(i).at(j).first;
      int end1 = self_rois_w.at(i).at(j).second;
      int length1 = end1-start1+1;
      for (int k=0; k!=self_rois_w.at(i+2).size();k++){
	int start2 = self_rois_w.at(i+2).at(k).first;
	int end2 = self_rois_w.at(i+2).at(k).second;
	int length2 = end2 - start2 + 1;
	if ( fabs(length2 - length1) < (length2 + length1) * asy){
	  int start3 = (start1+start2)/2.;
	  int end3 = (end1+end2)/2.;
	  if (start3 < end3 && start3 <= end1 && start3 <=end2 && end3 >= start1 && end3 >=start2){
	    // go through existing ones to make sure there is no overlap
	    int flag = 0; 
	    for (int i1 = 0; i1!=self_rois_w.at(i+1).size();i1++){
	      int max_start = start3;
	      if (self_rois_w.at(i+1).at(i1).first > max_start)
		max_start = self_rois_w.at(i+1).at(i1).first;
	      int min_end = end3;
	      if (self_rois_w.at(i+1).at(i1).second < min_end)
		min_end = self_rois_w.at(i+1).at(i1).second ;
	      if (max_start < min_end){
		flag = 1;
		break;
	      }
	    }
	    if (flag == 0)
	      self_rois_w.at(i+1).push_back(std::make_pair(start3,end3));
	  }
	}
      } 
    }
  }

}

// void WireCell2dToy::uBooNEDataROI::extend_ROI_others(int pad){
//   const int nbins = raw_fds.Get_Bins_Per_Frame();

//   for (int i=0;i!=others_rois_u.size();i++){
//     std::vector<std::pair<int,int>> temp_rois;
//     int temp_begin, temp_end;
//     for (int j=0;j!=others_rois_u.at(i).size();j++){
//       temp_begin = others_rois_u.at(i).at(j).first - pad;
//       if (temp_begin < 0 ) temp_begin = 0;
//       temp_end = others_rois_u.at(i).at(j).second + pad;
//       if (temp_end >= nbins) temp_end = nbins - 1;
//       // merge 
//       if (temp_rois.size() == 0){
// 	temp_rois.push_back(std::make_pair(temp_begin,temp_end));
//       }else{
// 	if (temp_begin < temp_rois.back().second){
// 	  if (temp_end > temp_rois.back().second){
// 	    temp_rois.back().second = temp_end;
// 	  }
// 	}else{
// 	  temp_rois.push_back(std::make_pair(temp_begin,temp_end));
// 	}
//       }
//     }
//     others_rois_u.at(i) = temp_rois;
//   }

//   for (int i=0;i!=others_rois_v.size();i++){
//     std::vector<std::pair<int,int>> temp_rois;
//     int temp_begin, temp_end;
//     for (int j=0;j!=others_rois_v.at(i).size();j++){
//       temp_begin = others_rois_v.at(i).at(j).first - pad;
//       if (temp_begin < 0 ) temp_begin = 0;
//       temp_end = others_rois_v.at(i).at(j).second + pad;
//       if (temp_end >= nbins) temp_end = nbins - 1;
//       // merge 
//       if (temp_rois.size() == 0){
// 	temp_rois.push_back(std::make_pair(temp_begin,temp_end));
//       }else{
// 	if (temp_begin < temp_rois.back().second){
// 	  if (temp_end > temp_rois.back().second){
// 	    temp_rois.back().second = temp_end;
// 	  }
// 	}else{
// 	  temp_rois.push_back(std::make_pair(temp_begin,temp_end));
// 	}
//       }
//     }
//     others_rois_v.at(i) = temp_rois;
//   }



// }


void WireCell2dToy::uBooNEDataROI::find_ROI_loose(int rebin){
  
  // some example thresholds
  Double_t factor = 3.5; // regular threshold
  Double_t max_th = 10000; // maximum threshold
  Double_t factor1 = 0.7; //special threshold
  Int_t short_length = 3; // short length
  Double_t fixed_threshold = 4000;



  // deconvolute with low-frequency filter 
  const int nbins = fds.Get_Bins_Per_Frame();

  // make a histogram for induction planes and fold it with a low-frequency stuff
  // if collection, just fill the histogram ... 
  TH1F *hresult = new TH1F("hresult","hresult",nbins,0,nbins);
  TH1F *hresult_filter = new TH1F("hresult_filter","hresult_filter",int(nbins/rebin),0,int(nbins/rebin));
  TF1 *filter_low = 0;
  
  TF1 *filter_low1 = new TF1("filter_low","1-exp(-pow(x/0.0025,2))");
  TF1 *filter_low2 = new TF1("filter_low","1-exp(-pow(x/0.02,2))");
  
  // load the data and do the convolution ... 
  const Frame& frame1 = fds.get();
  size_t ntraces = frame1.traces.size();
  double value_re[nbins];
  double value_im[nbins];
  for (int i=0;i!=nbins;i++){
    value_re[i] = 0;
    value_im[i] = 0;
  }

  for (int i=0;i!=ntraces;i++){
    const Trace& trace = frame1.traces[i];
    int tbin = trace.tbin;
    int chid = trace.chid;

    // if (lf_noisy_channels.find(chid)!=lf_noisy_channels.end()){
    //   std::cout << chid << " " << 1 << std::endl;
    // }else{
    //   std::cout << chid << " " << 0 << std::endl;
    // }
    
    //decide filter_low 
    filter_low = filter_low1;
    if (lf_noisy_channels.find(chid)!=lf_noisy_channels.end()){
      filter_low = filter_low2;
    }else if (lf_noisy_channels.find(chid+1)!=lf_noisy_channels.end()&&chid!=nwire_u-1){
      filter_low = filter_low2;
    }else if (lf_noisy_channels.find(chid-1)!=lf_noisy_channels.end()&&chid!=nwire_u){
      filter_low = filter_low2;
    }
    

    int nticks = trace.charge.size();
    hresult->Reset();
    hresult_filter->Reset();
    
    int dead_start = -1;
    int dead_end = -1;
    
    if (chid < nwire_u){
      if (umap.find(chid) != umap.end()){
	dead_start = umap[chid].first;
	dead_end = umap[chid].second;
      }
    }else if (chid < nwire_u + nwire_v){
      if (vmap.find(chid-nwire_u) != vmap.end()){
	
	dead_start = vmap[chid-nwire_u].first;
	dead_end = vmap[chid-nwire_v].second;
      }
    }else{
      if (wmap.find(chid-nwire_u-nwire_v) != wmap.end()){
	dead_start = wmap[chid-nwire_u-nwire_v].first;
	dead_end = wmap[chid-nwire_u-nwire_v].second;
      }
    }
    
    // if (chid == 7500) std::cout << 7500 << " " << dead_start << " " << dead_end << std::endl;

    for (int j=0;j!=nticks;j++){
      if (j < dead_start || j > dead_end)
	hresult->SetBinContent(j+1,trace.charge.at(j));
    }
    
    if (chid < nwire_u + nwire_v){
      // induction planes fold with low frequench stuff
      TH1 *hm = hresult->FFT(0,"MAG");
      TH1 *hp = hresult->FFT(0,"PH");
      TVirtualFFT *ifft2 = TVirtualFFT::FFT(1,&nticks,"C2R M K");

      for (int j=0;j!=nticks;j++){
	Double_t freq=0;
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
	if (j<dead_start || j > dead_end)
	  hresult->SetBinContent(j+1,fb->GetBinContent(j+1));
      }
      restore_baseline(hresult);
      
      for (int j=0;j!=int(nticks/rebin);j++){
	int sum = 0;
	for (int k=0;k!=rebin;k++){
	  int binno = j*rebin+k;
	  if (binno < dead_start || binno > dead_end)
	    sum += hresult->GetBinContent(binno+1);
	}
	hresult_filter->SetBinContent(j+1,sum);
      }

      delete fb;
      delete hm;
      delete hp;
      delete ifft2;
    }
    

    float th = cal_rms(hresult,chid) * rebin * factor;

    //    if (chid==746) std::cout << "a " << th << " " << max_th << std::endl;

    if (th > max_th) th = max_th;
    
    std::vector<std::pair <int,int> > ROIs_1;
    std::vector<Int_t> max_bins_1;
    int ntime = hresult_filter->GetNbinsX();
    
    for (Int_t j=1; j<ntime-1;j++){
      Double_t content = hresult_filter->GetBinContent(j+1);
      Double_t prev_content = hresult_filter->GetBinContent(j);
      Double_t next_content = hresult_filter->GetBinContent(j+2);
      Int_t flag_ROI = 0;
      Int_t begin=0;
      Int_t end=0;
      Int_t max_bin=0;
      if (content > th){
	begin = find_ROI_begin(hresult_filter,j, th*factor1) ;
	end = find_ROI_end(hresult_filter,j, th*factor1) ;
	max_bin = begin;
	for (Int_t k=begin;k<=end;k++){
	  if (hresult_filter->GetBinContent(k+1) > hresult_filter->GetBinContent(max_bin+1)){
	    max_bin = k;
	  }
	}
	flag_ROI = 1;
      }else{
	if (content > prev_content && content > next_content){
	  begin = find_ROI_begin(hresult_filter,j, prev_content);
	  end = find_ROI_end(hresult_filter,j, next_content );
	  max_bin = begin;
	  for (Int_t k=begin;k<=end;k++){
	    if (hresult_filter->GetBinContent(k+1) > hresult_filter->GetBinContent(max_bin+1)){
	      max_bin = k;
	    }
	  }
	  if (hresult_filter->GetBinContent(max_bin+1) - hresult_filter->GetBinContent(begin+1) +
	      hresult_filter->GetBinContent(max_bin+1) - hresult_filter->GetBinContent(end+1) > th*2){
	    flag_ROI = 1;
	  }
	  Int_t temp_begin = max_bin-short_length;
	  if (temp_begin < begin) temp_begin = begin;
	  Int_t temp_end = max_bin + short_length;
	  if (temp_end > end) temp_end = end;
	  if ((hresult_filter->GetBinContent(max_bin+1) - hresult_filter->GetBinContent(temp_begin+1) > th *factor1 &&
	       hresult_filter->GetBinContent(max_bin+1) - hresult_filter->GetBinContent(temp_end+1) > th*factor1) ){
	    flag_ROI = 1;
	  }
	  
	}
      }

       if (flag_ROI == 1){
      	if (ROIs_1.size()>0){
      	  if (begin <= ROIs_1.back().second){
      	    ROIs_1.back().second = end;
      	    if (hresult_filter->GetBinContent(max_bin+1) > hresult_filter->GetBinContent(max_bins_1.back()+1))
      	      max_bins_1.back() = max_bin;
      	  }else{
      	    ROIs_1.push_back(std::make_pair(begin,end));
      	    max_bins_1.push_back(max_bin);
      	  }
      	}else{
      	  ROIs_1.push_back(std::make_pair(begin,end));
      	  max_bins_1.push_back(max_bin);
      	}
	
      	if (end <hresult_filter->GetNbinsX()){
      	  j = end;
      	}else{
      	  j=hresult_filter->GetNbinsX();
      	}
      }
    }
    
    
    if (ROIs_1.size()==1){
    }else if (ROIs_1.size()>1){
      Int_t flag_repeat = 0;
      //  cout << "Xin1: " << ROIs_1.size() << endl;;
      while(flag_repeat){
    	flag_repeat = 1;
    	for (Int_t k=0;k<int(ROIs_1.size()-1);k++){
    	  Int_t begin = ROIs_1.at(k).first;
    	  Int_t end = ROIs_1.at(k+1).second;
	  
    	  Double_t begin_content = hresult_filter->GetBinContent(begin+1);
    	  Double_t end_content = hresult_filter->GetBinContent(end+1);
	  
    	  Int_t begin_1 = ROIs_1.at(k).second;
    	  Int_t end_1 = ROIs_1.at(k+1).first;	
	  
    	  Int_t flag_merge = 1;
    	  //Double_t sum1 = 0, sum2 = 0;
    	  for (Int_t j=begin_1; j<=end_1;j++){
    	    Double_t current_content = hresult_filter->GetBinContent(j+1);
    	    Double_t content = current_content - ((end_content - begin_content)*(j*1.-begin)/(end-begin*1.) + begin_content);
	    
    	    if (content < th*factor1){
    	      flag_merge = 0;
    	      break;
    	    }
    	    // sum1 += content;
    	    // sum2 ++;
    	    // cout << j << " " << content << endl;
    	  }
    	  // if (sum2 >0){
    	  //   if (sum1/sum2 < th*factor1) flag_merge = 0;
    	  // }
	  
    	  if (flag_merge == 1){
    	    ROIs_1.at(k).second = ROIs_1.at(k+1).second;
    	    ROIs_1.erase(ROIs_1.begin()+k+1);
    	    flag_repeat = 1;
    	    break;
    	  }
    	}
    	//	cout << "Xin2: " << ROIs_1.size() << endl;

      }
    }
    
    // scale back ... 
    for (int j = 0; j!=ROIs_1.size();j++){
      int begin = ROIs_1.at(j).first * rebin;
      int end = ROIs_1.at(j).second *rebin + (rebin-1);
      
      ROIs_1.at(j).first = begin;
      ROIs_1.at(j).second = end;

      //if (chid ==1195) std::cout << ROIs_1.at(j).first << " " << ROIs_1.at(j).second << std::endl;
    }


    if (chid < nwire_u){
      loose_rois_u.at(chid) = ROIs_1;
    }else if (chid < nwire_u + nwire_v){
      loose_rois_v.at(chid-nwire_u) = ROIs_1;
    }else{
      loose_rois_w.at(chid-nwire_u-nwire_v) = ROIs_1;
    }
  }

  
  delete hresult;
  delete hresult_filter;
  //delete filter_low;
  delete filter_low1;
  delete filter_low2;
}


void WireCell2dToy::uBooNEDataROI::extend_ROI_self(int pad){
  const int nbins = raw_fds.Get_Bins_Per_Frame();

  for (int i=0;i!=self_rois_u.size();i++){
    std::vector<std::pair<int,int>> temp_rois;
    int temp_begin=0, temp_end=0;
    for (int j=0;j!=self_rois_u.at(i).size();j++){
      temp_begin = self_rois_u.at(i).at(j).first - pad;
      if (temp_begin < 0 ) temp_begin = 0;
      temp_end = self_rois_u.at(i).at(j).second + pad;
      if (temp_end >= nbins) temp_end = nbins - 1;
      // merge 
      if (temp_rois.size() == 0){
	temp_rois.push_back(std::make_pair(temp_begin,temp_end));
      }else{
	if (temp_begin < temp_rois.back().second){
	  if (temp_end > temp_rois.back().second){
	    temp_rois.back().second = temp_end;
	  }
	}else{
	  temp_rois.push_back(std::make_pair(temp_begin,temp_end));
	}
      }
    }
    self_rois_u.at(i) = temp_rois;
  }

  for (int i=0;i!=self_rois_v.size();i++){
    std::vector<std::pair<int,int>> temp_rois;
    int temp_begin=0, temp_end=0;
    for (int j=0;j!=self_rois_v.at(i).size();j++){
      temp_begin = self_rois_v.at(i).at(j).first - pad;
      if (temp_begin < 0 ) temp_begin = 0;
      temp_end = self_rois_v.at(i).at(j).second + pad;
      if (temp_end >= nbins) temp_end = nbins - 1;
      // merge 
      if (temp_rois.size() == 0){
	temp_rois.push_back(std::make_pair(temp_begin,temp_end));
      }else{
	if (temp_begin < temp_rois.back().second){
	  if (temp_end > temp_rois.back().second){
	    temp_rois.back().second = temp_end;
	  }
	}else{
	  temp_rois.push_back(std::make_pair(temp_begin,temp_end));
	}
      }
    }
    self_rois_v.at(i) = temp_rois;
  }

  for (int i=0;i!=self_rois_w.size();i++){
    std::vector<std::pair<int,int>> temp_rois;
    int temp_begin=0, temp_end=0;
    for (int j=0;j!=self_rois_w.at(i).size();j++){
      temp_begin = self_rois_w.at(i).at(j).first - pad;
      if (temp_begin < 0 ) temp_begin = 0;
      temp_end = self_rois_w.at(i).at(j).second + pad;
      if (temp_end >= nbins) temp_end = nbins - 1;
      // merge 
      if (temp_rois.size() == 0){
	temp_rois.push_back(std::make_pair(temp_begin,temp_end));
      }else{
	if (temp_begin < temp_rois.back().second){
	  if (temp_end > temp_rois.back().second){
	    temp_rois.back().second = temp_end;
	  }
	}else{
	  temp_rois.push_back(std::make_pair(temp_begin,temp_end));
	}
      }
    }
    self_rois_w.at(i) = temp_rois;
  }


}


WireCell2dToy::uBooNEDataROI::~uBooNEDataROI()
{
  
}

// void WireCell2dToy::uBooNEDataROI::merge_ROIs(){

//   std::vector<std::vector<std::pair<int,int>>> others_rois_u1;
//   std::vector<std::vector<std::pair<int,int>>> others_rois_v1;
//   others_rois_u1.resize(nwire_u);
//   others_rois_v1.resize(nwire_v);
  
//   for (int i=0;i!=self_rois_u.size();i++){
//     others_rois_u1.at(i).insert(others_rois_u1.at(i).end(),self_rois_u.at(i).begin(),self_rois_u.at(i).end());
//     others_rois_u1.at(i).insert(others_rois_u1.at(i).end(),others_rois_u.at(i).begin(),others_rois_u.at(i).end());
//     std::sort(others_rois_u1.at(i).begin(),others_rois_u1.at(i).end());
//     for (int j=0;j!=others_rois_u1.at(i).size();j++){
//       if (combined_rois_u.at(i).size() == 0){
//    	combined_rois_u.at(i).push_back(others_rois_u1.at(i).at(j));
//       }else{
//    	if (others_rois_u1.at(i).at(j).first < combined_rois_u.at(i).back().second){
//    	  if (others_rois_u1.at(i).at(j).second > combined_rois_u.at(i).back().second)
//    	    combined_rois_u.at(i).back().second = others_rois_u1.at(i).at(j).second;
//    	}else{
//    	  combined_rois_u.at(i).push_back(others_rois_u1.at(i).at(j));
//    	}
//       }
//     }

//   }
  


//   for (int i=0;i!=self_rois_v.size();i++){
//     others_rois_v1.at(i).insert(others_rois_v1.at(i).end(),self_rois_v.at(i).begin(),self_rois_v.at(i).end());
//     others_rois_v1.at(i).insert(others_rois_v1.at(i).end(),others_rois_v.at(i).begin(),others_rois_v.at(i).end());
//     std::sort(others_rois_v1.at(i).begin(),others_rois_v1.at(i).end());
//     for (int j=0;j!=others_rois_v1.at(i).size();j++){
//       if (combined_rois_v.at(i).size() == 0){
//    	combined_rois_v.at(i).push_back(others_rois_v1.at(i).at(j));
//       }else{
//    	if (others_rois_v1.at(i).at(j).first < combined_rois_v.at(i).back().second){
//    	  if (others_rois_v1.at(i).at(j).second > combined_rois_v.at(i).back().second)
//    	    combined_rois_v.at(i).back().second = others_rois_v1.at(i).at(j).second;
//    	}else{
//    	  combined_rois_v.at(i).push_back(others_rois_v1.at(i).at(j));
//    	}
//       }
//     }
//   }
  

//   for (int i=0;i!=self_rois_w.size();i++){
//     combined_rois_w.at(i).insert(combined_rois_w.at(i).end(),self_rois_w.at(i).begin(),self_rois_w.at(i).end());
//   }
// }


// void WireCell2dToy::uBooNEDataROI::find_ROI_by_others(){
//   double u_pitch, v_pitch, w_pitch;
//   u_pitch = gds.pitch(kUwire);
//   v_pitch = gds.pitch(kVwire);
//   w_pitch = gds.pitch(kYwire);

//   int flag_dead = 0;
  
//   std::vector<float> udis,vdis,wdis;
//   for (int i=0; i!=self_rois_u.size();i++){
//     const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(0),i);
//     udis.push_back(gds.wire_dist(*wire));
//   }
//   for (int i=0; i!=self_rois_v.size();i++){
//     const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(1),i);
//     vdis.push_back(gds.wire_dist(*wire));
//   }
//   for (int i=0; i!=self_rois_w.size();i++){
//     const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(2),i);
//     wdis.push_back(gds.wire_dist(*wire));
//   }
  
//   float dis_u[3],dis_v[3],dis_w[3],
//     dis_puv[5],dis_puw[5],dis_pwv[5];

//   std::vector<std::vector<std::pair<int,int>>> others_rois_u1;
//   std::vector<std::vector<std::pair<int,int>>> others_rois_v1;
//   // std::vector<std::vector<std::pair<int,int>>> others_rois_w1;

//   others_rois_u1.resize(nwire_u);
//   others_rois_v1.resize(nwire_v);
//   //others_rois_w1.resize(nwire_w);


//   // // tiling the U and V
//   // for (int i = 0; i!= self_rois_u.size(); i++){
//   //   const GeomWire *wire1 = gds.by_planeindex(WirePlaneType_t(0),i);
//   //   for (int j = 0;j!=self_rois_v.size(); j++){
//   //     const GeomWire *wire2 = gds.by_planeindex(WirePlaneType_t(1),j);
//   //     GeomWireSelection wires;
//   //     int flag_wires = 0;
      
//   //     for (int i1 = 0; i1 != self_rois_u.at(i).size(); i1++){
//   // 	for (int j1 = 0; j1 != self_rois_v.at(j).size(); j1++){
	  
//   // 	  if (flag_wires == 1 && wires.size()==0) continue;

//   // 	  int start_1 = self_rois_u.at(i).at(i1).first;
//   // 	  int end_1 = self_rois_u.at(i).at(i1).second;
//   // 	  int start_2 = self_rois_v.at(j).at(j1).first;
//   // 	  int end_2 = self_rois_v.at(j).at(j1).second;
//   // 	  int start = start_1;
//   // 	  int end = end_1;
//   // 	  if (start_1 < start_2 ) start = start_2;
//   // 	  if (end_1 > end_2) end = end_2; 
	  
//   // 	  if ( end - start > 4){ // current cut ... 
//   // 	    // find wires ... 
//   // 	    if (flag_wires == 0){
//   // 	      flag_wires = 1;

//   // 	      dis_u[0] = udis.at(i) - u_pitch/2.;
//   // 	      dis_u[1] = dis_u[0] + u_pitch;
//   // 	      dis_u[2] = udis.at(i);
	      
//   // 	      dis_v[0] = vdis.at(j) - v_pitch/2.;
//   // 	      dis_v[1] = dis_v[0] + v_pitch;
//   // 	      dis_v[2] = vdis.at(j);

//   // 	      std::vector<Vector> puv(5);
	      
//   // 	      if(!gds.crossing_point(dis_u[2],dis_v[2],kUwire,kVwire, puv[4])) continue;
//   // 	      gds.crossing_point(dis_u[0],dis_v[0],kUwire,kVwire, puv[0]);
//   // 	      gds.crossing_point(dis_u[0],dis_v[1],kUwire,kVwire, puv[1]);
//   // 	      gds.crossing_point(dis_u[1],dis_v[1],kUwire,kVwire, puv[2]);
//   // 	      gds.crossing_point(dis_u[1],dis_v[0],kUwire,kVwire, puv[3]);
	      
//   // 	      puv[0] = 0.9*(puv[0]-puv[4])+puv[4];
//   // 	      puv[1] = 0.9*(puv[1]-puv[4])+puv[4];	    
//   // 	      puv[2] = 0.9*(puv[2]-puv[4])+puv[4];
//   // 	      puv[3] = 0.9*(puv[3]-puv[4])+puv[4];
	      
//   // 	      for (int a1=0;a1!=5;a1++){
//   // 		const GeomWire *n_wire = gds.closest(puv[a1],kYwire);
//   // 		if (n_wire == 0) continue;
//   // 		if (find(wires.begin(),wires.end(),n_wire) == wires.end())
//   // 		  wires.push_back(n_wire);
//   // 	      }
//   // 	    }
//   // 	    // find the index ... 
//   // 	    for (int k=0;k!=wires.size();k++){
//   // 	      others_rois_w1.at(wires.at(k)->index()).push_back(std::make_pair(start,end));
//   // 	    }
//   // 	  }
//   // 	}
//   //     }
//   //   }
//   // }
  
//   // for (int i=0;i!=others_rois_w1.size();i++){
//   //   if (others_rois_w1.at(i).size()==0) continue;
//   //   std::sort(others_rois_w1.at(i).begin(),others_rois_w1.at(i).end());
//   //   for (int j=0;j!=others_rois_w1.at(i).size();j++){
//   //     if (others_rois_w.at(i).size() == 0){
//   // 	others_rois_w.at(i).push_back(others_rois_w1.at(i).at(j));
//   //     }else{
//   // 	if (others_rois_w1.at(i).at(j).first <= others_rois_w.at(i).back().second){
//   // 	  if (others_rois_w1.at(i).at(j).second > others_rois_w.at(i).back().second)
//   // 	    others_rois_w.at(i).back().second = others_rois_w.at(i).back().second;
//   // 	}else{
//   // 	  others_rois_w.at(i).push_back(others_rois_w1.at(i).at(j));
//   // 	}
//   //     }
//   //   }
//   // }
  


//   // Tiling the U and W
  
//   for (int i = 0; i!= self_rois_u.size(); i++){
//     const GeomWire *wire1 = gds.by_planeindex(WirePlaneType_t(0),i);
//     for (int j = 0;j!=self_rois_w.size(); j++){
//       const GeomWire *wire2 = gds.by_planeindex(WirePlaneType_t(2),j);
//       GeomWireSelection wires;
//       int flag_wires = 0;
      
//       for (int i1 = 0; i1 != self_rois_u.at(i).size(); i1++){
// 	for (int j1 = 0; j1 != self_rois_w.at(j).size(); j1++){
	  
// 	  if (flag_wires == 1 && wires.size()==0) continue;
	  
// 	  int start_1 = self_rois_u.at(i).at(i1).first;
// 	  int end_1 = self_rois_u.at(i).at(i1).second;
// 	  int start_2 = self_rois_w.at(j).at(j1).first;
// 	  int end_2 = self_rois_w.at(j).at(j1).second;
// 	  int start = start_1;
// 	  int end = end_1;
	  
// 	  if (start_1 < start_2 ) start = start_2;
// 	  if (end_1 > end_2) end = end_2; 
	  
// 	  if ( end - start > 4){ // current cut ... 
// 	    // find wires ... 
// 	    if (flag_wires == 0){
// 	      flag_wires = 1;
	      
// 	      dis_u[0] = udis.at(i) - u_pitch/2.;
// 	      dis_u[1] = dis_u[0] + u_pitch;
// 	      dis_u[2] = udis.at(i);
	      
// 	      dis_w[0] = wdis.at(j) - w_pitch/2.;
// 	      dis_w[1] = dis_w[0] + w_pitch;
// 	      dis_w[2] = wdis.at(j);
	      
// 	      std::vector<Vector> puv(5);
	      
// 	      if(!gds.crossing_point(dis_u[2],dis_w[2],kUwire,kYwire, puv[4])) continue;
// 	      // gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, puv[0]);
// 	      // gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, puv[1]);
// 	      // gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, puv[2]);
// 	      // gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, puv[3]);
	      
// 	      // puv[0] = 0.9*(puv[0]-puv[4])+puv[4];
// 	      // puv[1] = 0.9*(puv[1]-puv[4])+puv[4];	    
// 	      // puv[2] = 0.9*(puv[2]-puv[4])+puv[4];
// 	      // puv[3] = 0.9*(puv[3]-puv[4])+puv[4];
	      
// 	      for (int a1=4;a1!=5;a1++){
// 		const GeomWire *n_wire = gds.closest(puv[a1],kVwire);
// 		if (n_wire == 0) continue;
// 		if (find(wires.begin(),wires.end(),n_wire) == wires.end())
// 		  wires.push_back(n_wire);
// 	      }
// 	    }
// 	    // find the index ... 
// 	    for (int k=0;k!=wires.size();k++){
// 	      others_rois_v1.at(wires.at(k)->index()).push_back(std::make_pair(start,end));
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
  
//   //  std::cout << others_rois_v1.at(3500-2400).size() << std::endl;


//   if (flag_dead == 1){
//     // Add U dead + W live
    
//     for (int i = 0; i!= self_rois_u.size(); i++){
//       if (umap.find(i) == umap.end()) continue;
//       const GeomWire *wire1 = gds.by_planeindex(WirePlaneType_t(0),i);
//       for (int j = 0;j!=self_rois_w.size(); j++){
// 	const GeomWire *wire2 = gds.by_planeindex(WirePlaneType_t(2),j);
// 	GeomWireSelection wires;
// 	int flag_wires = 0;
	
// 	{
// 	  for (int j1 = 0; j1 != self_rois_w.at(j).size(); j1++){
	    
// 	    if (flag_wires == 1 && wires.size()==0) continue;
	    
// 	    int start_1 = umap[i].first;
// 	    int end_1 = umap[i].second;
// 	    int start_2 = self_rois_w.at(j).at(j1).first;
// 	    int end_2 = self_rois_w.at(j).at(j1).second;
// 	    int start = start_1;
// 	    int end = end_1;
	    
// 	    if (start_1 < start_2 ) start = start_2;
// 	    if (end_1 > end_2) end = end_2; 
	    
// 	    if ( end - start > 4){ // current cut ... 
// 	      // find wires ... 
// 	      if (flag_wires == 0){
// 		flag_wires = 1;
		
// 		dis_u[0] = udis.at(i) - u_pitch/2.;
// 		dis_u[1] = dis_u[0] + u_pitch;
// 		dis_u[2] = udis.at(i);
		
// 		dis_w[0] = wdis.at(j) - w_pitch/2.;
// 		dis_w[1] = dis_w[0] + w_pitch;
// 		dis_w[2] = wdis.at(j);
		
// 		std::vector<Vector> puv(5);
		
// 		if(!gds.crossing_point(dis_u[2],dis_w[2],kUwire,kYwire, puv[4])) continue;
// 		gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, puv[0]);
// 		gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, puv[1]);
// 		gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, puv[2]);
// 		gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, puv[3]);
		
// 		puv[0] = 0.9*(puv[0]-puv[4])+puv[4];
// 		puv[1] = 0.9*(puv[1]-puv[4])+puv[4];	    
// 		puv[2] = 0.9*(puv[2]-puv[4])+puv[4];
// 		puv[3] = 0.9*(puv[3]-puv[4])+puv[4];
		
// 		for (int a1=0;a1!=5;a1++){
// 		  const GeomWire *n_wire = gds.closest(puv[a1],kVwire);
// 		  if (n_wire == 0) continue;
// 		  if (find(wires.begin(),wires.end(),n_wire) == wires.end())
// 		    wires.push_back(n_wire);
// 		}
// 	      }
// 	      // find the index ... 
// 	      for (int k=0;k!=wires.size();k++){
// 		others_rois_v1.at(wires.at(k)->index()).push_back(std::make_pair(start,end));
// 	      }
// 	    }
// 	  }
// 	}
//       }
//     }
    
    
//     // std::cout << others_rois_v1.at(3500-2400).size() << std::endl;
    
//     // Add W dead + U live
    
//     for (int i = 0; i!= self_rois_u.size(); i++){
//       const GeomWire *wire1 = gds.by_planeindex(WirePlaneType_t(0),i);
//       for (int j = 0;j!=self_rois_w.size(); j++){
// 	if (wmap.find(j) == wmap.end()) continue;
// 	const GeomWire *wire2 = gds.by_planeindex(WirePlaneType_t(2),j);
// 	GeomWireSelection wires;
// 	int flag_wires = 0;
	
// 	for (int i1 = 0; i1 != self_rois_u.at(i).size(); i1++){
// 	  {
	    
// 	    if (flag_wires == 1 && wires.size()==0) continue;
	    
// 	    int start_1 = self_rois_u.at(i).at(i1).first;
// 	    int end_1 = self_rois_u.at(i).at(i1).second;
// 	    int start_2 = wmap[j].first;
// 	    int end_2 = wmap[j].second;
// 	    int start = start_1;
// 	    int end = end_1;
	    
// 	    if (start_1 < start_2 ) start = start_2;
// 	    if (end_1 > end_2) end = end_2; 
	    
// 	    if ( end - start > 4){ // current cut ... 
// 	      // find wires ... 
// 	      if (flag_wires == 0){
// 		flag_wires = 1;
		
// 		dis_u[0] = udis.at(i) - u_pitch/2.;
// 		dis_u[1] = dis_u[0] + u_pitch;
// 		dis_u[2] = udis.at(i);
		
// 		dis_w[0] = wdis.at(j) - w_pitch/2.;
// 		dis_w[1] = dis_w[0] + w_pitch;
// 		dis_w[2] = wdis.at(j);
		
// 		std::vector<Vector> puv(5);
		
// 		if(!gds.crossing_point(dis_u[2],dis_w[2],kUwire,kYwire, puv[4])) continue;
// 		gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, puv[0]);
// 		gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, puv[1]);
// 		gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, puv[2]);
// 		gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, puv[3]);
		
// 		puv[0] = 0.9*(puv[0]-puv[4])+puv[4];
// 		puv[1] = 0.9*(puv[1]-puv[4])+puv[4];	    
// 		puv[2] = 0.9*(puv[2]-puv[4])+puv[4];
// 		puv[3] = 0.9*(puv[3]-puv[4])+puv[4];
		
// 		for (int a1=0;a1!=5;a1++){
// 		  const GeomWire *n_wire = gds.closest(puv[a1],kVwire);
// 		  if (n_wire == 0) continue;
// 		  if (find(wires.begin(),wires.end(),n_wire) == wires.end())
// 		    wires.push_back(n_wire);
// 		}
// 	      }
// 	      // find the index ... 
// 	      for (int k=0;k!=wires.size();k++){
// 		others_rois_v1.at(wires.at(k)->index()).push_back(std::make_pair(start,end));
// 	      }
// 	    }
// 	  }
// 	}
//       }
//     }
//   }

//    // std::cout << others_rois_v1.at(3500-2400).size() << std::endl;



//   for (int i=0;i!=others_rois_v1.size();i++){
//     if (others_rois_v1.at(i).size()==0) continue;
//     std::sort(others_rois_v1.at(i).begin(),others_rois_v1.at(i).end());
//     for (int j=0;j!=others_rois_v1.at(i).size();j++){
//       if (others_rois_v.at(i).size() == 0){
//       others_rois_v.at(i).push_back(others_rois_v1.at(i).at(j));
//       }else{
//      	if (others_rois_v1.at(i).at(j).first < others_rois_v.at(i).back().second){
//      	  if (others_rois_v1.at(i).at(j).second > others_rois_v.at(i).back().second)
//      	    others_rois_v.at(i).back().second = others_rois_v1.at(i).at(j).second;
//      	}else{
//      	  others_rois_v.at(i).push_back(others_rois_v1.at(i).at(j));
//      	}
//       }
//     }
//   }

//   // Tiling the V and W

//    for (int i = 0; i!= self_rois_v.size(); i++){
//     const GeomWire *wire1 = gds.by_planeindex(WirePlaneType_t(1),i);
//     for (int j = 0;j!=self_rois_w.size(); j++){
//       const GeomWire *wire2 = gds.by_planeindex(WirePlaneType_t(2),j);
//       GeomWireSelection wires;
//       int flag_wires = 0;
      
//       for (int i1 = 0; i1 != self_rois_v.at(i).size(); i1++){
// 	for (int j1 = 0; j1 != self_rois_w.at(j).size(); j1++){
	  
// 	  if (flag_wires == 1 && wires.size()==0) continue;

// 	  int start_1 = self_rois_v.at(i).at(i1).first;
// 	  int end_1 = self_rois_v.at(i).at(i1).second;
// 	  int start_2 = self_rois_w.at(j).at(j1).first;
// 	  int end_2 = self_rois_w.at(j).at(j1).second;
// 	  int start = start_1;
// 	  int end = end_1;
// 	  if (start_1 < start_2 ) start = start_2;
// 	  if (end_1 > end_2) end = end_2; 
	  
// 	  if ( end - start > 4){ // current cut ... 
// 	    // find wires ... 
// 	    if (flag_wires == 0){
// 	      flag_wires = 1;

// 	      dis_v[0] = vdis.at(i) - v_pitch/2.;
// 	      dis_v[1] = dis_v[0] + v_pitch;
// 	      dis_v[2] = vdis.at(i);
	      
// 	      dis_w[0] = wdis.at(j) - w_pitch/2.;
// 	      dis_w[1] = dis_w[0] + w_pitch;
// 	      dis_w[2] = wdis.at(j);

// 	      std::vector<Vector> puv(5);
	      
// 	      if(!gds.crossing_point(dis_v[2],dis_w[2],kVwire,kYwire, puv[4])) continue;
// 	      // gds.crossing_point(dis_v[0],dis_w[0],kVwire,kYwire, puv[0]);
// 	      // gds.crossing_point(dis_v[0],dis_w[1],kVwire,kYwire, puv[1]);
// 	      // gds.crossing_point(dis_v[1],dis_w[1],kVwire,kYwire, puv[2]);
// 	      // gds.crossing_point(dis_v[1],dis_w[0],kVwire,kYwire, puv[3]);
	      
// 	      // puv[0] = 0.9*(puv[0]-puv[4])+puv[4];
// 	      // puv[1] = 0.9*(puv[1]-puv[4])+puv[4];	    
// 	      // puv[2] = 0.9*(puv[2]-puv[4])+puv[4];
// 	      // puv[3] = 0.9*(puv[3]-puv[4])+puv[4];
	      
// 	      for (int a1=4;a1!=5;a1++){
// 		const GeomWire *n_wire = gds.closest(puv[a1],kUwire);
// 		if (n_wire == 0) continue;
// 		if (find(wires.begin(),wires.end(),n_wire) == wires.end())
// 		  wires.push_back(n_wire);
// 	      }
// 	    }
// 	    // find the index ... 
// 	    for (int k=0;k!=wires.size();k++){
// 	      others_rois_u1.at(wires.at(k)->index()).push_back(std::make_pair(start,end));
// 	    }
// 	  }
// 	}
//       }
//     }
//   }


//    if (flag_dead == 1){
//      // V dead
     
//      for (int i = 0; i!= self_rois_v.size(); i++){
//        if (vmap.find(i) == vmap.end()) continue;
//        const GeomWire *wire1 = gds.by_planeindex(WirePlaneType_t(1),i);
//        for (int j = 0;j!=self_rois_w.size(); j++){
// 	 const GeomWire *wire2 = gds.by_planeindex(WirePlaneType_t(2),j);
// 	 GeomWireSelection wires;
// 	 int flag_wires = 0;
	 
// 	 {
// 	   for (int j1 = 0; j1 != self_rois_w.at(j).size(); j1++){
	     
// 	     if (flag_wires == 1 && wires.size()==0) continue;
	     
// 	     int start_1 = vmap[i].first;
// 	     int end_1 = vmap[i].second;
// 	     int start_2 = self_rois_w.at(j).at(j1).first;
// 	     int end_2 = self_rois_w.at(j).at(j1).second;
// 	     int start = start_1;
// 	     int end = end_1;
// 	     if (start_1 < start_2 ) start = start_2;
// 	     if (end_1 > end_2) end = end_2; 
	     
// 	     if ( end - start > 4){ // current cut ... 
// 	       // find wires ... 
// 	       if (flag_wires == 0){
// 		 flag_wires = 1;
		 
// 		 dis_v[0] = vdis.at(i) - v_pitch/2.;
// 		 dis_v[1] = dis_v[0] + v_pitch;
// 		 dis_v[2] = vdis.at(i);
		 
// 		 dis_w[0] = wdis.at(j) - w_pitch/2.;
// 		 dis_w[1] = dis_w[0] + w_pitch;
// 		 dis_w[2] = wdis.at(j);
		 
// 		 std::vector<Vector> puv(5);
		 
// 		 if(!gds.crossing_point(dis_v[2],dis_w[2],kVwire,kYwire, puv[4])) continue;
// 		 gds.crossing_point(dis_v[0],dis_w[0],kVwire,kYwire, puv[0]);
// 		 gds.crossing_point(dis_v[0],dis_w[1],kVwire,kYwire, puv[1]);
// 		 gds.crossing_point(dis_v[1],dis_w[1],kVwire,kYwire, puv[2]);
// 		 gds.crossing_point(dis_v[1],dis_w[0],kVwire,kYwire, puv[3]);
		 
// 		 puv[0] = 0.9*(puv[0]-puv[4])+puv[4];
// 		 puv[1] = 0.9*(puv[1]-puv[4])+puv[4];	    
// 		 puv[2] = 0.9*(puv[2]-puv[4])+puv[4];
// 		 puv[3] = 0.9*(puv[3]-puv[4])+puv[4];
		 
// 		 for (int a1=0;a1!=5;a1++){
// 		   const GeomWire *n_wire = gds.closest(puv[a1],kUwire);
// 		   if (n_wire == 0) continue;
// 		   if (find(wires.begin(),wires.end(),n_wire) == wires.end())
// 		     wires.push_back(n_wire);
// 		 }
// 	       }
// 	       // find the index ... 
// 	       for (int k=0;k!=wires.size();k++){
// 		 others_rois_u1.at(wires.at(k)->index()).push_back(std::make_pair(start,end));
// 	       }
// 	     }
// 	   }
// 	 }
//        }
//      }
     
//      // W dead
     
//      for (int i = 0; i!= self_rois_v.size(); i++){
//        const GeomWire *wire1 = gds.by_planeindex(WirePlaneType_t(1),i);
//        for (int j = 0;j!=self_rois_w.size(); j++){
// 	 if (wmap.find(j) == wmap.end()) continue;
// 	 const GeomWire *wire2 = gds.by_planeindex(WirePlaneType_t(2),j);
// 	 GeomWireSelection wires;
// 	 int flag_wires = 0;
	 
// 	 for (int i1 = 0; i1 != self_rois_v.at(i).size(); i1++){
// 	   {
	     
// 	     if (flag_wires == 1 && wires.size()==0) continue;
	     
// 	     int start_1 = self_rois_v.at(i).at(i1).first;
// 	     int end_1 = self_rois_v.at(i).at(i1).second;
// 	     int start_2 = wmap[j].first;
// 	     int end_2 = wmap[j].second;
// 	     int start = start_1;
// 	     int end = end_1;
// 	     if (start_1 < start_2 ) start = start_2;
// 	     if (end_1 > end_2) end = end_2; 
	     
// 	     if ( end - start > 4){ // current cut ... 
// 	       // find wires ... 
// 	       if (flag_wires == 0){
// 		 flag_wires = 1;
		 
// 		 dis_v[0] = vdis.at(i) - v_pitch/2.;
// 		 dis_v[1] = dis_v[0] + v_pitch;
// 		 dis_v[2] = vdis.at(i);
		 
// 		 dis_w[0] = wdis.at(j) - w_pitch/2.;
// 		 dis_w[1] = dis_w[0] + w_pitch;
// 		 dis_w[2] = wdis.at(j);
		 
// 		 std::vector<Vector> puv(5);
		 
// 		 if(!gds.crossing_point(dis_v[2],dis_w[2],kVwire,kYwire, puv[4])) continue;
// 		 gds.crossing_point(dis_v[0],dis_w[0],kVwire,kYwire, puv[0]);
// 		 gds.crossing_point(dis_v[0],dis_w[1],kVwire,kYwire, puv[1]);
// 		 gds.crossing_point(dis_v[1],dis_w[1],kVwire,kYwire, puv[2]);
// 		 gds.crossing_point(dis_v[1],dis_w[0],kVwire,kYwire, puv[3]);
		 
// 		 puv[0] = 0.9*(puv[0]-puv[4])+puv[4];
// 		 puv[1] = 0.9*(puv[1]-puv[4])+puv[4];	    
// 		 puv[2] = 0.9*(puv[2]-puv[4])+puv[4];
// 		 puv[3] = 0.9*(puv[3]-puv[4])+puv[4];
		 
// 		 for (int a1=0;a1!=5;a1++){
// 		   const GeomWire *n_wire = gds.closest(puv[a1],kUwire);
// 		   if (n_wire == 0) continue;
// 		   if (find(wires.begin(),wires.end(),n_wire) == wires.end())
// 		     wires.push_back(n_wire);
// 		 }
// 	       }
// 	       // find the index ... 
// 	       for (int k=0;k!=wires.size();k++){
// 		 others_rois_u1.at(wires.at(k)->index()).push_back(std::make_pair(start,end));
// 	       }
// 	     }
// 	   }
// 	 }
//        }
//      }
//    }



//   for (int i=0;i!=others_rois_u1.size();i++){
//     if (others_rois_u1.at(i).size()==0) continue;
//     std::sort(others_rois_u1.at(i).begin(),others_rois_u1.at(i).end());
//     for (int j=0;j!=others_rois_u1.at(i).size();j++){
//       if (others_rois_u.at(i).size() == 0){
//       others_rois_u.at(i).push_back(others_rois_u1.at(i).at(j));
//       }else{
//      	if (others_rois_u1.at(i).at(j).first < others_rois_u.at(i).back().second){
//      	  if (others_rois_u1.at(i).at(j).second > others_rois_u.at(i).back().second)
//      	    others_rois_u.at(i).back().second = others_rois_u1.at(i).at(j).second;
//      	}else{
//      	  others_rois_u.at(i).push_back(others_rois_u1.at(i).at(j));
//      	}
//       }
//     }
//   }

  
// }


// void WireCell2dToy::uBooNEDataROI::find_ROI_by_raw_itself(int th_factor , int pad ){
//   const int nbins = raw_fds.Get_Bins_Per_Frame();
//   TH1F *hresult = new TH1F("hresult","hresult",nbins,0,nbins);
    
//   const Frame& frame1 = raw_fds.get();
//   size_t ntraces = frame1.traces.size();
  

//   for (int i=0;i!=ntraces;i++){
//     const Trace& trace = frame1.traces[i];
//     int tbin = trace.tbin;
//     int chid = trace.chid;
//     int nticks = trace.charge.size();
//     hresult->Reset();
        
//     int dead_start = -1;
//     int dead_end = -1;



//     if (chid >= nwire_u + nwire_v) continue; // no need for the collection plane ... 

//     if (chid < nwire_u){
//       if (umap.find(chid) != umap.end()){
	
// 	dead_start = umap[chid].first;
// 	dead_end = umap[chid].second;
//       }
//     }else if (chid < nwire_u + nwire_v){
//       if (vmap.find(chid-nwire_u) != vmap.end()){
	
// 	dead_start = vmap[chid-nwire_u].first;
// 	dead_end = vmap[chid-nwire_v].second;
//       }
//     }else{
//       if (wmap.find(chid-nwire_u-nwire_v) != wmap.end()){
// 	dead_start = wmap[chid-nwire_u-nwire_v].first;
// 	dead_end = wmap[chid-nwire_u-nwire_v].second;
//       }
//     }
    
//     //if (chid == 7500) std::cout << 7500 << " " << dead_start << " " << dead_end << std::endl;

//     for (int j=0;j!=nticks;j++){
//       if (j < dead_start || j > dead_end)
// 	hresult->SetBinContent(j+1,trace.charge.at(j));
//     }
    
//     //std::cout << i << " " << ntraces << " " << nticks << std::endl;

//     float th = cal_rms(hresult,chid);
//     float threshold = th_factor * th + 1e-6;
    
//     int roi_begin = -1;
//     int roi_end = -1;
//     std::vector<std::pair<int,int>> temp_rois;

//     // search the things above threshold (positive) add pad after it 
//     // search the things below -threshold (negative) add pad before it
//     for (int j=0;j<hresult->GetNbinsX()-1;j++){
//       double content = hresult->GetBinContent(j+1);
//       if (content > threshold && chid >= nwire_u){ // not u plane, look at the positive
// 	roi_begin = j;
// 	roi_end = j;
// 	for (int k=j+1;k<hresult->GetNbinsX();k++){
// 	  if (hresult->GetBinContent(k+1) > threshold){
// 	    roi_end = k;
// 	  }else{
// 	    break;
// 	  }
// 	}
// 	if (roi_end - roi_begin >=2){
// 	  int temp_roi_end = roi_end + pad;
// 	  if (temp_roi_end >= hresult->GetNbinsX())
// 	    temp_roi_end = hresult->GetNbinsX()-1;
// 	  temp_rois.push_back(std::make_pair(roi_begin,temp_roi_end));
// 	}
// 	j = roi_end + 1;
//       }

//       if (content < -threshold){
// 	roi_begin = j;
// 	roi_end = j;
// 	for (int k=j+1;k<hresult->GetNbinsX();k++){
// 	  if (hresult->GetBinContent(k+1) < -threshold){
// 	    roi_end = k;
// 	  }else{
// 	    break;
// 	  }
// 	}
// 	if (roi_end - roi_begin >=2){
// 	  int temp_roi_begin;
// 	  if (chid <nwire_u){
// 	    temp_roi_begin = roi_begin - pad;
// 	  }else{
// 	    temp_roi_begin = roi_begin - pad * 2;
// 	  }

// 	  if (temp_roi_begin < 0)
// 	    temp_roi_begin = 0;
// 	  temp_rois.push_back(std::make_pair(temp_roi_begin,roi_end));
// 	}
// 	j = roi_end + 1;
//       }
      
//     }

//     // load the ROIs from the self in
//     if (chid < nwire_u){
//       temp_rois.insert(temp_rois.end(),self_rois_u.at(chid).begin(), self_rois_u.at(chid).end());
//     }else if (chid < nwire_u + nwire_v){
//       temp_rois.insert(temp_rois.end(),self_rois_v.at(chid - nwire_u).begin(), self_rois_v.at(chid - nwire_u).end());
//     }
//     // sort them ...
//     std::sort(temp_rois.begin(),temp_rois.end());
//     // merge them ...
//     std::vector<std::pair<int,int>> temp_rois1 ; // to store 
    
//     for (int j=0;j!=temp_rois.size();j++){
//       if (temp_rois1.size() == 0){
//    	temp_rois1.push_back(temp_rois.at(j));
//       }else{
//    	if (temp_rois.at(j).first < temp_rois1.back().second){
//    	  if (temp_rois.at(j).second > temp_rois1.back().second)
//    	    temp_rois1.back().second = temp_rois.at(j).second;
//    	}else{
//    	  temp_rois1.push_back(temp_rois.at(j));
//    	}
//       }
//     }
//     // copy them back
//     if (chid < nwire_u){
//       self_rois_u.at(chid) = temp_rois1;
//     }else if (chid < nwire_u + nwire_v){
//       self_rois_v.at(chid - nwire_u) = temp_rois1;
//     }

//   }

//   delete hresult;
  
// }


void WireCell2dToy::uBooNEDataROI::find_ROI_by_decon_itself(int th_factor_ind, int th_factor_col , int pad ){
  
  const int nbins = fds.Get_Bins_Per_Frame();

  // make a histogram for induction planes and fold it with a low-frequency stuff
  // if collection, just fill the histogram ... 
  TH1F *hresult = new TH1F("hresult","hresult",nbins,0,nbins);
  TH1F *hresult_filter = new TH1F("hresult_filter","hresult_filter",nbins,0,nbins);
  TH1F *hresult_filter_tight = new TH1F("hresult_filter_tight","hresult_filter_tight",nbins,0,nbins);
  TH1F *hresult_roi = new TH1F("hresult_roi","hresult_roi",nbins,0,nbins);

  
  TF1 *filter_low = new TF1("filter_low","1-exp(-pow(x/0.02,2))");
  TF1 *filter_low_tight = new TF1("filter_low","1-exp(-pow(x/0.1,2))");
  int filter_pad = 60; //ticks

  // load the data and do the convolution ... 
  const Frame& frame1 = fds.get();
  size_t ntraces = frame1.traces.size();
  double value_re[nbins];
  double value_im[nbins];
  for (int i=0;i!=nbins;i++){
    value_re[i] = 0;
    value_im[i] = 0;
  }
  
  for (int i=0;i!=ntraces;i++){
    const Trace& trace = frame1.traces[i];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nticks = trace.charge.size();
    hresult->Reset();
    hresult_filter->Reset();
    hresult_filter_tight->Reset();
    hresult_roi->Reset();
    
    int dead_start = -1;
    int dead_end = -1;
    
    if (chid < nwire_u){
      if (umap.find(chid) != umap.end()){
	dead_start = umap[chid].first;
	dead_end = umap[chid].second;
      }
    }else if (chid < nwire_u + nwire_v){
      if (vmap.find(chid-nwire_u) != vmap.end()){
	
	dead_start = vmap[chid-nwire_u].first;
	dead_end = vmap[chid-nwire_v].second;
      }
    }else{
      if (wmap.find(chid-nwire_u-nwire_v) != wmap.end()){
	dead_start = wmap[chid-nwire_u-nwire_v].first;
	dead_end = wmap[chid-nwire_u-nwire_v].second;
      }
    }
    
    // if (chid == 7500) std::cout << 7500 << " " << dead_start << " " << dead_end << std::endl;

    for (int j=0;j!=nticks;j++){
      if (j < dead_start || j > dead_end)
	hresult->SetBinContent(j+1,trace.charge.at(j));
    }
    
    if (chid < nwire_u + nwire_v){
      // induction planes fold with low frequench stuff
      TH1 *hm = hresult->FFT(0,"MAG");
      TH1 *hp = hresult->FFT(0,"PH");
      TVirtualFFT *ifft2 = TVirtualFFT::FFT(1,&nticks,"C2R M K");

      for (int j=0;j!=nticks;j++){
	Double_t freq=0;
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
	if (j < dead_start || j > dead_end)
	  hresult_filter->SetBinContent(j+1,fb->GetBinContent(j+1));
      }

      delete fb;


      for (int j=0;j!=nticks;j++){
	Double_t freq=0;
	if (j < nticks/2.){
	  freq = j/(1.*nticks)*2.;
	}else{
	  freq = (nticks - j)/(1.*nticks)*2.;
	}
	
	value_re[j] = hm->GetBinContent(j+1) * cos(hp->GetBinContent(j+1)) / nticks * filter_low_tight->Eval(freq);
	value_im[j] = hm->GetBinContent(j+1) * sin(hp->GetBinContent(j+1)) / nticks * filter_low_tight->Eval(freq);
      }
      ifft2->SetPointsComplex(value_re,value_im);
      ifft2->Transform();
      fb = 0;
      fb = TH1::TransformHisto(ifft2,fb,"Re");
      
      for (int j=0;j!=nticks;j++){
	if (j < dead_start || j > dead_end)
	  hresult_filter_tight->SetBinContent(j+1,fb->GetBinContent(j+1));
      }
      delete fb;

      


      delete hm;
      delete hp;
      delete ifft2;
    }else{
      for (int j=0;j!=nticks;j++){
	if (j < dead_start || j > dead_end)
	  hresult_filter->SetBinContent(j+1,hresult->GetBinContent(j+1));
      }
    }
    restore_baseline(hresult_filter);
    restore_baseline(hresult_filter_tight);
    //std::cout << chid << " " << cal_rms(hresult,chid) << std::endl;
    float th = cal_rms(hresult_filter,chid);
    float th_tight = cal_rms(hresult_filter_tight,chid);
    float threshold=0;// = th_factor * th + 1;
    float threshold_tight = 0;
    //int pad = 5;

    if (chid < nwire_u){
      uplane_rms.at(chid) = th;
      threshold = th_factor_ind * th + 1;
      threshold_tight = th_factor_ind * th_tight + 1;
    }else if (chid < nwire_u + nwire_v){
      vplane_rms.at(chid-nwire_u) = th;
      threshold = th_factor_ind * th + 1;
      threshold_tight = th_factor_ind * th_tight + 1;
    }else{
      wplane_rms.at(chid-nwire_u-nwire_v) = th;
      threshold = th_factor_col * th + 1;
      threshold_tight = th_factor_col * th_tight + 1;
    }

    
    int roi_begin=-1;
    int roi_end=-1;
    
    std::vector<std::pair<int,int>> temp_rois;
    // now find ROI, above five sigma, and pad with +- six time ticks
    for (int j=0;j<hresult_filter->GetNbinsX()-1;j++){
      double content = hresult_filter->GetBinContent(j+1);
      double content_tight = hresult_filter_tight->GetBinContent(j+1);

      if (content > threshold || 
	  (content_tight > threshold )){
	roi_begin = j;
	roi_end = j;
	for (int k=j+1;k<hresult_filter->GetNbinsX();k++){
	  if (hresult_filter->GetBinContent(k+1) > threshold ||
	      (hresult_filter_tight->GetBinContent(k+1) > threshold)){
	    roi_end = k;
	  }else{
	    break;
	  }
	}
	int temp_roi_begin = roi_begin - pad ; // filter_pad;
	if (temp_roi_begin <0 ) temp_roi_begin = 0;
	int temp_roi_end = roi_end + pad; // filter_pad;
	if (temp_roi_end >hresult_filter->GetNbinsX()-1) temp_roi_end = hresult_filter->GetNbinsX()-1;

	//if (chid == 1151) std::cout << temp_roi_begin << " " << temp_roi_end << std::endl;

	
	if (temp_rois.size() == 0){
	  temp_rois.push_back(std::make_pair(temp_roi_begin,temp_roi_end));
	}else{
	  if (temp_roi_begin <= temp_rois.back().second){
	    temp_rois.back().second = temp_roi_end;
	  }else{
	    temp_rois.push_back(std::make_pair(temp_roi_begin,temp_roi_end));
	  }
	}
	j = roi_end + 1;
      }
    }

    // // fill the new roi histogram
    // for (int j=0;j!=temp_rois.size();j++){
    //   int start  = temp_rois.at(j).first;
    //   int end = temp_rois.at(j).second;
    //   float content_start = hresult->GetBinContent(start+1);
    //   float content_end = hresult->GetBinContent(end+1);
      
    //   for (int k=start;k<=end;k++){
    // 	float content = (hresult->GetBinContent(k+1) - content_start)*(k-start)*1.0/(end-start) + content_start;
    // 	hresult_roi->SetBinContent(k+1,hresult->GetBinContent(k+1) - content);
    //   }
    // }

    
    // roi_begin=-1;
    // roi_end=-1;
    // temp_rois.clear();
    
    // // now find ROI, above five sigma, and pad with +- six time ticks
    // for (int j=0;j<hresult_roi->GetNbinsX()-1;j++){
    //   double content = hresult_roi->GetBinContent(j+1);
    //   if (content > threshold){
    // 	roi_begin = j;
    // 	roi_end = j;
    // 	for (int k=j+1;k<hresult_roi->GetNbinsX();k++){
    // 	  if (hresult_roi->GetBinContent(k+1) > threshold){
    // 	    roi_end = k;
    // 	  }else{
    // 	    break;
    // 	  }
    // 	}
    // 	int temp_roi_begin = roi_begin - pad;
    // 	if (temp_roi_begin <0 ) temp_roi_begin = 0;
    // 	int temp_roi_end = roi_end + pad;
    // 	if (temp_roi_end >hresult_roi->GetNbinsX()-1) temp_roi_end = hresult_roi->GetNbinsX()-1;
	
    // 	if (temp_rois.size() == 0){
    // 	  temp_rois.push_back(std::make_pair(temp_roi_begin,temp_roi_end));
    // 	}else{
    // 	  if (temp_roi_begin <= temp_rois.back().second){
    // 	    temp_rois.back().second = temp_roi_end;
    // 	  }else{
    // 	    temp_rois.push_back(std::make_pair(temp_roi_begin,temp_roi_end));
    // 	  }
    // 	}
    // 	j = roi_end + 1;
    //   }
    // }

    


    if (chid < nwire_u){
      self_rois_u.at(chid) = temp_rois;
    }else if (chid < nwire_u + nwire_v){
      self_rois_v.at(chid-nwire_u) = temp_rois;
    }else{
      self_rois_w.at(chid-nwire_u-nwire_v) = temp_rois;
    }
    //    std::cout << chid << " " << temp_rois.size() << std::endl;
  }

  
  delete hresult;
  delete hresult_filter;
  delete hresult_filter_tight;
  delete hresult_roi;
  delete filter_low;
  delete filter_low_tight;
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
    
    //if (chid == 746) std::cout << min1 << " " << max1 << " " << start << " " << end << std::endl;

    TH1F *h6 = new TH1F("h6","h6",int(max1-min1+1),min1,max1+1);
    for (int i=0;i!=htemp->GetNbinsX();i++){
      if (i < start || i > end){
    	h6->Fill(int(htemp->GetBinContent(i+1)));
      }
    }
    if (h6->GetSum()>0){
      //calculate 0.16, 0.84 percentile ...  
      double xq=0;
      xq = 0.16;
      double par[2];
      h6->GetQuantiles(1,&par[0],&xq);
      xq = 0.84;
      h6->GetQuantiles(1,&par[1],&xq);
      rms = (par[1]-par[0])/2.;
      
      //      if (chid == 746) std::cout << rms << " " << (par[0]+par[1])/2.<< std::endl;
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
    
    //    if (chid == 746) std::cout << rms1 << std::endl;
    return rms1;
}
