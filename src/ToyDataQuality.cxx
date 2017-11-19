#include "WireCell2dToy/ToyDataQuality.h"

bool WireCell2dToy::Noisy_Event_ID(TH2F *hu_decon, TH2F *hv_decon, TH2F *hw_decon, TH1F *hu_threshold, TH1F *hv_threshold, TH1F *hw_threshold, WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map, TH2F *hu_decon_g, TH2F *hv_decon_g, TH2F *hw_decon_g){

  int nwire_u = hu_decon->GetNbinsX();
  int nwire_v = hu_decon->GetNbinsX();
  int nwire_w = hu_decon->GetNbinsX();
  
  // filter noisy events ... 
  int n_time_slice = hu_decon->GetNbinsY();
  int length_cut = 10;
  int n_rebin = 1;

  std::vector<std::tuple<int,int,float>> summary_u;
  std::vector<std::tuple<int,int,float>> summary_v;
  std::vector<std::tuple<int,int,float>> summary_w;

  
  // ID the events
  for (int i=0;i<n_time_slice/n_rebin;i++){
    bool flag_roi = false;
    int start;
    int end;
    std::vector<std::pair<int,int>> ROI_u;
    std::vector<std::pair<int,int>> ROI_v;
    std::vector<std::pair<int,int>> ROI_w;
    // U plane
    for (int j = 0; j!= nwire_u; j++){
      float content = 0;
      for (int k=0;k!=n_rebin;k++){
	content += hu_decon->GetBinContent(j+1,i*n_rebin+1+k);
      }
      float rms = hu_threshold->GetBinContent(j+1);
      if (!flag_roi){
	if (content > rms){
	  start = j;
	  flag_roi = true;
	}
      }else{
	if (content <=rms || j==nwire_u-1){
	  end = j;
	  flag_roi = false;
	  if (end-start>length_cut)
	    ROI_u.push_back(std::make_pair(start,end));
	  //	    std::cout << i << " " << start << " " << end << std::endl;
	}
      }
    }
    // V plane
    for (int j = 0; j!= nwire_v; j++){
      float content =  0;
      for (int k=0;k!=n_rebin;k++){
	content += hv_decon->GetBinContent(j+1,i*n_rebin+1+k);
      }
      float rms = hv_threshold->GetBinContent(j+1);
      if (!flag_roi){
	if (content > rms){
	  start = j;
	  flag_roi = true;
	}
      }else{
	if (content <=rms || j==nwire_v-1){
	  end = j;
	  flag_roi = false;
	  if (end-start>length_cut)
	    ROI_v.push_back(std::make_pair(start,end));
	  //	    std::cout << i << " " << start << " " << end << std::endl;
	}
      }
    }
    // W plane
    for (int j = 0; j!= nwire_w; j++){
      float content  = 0;
      for (int k=0;k!=n_rebin;k++){
	content += hw_decon->GetBinContent(j+1,i*n_rebin+1+k);
      }
      float rms = hw_threshold->GetBinContent(j+1);
      if (!flag_roi){
	if (content > rms){
	  start = j;
	  flag_roi = true;
	}
      }else{
	if (content <=rms || j==nwire_w-1){
	  end = j;
	  flag_roi = false;
	  if (end-start>length_cut)
	    ROI_w.push_back(std::make_pair(start,end));
	  //	    std::cout << i << " " << start << " " << end << std::endl;
	}
      }
    }
    float nfired_u = 0;
    float nfired_v = 0;
    float nfired_w = 0;
    int ncover_u = 0;
    int ncover_v = 0;
    int ncover_w = 0;

    if (ROI_u.size()>0){
      for (int j=0;j!=ROI_u.size();j++){
	nfired_u += ROI_u.at(j).second - ROI_u.at(j).first;
      }
      ncover_u = ROI_u.back().second - ROI_u.front().first;
      nfired_u/=ncover_u;
    }
    
    if (ROI_v.size()>0){
      for (int j=0;j!=ROI_v.size();j++){
	nfired_v += ROI_v.at(j).second - ROI_v.at(j).first;
      }
      ncover_v = ROI_v.back().second - ROI_v.front().first;
      nfired_v/=ncover_v;
    }
    if (ROI_w.size()>0){
      for (int j=0;j!=ROI_w.size();j++){
	nfired_w += ROI_w.at(j).second - ROI_w.at(j).first;
      }
      ncover_w = ROI_w.back().second - ROI_w.front().first;
      nfired_w/=ncover_w;
    }
    
    if (ncover_u * nfired_u > 100){
      summary_u.push_back(std::make_tuple(i,ncover_u,nfired_u));
    }
    if (ncover_v * nfired_v > 100){
      summary_v.push_back(std::make_tuple(i,ncover_v,nfired_v));
    }
    if (ncover_w * nfired_w > 100){
      summary_w.push_back(std::make_tuple(i,ncover_w,nfired_w));
    }
    
    //       if (ROI_u.size()>0 && ROI_v.size() > 0 && ROI_w.size() > 0)
    //std::cout << i*n_rebin * nrebin << " " << nfired_u * ncover_u << " " << ncover_u
    // 		<< " " << nfired_v * ncover_v << " " << ncover_v
    //		<< " " << nfired_w * ncover_w << " " << ncover_w << std::endl;
  }

  bool flag_u = false;
  int prev_time = -1;
  int n_cover = 0;
  int n_fire = 0;
  int start_time =0;
  int end_time = 0;
  
  for (size_t i=0;i!=summary_u.size();i++){
    int time = std::get<0>(summary_u.at(i));
    int ncover = std::get<1>(summary_u.at(i));
    float percentage = std::get<2>(summary_u.at(i));
    
    if (time > prev_time+3){
      if (n_cover >=6 && n_fire>12 || n_cover >6 && n_fire > 8){
	end_time = prev_time;
	std::cout << "U: " << n_cover << " " << n_fire << " " << start_time << " " << end_time << std::endl;
	flag_u = true;
      }
      
      n_cover = 0;
      n_fire = 0;
      start_time = time;
      end_time = time;
    }
    // std::cout << time << " " << ncover << " " << ncover*percentage << std::endl;
    
    if (ncover>1200. && (percentage > 0.25 || ncover *percentage > 200))
      n_cover ++;
    if (ncover*percentage>150)
      n_fire ++;
    prev_time = time;
  }
  
  if (n_cover >=6 && n_fire>12|| n_cover >6 && n_fire > 8){
    end_time = prev_time;
    flag_u = true;
    std::cout << "U: " << n_cover << " " << n_fire << " " << start_time << " " << end_time << std::endl;
  }

  bool flag_v = false;
  prev_time = -1;
  n_cover = 0;
  n_fire = 0;
  start_time =0;
  end_time = 0;
  
  for (size_t i=0;i!=summary_v.size();i++){
    int time = std::get<0>(summary_v.at(i));
    int ncover = std::get<1>(summary_v.at(i));
    float percentage = std::get<2>(summary_v.at(i));
    
    if (time > prev_time+3){
      //std::cout << n_cover << " " << n_fire << std::endl;
      if (n_cover >=5 && n_fire>=7){
	end_time = prev_time;
	flag_v = true;
	std::cout << "V: " << n_cover << " " << n_fire << " " << start_time << " " << prev_time << std::endl;
      }
      
      n_cover = 0;
      n_fire = 0;
      start_time = time;
      end_time = time;
    }
    //    std::cout << time << " " << ncover << " " << ncover*percentage << std::endl;
    
    if (ncover>1200. && (percentage > 0.25 || ncover *percentage > 200))
      n_cover ++;
    if (ncover*percentage>150)
      n_fire ++;
    prev_time = time;
  }
  if (n_cover >=5 && n_fire>=7){
    end_time = prev_time;
    flag_v = true;
    std::cout << "V: " << n_cover << " " << n_fire << " " << start_time << " " << end_time << std::endl;
  }

  bool flag_w = false;
  prev_time = -1;
  n_cover = 0;
  n_fire = 0;
  start_time = 0;
  end_time = 0;
  
  for (size_t i=0;i!=summary_w.size();i++){
    int time = std::get<0>(summary_w.at(i));
    int ncover = std::get<1>(summary_w.at(i));
    float percentage = std::get<2>(summary_w.at(i));
    
    if (time > prev_time+3){
      //std::cout << n_cover << " " << n_fire << std::endl;
      if (n_cover >=2 && n_fire>=3){
	end_time = prev_time;
	flag_w = true;
	std::cout << "W: " << n_cover << " " << n_fire << " " << start_time << " " << end_time << std::endl;
      }
      
      n_cover = 0;
      n_fire = 0;
      start_time = time;
      end_time = time;
    }
    //   std::cout << time << " " << ncover << " " << ncover*percentage << std::endl;
    
    if (ncover>1800. && (percentage > 0.25 || ncover *percentage > 200))
      n_cover ++;
    if (ncover*percentage>150)
      n_fire ++;
    prev_time = time;
  }
  if (n_cover >=2 && n_fire>=3){
    end_time = prev_time;
    flag_w = true;
    std::cout << "W: " << n_cover << " " << n_fire << " " << start_time << " " << end_time << std::endl;
  }
  
  std::cout << "Xin: " << " " << flag_u << " " << flag_v << " " << flag_w << std::endl;
  
  //  std::cout << n_cover << " " << n_fire << std::endl;
  // ID the region ...
}
