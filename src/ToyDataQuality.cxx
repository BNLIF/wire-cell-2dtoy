#include "WireCell2dToy/ToyDataQuality.h"

void WireCell2dToy::Organize_Dead_Channels(WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map, int nbins, int n_div){
  std::vector<int> boundaries;
  for (int i=0;i!=n_div+1;i++){
    boundaries.push_back(int(i*nbins/n_div));
    //std::cout << boundaries.back() << std::endl;
  }
  
  for (auto it = uplane_map.begin(); it!=uplane_map.end(); it++){
    int start = it->second.first;
    int end = it->second.second;
    for (int i=0;i!=n_div;i++){
      if (start >= boundaries.at(i) && start < boundaries.at(i+1))
	it->second.first = boundaries.at(i);
      if (end > boundaries.at(i) && end <= boundaries.at(i+1))
	it->second.second = boundaries.at(i+1);
    }
  }
  for (auto it = vplane_map.begin(); it!=vplane_map.end(); it++){
    int start = it->second.first;
    int end = it->second.second;
    for (int i=0;i!=n_div;i++){
      if (start >= boundaries.at(i) && start < boundaries.at(i+1))
	it->second.first = boundaries.at(i);
      if (end > boundaries.at(i) && end <= boundaries.at(i+1))
	it->second.second = boundaries.at(i+1);
    }
  }
  for (auto it = wplane_map.begin(); it!=wplane_map.end(); it++){
    int start = it->second.first;
    int end = it->second.second;
    for (int i=0;i!=n_div;i++){
      if (start >= boundaries.at(i) && start < boundaries.at(i+1))
	it->second.first = boundaries.at(i);
      if (end > boundaries.at(i) && end <= boundaries.at(i+1))
	it->second.second = boundaries.at(i+1);
    }
  }

  
  
  // for (auto it = vplane_map.begin(); it!=vplane_map.end(); it++){
  //   int start = it->second.first;
  //   int end = it->second.second;
  //   std::cout << start << " " << end << std::endl;
  // }
  
  
  
}

int WireCell2dToy::Noisy_Event_ID(TH2F *hu_decon, TH2F *hv_decon, TH2F *hw_decon, std::vector<float>& uplane_rms, std::vector<float>& vplane_rms, std::vector<float>& wplane_rms, WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map, TH2F *hu_decon_g, TH2F *hv_decon_g, TH2F *hw_decon_g, int nrebin, TH2F *hv_raw, bool flag_corr){

  int nwire_u = hu_decon->GetNbinsX();
  int nwire_v = hv_decon->GetNbinsX();
  int nwire_w = hw_decon->GetNbinsX();
  
  // filter noisy events ... 
  int n_time_slice = hu_decon->GetNbinsY();
  int length_cut = 3;
  int n_rebin = 1;
  int time_cut = 3;

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
      float rms = uplane_rms.at(j);
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
      float rms = vplane_rms.at(j);
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
      float rms = wplane_rms.at(j);
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
  int acc_cover = 0;
  int acc_total = 0;
  int acc_fire = 0;

  std::vector<std::pair<int,int>> noisy_u;
  std::vector<std::pair<int,int>> noisy_v;
  std::vector<std::pair<int,int>> noisy_w;
  
  
  for (size_t i=0;i!=summary_u.size();i++){
    int time = std::get<0>(summary_u.at(i));
    int ncover = std::get<1>(summary_u.at(i));
    float percentage = std::get<2>(summary_u.at(i));
    
    if (time > prev_time+time_cut){
      //std::cout << n_cover << " " << n_fire << " " << acc_fire*1./acc_cover << " " << acc_fire * 1./acc_total << " " << std::endl;
      if (n_cover >=12 && n_fire>=14 || (n_cover>=6 && n_fire>=6 &&acc_fire > 0.22 * acc_total)){
	end_time = prev_time;
	std::cout << "U: " << n_cover << " " << n_fire << " " << start_time << " " << end_time << std::endl;
	noisy_u.push_back(std::make_pair(start_time,end_time));
	flag_u = true;
      }
      
      n_cover = 0;
      n_fire = 0;
      acc_cover = 0;
      acc_total = 0;
      acc_fire = 0;
      start_time = time;
      end_time = time;
    }
    //std::cout << time << " " << ncover << " " << ncover*percentage << std::endl;

    
    
    if (ncover>1200. && (percentage > 0.25 || ncover *percentage > 300))
      n_cover ++;
    if (ncover*percentage>150){
      n_fire ++;
      acc_total += nwire_u;
      acc_cover += ncover;
      acc_fire  += ncover*percentage;
    }
    prev_time = time;
  }
  //std::cout << n_cover << " " << n_fire << " " << acc_fire*1./acc_cover << " " << acc_fire * 1./acc_total << " " << std::endl;
  if (n_cover >=12 && n_fire>=14 || (n_cover>=6 && n_fire>=6 &&acc_fire > 0.22 * acc_total)){
    end_time = prev_time;
    flag_u = true;
    std::cout << "U: " << n_cover << " " << n_fire << " " << start_time << " " << end_time << std::endl;
    noisy_u.push_back(std::make_pair(start_time,end_time));
  }

  bool flag_v = false;
  prev_time = -1;
  n_cover = 0;
  n_fire = 0;
  start_time =0;
  end_time = 0;
  acc_cover = 0;
  acc_total = 0;
  acc_fire = 0;
  
  for (size_t i=0;i!=summary_v.size();i++){
    int time = std::get<0>(summary_v.at(i));
    int ncover = std::get<1>(summary_v.at(i));
    float percentage = std::get<2>(summary_v.at(i));
    
    if (time > prev_time+time_cut){
      // std::cout << n_cover << " " << n_fire << std::endl;
      if (n_cover >=12 && n_fire>=14 || (n_cover>=6 && n_fire>=6 &&acc_fire > 0.22 * acc_total)){
	end_time = prev_time;
	flag_v = true;
	std::cout << "V: " << n_cover << " " << n_fire << " " << start_time << " " << prev_time << std::endl;
	noisy_v.push_back(std::make_pair(start_time,end_time));
      }
      
      n_cover = 0;
      n_fire = 0;
      acc_cover = 0;
      acc_total = 0;
      acc_fire = 0;
      start_time = time;
      end_time = time;
    }
    //std::cout << time << " " << ncover << " " << ncover*percentage << std::endl;
    
    if (ncover>1200. && (percentage > 0.25 || ncover *percentage > 300))
      n_cover ++;
    if (ncover*percentage>150){
      n_fire ++;
      acc_total += nwire_u;
      acc_cover += ncover;
      acc_fire  += ncover*percentage;
    }
    prev_time = time;
  }
  //std::cout << n_cover << " " << n_fire << std::endl;
  if (n_cover >=12 && n_fire>=14 || (n_cover>=6 && n_fire>=6 &&acc_fire > 0.22 * acc_total)){
    end_time = prev_time;
    flag_v = true;
    std::cout << "V: " << n_cover << " " << n_fire << " " << start_time << " " << end_time << std::endl;
    noisy_v.push_back(std::make_pair(start_time,end_time));
  }

  bool flag_w = false;
  prev_time = -1;
  n_cover = 0;
  n_fire = 0;
  start_time = 0;
  end_time = 0;
  acc_cover = 0;
  acc_total = 0;
  acc_fire = 0;
  
  for (size_t i=0;i!=summary_w.size();i++){
    int time = std::get<0>(summary_w.at(i));
    int ncover = std::get<1>(summary_w.at(i));
    float percentage = std::get<2>(summary_w.at(i));
    
    if (time > prev_time+time_cut){
      //std::cout << n_cover << " " << n_fire << std::endl;
      if (n_cover >=12 && n_fire>=14 || (n_cover>=6 && n_fire>=6 &&acc_fire > 0.22 * acc_total)){
	end_time = prev_time;
	flag_w = true;
	std::cout << "W: " << n_cover << " " << n_fire << " " << start_time << " " << end_time << std::endl;
	noisy_w.push_back(std::make_pair(start_time,end_time));
      }
      
      n_cover = 0;
      n_fire = 0;
      acc_cover = 0;
      acc_total = 0;
      acc_fire = 0;
      start_time = time;
      end_time = time;
    }
    //   std::cout << time << " " << ncover << " " << ncover*percentage << std::endl;
    
    if (ncover>1800. && (percentage > 0.2 || ncover *percentage > 360))
      n_cover ++;
    if (ncover*percentage>180){
      n_fire ++;
      acc_total += nwire_u;
      acc_cover += ncover;
      acc_fire  += ncover*percentage;
    }
    prev_time = time;
  }
  if (n_cover >=12 && n_fire>=14 || (n_cover>=6 && n_fire>=6 &&acc_fire > 0.22 * acc_total)){
    end_time = prev_time;
    flag_w = true;
    std::cout << "W: " << n_cover << " " << n_fire << " " << start_time << " " << end_time << std::endl;
    noisy_w.push_back(std::make_pair(start_time,end_time));
  }

  
  std::cout << "Xin: " << " " << flag_u << " " << flag_v << " " << flag_w << std::endl;

  
  int min_time = 3180;
  int max_time = 7870;
  bool flag_active = false;
    
  if (flag_corr){
    // ID the region ...
    for (auto it = noisy_u.begin(); it!= noisy_u.end(); it++){
      int start_time = (it->first-1) * nrebin+1;
      int end_time = it->second * nrebin+1;

      if (end_time > min_time && start_time < max_time)
	flag_active = true;
      
      for (int i=0;i!=nwire_u;i++){
	for (int j=it->first; j!=it->second+1;j++){
	  hu_decon->SetBinContent(i+1, j+1, 0);
	  hu_decon_g->SetBinContent(i+1, j+1, 0);
	}
	if (uplane_map.find(i)==uplane_map.end()){
	  uplane_map[i] = std::make_pair(start_time,end_time);
	}else{
	  uplane_map[i] = std::make_pair(std::min(uplane_map[i].first,start_time),std::max(uplane_map[i].second,end_time));
	}
      }
      for (int i=0;i!=nwire_v;i++){
	for (int j=it->first; j!=it->second+1;j++){
	  hv_decon->SetBinContent(i+1, j+1, 0);
	  hv_decon_g->SetBinContent(i+1, j+1, 0);
	  for (int k=0;k!=nrebin;k++){
	    hv_raw->SetBinContent(i+1,nrebin*j+k+1,0);
	  }
	}
	if (vplane_map.find(i)==vplane_map.end()){
	  vplane_map[i] = std::make_pair(start_time,end_time);
	}else{
	  vplane_map[i] = std::make_pair(std::min(vplane_map[i].first,start_time),std::max(vplane_map[i].second,end_time));
	}
      }
      for (int i=0;i!=nwire_w;i++){
	for (int j=it->first; j!=it->second+1;j++){
	  hw_decon->SetBinContent(i+1, j+1, 0);
	  hw_decon_g->SetBinContent(i+1, j+1, 0);
	}
	if (wplane_map.find(i)==wplane_map.end()){
	  wplane_map[i] = std::make_pair(start_time,end_time);
	}else{
	  wplane_map[i] = std::make_pair(std::min(wplane_map[i].first,start_time),std::max(wplane_map[i].second,end_time));
	}
      }
    }

    for (auto it = noisy_v.begin(); it!=noisy_v.end(); it++){
      // bool flag_veto = false;
      // for (auto it1 = noisy_u.begin(); it1!=noisy_u.end(); it1++){
      // 	if ( it1->second >= it->first && it->second >= it1->first){
      // 	  flag_veto = true;
      // 	  break;
      // 	}
      // }
      // if (flag_veto){
      int start_time = (it->first-1) * nrebin+1;
      int end_time = it->second * nrebin+1;

      if (end_time > min_time && start_time < max_time)
	flag_active = true;
      
      for (int i=0;i!=nwire_u;i++){
	for (int j=it->first; j!=it->second+1;j++){
	  hu_decon->SetBinContent(i+1, j+1, 0);
	  hu_decon_g->SetBinContent(i+1, j+1, 0);
	}
	if (uplane_map.find(i)==uplane_map.end()){
	  uplane_map[i] = std::make_pair(start_time,end_time);
	}else{
	  uplane_map[i] = std::make_pair(std::min(uplane_map[i].first,start_time),std::max(uplane_map[i].second,end_time));
	}
      }
      for (int i=0;i!=nwire_v;i++){
	for (int j=it->first; j!=it->second+1;j++){
	  hv_decon->SetBinContent(i+1, j+1, 0);
	  hv_decon_g->SetBinContent(i+1, j+1, 0);
	  for (int k=0;k!=nrebin;k++){
	    hv_raw->SetBinContent(i+1,nrebin*j+k+1,0);
	  }
	}
	if (vplane_map.find(i)==vplane_map.end()){
	  vplane_map[i] = std::make_pair(start_time,end_time);
	}else{
	  vplane_map[i] = std::make_pair(std::min(vplane_map[i].first,start_time),std::max(vplane_map[i].second,end_time));
	}
      }
      for (int i=0;i!=nwire_w;i++){
	for (int j=it->first; j!=it->second+1;j++){
	  hw_decon->SetBinContent(i+1, j+1, 0);
	  hw_decon_g->SetBinContent(i+1, j+1, 0);
	}
	if (wplane_map.find(i)==wplane_map.end()){
	  wplane_map[i] = std::make_pair(start_time,end_time);
	}else{
	  wplane_map[i] = std::make_pair(std::min(wplane_map[i].first,start_time),std::max(wplane_map[i].second,end_time));
	}
      }
	//}
    }

    for (auto it = noisy_w.begin(); it!=noisy_w.end(); it++){
      // bool flag_veto1 = false;
      // bool flag_veto2 = false;
      // for (auto it1 = noisy_u.begin(); it1!=noisy_u.end(); it1++){
      // 	if ( it1->second >= it->first && it->second >= it1->first){
      // 	  flag_veto1 = true;
      // 	  break;
      // 	}
      // }
      // for (auto it1 = noisy_v.begin(); it1!=noisy_v.end(); it1++){
      // 	if ( it1->second >= it->first && it->second >= it1->first){
      // 	  flag_veto2 = true;
      // 	  break;
      // 	}
      // }
      // if (flag_veto1 && flag_veto2){
      int start_time = (it->first-1) * nrebin+1;
      int end_time = it->second * nrebin+1;

      if (end_time > min_time && start_time < max_time)
	flag_active = true;
      
      for (int i=0;i!=nwire_u;i++){
	for (int j=it->first; j!=it->second+1;j++){
	  hu_decon->SetBinContent(i+1, j+1, 0);
	  hu_decon_g->SetBinContent(i+1, j+1, 0);
	}
	if (uplane_map.find(i)==uplane_map.end()){
	  uplane_map[i] = std::make_pair(start_time,end_time);
	}else{
	  uplane_map[i] = std::make_pair(std::min(uplane_map[i].first,start_time),std::max(uplane_map[i].second,end_time));
	}
      }
      for (int i=0;i!=nwire_v;i++){
	for (int j=it->first; j!=it->second+1;j++){
	  hv_decon->SetBinContent(i+1, j+1, 0);
	  hv_decon_g->SetBinContent(i+1, j+1, 0);
	  for (int k=0;k!=nrebin;k++){
	    hv_raw->SetBinContent(i+1,nrebin*j+k+1,0);
	  }
	}
	if (vplane_map.find(i)==vplane_map.end()){
	  vplane_map[i] = std::make_pair(start_time,end_time);
	}else{
	  vplane_map[i] = std::make_pair(std::min(vplane_map[i].first,start_time),std::max(vplane_map[i].second,end_time));
	}
      }
      for (int i=0;i!=nwire_w;i++){
	for (int j=it->first; j!=it->second+1;j++){
	  hw_decon->SetBinContent(i+1, j+1, 0);
	  hw_decon_g->SetBinContent(i+1, j+1, 0);
	}
	if (wplane_map.find(i)==wplane_map.end()){
	  wplane_map[i] = std::make_pair(start_time,end_time);
	}else{
	  wplane_map[i] = std::make_pair(std::min(wplane_map[i].first,start_time),std::max(wplane_map[i].second,end_time));
	}
      }
    }
    
    //}
  }

  if (flag_u||flag_v||flag_w){
    if (flag_active){
      return 2;
    }else{
      return 1;
    }
  }else{
    return 0;
  }
  // std::cout << n_cover << " " << n_fire << std::endl;
  
}
