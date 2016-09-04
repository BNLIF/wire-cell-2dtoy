#include "WireCell2dToy/uBooNE_Data_After_ROI.h"
#include "TSpectrum.h"


using namespace WireCell;


WireCell2dToy::uBooNEDataAfterROI::uBooNEDataAfterROI(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, WireCell2dToy::uBooNEDataROI& rois, int rebin)
  : fds(fds)
  , gds(gds)
  , rois(rois)
  , rebin(rebin)
{
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  // rois_u_tight.resize(nwire_u);
  // rois_v_tight.resize(nwire_v);
  // rois_w_tight.resize(nwire_w);

  // rois_u_loose.resize(nwire_u);
  // rois_v_loose.resize(nwire_v);
  // rois_w_loose.resize(nwire_w);
  
  for (Int_t i=0;i!=nwire_u;i++){
    SignalROIList temp_rois;
    rois_u_tight.push_back(temp_rois);
  }
  for (Int_t i=0;i!=nwire_v;i++){
    SignalROIList temp_rois;
    rois_v_tight.push_back(temp_rois);
  }
  for (Int_t i=0;i!=nwire_w;i++){
    SignalROIList temp_rois;
    rois_w_tight.push_back(temp_rois);
  }
  
  for (Int_t i=0;i!=nwire_u;i++){
    SignalROIList temp_rois;
    rois_u_loose.push_back(temp_rois);
  }
  for (Int_t i=0;i!=nwire_v;i++){
    SignalROIList temp_rois;
    rois_v_loose.push_back(temp_rois);
  }
  for (Int_t i=0;i!=nwire_w;i++){
    SignalROIList temp_rois;
    rois_w_loose.push_back(temp_rois);
  }
  

  bins_per_frame = int(fds.Get_Bins_Per_Frame()/rebin);
}


WireCell2dToy::uBooNEDataAfterROI::~uBooNEDataAfterROI(){
  
}

int WireCell2dToy::uBooNEDataAfterROI::size() const{
  return 1;
}

void WireCell2dToy::uBooNEDataAfterROI::Clear(){
  for (int i=0;i!=nwire_u;i++){
    for (auto it = rois_u_tight.at(i).begin(); it!=rois_u_tight.at(i).end();it++){
      delete *it;
    }
    rois_u_tight.at(i).clear();
    for (auto it = rois_u_loose.at(i).begin(); it!=rois_u_loose.at(i).end();it++){
      delete *it;
    }
    rois_u_loose.at(i).clear();
  }

  for (int i=0;i!=nwire_v;i++){
    for (auto it = rois_v_tight.at(i).begin(); it!=rois_v_tight.at(i).end();it++){
      delete *it;
    }
    rois_v_tight.at(i).clear();
    for (auto it = rois_v_loose.at(i).begin(); it!=rois_v_loose.at(i).end();it++){
      delete *it;
    }
    rois_v_loose.at(i).clear();
  }

  for (int i=0;i!=nwire_w;i++){
    for (auto it = rois_w_tight.at(i).begin(); it!=rois_w_tight.at(i).end();it++){
      delete *it;
    }
    rois_w_tight.at(i).clear();
    for (auto it = rois_w_loose.at(i).begin(); it!=rois_w_loose.at(i).end();it++){
      delete *it;
    }
    rois_w_loose.at(i).clear();
  }
   
   front_rois.clear();
   back_rois.clear();
   contained_rois.clear();
   


}

void WireCell2dToy::uBooNEDataAfterROI::CleanUpROIs(){
  // clean up ROIs
  std::map<SignalROI*, int> ROIsaved_map;
  //int counter = 0;
  for (int i=0;i!=rois_u_loose.size();i++){
    // counter += rois_u_loose.at(i).size();
    for (auto it = rois_u_loose.at(i).begin(); it!= rois_u_loose.at(i).end();it++){
      SignalROI *roi = *it;
      if (ROIsaved_map.find(roi)==ROIsaved_map.end()){
     	if (contained_rois.find(roi) != contained_rois.end()){
	  // contain good stuff
	  SignalROISelection temp_rois;
	  temp_rois.push_back(roi);
	  ROIsaved_map[roi] = 1;
	  
	  while(temp_rois.size()){
	    SignalROI *temp_roi = temp_rois.back();
	    temp_rois.pop_back();
	    // save all its neighbour into a temporary holder
	    if (front_rois.find(temp_roi)!=front_rois.end()){
	      for (auto it1 = front_rois[temp_roi].begin();it1!=front_rois[temp_roi].end();it1++){
		if (ROIsaved_map.find(*it1)==ROIsaved_map.end()){
		  temp_rois.push_back(*it1);
		  ROIsaved_map[*it1] = 1;
		}
	      }
	    }
	    if (back_rois.find(temp_roi)!=back_rois.end()){
	      for (auto it1 = back_rois[temp_roi].begin();it1!=back_rois[temp_roi].end();it1++){
		if (ROIsaved_map.find(*it1)==ROIsaved_map.end()){
		  temp_rois.push_back(*it1);
		  ROIsaved_map[*it1] = 1;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  //remove the bad ones ...
  //int counter2 = 0;
  for (int i=0;i!=rois_u_loose.size();i++){
    SignalROISelection to_be_removed;
    for (auto it = rois_u_loose.at(i).begin(); it!= rois_u_loose.at(i).end();it++){
      SignalROI *roi = *it;
      if (ROIsaved_map.find(roi) == ROIsaved_map.end()){
	//	counter2 ++;
	to_be_removed.push_back(roi);
	//it = rois_u_loose.at(i).erase(it);
	// check contained map
	if (contained_rois.find(roi)!= contained_rois.end()){
	  std::cout << "Wrong! " << std::endl;
	}
	// check front map
	if (front_rois.find(roi)!=front_rois.end()){
	  for (auto it1 = front_rois[roi].begin(); it1 != front_rois[roi].end(); it1++){
	    auto it2 = find(back_rois[*it1].begin(),back_rois[*it1].end(),roi);
	    back_rois[*it1].erase(it2);
	  }
	  front_rois.erase(roi);
	}
	// check back map
	if (back_rois.find(roi)!=back_rois.end()){
	  for (auto it1 = back_rois[roi].begin(); it1!=back_rois[roi].end(); it1++){
	    auto it2 = find(front_rois[*it1].begin(),front_rois[*it1].end(),roi);
	    front_rois[*it1].erase(it2);
	  }
	  back_rois.erase(roi);
	}
      }
    }

    for (auto it = to_be_removed.begin(); it!= to_be_removed.end(); it++){
      auto it1 = find(rois_u_loose.at(i).begin(), rois_u_loose.at(i).end(),*it);
      rois_u_loose.at(i).erase(it1);
    }
  }


  // int counter = 0;
  for (int i=0;i!=rois_v_loose.size();i++){
    //  counter += rois_v_loose.at(i).size();
    for (auto it = rois_v_loose.at(i).begin(); it!= rois_v_loose.at(i).end();it++){
      SignalROI *roi = *it;
      if (ROIsaved_map.find(roi)==ROIsaved_map.end()){
     	if (contained_rois.find(roi) != contained_rois.end()){
	  // contain good stuff
	  SignalROISelection temp_rois;
	  temp_rois.push_back(roi);
	  ROIsaved_map[roi] = 1;
	  
	  while(temp_rois.size()){
	    SignalROI *temp_roi = temp_rois.back();
	    temp_rois.pop_back();
	    // save all its neighbour into a temporary holder
	    if (front_rois.find(temp_roi)!=front_rois.end()){
	      for (auto it1 = front_rois[temp_roi].begin();it1!=front_rois[temp_roi].end();it1++){
		if (ROIsaved_map.find(*it1)==ROIsaved_map.end()){
		  temp_rois.push_back(*it1);
		  ROIsaved_map[*it1] = 1;
		}
	      }
	    }
	    if (back_rois.find(temp_roi)!=back_rois.end()){
	      for (auto it1 = back_rois[temp_roi].begin();it1!=back_rois[temp_roi].end();it1++){
		if (ROIsaved_map.find(*it1)==ROIsaved_map.end()){
		  temp_rois.push_back(*it1);
		  ROIsaved_map[*it1] = 1;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  //remove the bad ones ...
  //int counter2 = 0;
  for (int i=0;i!=rois_v_loose.size();i++){
    SignalROISelection to_be_removed;
    for (auto it = rois_v_loose.at(i).begin(); it!= rois_v_loose.at(i).end();it++){
      SignalROI *roi = *it;
      if (ROIsaved_map.find(roi) == ROIsaved_map.end()){
	//	counter2 ++;
	to_be_removed.push_back(roi);
	//it = rois_v_loose.at(i).erase(it);
	// check contained map
	if (contained_rois.find(roi)!= contained_rois.end()){
	  std::cout << "Wrong! " << std::endl;
	}
	// check front map
	if (front_rois.find(roi)!=front_rois.end()){
	  for (auto it1 = front_rois[roi].begin(); it1 != front_rois[roi].end(); it1++){
	    auto it2 = find(back_rois[*it1].begin(),back_rois[*it1].end(),roi);
	    back_rois[*it1].erase(it2);
	  }
	  front_rois.erase(roi);
	}
	// check back map
	if (back_rois.find(roi)!=back_rois.end()){
	  for (auto it1 = back_rois[roi].begin(); it1!=back_rois[roi].end(); it1++){
	    auto it2 = find(front_rois[*it1].begin(),front_rois[*it1].end(),roi);
	    front_rois[*it1].erase(it2);
	  }
	  back_rois.erase(roi);
	}
      }
    }

    for (auto it = to_be_removed.begin(); it!= to_be_removed.end(); it++){
      auto it1 = find(rois_v_loose.at(i).begin(), rois_v_loose.at(i).end(),*it);
      rois_v_loose.at(i).erase(it1);
    }
  }




  // int counter1 = 0;
  // for (int i=0;i!=rois_v_loose.size();i++){
  //   counter1+=rois_v_loose.at(i).size();
  // }
  
  // std::cout << counter << " " << ROIsaved_map.size() << " " << counter1 << " " << counter2 << std::endl;


}

void WireCell2dToy::uBooNEDataAfterROI::BreakROI(SignalROI* roi, float rms){
  // main algorithm 
  int start_bin = roi->get_start_bin();
  int end_bin = roi->get_end_bin();
  TH1F *htemp = new TH1F("htemp","htemp",end_bin-start_bin+1,start_bin,end_bin+1);
  std::vector<float>& contents = roi->get_contents();
  for (Int_t i=0;i!=htemp->GetNbinsX();i++){
    htemp->SetBinContent(i+1,contents.at(i));
  }
  TSpectrum *s = new TSpectrum(100);
  Int_t nfound = s->Search(htemp,2,"nobackground new",0.1);

  Int_t flag_single_peak = 0;

  float th_peak = 3.0;
  float tail_th = 2.5;
  float th_1d = 3.0;
    
  std::set<int> saved_boundaries;

  if (nfound > 1){
    Int_t npeaks = s->GetNPeaks();
    Double_t *peak_pos = s->GetPositionX();
    Double_t *peak_height = s->GetPositionY();
    
    int order_peak_pos[105];
    int npeaks_threshold = 0;
    for (Int_t j=0;j!=npeaks;j++){
      order_peak_pos[j] = *(peak_pos+j);
      if (*(peak_height+j)>th_peak*rms){
	npeaks_threshold ++;
      }
    }

    if (npeaks_threshold >1){
      std::sort(order_peak_pos,order_peak_pos + npeaks);
      
      float valley_pos[100];
      valley_pos[0] = start_bin;
      
      for (Int_t j=0;j!=npeaks-1;j++){
	Float_t min = 1e9;
	for (Int_t k = order_peak_pos[j]-start_bin; k< order_peak_pos[j+1]-start_bin;k++){
	  
	  if (htemp->GetBinContent(k+1) < min){
	    min = htemp->GetBinContent(k+1);
	    valley_pos[j+1] = k+start_bin;
	  }
	}
	//	std::cout << valley_pos[j+1] << std::endl;
	//std::cout << *(peak_pos+j) << std::endl;
      }
      
      valley_pos[npeaks] = end_bin;
      
     
      
      for (Int_t j=0;j!=npeaks;j++){
	int flag = 0;
	Int_t start_pos = valley_pos[j];
	Double_t start_content = htemp->GetBinContent(valley_pos[j]-start_bin+1);
	Int_t end_pos = valley_pos[j+1];
	Double_t end_content = htemp->GetBinContent(valley_pos[j+1]-start_bin+1);
	Int_t Peak_pos = order_peak_pos[j];
	Double_t peak_content = htemp->GetBinContent(order_peak_pos[j]-start_bin+1);
      
	if ((start_content > rms || end_content > rms) &&
	    start_content < peak_content /2. && end_content < peak_content/2.) flag = 1;
	
	// deal with the small peaks ... 
	if (peak_content > rms && peak_content < rms * 1.5) flag =1;
	
	if (flag==1){
	  saved_boundaries.insert(start_pos);
	  saved_boundaries.insert(end_pos);
	  //	  std::cout << start_pos << " " << start_content << " " << end_pos << 
	  //  " " << end_content << std::endl;
	}
      }
      

      TH1F *htemp1 = (TH1F*)htemp->Clone("htemp1");
      for (Int_t j=0;j!=npeaks;j++){
	int flag = 0;
	Int_t start_pos = valley_pos[j];
	Double_t start_content = htemp1->GetBinContent(valley_pos[j]-start_bin+1);
	Int_t end_pos = valley_pos[j+1];
	Double_t end_content = htemp1->GetBinContent(valley_pos[j+1]-start_bin+1);
	Int_t Peak_pos = order_peak_pos[j];
	Double_t peak_content = htemp1->GetBinContent(order_peak_pos[j]-start_bin+1);
	
	if (saved_boundaries.find(start_pos) != saved_boundaries.end() ||
	    saved_boundaries.find(end_pos) != saved_boundaries.end()){
	  
	  //std::cout << j << std::endl;
	  
	  float max = 0;
	  float ave = 0;
	  float ave1 = 0;
	  
	  for (Int_t k = start_pos; k!=end_pos+1;k++){
	    Double_t temp_content = htemp1->GetBinContent(k-start_bin+1) - (start_content + (end_content-start_content) * (k-start_pos) / (end_pos - start_pos));
	    if (temp_content > rms){
	      ave += temp_content;
	      ave1 += 1;
	    }
	    if (temp_content > max)
	      max = temp_content;
	    //	      htemp->SetBinContent(k-begin+1,temp_content);
	  }
	  //std::cout << "abc: " << start_pos << " " << end_pos << " " << max << " " << ave << " " << ave1 << std::endl;
	  
	  if (max *ave1 > tail_th*ave){
	    //std::cout << "abc: " << start_pos << " " << end_pos << " " << max << " " << ave << " " << ave1 << std::endl;
	    for (Int_t k = start_pos; k!=end_pos+1;k++){
	      Double_t temp_content = htemp1->GetBinContent(k-start_bin+1) - (start_content + (end_content-start_content) * (k-start_pos) / (end_pos - start_pos));
	      if (temp_content > th_1d*rms){
		htemp->SetBinContent(k-start_bin+1,temp_content);
	      }else{
		htemp->SetBinContent(k-start_bin+1,0);
	      }
	    }
	  }else{
	    for (Int_t k = start_pos; k!=end_pos+1;k++){
	      Double_t temp_content = htemp1->GetBinContent(k-start_bin+1) - (start_content + (end_content-start_content) * (k-start_pos) / (end_pos - start_pos));
	      htemp->SetBinContent(k-start_bin+1,temp_content);
	    }
	  }
	  
	     
	}
      }
      delete htemp1;

    } else{
      flag_single_peak = 1;
    }
  } else{
    flag_single_peak = 1;
  }

  if (flag_single_peak == 1){
    // find maximum and average
    float max = 0;
    float ave = 0;
    float ave1 = 0;
    for (Int_t j=0;j!=htemp->GetNbinsX();j++){
      if (htemp->GetBinContent(j+1) > rms){
  	ave += htemp->GetBinContent(j+1);
  	ave1 += 1;
      }
      if (htemp->GetBinContent(j+1) > max){
  	max = htemp->GetBinContent(j+1);
      }
    }
    
    if (max *ave1 > tail_th*ave){
      for (Int_t j=0;j!=htemp->GetNbinsX();j++){
  	if (htemp->GetBinContent(j+1) < th_1d*rms){
  	  htemp->SetBinContent(j+1,0);
  	}
      }
    }
    
  }
  
   for (Int_t i=0;i!=htemp->GetNbinsX();i++){
     //htemp->SetBinContent(i+1,contents.at(i));
     contents.at(i) = htemp->GetBinContent(i+1);
   }
  

  delete s;
  delete htemp;
}

void WireCell2dToy::uBooNEDataAfterROI::BreakROIs(){
  // get RMS value, and put in 
  std::vector<float>& rms_u = rois.get_uplane_rms();
  std::vector<float>& rms_v = rois.get_vplane_rms();

  for (int i=0;i!=rois_u_loose.size();i++){
    for (auto it = rois_u_loose.at(i).begin(); it!= rois_u_loose.at(i).end(); it++){
      BreakROI(*it,rms_u.at(i));
    }
  }
  
  for (int i=0;i!=rois_v_loose.size();i++){
    for (auto it = rois_v_loose.at(i).begin(); it!= rois_v_loose.at(i).end(); it++){
      BreakROI(*it,rms_v.at(i));
    }
  }
}

void WireCell2dToy::uBooNEDataAfterROI::ShrinkROIs(){
  
}


int WireCell2dToy::uBooNEDataAfterROI::jump(int frame_number){
  Clear();
  if (frame.index == frame_number) {
    return frame_number;
  }
  frame.clear();

  fds.jump(frame_number);
  const Frame& frame1 = fds.get();
  size_t ntraces = frame1.traces.size();
  int nticks = fds.Get_Bins_Per_Frame();

  TH1F *htemp_signal = new TH1F("htemp_signal","htemp_signal",nticks,0,nticks);
  
  for (size_t ind =0; ind<ntraces; ++ind) {
    htemp_signal->Reset();
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();

    //copy signal
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp_signal->SetBinContent(tt,trace.charge.at(i));
    }
    
    std::vector<std::pair<int,int>>& uboone_rois = rois.get_self_rois(chid);
    for (int i=0;i!=uboone_rois.size();i++){
      //std::cout << nbins << " " << i << " " << uboone_rois.size() << " " << chid << " " << uboone_rois.at(i).first << " " << uboone_rois.at(i).second << std::endl;
      SignalROI *tight_roi;
      if (chid < nwire_u){
	tight_roi = new SignalROI(WirePlaneType_t(0), chid,uboone_rois.at(i).first,uboone_rois.at(i).second, htemp_signal);
      }else if (chid < nwire_u+nwire_v){
	tight_roi = new SignalROI(WirePlaneType_t(1), chid,uboone_rois.at(i).first,uboone_rois.at(i).second, htemp_signal);
      }else{
	tight_roi = new SignalROI(WirePlaneType_t(2), chid,uboone_rois.at(i).first,uboone_rois.at(i).second, htemp_signal);
      }
      // std::cout << i << std::endl;

      if (chid < nwire_u){
	rois_u_tight[chid].push_back(tight_roi);
       	if (chid>0){
       	  //form connectivity map
       	  for (auto it = rois_u_tight[chid-1].begin();it!=rois_u_tight[chid-1].end();it++){
       	    SignalROI *prev_roi = *it;
      	    if (tight_roi->overlap(prev_roi)){
      	      if (front_rois.find(prev_roi) == front_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(tight_roi);
      		front_rois[prev_roi] = temp_rois;
      	      }else{
      		front_rois[prev_roi].push_back(tight_roi);
      	      }
      	      if (back_rois.find(tight_roi) == back_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(prev_roi);
      		back_rois[tight_roi] = temp_rois;
      	      }else{
      		back_rois[tight_roi].push_back(prev_roi);
      	      }
      	    }
       	  }
     	}

	if (chid < nwire_u-1){
	  // add the code for the next one to be completed
	  for (auto it = rois_u_tight[chid+1].begin();it!=rois_u_tight[chid+1].end();it++){
	    SignalROI *next_roi = *it;
	    if (tight_roi->overlap(next_roi)){
	      if (back_rois.find(next_roi) == back_rois.end()){
		SignalROISelection temp_rois;
		temp_rois.push_back(tight_roi);
		back_rois[next_roi] = temp_rois;
	      }else{
		back_rois[next_roi].push_back(tight_roi);
	      }
	      
	      if (front_rois.find(tight_roi) == front_rois.end()){
		SignalROISelection temp_rois;
		temp_rois.push_back(next_roi);
		front_rois[tight_roi] = temp_rois;
	      }else{
		front_rois[tight_roi].push_back(next_roi);
	      }
	    }
	  }
	}
	
	
      }else if (chid < nwire_u + nwire_v){
       	rois_v_tight[chid-nwire_u].push_back(tight_roi);

      	if (chid>nwire_u){
      	  //form connectivity map
      	  for (auto it = rois_v_tight[chid - nwire_u-1].begin();it!=rois_v_tight[chid-nwire_u-1].end();it++){
      	    SignalROI *prev_roi = *it;
      	    if (tight_roi->overlap(prev_roi)){
      	      if (front_rois.find(prev_roi) == front_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(tight_roi);
      		front_rois[prev_roi] = temp_rois;
      	      }else{
      		front_rois[prev_roi].push_back(tight_roi);
      	      }
      	      if (back_rois.find(tight_roi) == back_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(prev_roi);
      		back_rois[tight_roi] = temp_rois;
      	      }else{
      		back_rois[tight_roi].push_back(prev_roi);
      	      }
      	    }
      	  }
      	}

	if (chid<nwire_u+nwire_v-1){
      	  //form connectivity map
      	  for (auto it = rois_v_tight[chid - nwire_u+1].begin();it!=rois_v_tight[chid-nwire_u+1].end();it++){
      	    SignalROI *next_roi = *it;
      	    if (tight_roi->overlap(next_roi)){
      	      if (back_rois.find(next_roi) == back_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(tight_roi);
      		back_rois[next_roi] = temp_rois;
      	      }else{
      		back_rois[next_roi].push_back(tight_roi);
      	      }
      	      if (front_rois.find(tight_roi) == front_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(next_roi);
      		front_rois[tight_roi] = temp_rois;
      	      }else{
      		front_rois[tight_roi].push_back(next_roi);
      	      }
      	    }
      	  }
      	}

	

      }else{
       	rois_w_tight[chid-nwire_u-nwire_v].push_back(tight_roi);

      	if (chid>nwire_u+nwire_v){
      	  //form connectivity map
      	  for (auto it = rois_w_tight[chid-nwire_u-nwire_v-1].begin();it!=rois_w_tight[chid-nwire_u-nwire_v-1].end();it++){
      	    SignalROI *prev_roi = *it;
      	    if (tight_roi->overlap(prev_roi)){
      	      if (front_rois.find(prev_roi) == front_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(tight_roi);
      		front_rois[prev_roi] = temp_rois;
      	      }else{
      		front_rois[prev_roi].push_back(tight_roi);
      	      }
      	      if (back_rois.find(tight_roi) == back_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(prev_roi);
      		back_rois[tight_roi] = temp_rois;
      	      }else{
      		back_rois[tight_roi].push_back(prev_roi);
      	      }
      	    }
      	  }
      	}


	if (chid<nwire_u+nwire_v+nwire_w-1){
      	  //form connectivity map
      	  for (auto it = rois_w_tight[chid-nwire_u-nwire_v+1].begin();it!=rois_w_tight[chid-nwire_u-nwire_v+1].end();it++){
      	    SignalROI *next_roi = *it;
      	    if (tight_roi->overlap(next_roi)){
      	      if (back_rois.find(next_roi) == back_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(tight_roi);
      		back_rois[next_roi] = temp_rois;
      	      }else{
      		back_rois[next_roi].push_back(tight_roi);
      	      }
      	      if (front_rois.find(tight_roi) == front_rois.end()){
      		SignalROISelection temp_rois;
      		temp_rois.push_back(next_roi);
      		front_rois[tight_roi] = temp_rois;
      	      }else{
      		front_rois[tight_roi].push_back(next_roi);
      	      }
      	    }
      	  }
      	}


      }
    }
    
    uboone_rois = rois.get_loose_rois(chid);
    for (int i = 0; i!=uboone_rois.size();i++){
      SignalROI *loose_roi;
      if (chid < nwire_u){
	loose_roi = new SignalROI(WirePlaneType_t(0),chid,uboone_rois.at(i).first,uboone_rois.at(i).second,htemp_signal);
      }else if (chid < nwire_u + nwire_v){
	loose_roi = new SignalROI(WirePlaneType_t(1),chid,uboone_rois.at(i).first,uboone_rois.at(i).second,htemp_signal);
      }else{
	loose_roi = new SignalROI(WirePlaneType_t(2),chid,uboone_rois.at(i).first,uboone_rois.at(i).second,htemp_signal);
      }
      
      if (chid < nwire_u){
    	rois_u_loose[chid].push_back(loose_roi);

    	if (chid>0){
    	  //form connectivity map
    	  for (auto it=rois_u_loose[chid-1].begin();it!=rois_u_loose[chid-1].end();it++){
    	    SignalROI *prev_roi = *it;
    	    if (loose_roi->overlap(prev_roi)){
    	      if (front_rois.find(prev_roi) == front_rois.end()){
    		SignalROISelection temp_rois;
    		temp_rois.push_back(loose_roi);
    		front_rois[prev_roi] = temp_rois;
    	      }else{
    		front_rois[prev_roi].push_back(loose_roi);
    	      }
    	      if (back_rois.find(loose_roi) == back_rois.end()){
    		SignalROISelection temp_rois;
    		temp_rois.push_back(prev_roi);
    		back_rois[loose_roi] = temp_rois;
    	      }else{
    		back_rois[loose_roi].push_back(prev_roi);
    	      }
    	    }
    	  }
    	}

	if (chid<nwire_u-1){
    	  //form connectivity map
    	  for (auto it=rois_u_loose[chid+1].begin();it!=rois_u_loose[chid+1].end();it++){
    	    SignalROI *next_roi = *it;
    	    if (loose_roi->overlap(next_roi)){
    	      if (back_rois.find(next_roi) == back_rois.end()){
    		SignalROISelection temp_rois;
    		temp_rois.push_back(loose_roi);
    		back_rois[next_roi] = temp_rois;
    	      }else{
    		back_rois[next_roi].push_back(loose_roi);
    	      }
    	      if (front_rois.find(loose_roi) == front_rois.end()){
    		SignalROISelection temp_rois;
    		temp_rois.push_back(next_roi);
    		front_rois[loose_roi] = temp_rois;
    	      }else{
    		front_rois[loose_roi].push_back(next_roi);
    	      }
    	    }
    	  }
    	}
	

    	//form contained map ... 
    	for (auto it=rois_u_tight[chid].begin();it!=rois_u_tight[chid].end();it++){
    	  SignalROI *tight_roi = *it;
    	  if (tight_roi->overlap(loose_roi)){
    	    if (contained_rois.find(loose_roi)==contained_rois.end()){
    	      SignalROISelection temp_rois;
    	      temp_rois.push_back(tight_roi);
    	      contained_rois[loose_roi] = temp_rois;
    	    }else{
    	      contained_rois[loose_roi].push_back(tight_roi);
    	    }
    	  }
    	}

      }else if (chid < nwire_u + nwire_v){
    	rois_v_loose[chid-nwire_u].push_back(loose_roi);

    	if (chid>nwire_u){
    	  //form connectivity map
    	  for (auto it = rois_v_loose[chid-nwire_u-1].begin();it!=rois_v_loose[chid-nwire_u-1].end();it++){
    	    SignalROI *prev_roi = *it;
    	    if (loose_roi->overlap(prev_roi)){
    	      if (front_rois.find(prev_roi) == front_rois.end()){
    		SignalROISelection temp_rois;
    		temp_rois.push_back(loose_roi);
    		front_rois[prev_roi] = temp_rois;
    	      }else{
    		front_rois[prev_roi].push_back(loose_roi);
    	      }
    	      if (back_rois.find(loose_roi) == back_rois.end()){
    		SignalROISelection temp_rois;
    		temp_rois.push_back(prev_roi);
    		back_rois[loose_roi] = temp_rois;
    	      }else{
    		back_rois[loose_roi].push_back(prev_roi);
    	      }
    	    }
    	  }
    	}

	if (chid<nwire_u+nwire_v-1){
    	  //form connectivity map
    	  for (auto it = rois_v_loose[chid-nwire_u+1].begin();it!=rois_v_loose[chid-nwire_u+1].end();it++){
    	    SignalROI *next_roi = *it;
    	    if (loose_roi->overlap(next_roi)){
    	      if (back_rois.find(next_roi) == back_rois.end()){
    		SignalROISelection temp_rois;
    		temp_rois.push_back(loose_roi);
    		back_rois[next_roi] = temp_rois;
    	      }else{
    		back_rois[next_roi].push_back(loose_roi);
    	      }
    	      if (front_rois.find(loose_roi) == front_rois.end()){
    		SignalROISelection temp_rois;
    		temp_rois.push_back(next_roi);
    		front_rois[loose_roi] = temp_rois;
    	      }else{
    		front_rois[loose_roi].push_back(next_roi);
    	      }
    	    }
    	  }
    	}

    	//form contained map ... 
    	for (auto it = rois_v_tight[chid-nwire_u].begin();it!=rois_v_tight[chid-nwire_u].end();it++){
    	  SignalROI *tight_roi = *it;
    	  if (tight_roi->overlap(loose_roi)){
    	    if (contained_rois.find(loose_roi)==contained_rois.end()){
    	      SignalROISelection temp_rois;
    	      temp_rois.push_back(tight_roi);
    	      contained_rois[loose_roi] = temp_rois;
    	    }else{
    	      contained_rois[loose_roi].push_back(tight_roi);
    	    }
    	  }
    	}

      }else{
    	rois_w_loose[chid-nwire_u-nwire_v].push_back(loose_roi);

    	if (chid>nwire_u+nwire_v){
    	  //form connectivity map
    	  for (auto it=rois_w_loose[chid-nwire_u-nwire_v-1].begin();it!=rois_w_loose[chid-nwire_u-nwire_v-1].end();it++){
    	    SignalROI *prev_roi = *it;
    	    if (loose_roi->overlap(prev_roi)){
    	      if (front_rois.find(prev_roi) == front_rois.end()){
    		SignalROISelection temp_rois;
    		temp_rois.push_back(loose_roi);
    		front_rois[prev_roi] = temp_rois;
    	      }else{
    		front_rois[prev_roi].push_back(loose_roi);
    	      }
    	      if (back_rois.find(loose_roi) == back_rois.end()){
    		SignalROISelection temp_rois;
    		temp_rois.push_back(prev_roi);
    		back_rois[loose_roi] = temp_rois;
    	      }else{
    		back_rois[loose_roi].push_back(prev_roi);
    	      }
    	    }
    	  }
    	}

	if (chid<nwire_u+nwire_v+nwire_w-1){
    	  //form connectivity map
    	  for (auto it=rois_w_loose[chid-nwire_u-nwire_v+1].begin();it!=rois_w_loose[chid-nwire_u-nwire_v+1].end();it++){
    	    SignalROI *next_roi = *it;
    	    if (loose_roi->overlap(next_roi)){
    	      if (back_rois.find(next_roi) == back_rois.end()){
    		SignalROISelection temp_rois;
    		temp_rois.push_back(loose_roi);
    		back_rois[next_roi] = temp_rois;
    	      }else{
    		back_rois[next_roi].push_back(loose_roi);
    	      }
    	      if (front_rois.find(loose_roi) == front_rois.end()){
    		SignalROISelection temp_rois;
    		temp_rois.push_back(next_roi);
    		front_rois[loose_roi] = temp_rois;
    	      }else{
    		front_rois[loose_roi].push_back(next_roi);
    	      }
    	    }
    	  }
    	}
	
    	//form contained map ... 
    	for (auto it=rois_w_tight[chid-nwire_u-nwire_v].begin();it!=rois_w_tight[chid-nwire_u-nwire_v].end();it++){
    	  SignalROI *tight_roi = *it;
    	  if (tight_roi->overlap(loose_roi)){
    	    if (contained_rois.find(loose_roi)==contained_rois.end()){
    	      SignalROISelection temp_rois;
    	      temp_rois.push_back(tight_roi);
    	      contained_rois[loose_roi] = temp_rois;
    	    }else{
    	      contained_rois[loose_roi].push_back(tight_roi);
    	    }
    	  }
	}
	
      }
    }
  }

  std::cout << rois_u_tight.size() << " " << rois_v_tight.size() << " " << rois_w_tight.size() << " " << rois_u_loose.size() << " " << rois_v_loose.size() << " " << rois_w_loose.size() << " " << front_rois.size() << " " << back_rois.size() << " " << contained_rois.size() << std::endl;
  CleanUpROIs();
  BreakROIs();
  

  // load results back into the data ... 
  
  for (size_t ind=0; ind<ntraces; ++ind) {
    htemp_signal->Reset();
    
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();

    // load results back into the histogram
    if (chid < nwire_u){
      for (auto it = rois_u_loose.at(chid).begin(); it!= rois_u_loose.at(chid).end();it++){
	SignalROI *roi =  *it;
	std::vector<float>& contents = roi->get_contents();
	int start_bin = roi->get_start_bin();
	for (int i=0;i!=contents.size();i++){
	  htemp_signal->SetBinContent(start_bin+1+i,contents.at(i));
	}
      }
    }else if (chid < nwire_u + nwire_v){
      for (auto it = rois_v_loose.at(chid-nwire_u).begin(); it!= rois_v_loose.at(chid-nwire_u).end();it++){
      	SignalROI *roi =  *it;
      	std::vector<float>& contents = roi->get_contents();
      	int start_bin = roi->get_start_bin();
      	for (int i=0;i!=contents.size();i++){
      	  htemp_signal->SetBinContent(start_bin+1+i,contents.at(i));
      	}
      }
    }else{
      for (auto it = rois_w_tight.at(chid-nwire_u-nwire_v).begin(); it!= rois_w_tight.at(chid-nwire_u-nwire_v).end();it++){
      	SignalROI *roi =  *it;
      	std::vector<float>& contents = roi->get_contents();
      	int start_bin = roi->get_start_bin();
      	for (int i=0;i!=contents.size();i++){
      	  htemp_signal->SetBinContent(start_bin+1+i,contents.at(i));
      	}
      }
    }
        
    // save into frames ... 
    WireCell::Trace trace1;
    trace1.chid = chid;
    trace1.tbin = 0;
    trace1.charge.resize(bins_per_frame, 0.0);
    
    for (int i=0;i!=bins_per_frame;i++){
      float sum = 0;
      for (int j=0;j!=rebin;j++){
	sum += htemp_signal->GetBinContent(i*rebin+j+1);
      }
      trace1.charge.at(i) = sum;
    }
    frame.traces.push_back(trace1);
  }
  
  delete htemp_signal;

  frame.index = frame_number;
  return frame.index;

}
