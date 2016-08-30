#include "WireCell2dToy/uBooNE_Data_After_ROI.h"

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
  TH1F *htemp_flag = new TH1F("htemp_flag","htemp_flag",nticks,0,nticks);
  
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

  


  std::vector<std::pair<int,int>>& uboone_rois = rois.get_loose_rois(0);

  
  
  for (size_t ind=0; ind<ntraces; ++ind) {
    htemp_signal->Reset();
    htemp_flag->Reset();
    
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();

    //copy signal
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp_signal->SetBinContent(tt,trace.charge.at(i));
    }
    
    //set ROIs
    uboone_rois = rois.get_loose_rois(chid);
   

    // std::cout << htemp_signal->GetBinContent(100) << " " << uboone_rois.size() << std::endl;
    
    for (int i=0;i!=uboone_rois.size();i++){
      for (int j=uboone_rois.at(i).first; j<=uboone_rois.at(i).second; j++){
	htemp_flag->SetBinContent(j+1,1);
      }
    }
    //remove outside ROIs
    for (int i=0;i!=htemp_signal->GetNbinsX();i++){
      if (htemp_flag->GetBinContent(i+1)==0)
	htemp_signal->SetBinContent(i+1,0);
    }

    //apply adaptive baseline ... 
    for (int i=0;i!=uboone_rois.size();i++){
      int start_bin = uboone_rois.at(i).first;
      int end_bin = uboone_rois.at(i).second;
      float start_content = htemp_signal->GetBinContent(start_bin+1);
      float end_content = htemp_signal->GetBinContent(end_bin+1);

      for (int j=start_bin;j<=end_bin;j++){
	float content = start_content + (j-start_bin) * 1.0 /(end_bin - start_bin) * (end_content-start_content);
	htemp_signal->SetBinContent(j+1,htemp_signal->GetBinContent(j+1) - content);
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
	//sum += htemp_flag->GetBinContent(i*rebin+j+1);
      }
      trace1.charge.at(i) = sum;
      //      std::cout << rebin << " " << i << " " << sum << std::endl;
    }
    frame.traces.push_back(trace1);
  }
  
  delete htemp_signal;
  delete htemp_flag;

  frame.index = frame_number;
  return frame.index;

}
