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

void WireCell2dToy::uBooNEDataAfterROI::CleanUpInductionROIs(){
  // deal with loose ROIs
  // focus on the isolated ones first
  float threshold = 1500;
  std::list<SignalROI*> Bad_ROIs;
  for (int i=0;i!=nwire_u;i++){
    for (auto it = rois_u_loose.at(i).begin();it!=rois_u_loose.at(i).end();it++){
      SignalROI* roi = *it;
      if (front_rois.find(roi)==front_rois.end() && back_rois.find(roi)==back_rois.end()){
	if (roi->get_above_threshold(threshold).size()==0)
	  Bad_ROIs.push_back(roi);
      }
    }
  }
  for (auto it = Bad_ROIs.begin(); it!=Bad_ROIs.end(); it ++){
    SignalROI* roi = *it;
    int chid = roi->get_chid();
    auto it1 = find(rois_u_loose.at(chid).begin(), rois_u_loose.at(chid).end(),roi);
    if (it1 != rois_u_loose.at(chid).end())
      rois_u_loose.at(chid).erase(it1);
    delete roi;
  }
  Bad_ROIs.clear();
  for (int i=0;i!=nwire_v;i++){
    for (auto it = rois_v_loose.at(i).begin();it!=rois_v_loose.at(i).end();it++){
      SignalROI* roi = *it;
      if (front_rois.find(roi)==front_rois.end() && back_rois.find(roi)==back_rois.end()){
	if (roi->get_above_threshold(threshold).size()==0)
	  Bad_ROIs.push_back(roi);
      }
    }
  }
  for (auto it = Bad_ROIs.begin(); it!=Bad_ROIs.end(); it ++){
    SignalROI* roi = *it;
    int chid = roi->get_chid()-nwire_u;
    auto it1 = find(rois_v_loose.at(chid).begin(), rois_v_loose.at(chid).end(),roi);
    if (it1 != rois_v_loose.at(chid).end())
      rois_v_loose.at(chid).erase(it1);
    delete roi;
  }


  threshold = 1200;
  std::set<SignalROI*> Good_ROIs;
  for (int i=0;i!=nwire_u;i++){
    for (auto it = rois_u_loose.at(i).begin();it!=rois_u_loose.at(i).end();it++){
      SignalROI* roi = *it;
      if (roi->get_above_threshold(threshold).size()!=0)
	Good_ROIs.insert(roi);
    }
  }
  Bad_ROIs.clear();
  for (int i=0;i!=nwire_u;i++){
    for (auto it = rois_u_loose.at(i).begin();it!=rois_u_loose.at(i).end();it++){
      SignalROI* roi = *it;
      
      if (Good_ROIs.find(roi)!=Good_ROIs.end()) continue;
      if (front_rois.find(roi)!=front_rois.end()){
	SignalROISelection next_rois = front_rois[roi];
	int flag_qx = 0;
	for (int i=0;i!=next_rois.size();i++){
	  SignalROI* roi1 = next_rois.at(i);
	  if (Good_ROIs.find(roi1)!=Good_ROIs.end()) {
	    flag_qx = 1;
	    continue;
	  }
	}
	if (flag_qx == 1) continue;
      }
      
      if (back_rois.find(roi)!=back_rois.end()){
	SignalROISelection next_rois = back_rois[roi];
	int flag_qx = 0;
	for (int i=0;i!=next_rois.size();i++){
	  SignalROI* roi1 = next_rois.at(i);
	  if (Good_ROIs.find(roi1)!=Good_ROIs.end()) {
	    flag_qx = 1;
	    continue;
	  }
	}
	if (flag_qx == 1) continue;
      }
      
      Bad_ROIs.push_back(roi);
    }
  }
  for (auto it = Bad_ROIs.begin(); it!=Bad_ROIs.end(); it ++){
    SignalROI* roi = *it;
    int chid = roi->get_chid();
    //std::cout << chid << std::endl;
    if (front_rois.find(roi)!=front_rois.end()){
      SignalROISelection next_rois = front_rois[roi];
      for (int i=0;i!=next_rois.size();i++){
   	//unlink the current roi
   	unlink(roi,next_rois.at(i));
      }
      front_rois.erase(roi);
    }
    
    if (back_rois.find(roi)!=back_rois.end()){
      SignalROISelection next_rois = back_rois[roi];
      for (int i=0;i!=next_rois.size();i++){
   	//unlink the current roi
   	unlink(roi,next_rois.at(i));
      }
      back_rois.erase(roi);
    }
    auto it1 = find(rois_u_loose.at(chid).begin(), rois_u_loose.at(chid).end(),roi);
    if (it1 != rois_u_loose.at(chid).end())
      rois_u_loose.at(chid).erase(it1);
    
    delete roi;
  }


  Good_ROIs.clear();
  for (int i=0;i!=nwire_v;i++){
    for (auto it = rois_v_loose.at(i).begin();it!=rois_v_loose.at(i).end();it++){
      SignalROI* roi = *it;
      if (roi->get_above_threshold(threshold).size()!=0)
	Good_ROIs.insert(roi);
    }
  }
  Bad_ROIs.clear();
  for (int i=0;i!=nwire_v;i++){
    for (auto it = rois_v_loose.at(i).begin();it!=rois_v_loose.at(i).end();it++){
      SignalROI* roi = *it;
      
      if (Good_ROIs.find(roi)!=Good_ROIs.end()) continue;
      if (front_rois.find(roi)!=front_rois.end()){
	SignalROISelection next_rois = front_rois[roi];
	int flag_qx = 0;
	for (int i=0;i!=next_rois.size();i++){
	  SignalROI* roi1 = next_rois.at(i);
	  if (Good_ROIs.find(roi1)!=Good_ROIs.end()) {
	    flag_qx = 1;
	    continue;
	  }
	}
	if (flag_qx == 1) continue;
      }
      
      if (back_rois.find(roi)!=back_rois.end()){
	SignalROISelection next_rois = back_rois[roi];
	int flag_qx = 0;
	for (int i=0;i!=next_rois.size();i++){
	  SignalROI* roi1 = next_rois.at(i);
	  if (Good_ROIs.find(roi1)!=Good_ROIs.end()) {
	    flag_qx = 1;
	    continue;
	  }
	}
	if (flag_qx == 1) continue;
      }
      
      Bad_ROIs.push_back(roi);
    }
  }
  for (auto it = Bad_ROIs.begin(); it!=Bad_ROIs.end(); it ++){
    SignalROI* roi = *it;
    int chid = roi->get_chid()-nwire_u;
    //std::cout << chid << std::endl;
    if (front_rois.find(roi)!=front_rois.end()){
      SignalROISelection next_rois = front_rois[roi];
      for (int i=0;i!=next_rois.size();i++){
   	//unlink the current roi
   	unlink(roi,next_rois.at(i));
      }
      front_rois.erase(roi);
    }
    
    if (back_rois.find(roi)!=back_rois.end()){
      SignalROISelection next_rois = back_rois[roi];
      for (int i=0;i!=next_rois.size();i++){
   	//unlink the current roi
   	unlink(roi,next_rois.at(i));
      }
      back_rois.erase(roi);
    }
    auto it1 = find(rois_v_loose.at(chid).begin(), rois_v_loose.at(chid).end(),roi);
    if (it1 != rois_v_loose.at(chid).end())
      rois_v_loose.at(chid).erase(it1);
    
    delete roi;
  }

}


void WireCell2dToy::uBooNEDataAfterROI::CleanUpCollectionROIs(){
  // deal with tight ROIs, 
  // scan with all the tight ROIs to look for peaks above certain threshold, put in a temporary set
  float threshold = 1200; //electrons, about 1/2 of MIP per tick ...
  std::set<SignalROI*> Good_ROIs;
  for (int i=0;i!=nwire_w;i++){
    for (auto it = rois_w_tight.at(i).begin();it!=rois_w_tight.at(i).end();it++){
      SignalROI* roi = *it;
      if (roi->get_above_threshold(threshold).size()!=0)
	Good_ROIs.insert(roi);
    }
  }
  // for a particular ROI if it is not in, or it is not connected with one in the temporary map, then remove it
  std::list<SignalROI*> Bad_ROIs;
  for (int i=0;i!=nwire_w;i++){
    for (auto it = rois_w_tight.at(i).begin();it!=rois_w_tight.at(i).end();it++){
      SignalROI* roi = *it;
      
      if (Good_ROIs.find(roi)!=Good_ROIs.end()) continue;
      if (front_rois.find(roi)!=front_rois.end()){
	SignalROISelection next_rois = front_rois[roi];
	int flag_qx = 0;
	for (int i=0;i!=next_rois.size();i++){
	  SignalROI* roi1 = next_rois.at(i);
	  if (Good_ROIs.find(roi1)!=Good_ROIs.end()) {
	    flag_qx = 1;
	    continue;
	  }
	}
	if (flag_qx == 1) continue;
      }
      
      if (back_rois.find(roi)!=back_rois.end()){
	SignalROISelection next_rois = back_rois[roi];
	int flag_qx = 0;
	for (int i=0;i!=next_rois.size();i++){
	  SignalROI* roi1 = next_rois.at(i);
	  if (Good_ROIs.find(roi1)!=Good_ROIs.end()) {
	    flag_qx = 1;
	    continue;
	  }
	}
	if (flag_qx == 1) continue;
      }
      
      Bad_ROIs.push_back(roi);
    }
  }
  
  // remove the ROI and then update the map
  
  for (auto it = Bad_ROIs.begin(); it!=Bad_ROIs.end(); it ++){
    SignalROI* roi = *it;
    int chid = roi->get_chid()-nwire_u-nwire_v;
    //std::cout << chid << std::endl;
    if (front_rois.find(roi)!=front_rois.end()){
      SignalROISelection next_rois = front_rois[roi];
      for (int i=0;i!=next_rois.size();i++){
   	//unlink the current roi
   	unlink(roi,next_rois.at(i));
      }
      front_rois.erase(roi);
    }
  
    if (back_rois.find(roi)!=back_rois.end()){
      SignalROISelection next_rois = back_rois[roi];
      for (int i=0;i!=next_rois.size();i++){
   	//unlink the current roi
   	unlink(roi,next_rois.at(i));
      }
      back_rois.erase(roi);
    }
    auto it1 = find(rois_w_tight.at(chid).begin(), rois_w_tight.at(chid).end(),roi);
    if (it1 != rois_w_tight.at(chid).end())
      rois_w_tight.at(chid).erase(it1);
    
    delete roi;
  }
  
}


void WireCell2dToy::uBooNEDataAfterROI::BreakROI1(SignalROI* roi){
  int start_bin = roi->get_start_bin();
  int end_bin = roi->get_end_bin();
  
  if (start_bin <0 || end_bin < 0) return;

  TH1F *htemp = new TH1F("htemp","htemp",end_bin-start_bin+1,start_bin,end_bin+1);
  std::vector<float>& contents = roi->get_contents();
  for (Int_t i=0;i!=htemp->GetNbinsX();i++){
    htemp->SetBinContent(i+1,contents.at(i));
  }

   // now create many new ROIs
   std::vector<int> bins;
   for (int i=0;i!=htemp->GetNbinsX();i++){
     if (fabs(htemp->GetBinContent(i+1))<1e-3)
       bins.push_back(i+start_bin);
   }
   int chid = roi->get_chid();
   WirePlaneType_t plane = roi->get_plane();
   SignalROISelection new_rois;

   // if (chid == 1274)
   //   std::cout << "BreakROI1: " << chid << " " << roi->get_start_bin() << " " << roi->get_end_bin() << " " << bins.size()  << " " << htemp->GetBinContent(1) << " " << htemp->GetBinContent(end_bin-start_bin+1) << std::endl;
   
   TH1F *h1 = new TH1F("h1","h1",end_bin+1,0,end_bin+1);
   for (int i=0;i!=bins.size()-1;i++){
     int start_bin1 = bins.at(i);
     int end_bin1 = bins.at(i+1);
     // if (chid == 1274)
     //   std::cout << start_bin1 << " " << end_bin1 << std::endl;
     h1->Reset();
     for (int j=start_bin1;j<=end_bin1;j++){
       h1->SetBinContent(j+1,htemp->GetBinContent(j-start_bin+1));
     }
     if (start_bin1 >=0 && end_bin1 >start_bin1){
       SignalROI *sub_roi = new SignalROI(plane,chid,start_bin1,end_bin1,h1);
       new_rois.push_back(sub_roi);
     }
   }

   // update the list 
   if (chid < nwire_u){
     auto it = find(rois_u_loose.at(chid).begin(),rois_u_loose.at(chid).end(),roi);
     rois_u_loose.at(chid).erase(it);
     for (int i=0;i!=new_rois.size();i++){
       rois_u_loose.at(chid).push_back(new_rois.at(i));
     }
   }else if (chid < nwire_u+nwire_v){
     auto it = find(rois_v_loose.at(chid-nwire_u).begin(),rois_v_loose.at(chid-nwire_u).end(),roi);
     rois_v_loose.at(chid-nwire_u).erase(it);
     for (int i=0;i!=new_rois.size();i++){
       rois_v_loose.at(chid-nwire_u).push_back(new_rois.at(i));
     }
   }else{
     auto it = find(rois_w_loose.at(chid-nwire_u-nwire_v).begin(),rois_w_loose.at(chid-nwire_u-nwire_v).end(),roi);
     rois_w_loose.at(chid-nwire_u-nwire_v).erase(it);
     for (int i=0;i!=new_rois.size();i++){
       rois_w_loose.at(chid-nwire_u-nwire_v).push_back(new_rois.at(i));
     } 
   }

   // update all the maps 
   // update front map
   if (front_rois.find(roi)!=front_rois.end()){
     SignalROISelection next_rois = front_rois[roi];
     for (int i=0;i!=next_rois.size();i++){
       //unlink the current roi
       unlink(roi,next_rois.at(i));
       //loop new rois and link them
       for (int j=0;j!=new_rois.size();j++){
   	 if (new_rois.at(j)->overlap(next_rois.at(i)))
   	   link(new_rois.at(j),next_rois.at(i));
       }
     }
     front_rois.erase(roi);
   }
   // update back map
   if (back_rois.find(roi)!=back_rois.end()){
     SignalROISelection prev_rois = back_rois[roi];
     for (int i=0;i!=prev_rois.size();i++){
       // unlink the current roi
       unlink(prev_rois.at(i),roi);
       // loop new rois and link them
       for (int j=0;j!=new_rois.size();j++){
   	 if (new_rois.at(j)->overlap(prev_rois.at(i)))
   	   link(prev_rois.at(i),new_rois.at(j));
       }
     }
     back_rois.erase(roi);
   }

   // update contained map 
   if (contained_rois.find(roi)!=contained_rois.end()){
     SignalROISelection tight_rois = contained_rois[roi];
     for (int i=0;i!=tight_rois.size();i++){
       for (int j=0;j!=new_rois.size();j++){
       	 if (new_rois.at(j)->overlap(tight_rois.at(i))){
   	   if (contained_rois.find(new_rois.at(j)) == contained_rois.end()){
   	     SignalROISelection temp_rois;
   	     temp_rois.push_back(tight_rois.at(i));
   	     contained_rois[new_rois.at(j)] = temp_rois;
   	   }else{
   	     contained_rois[new_rois.at(j)].push_back(tight_rois.at(i));
   	   }
   	 }
       }
     }
     contained_rois.erase(roi);
   }
   
   // delete the old ROI
   delete roi;
   delete h1;
   delete htemp;

}

void WireCell2dToy::uBooNEDataAfterROI::BreakROI(SignalROI* roi, float rms){
  // main algorithm 
  int start_bin = roi->get_start_bin();
  int end_bin = roi->get_end_bin();
  
  if (start_bin <0 || end_bin <0 ) return;
  // if (roi->get_chid()==1240){
  //   std::cout << "xin: " << start_bin << " " << end_bin << std::endl;
  // }

  TH1F *htemp = new TH1F("htemp","htemp",end_bin-start_bin+1,start_bin,end_bin+1);
  std::vector<float>& contents = roi->get_contents();
  for (Int_t i=0;i!=htemp->GetNbinsX();i++){
    htemp->SetBinContent(i+1,contents.at(i));
  }
  TSpectrum *s = new TSpectrum(200);
  Int_t nfound = s->Search(htemp,2,"nobackground new",0.1);
  float th_peak = 3.0;
  float sep_peak = 6.0;
  
  std::set<int> saved_boundaries;

  if (nfound > 1){
    Int_t npeaks = s->GetNPeaks();
    Double_t *peak_pos = s->GetPositionX();
    Double_t *peak_height = s->GetPositionY();
    int order_peak_pos[205];
    int npeaks_threshold = 0;
    for (Int_t j=0;j!=npeaks;j++){
      order_peak_pos[j] = *(peak_pos+j);
      if (*(peak_height+j)>th_peak*rms){
	npeaks_threshold ++;
      }
    }
    if (npeaks_threshold >1){
      std::sort(order_peak_pos,order_peak_pos + npeaks);
      float valley_pos[205];
      valley_pos[0] = start_bin;

      // find the first real valley
      float min = 1e9;
      for (int k=0; k< order_peak_pos[0]-start_bin;k++){
       	if (htemp->GetBinContent(k+1) < min){
       	  min = htemp->GetBinContent(k+1);
	  valley_pos[0] = k+start_bin;
       	}
      }
      if (valley_pos[0] != start_bin){
	for (int j=npeaks-1;j>=0;j--){
	  order_peak_pos[j+1] = order_peak_pos[j];
	}
	npeaks ++;
	order_peak_pos[0] = start_bin;
	for (int j=start_bin; j!=valley_pos[0];j++){
	  if (htemp->GetBinContent(j-start_bin+1) > htemp->GetBinContent(order_peak_pos[0]-start_bin+1))
	    order_peak_pos[0] = j;
	}
	valley_pos[0] = start_bin;
      }

      for (Int_t j=0;j!=npeaks-1;j++){
	Float_t min = 1e9;
	valley_pos[j+1] = order_peak_pos[j];
	for (Int_t k = order_peak_pos[j]-start_bin; k< order_peak_pos[j+1]-start_bin;k++){
	  if (htemp->GetBinContent(k+1) < min){
	    min = htemp->GetBinContent(k+1);
	    valley_pos[j+1] = k+start_bin;
	  }
	}
      }

      //find the end ... 
      valley_pos[npeaks] = end_bin;
      min = 1e9;
      for (int k=order_peak_pos[npeaks-1]-start_bin; k<= end_bin-start_bin;k++){
      	if (htemp->GetBinContent(k+1) < min){
      	  min = htemp->GetBinContent(k+1);
      	  valley_pos[npeaks] = k+start_bin;
      	}
      }
      if (valley_pos[npeaks]!=end_bin){
	npeaks ++;
	valley_pos[npeaks] = end_bin;
	order_peak_pos[npeaks-1] = end_bin;
	for (int j=valley_pos[npeaks-1];j!=valley_pos[npeaks];j++){
	  if (htemp->GetBinContent(j-start_bin+1) > htemp->GetBinContent(order_peak_pos[npeaks-1] -start_bin+1))
	    order_peak_pos[npeaks-1] = j;
	}
      }
      
      
      // if (roi->get_chid() == 1240 && roi->get_plane() == WirePlaneType_t(0)){
      // 	for (int j=0;j!=npeaks;j++){
      // 	  std::cout << valley_pos[j] << " " << htemp->GetBinContent(valley_pos[j]-start_bin+1)<< " " << order_peak_pos[j] << " " << htemp->GetBinContent( order_peak_pos[j]-start_bin+1) << " " << valley_pos[j+1] << " " << htemp->GetBinContent(valley_pos[j+1] - start_bin+1)<< std::endl ;
      // 	}
      // }
	

      
      // need to organize the peaks and valleys ... 
      float valley_pos1[205];
      float peak_pos1[205];
      int npeaks1 = 0;
      // fill in the first valley;
      valley_pos1[0] = valley_pos[0];
      for (int j=0;j<npeaks;j++){
	if (npeaks1 >0){
	  // find the lowest valley except the first peak, except the first one
	  if (htemp->GetBinContent(valley_pos[j]-start_bin+1)< htemp->GetBinContent(valley_pos1[npeaks1]-start_bin+1)){
	    valley_pos1[npeaks1] = valley_pos[j];
	  }
	}

	// find the next peak
	if ( htemp->GetBinContent(order_peak_pos[j]-start_bin+1) - htemp->GetBinContent(valley_pos1[npeaks1]-start_bin+1) > rms * sep_peak){
	  peak_pos1[npeaks1] = order_peak_pos[j] ;
	  npeaks1 ++;
	  int flag1 = 0;

	  for (int k=j+1;k!=npeaks+1;k++){
	    // find the highest peak before ... 
	    if (k<=npeaks){
	      if (htemp->GetBinContent(order_peak_pos[k-1]-start_bin+1) > htemp->GetBinContent(peak_pos1[npeaks1-1]-start_bin+1))
		peak_pos1[npeaks1-1] = order_peak_pos[k-1];
	    }

	    if (htemp->GetBinContent(peak_pos1[npeaks1-1]-start_bin+1) - htemp->GetBinContent(valley_pos[k]-start_bin+1) > rms * sep_peak){
	      valley_pos1[npeaks1] = valley_pos[k];
	      j = k-1;
	      flag1 = 1;
	      break;
	    }
	    // find the next valley
	  }
	  if (flag1 == 0){
	    valley_pos1[npeaks1] = valley_pos[npeaks];
	    j = npeaks;
	  }

	  // if (roi->get_chid() == 1240 && roi->get_plane() == WirePlaneType_t(0))
	  //   std::cout << "c: " << npeaks << " " << valley_pos1[npeaks1-1] << " " << peak_pos1[npeaks1-1] << " " << valley_pos1[npeaks1] << " " << rms * sep_peak << std::endl;
	}
      }
      // fill the last valley
      valley_pos1[npeaks1] = valley_pos[npeaks];
      //std::cout << roi->get_plane() << " " << roi->get_chid() << " " << npeaks << " " << npeaks1 << " " ;
      //std::cout << start_bin << " " << end_bin << std::endl;


      // if (roi->get_chid() == 1240 && roi->get_plane() == WirePlaneType_t(0)){
      // 	for (int j=0;j!=npeaks1;j++){
      // 	  std::cout << valley_pos1[j] << " " << peak_pos1[j] << " " <<  valley_pos1[j+1] << std::endl ;
      // 	}
      // }


      for (Int_t j=0;j!=npeaks1;j++){
	Int_t start_pos = valley_pos1[j];
	Int_t end_pos = valley_pos1[j+1];
	//std::cout << "b " << start_pos << " " << end_pos << std::endl;
	saved_boundaries.insert(start_pos);
	saved_boundaries.insert(end_pos);
      }
            
      TH1F *htemp1 = (TH1F*)htemp->Clone("htemp1");
      htemp->Reset();
      for (Int_t j=0;j!=npeaks1;j++){
	int flag = 0;
	Int_t start_pos = valley_pos1[j];
	Double_t start_content = htemp1->GetBinContent(valley_pos1[j]-start_bin+1);
	Int_t end_pos = valley_pos1[j+1];
	Double_t end_content = htemp1->GetBinContent(valley_pos1[j+1]-start_bin+1);
	
	//	std::cout << "a " << start_pos << " " << end_pos << std::endl;

	if (saved_boundaries.find(start_pos) != saved_boundaries.end() ||
	    saved_boundaries.find(end_pos) != saved_boundaries.end()){
	  // if (roi->get_chid() == 1240 && roi->get_plane() == WirePlaneType_t(0))
	  //   std::cout << "d: " << start_pos << " " << end_pos << std::endl;

	  for (Int_t k = start_pos; k!=end_pos+1;k++){
	    Double_t temp_content = htemp1->GetBinContent(k-start_bin+1) - (start_content + (end_content-start_content) * (k-start_pos) / (end_pos - start_pos));
	    htemp->SetBinContent(k-start_bin+1,temp_content);
	  }
	}
      }
      delete htemp1;
    }
  }
  
  
  // if (roi->get_chid() == 1308)
  //   std::cout << "Break:  "  << roi->get_chid() << " " << start_bin << " " << end_bin << std::endl;
  
  for (int qx = 0; qx!=2; qx++){
    

    // Now we should go through the system again and re-adjust the content
    std::vector<std::pair<int,int>> bins;
    for (int i=0;i<htemp->GetNbinsX();i++){
      if (htemp->GetBinContent(i+1) < 3*rms){
	int start = i;
	int end = i;
	for (int j=i+1;j<htemp->GetNbinsX();j++){
	  if (htemp->GetBinContent(j+1) < 3*rms){
	    end = j;
	  }else{
	    break;
	  }
	}
	bins.push_back(std::make_pair(start,end));
	// if (roi->get_chid() == 1308)
	//   std::cout << qx << " " <<  start+start_bin << " " << end + start_bin << " " << 3*rms << " " << htemp->GetBinContent(6529-start_bin+1) << std::endl;
	i = end;
      }
    }


    std::vector<int> saved_b;
    for (int i=0;i!=bins.size();i++){
      int start = bins.at(i).first;
      int end = bins.at(i).second;
      // find minimum or zero
      float min = 1e9;
      int bin_min = start;
      for (int j=start;j<=end;j++){
	if (fabs(htemp->GetBinContent(j+1))< 1e-3){
	  bin_min = j;
	  break;
	}else{
	  if (htemp->GetBinContent(j+1) < min){
	    min = htemp->GetBinContent(j+1);
	    bin_min = j;
	  }
	}
      }
      
      // if (roi->get_chid() == 1308)
      // 	std::cout << bin_min+start_bin << std::endl;
      
      saved_b.push_back(bin_min);
    }
    
    if (saved_b.size() >=0){
      TH1F *htemp1 = (TH1F*)htemp->Clone("htemp1");
      htemp->Reset();
      // std::cout << saved_b.size() << " " << bins.size() << " " << htemp->GetNbinsX() << std::endl;
      for (Int_t j=0;j<saved_b.size()-1;j++){
	int flag = 0;
	int start_pos = saved_b.at(j);
	float start_content = htemp1->GetBinContent(start_pos+1);
	int end_pos = saved_b.at(j+1);
	float end_content = htemp1->GetBinContent(end_pos+1);
	
	for (Int_t k = start_pos; k!=end_pos+1;k++){
	  Double_t temp_content = htemp1->GetBinContent(k+1) - (start_content + (end_content-start_content) * (k-start_pos) / (end_pos - start_pos));
	  htemp->SetBinContent(k+1,temp_content);
	}
      }
      delete htemp1;
    }
  }

  // get back to the original content 
  for (Int_t i=0;i!=htemp->GetNbinsX();i++){
    contents.at(i) = htemp->GetBinContent(i+1);
  }
  
  



   delete s;
   delete htemp;
   

}

void WireCell2dToy::uBooNEDataAfterROI::unlink(SignalROI* prev_roi, SignalROI* next_roi){
  if (front_rois.find(prev_roi)!=front_rois.end()){
    SignalROISelection& temp_rois = front_rois[prev_roi];
    auto it = find(temp_rois.begin(),temp_rois.end(),next_roi);
    if (it != temp_rois.end())
      temp_rois.erase(it);
  }
  if (back_rois.find(next_roi)!=back_rois.end()){
    SignalROISelection& temp_rois = back_rois[next_roi];
    auto it = find(temp_rois.begin(),temp_rois.end(),prev_roi);
    if (it != temp_rois.end())
      temp_rois.erase(it);
  }
}

void WireCell2dToy::uBooNEDataAfterROI::link(SignalROI* prev_roi, SignalROI* next_roi){
  if (front_rois.find(prev_roi)!=front_rois.end()){
    SignalROISelection& temp_rois = front_rois[prev_roi];
    auto it = find(temp_rois.begin(),temp_rois.end(),next_roi);
    if (it == temp_rois.end())
      temp_rois.push_back(next_roi);
  }else{
    SignalROISelection temp_rois;
    temp_rois.push_back(next_roi);
    front_rois[prev_roi] = temp_rois;
  }

  if (back_rois.find(next_roi)!=back_rois.end()){
    SignalROISelection& temp_rois = back_rois[next_roi];
    auto it = find(temp_rois.begin(),temp_rois.end(),prev_roi);
    if (it == temp_rois.end())
      temp_rois.push_back(prev_roi);
  }else{
    SignalROISelection temp_rois;
    temp_rois.push_back(prev_roi);
    back_rois[next_roi] = temp_rois;
  }
}


void WireCell2dToy::uBooNEDataAfterROI::BreakROIs(){
  // get RMS value, and put in 
  std::vector<float>& rms_u = rois.get_uplane_rms();
  std::vector<float>& rms_v = rois.get_vplane_rms();

  SignalROISelection all_rois;

  for (int i=0;i!=rois_u_loose.size();i++){
    for (auto it = rois_u_loose.at(i).begin(); it!= rois_u_loose.at(i).end(); it++){
      BreakROI(*it,rms_u.at(i));
      all_rois.push_back(*it);
    }
  }
  
  for (int i=0;i!=rois_v_loose.size();i++){
    for (auto it = rois_v_loose.at(i).begin(); it!= rois_v_loose.at(i).end(); it++){
      BreakROI(*it,rms_v.at(i));
      all_rois.push_back(*it);
    }
  }



  for (int i=0;i!=all_rois.size();i++){
    BreakROI1(all_rois.at(i));
  }

  // int num_tight[3]={0,0,0};
  // int num_loose[3]={0,0,0};
  // for (int i=0;i!=nwire_u;i++){
  //   num_tight[0] += rois_u_tight.at(i).size();
  //   num_loose[0] += rois_u_loose.at(i).size();
  // }
  // for (int i=0;i!=nwire_v;i++){
  //   num_tight[1] += rois_v_tight.at(i).size();
  //   num_loose[1] += rois_v_loose.at(i).size();
  // }
  // for (int i=0;i!=nwire_w;i++){
  //   num_tight[2] += rois_w_tight.at(i).size();
  //   num_loose[2] += rois_w_loose.at(i).size();
  // }

  // std::cout << num_tight[0] << " " << num_tight[1] << " " << num_tight[2] << " " << num_loose[0] << " " << num_loose[1] << " " << num_loose[2] << " " << front_rois.size() << " " << back_rois.size() << " " << contained_rois.size() << std::endl;
}

void WireCell2dToy::uBooNEDataAfterROI::ShrinkROI(SignalROI *roi){
  
  // get tight ROI as a inner boundary
  // get the nearby ROIs with threshold as some sort of boundary 
  int start_bin = roi->get_start_bin();
  int end_bin = roi->get_end_bin();

  if (start_bin <0 || end_bin <0) return;
  
  int chid = roi->get_chid();
  WirePlaneType_t plane = roi->get_plane();
  std::vector<float>& contents = roi->get_contents();
  
  float threshold1;
  if (chid < nwire_u){
    threshold1 = rois.get_uplane_rms().at(chid) * 3.0;
  }else if (chid<nwire_u+nwire_v){
    threshold1 = rois.get_vplane_rms().at(chid-nwire_u)*3.0;
  }
  
  int channel_save = 1240;
  int print_flag = 0;

  // std::cout << "check tight ROIs " << std::endl;
  // use to save contents 
  TH1F *htemp = new TH1F("htemp","htemp",end_bin-start_bin+1,start_bin,end_bin+1);
  // check tight ROIs
  if (contained_rois.find(roi)!=contained_rois.end()){
    for (auto it = contained_rois[roi].begin();it!=contained_rois[roi].end();it++){
      SignalROI *tight_roi = *it;
      int start_bin1 = tight_roi->get_start_bin();
      int end_bin1 = tight_roi->get_end_bin();
      
      if (chid == channel_save && print_flag)
	std::cout << "Tight "  " " << start_bin1 << " " << end_bin1 << std::endl;

      for (int i=start_bin1;i<=end_bin1;i++){
	if (i-start_bin >=0 && i-start_bin <=htemp->GetNbinsX()){
	  htemp->SetBinContent(i-start_bin+1,1);
	}
      }
    }
  }

  // std::cout << "check front ROIs " << std::endl;

  //check front ROIs
  if (front_rois.find(roi)!=front_rois.end()){
    for (auto it=front_rois[roi].begin();it!=front_rois[roi].end();it++){
      SignalROI *next_roi = *it;
      int start_bin1 = next_roi->get_start_bin();
      int chid1 = next_roi->get_chid();
      WirePlaneType_t plane1 = next_roi->get_plane();
      float threshold;
      if (chid1 < nwire_u){
	threshold = rois.get_uplane_rms().at(chid1) * 3.0;
      }else if (chid1<nwire_u+nwire_v){
	threshold = rois.get_vplane_rms().at(chid1-nwire_u)*3.0;
      }
      std::vector<std::pair<int,int>> contents_above_threshold = next_roi->get_above_threshold(threshold);
      for (int i=0;i!=contents_above_threshold.size();i++){
	if (chid == channel_save && print_flag)
	  std::cout << "Front " << chid1 << " " << start_bin1 + contents_above_threshold.at(i).first << " " << start_bin1 + contents_above_threshold.at(i).second << std::endl;

	for (int j=contents_above_threshold.at(i).first;j<=contents_above_threshold.at(i).second;j++){
	  if (j+start_bin1-start_bin >=0 && j+start_bin1-start_bin <htemp->GetNbinsX()){
	    if (contents.at(j+start_bin1-start_bin) > threshold1) htemp->SetBinContent(j+start_bin1-start_bin+1,1);
	  }
	}
      }
    }
  }
  
  //std::cout << "check back ROIs " << std::endl;

  //check back ROIs
  if (back_rois.find(roi)!=back_rois.end()){
    for (auto it=back_rois[roi].begin();it!=back_rois[roi].end();it++){
      SignalROI *prev_roi = *it;
      int start_bin1 = prev_roi->get_start_bin();
      int chid1 = prev_roi->get_chid();
      WirePlaneType_t plane1 = prev_roi->get_plane();
      float threshold;
      if (chid1 < nwire_u){
	threshold = rois.get_uplane_rms().at(chid1) * 3.0;
      }else if (chid1<nwire_u+nwire_v){
	threshold = rois.get_vplane_rms().at(chid1-nwire_u)*3.0;
      }
      std::vector<std::pair<int,int>> contents_above_threshold = prev_roi->get_above_threshold(threshold);
      for (int i=0;i!=contents_above_threshold.size();i++){
	 if (chid == channel_save && print_flag)
	   std::cout << "Back " << chid1 << " " << start_bin1 + contents_above_threshold.at(i).first << " " << start_bin1 + contents_above_threshold.at(i).second << std::endl;

	for (int j=contents_above_threshold.at(i).first;j<=contents_above_threshold.at(i).second;j++){
	  if (j+start_bin1-start_bin >=0 && j+start_bin1-start_bin <htemp->GetNbinsX()){
	    if (contents.at(j+start_bin1-start_bin) > threshold1)  htemp->SetBinContent(j+start_bin1-start_bin+1,1);
	  }
	}
      }
    }
  }

  //std::cout << "check contents " << std::endl;

  // // consider the 1/2 of the peak as threshold;
  // float max = 0;
  // for (int i=0;i!=contents.size();i++){
  //   if (contents.at(i) > max)
  //     max = contents.at(i);
  // }
  // for (int i=0;i!=contents.size();i++){
  //   // if (contents.at(i) > max/2. && contents.at(i) > threshold1*2 ) htemp->SetBinContent(i+1,1);
  // }
  
  // get the first bin, and last bin, add pad
  int pad = 5;
  int new_start_bin=start_bin;
  int new_end_bin=end_bin;
  for (int i=0;i!=htemp->GetNbinsX();i++){
    if (htemp->GetBinContent(i+1) >0){
      new_start_bin = i + start_bin;
      break;
    }
  }
  for (int i = htemp->GetNbinsX()-1; i>=0;i--){
    if (htemp->GetBinContent(i+1)>0){
      new_end_bin = i + start_bin;
      break;
    }
  }
  new_start_bin -= pad;
  new_end_bin += pad;
  if (new_start_bin < start_bin) new_start_bin = start_bin;
  if (new_end_bin > end_bin) new_end_bin = end_bin;
  
  if (chid == channel_save && print_flag)
    std::cout << "check contents " << " " << start_bin << " " << end_bin << " " << new_start_bin << " " << new_end_bin << std::endl;
  

  // create a new ROI
  TH1F *h1 = new TH1F("h1","h1",end_bin+1,0,end_bin+1);
  for (int i=new_start_bin; i<=new_end_bin;i++){
    h1->SetBinContent(i+1,contents.at(i-start_bin));
  }
  
  SignalROISelection new_rois;
  if (new_start_bin >=0 && new_end_bin > new_start_bin){
    SignalROI *new_roi = new SignalROI(plane,chid,new_start_bin,new_end_bin,h1);
    new_rois.push_back(new_roi);
  }

  // std::cout << "update maps " << std::endl;
  
   // update the list 
   if (chid < nwire_u){
     auto it = find(rois_u_loose.at(chid).begin(),rois_u_loose.at(chid).end(),roi);
     rois_u_loose.at(chid).erase(it);
     for (int i=0;i!=new_rois.size();i++){
       rois_u_loose.at(chid).push_back(new_rois.at(i));
     }
   }else if (chid < nwire_u+nwire_v){
     auto it = find(rois_v_loose.at(chid-nwire_u).begin(),rois_v_loose.at(chid-nwire_u).end(),roi);
     rois_v_loose.at(chid-nwire_u).erase(it);
     for (int i=0;i!=new_rois.size();i++){
       rois_v_loose.at(chid-nwire_u).push_back(new_rois.at(i));
     }
   }else{
     auto it = find(rois_w_loose.at(chid-nwire_u-nwire_v).begin(),rois_w_loose.at(chid-nwire_u-nwire_v).end(),roi);
     rois_w_loose.at(chid-nwire_u-nwire_v).erase(it);
     for (int i=0;i!=new_rois.size();i++){
       rois_w_loose.at(chid-nwire_u-nwire_v).push_back(new_rois.at(i));
     } 
   }

   // update all the maps 
   // update front map
   if (front_rois.find(roi)!=front_rois.end()){
     SignalROISelection next_rois = front_rois[roi];
     for (int i=0;i!=next_rois.size();i++){
       //unlink the current roi
       unlink(roi,next_rois.at(i));
       //loop new rois and link them
       for (int j=0;j!=new_rois.size();j++){
	 if (new_rois.at(j)->overlap(next_rois.at(i)))
	   link(new_rois.at(j),next_rois.at(i));
       }
     }
     front_rois.erase(roi);
   }
   // update back map
   if (back_rois.find(roi)!=back_rois.end()){
     SignalROISelection prev_rois = back_rois[roi];
     for (int i=0;i!=prev_rois.size();i++){
       // unlink the current roi
       unlink(prev_rois.at(i),roi);
       // loop new rois and link them
       for (int j=0;j!=new_rois.size();j++){
	 if (new_rois.at(j)->overlap(prev_rois.at(i)))
	   link(prev_rois.at(i),new_rois.at(j));
       }
     }
     back_rois.erase(roi);
   }

   // update contained map 
   if (contained_rois.find(roi)!=contained_rois.end()){
     SignalROISelection tight_rois = contained_rois[roi];
     for (int i=0;i!=tight_rois.size();i++){
       for (int j=0;j!=new_rois.size();j++){
	 if (new_rois.at(j)->overlap(tight_rois.at(i))){
	   if (contained_rois.find(new_rois.at(j)) == contained_rois.end()){
	     SignalROISelection temp_rois;
	     temp_rois.push_back(tight_rois.at(i));
	     contained_rois[new_rois.at(j)] = temp_rois;
	   }else{
	     contained_rois[new_rois.at(j)].push_back(tight_rois.at(i));
	   }
	 }
       }
     }
     contained_rois.erase(roi);
   }
   
   // delete the old ROI
   delete roi;

  delete htemp;
  delete h1;
  
}

void WireCell2dToy::uBooNEDataAfterROI::ShrinkROIs(){
  // collect all ROIs
  SignalROISelection all_rois;
  for (int i=0;i!=rois_u_loose.size();i++){
    for (auto it = rois_u_loose.at(i).begin(); it!= rois_u_loose.at(i).end(); it++){
      all_rois.push_back(*it);
    }
  }  
  for (int i=0;i!=rois_v_loose.size();i++){
    for (auto it = rois_v_loose.at(i).begin(); it!= rois_v_loose.at(i).end(); it++){
      all_rois.push_back(*it);
    }
  }
  for (int i=0;i!=all_rois.size();i++){
    ShrinkROI(all_rois.at(i));
  }
}

void WireCell2dToy::uBooNEDataAfterROI::generate_merge_ROIs(){
  // find tight ROIs not contained by the loose ROIs
  for (int i = 0;i!=nwire_u;i++){
    std::map<SignalROI*,int> covered_tight_rois;
    for (auto it = rois_u_loose.at(i).begin();it!=rois_u_loose.at(i).end();it++){
      SignalROI *roi = *it;
      if (contained_rois.find(roi) != contained_rois.end()){
	for (auto it1 = contained_rois[roi].begin(); it1!= contained_rois[roi].end(); it1++){
	  if (covered_tight_rois.find(*it1)==covered_tight_rois.end()){
	    covered_tight_rois[*it1]  =1;
	  }
	}
      }
    }
    SignalROISelection saved_rois;
    for (auto it = rois_u_tight.at(i).begin();it!=rois_u_tight.at(i).end();it++){
      SignalROI *roi = *it;
      if (covered_tight_rois.find(roi) == covered_tight_rois.end()){
	saved_rois.push_back(roi);
      }
    }
    // if (i == 1212)
    //   std::cout << saved_rois.size() << " " << saved_rois.at(0)->get_start_bin() << " " << saved_rois.at(0)->get_end_bin() << std::endl;
    
    for (auto it = saved_rois.begin(); it!=saved_rois.end();it++){
      SignalROI *roi = *it;
      // Duplicate them 
      SignalROI *loose_roi = new SignalROI(roi);
      
      rois_u_loose.at(i).push_back(loose_roi);

      // update all the maps     
      // contained
      SignalROISelection temp_rois;
      temp_rois.push_back(roi);
      contained_rois[loose_roi] = temp_rois;
      // front map loose ROI
      if (i < nwire_u-1){
	for (auto it1 = rois_u_loose.at(i+1).begin(); it1!=rois_u_loose.at(i+1).end(); it1++){
	  SignalROI *next_roi = *it1;
	  
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
      // back map loose ROI
      if (i > 0){
	for (auto it1 = rois_u_loose.at(i-1).begin(); it1!=rois_u_loose.at(i-1).end(); it1++){
	  SignalROI *prev_roi = *it1;
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
    }
  }
 




  for (int i = 0;i!=nwire_v;i++){
    std::map<SignalROI*,int> covered_tight_rois;
    for (auto it = rois_v_loose.at(i).begin();it!=rois_v_loose.at(i).end();it++){
      SignalROI *roi = *it;
      if (contained_rois.find(roi) != contained_rois.end()){
	for (auto it1 = contained_rois[roi].begin(); it1!= contained_rois[roi].end(); it1++){
	  if (covered_tight_rois.find(*it1)==covered_tight_rois.end()){
	    covered_tight_rois[*it1]  =1;
	  }
	}
      }
    }
    SignalROISelection saved_rois;
    for (auto it = rois_v_tight.at(i).begin();it!=rois_v_tight.at(i).end();it++){
      SignalROI *roi = *it;
      if (covered_tight_rois.find(roi) == covered_tight_rois.end()){
	saved_rois.push_back(roi);
      }
    }
    //if (i == 3885-2400)
    //  std::cout << saved_rois.size() << std::endl;
    //   std::cout << saved_rois.size() << " " << saved_rois.at(0)->get_start_bin() << " " << saved_rois.at(0)->get_end_bin() << std::endl;
    
    for (auto it = saved_rois.begin(); it!=saved_rois.end();it++){
      SignalROI *roi = *it;
      // Duplicate them 
      SignalROI *loose_roi = new SignalROI(roi);
      
      rois_v_loose.at(i).push_back(loose_roi);

      // update all the maps     
      // contained
      SignalROISelection temp_rois;
      temp_rois.push_back(roi);
      contained_rois[loose_roi] = temp_rois;
      // front map loose ROI
      if (i < nwire_v-1){
	for (auto it1 = rois_v_loose.at(i+1).begin(); it1!=rois_v_loose.at(i+1).end(); it1++){
	  SignalROI *next_roi = *it1;
	  
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
      // back map loose ROI
      if (i > 0){
	for (auto it1 = rois_v_loose.at(i-1).begin(); it1!=rois_v_loose.at(i-1).end(); it1++){
	  SignalROI *prev_roi = *it1;
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
    }
  }
  
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
      float threshold;
      if (chid < nwire_u){
	tight_roi = new SignalROI(WirePlaneType_t(0), chid,uboone_rois.at(i).first,uboone_rois.at(i).second, htemp_signal);
	threshold = rois.get_uplane_rms().at(chid) * 3.0;
      }else if (chid < nwire_u+nwire_v){
	tight_roi = new SignalROI(WirePlaneType_t(1), chid,uboone_rois.at(i).first,uboone_rois.at(i).second, htemp_signal);
	threshold = rois.get_vplane_rms().at(chid-nwire_u) * 3.0;
      }else{
	tight_roi = new SignalROI(WirePlaneType_t(2), chid,uboone_rois.at(i).first,uboone_rois.at(i).second, htemp_signal);
	threshold = rois.get_wplane_rms().at(chid-nwire_u-nwire_v) * 3.0;
      }
      // std::cout << i << std::endl;
      
      
      // judge ...
      if (tight_roi->get_above_threshold(threshold).size()==0) {
	delete tight_roi;
	continue;
      }

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
      float threshold;
      if (chid < nwire_u){
	loose_roi = new SignalROI(WirePlaneType_t(0),chid,uboone_rois.at(i).first,uboone_rois.at(i).second,htemp_signal);
	threshold = rois.get_uplane_rms().at(chid) * 3.0;
      }else if (chid < nwire_u + nwire_v){
	loose_roi = new SignalROI(WirePlaneType_t(1),chid,uboone_rois.at(i).first,uboone_rois.at(i).second,htemp_signal);
	threshold = rois.get_vplane_rms().at(chid - nwire_u) * 3.0;
      }else{
	loose_roi = new SignalROI(WirePlaneType_t(2),chid,uboone_rois.at(i).first,uboone_rois.at(i).second,htemp_signal);
	threshold = rois.get_wplane_rms().at(chid-nwire_u-nwire_v) * 3.0;
      }

      if (loose_roi->get_above_threshold(threshold).size()==0) {
       	delete loose_roi;
       	continue;
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


  
  //std::cout << rois_u_tight.size() << " " << rois_v_tight.size() << " " << rois_w_tight.size() << " " << rois_u_loose.size() << " " << rois_v_loose.size() << " " << rois_w_loose.size() << " " << front_rois.size() << " " << back_rois.size() << " " << contained_rois.size() << std::endl;
  std::cout << "Clean up loose ROIs" << std::endl;
  CleanUpROIs();
  std::cout << "Generate more loose ROIs from isolated good tight ROIs" << std::endl;
  generate_merge_ROIs();

  for (int qx = 0; qx!=2; qx++){
    std::cout << "Break loose ROIs" << std::endl;
    BreakROIs();
    std::cout << "Clean up ROIs 2nd time" << std::endl;
    CheckROIs();
    CleanUpROIs();
  }
  
  std::cout << "Shrink ROIs" << std::endl;
  ShrinkROIs();
  std::cout << "Clean up ROIs 3rd time" << std::endl;
  CheckROIs();
  CleanUpROIs();


  // Further reduce fake hits
  std::cout << "Remove fake hits " << std::endl;
  CleanUpCollectionROIs();
  CleanUpInductionROIs();

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
	//for (auto it = rois_u_tight.at(chid).begin(); it!= rois_u_tight.at(chid).end();it++){
	SignalROI *roi =  *it;
	std::vector<float>& contents = roi->get_contents();
	int start_bin = roi->get_start_bin();
	for (int i=0;i!=contents.size();i++){
	  htemp_signal->SetBinContent(start_bin+1+i,contents.at(i));
	}
      }
    }else if (chid < nwire_u + nwire_v){
      for (auto it = rois_v_loose.at(chid-nwire_u).begin(); it!= rois_v_loose.at(chid-nwire_u).end();it++){
      //for (auto it = rois_v_tight.at(chid-nwire_u).begin(); it!= rois_v_tight.at(chid-nwire_u).end();it++){
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


void WireCell2dToy::uBooNEDataAfterROI::CheckROIs(){
  std::vector<float>& rms_u = rois.get_uplane_rms();
  std::vector<float>& rms_v = rois.get_vplane_rms();
  
  // for (int i=0;i!=rois_u_loose.size();i++){
  //   for (auto it = rois_u_loose.at(i).begin(); it!= rois_u_loose.at(i).end();it++){
  //     SignalROI *roi = *it;
  //     int chid = roi->get_chid();
  //     if (chid != i) 
  // 	std::cout << roi << std::endl;
      
  //     if (front_rois.find(roi)!=front_rois.end()){
  //  	for (auto it1 = front_rois[roi].begin();it1!=front_rois[roi].end();it1++){
  // 	  SignalROI *roi1 = *it1;
  // 	  int chid1 = roi1->get_chid();
  // 	  if (chid1!=i+1)
  // 	    std::cout << roi1 << " " << chid << " " << chid1 << std::endl;
  // 	}
  //     }
  //   }
  // }
  
  //loop over u loose
  for (int i=0;i!=rois_u_loose.size();i++){
    for (auto it = rois_u_loose.at(i).begin(); it!= rois_u_loose.at(i).end();it++){
      SignalROI *roi = *it;
      int chid = roi->get_chid();
      float th;
      th = 3*rms_u.at(chid);

      if (front_rois.find(roi)!=front_rois.end()){
	SignalROISelection temp_rois;
  	for (auto it1 = front_rois[roi].begin();it1!=front_rois[roi].end();it1++){
  	  SignalROI *roi1 = *it1;
  	  int chid1 = roi1->get_chid();
  	  //std::cout << "F " << i << " " << rois_u_loose.size() << " " << roi1 << " " << chid << " " << chid1 << std::endl;
  	  float th1;
  	  th1 = 3*rms_u.at(chid1);
  	  if (roi->overlap(roi1,th,th1)){
  	  }else{
	    temp_rois.push_back(roi1);
  	    //unlink(roi,roi1);
  	  }
  	}
	for (auto it2 = temp_rois.begin(); it2!= temp_rois.end();it2++){
	  unlink(roi,*it2);
	}
      }

      if (back_rois.find(roi)!=back_rois.end()){
	SignalROISelection temp_rois;
  	for (auto it1 = back_rois[roi].begin();it1!=back_rois[roi].end();it1++){
  	  SignalROI *roi1 = *it1;
  	  int chid1 = roi1->get_chid();
  	  //std::cout << "B " << roi1 << " " << chid << " " << chid1 << std::endl;
  	  float th1;
  	  th1 = 3*rms_u.at(chid1);
  	  if (roi->overlap(roi1,th,th1)){
  	  }else{
	    temp_rois.push_back(roi1);
  	    //unlink(roi,roi1);
  	  }
  	}
	for (auto it2 = temp_rois.begin(); it2!= temp_rois.end();it2++){
	  unlink(roi,*it2);
	}
      }

    }
  }


  //loop over v loose
  for (int i=0;i!=rois_v_loose.size();i++){
    for (auto it = rois_v_loose.at(i).begin(); it!= rois_v_loose.at(i).end();it++){
      SignalROI *roi = *it;
      int chid = roi->get_chid()-nwire_u;
      float th;
      th = 3*rms_v.at(chid);
      if (front_rois.find(roi)!=front_rois.end()){
	SignalROISelection temp_rois;
  	for (auto it1 = front_rois[roi].begin();it1!=front_rois[roi].end();it1++){
  	  SignalROI *roi1 = *it1;
  	  int chid1 = roi1->get_chid()-nwire_u;
  	  float th1;
  	  th1 = 3*rms_v.at(chid1);
  	  if (roi->overlap(roi1,th,th1)){
  	  }else{
	    temp_rois.push_back(roi1);
	    // unlink(roi,roi1);
  	  }
  	}
	for (auto it2 = temp_rois.begin(); it2!= temp_rois.end();it2++){
	  unlink(roi,*it2);
	}
      }
    
      if (back_rois.find(roi)!=back_rois.end()){
	SignalROISelection temp_rois;
  	for (auto it1 = back_rois[roi].begin();it1!=back_rois[roi].end();it1++){
  	  SignalROI *roi1 = *it1;
  	  int chid1 = roi1->get_chid()-nwire_u;
  	  float th1;
  	  th1 = 3*rms_v.at(chid1);
  	  if (roi->overlap(roi1,th,th1)){
  	  }else{
	    temp_rois.push_back(roi1);
	    // unlink(roi,roi1);
  	  }
  	}
	for (auto it2 = temp_rois.begin(); it2!= temp_rois.end();it2++){
	  unlink(roi,*it2);
	}
      }
    }
  }

}
