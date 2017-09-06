#include "WireCell2dToy/uBooNE_L1SP.h"

#include "WireCellRess/LassoModel.h"
#include "WireCellRess/ElasticNetModel.h"
#include <Eigen/Dense>

using namespace Eigen;
using namespace WireCell;

WireCell2dToy::uBooNE_L1SP::uBooNE_L1SP(TH2F *hv_raw, TH2F *hv_decon, int nrebin)
  : hv_raw(hv_raw)
  , hv_decon(hv_decon)
  , nrebin(nrebin)
{
}

WireCell2dToy::uBooNE_L1SP::~uBooNE_L1SP(){ 
}

void WireCell2dToy::uBooNE_L1SP::AddWires(int time_slice, GeomWireSelection& wires){
  if (wires.size() >0){
    for (auto it = wires.begin(); it!=wires.end(); it++){
      int wire_index = (*it)->index();
      if (init_map.find(wire_index)!=init_map.end()){
	init_map[wire_index].push_back(time_slice);
      }else{
	std::vector<int> times;
	times.push_back(time_slice);
	init_map[wire_index] = times;
      }
    }
  }
}

void WireCell2dToy::uBooNE_L1SP::Form_rois(int pad){
  
  for (auto it = init_map.begin(); it!=init_map.end(); it++){
    int wire_index = it->first;
    std::vector<int> time_slices = it->second;
    std::sort(time_slices.begin(), time_slices.end());

    std::vector<std::pair<int,int>> rois;
    std::vector<std::pair<int,int>> rois_save;
    
    rois.push_back(std::make_pair(time_slices.front(),time_slices.front()));
    for (size_t i=1; i<time_slices.size();i++){
      if (time_slices.at(i) - rois.back().second <= pad*2){
	rois.back().second = time_slices.at(i);
      }else{
	rois.push_back(std::make_pair(time_slices.at(i),time_slices.at(i)));
      }
    }
    
    // extend the rois to both side according to the bin content
    for (auto it = rois.begin(); it!= rois.end();  it++){
      int start_bin = it->first;
      int end_bin = it->second;
      //std::cout << start_bin << " " << end_bin << " " ;
      while(hv_decon->GetBinContent(wire_index+1,start_bin)>0){
	start_bin --;
	if (start_bin<=0) break;
      }
      while(hv_decon->GetBinContent(wire_index+1,end_bin)>0){
	end_bin ++;
	if (end_bin >= hv_decon->GetNbinsX()-1) break;
      }

      // add padding ...
      start_bin = start_bin - pad;
      while(hv_decon->GetBinContent(wire_index+1,start_bin)>0){
	start_bin -=pad;
	if (start_bin<=0) break;
      }
      end_bin = end_bin + pad;
      while(hv_decon->GetBinContent(wire_index+1,end_bin)>0){
	end_bin +=pad;
	if (end_bin >= hv_decon->GetNbinsX()-1) break;
      }
      
      if (start_bin <0) start_bin = 0;
      if (end_bin>hv_decon->GetNbinsX()-1) end_bin = hv_decon->GetNbinsX()-1;
      it->first = start_bin;
      it->second = end_bin;
      //std::cout << start_bin << " " << end_bin << std::endl;
    }
    
    std::cout << wire_index << " " << rois.size() << " ";    
    // merge them ...
    for (auto it = rois.begin(); it!= rois.end();  it++){
      if (rois_save.size()==0){
	rois_save.push_back(*it);
      }else if (it->first <= rois_save.back().second){
	rois_save.back().second = it->second;
      }else{
	rois_save.push_back(*it);
      }
    }
    
    // std::cout << wire_index << " " << rois_save.size() << std::endl;
    // std::cout << wire_index << " " << time_slices.size() << " " << time_slices.front() << " " << time_slices.back() << " " << hv_decon->GetBinContent(wire_index+1,time_slices.front()+1)<< std::endl;
  }
}
