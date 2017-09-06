#include "WireCell2dToy/uBooNE_L1SP.h"

#include "WireCellRess/LassoModel.h"
#include "WireCellRess/ElasticNetModel.h"
#include <Eigen/Dense>

using namespace Eigen;
using namespace WireCell;

WireCell2dToy::uBooNE_L1SP::uBooNE_L1SP(TH2F *hv_raw)
  : hv_raw(hv_raw)
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

void WireCell2dToy::uBooNE_L1SP::Form_rois(){
  
}
