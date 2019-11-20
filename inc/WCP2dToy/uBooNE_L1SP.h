#ifndef WIRECELL_UBOONE_L1SP_H
#define WIRECELL_UBOONE_L1SP_H

#include "WCPData/GeomWire.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TF1.h"


namespace WCP2dToy{
  class uBooNE_L1SP {
  public:
    uBooNE_L1SP(TH2F *hv_raw, TH2F* hv_decon, TH2F *hv_decon_g, int nrebin, double time_offset = 0);
    ~uBooNE_L1SP();

    void AddWires(int time_slice, WCP::GeomWireSelection& wires);

    void AddWireTime_Raw();
    
    void Form_rois(int pad = 2);
    void L1_fit(int wire_index, int start_tick, int end_tick);
    
    std::set<int>& get_time_slice_set(){return time_slice_set;};
    
  protected:
    TH2F *hv_raw;
    TH2F *hv_decon;
    TH2F *hv_decon_g;
    TGraph *gv, *gw;
    TF1 *filter_g;

    std::set<int> time_slice_set;
    std::set<int> init_time_slice_set;

    double time_offset;
    
    int nrebin;
    // wire index --> time
    std::map<int,std::set<int>> init_map;
  };
}

#endif
