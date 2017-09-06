#ifndef WIRECELL_UBOONE_L1SP_H
#define WIRECELL_UBOONE_L1SP_H

#include "WireCellData/GeomWire.h"
#include "TH2F.h"



namespace WireCell2dToy{
  class uBooNE_L1SP {
  public:
    uBooNE_L1SP(TH2F *hv_raw);
    ~uBooNE_L1SP();

    void AddWires(int time_slice, WireCell::GeomWireSelection& wires);
    void Form_rois();
  protected:
    TH2F *hv_raw;
    // wire index --> time
    std::map<int,std::vector<int>> init_map;
  };
}

#endif
