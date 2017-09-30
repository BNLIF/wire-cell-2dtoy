#ifndef WIRECELL_UBOONE_LIGHT_RECO_H
#define WIRECELL_UBOONE_LIGHT_RECO_H

#include "WireCellData/COphit.h"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include <vector>

namespace WireCell2dToy{
  class uBooNE_light_reco{
  public:
    uBooNE_light_reco(const char* root_file);
    ~uBooNE_light_reco();

    void load_event(int eve_num);
    
  protected:
    TFile *file;
    TTree *T;
    WireCell::COphitSelection op_hits;
  };
}

#endif
