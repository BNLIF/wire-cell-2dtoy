#ifndef WIRECELL_UBOONE_LIGHT_RECO_H
#define WIRECELL_UBOONE_LIGHT_RECO_H


#include "WireCellData/Opflash.h"

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
    TH1F* get_raw_hist(int ch){return hraw[ch];};
    TH1F* get_decon_hist(int ch){return hdecon[ch];};
    
  protected:
    void Process_beam_wfs();
    std::pair<double,double> cal_mean_rms(TH1 *hist, int nbin=1500);
    
    TFile *file;
    TTree *T;
    WireCell::OpflashSelection flashes;
    // WireCell::COphitSelection op_hits;

    double gain[32];
    TH1F **hraw; // raw
    TH1F **hdecon; // deconvolution and rebin ...
    
    
    
  };
}

#endif
