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

    void load_event_raw(int eve_num); // raw data from swizzler
    void load_event(int eve_num); // from 'saturation' producer w/in uboonecode
    TH1F* get_raw_hist(int ch){return hraw[ch];};
    TH1F* get_decon_hist(int ch){return hdecon[ch];};
    TH1F* get_l1_hist(int ch){return hl1[ch];};

    TH1F* get_totPE(){return h_totPE;};
    TH1F* get_mult(){return h_mult;};
    TH1F* get_l1_mult(){return h_l1_mult;};
    TH1F* get_l1_totPE(){return h_l1_totPE;};
    
    WireCell::OpflashSelection& get_flashes(){return flashes;};
    WireCell::OpflashSelection& get_cosmic_flashes(){return cosmic_flashes;};
    WireCell::OpflashSelection& get_beam_flashes(){return beam_flashes;};
    
  protected:
    void Process_beam_wfs();
    void sort_flashes();
    std::pair<double,double> cal_mean_rms(TH1 *hist, int nbin=1500);
    void mergeRaw(TH1S *h_hg, short hg_chan, double hg_timestamp,
		  Int_t nentries_lg, TH1S *h_lg, std::vector<short> *lg_chan, std::vector<double> *lg_timestamp,
		  bool beamDisc = true, std::vector<int> OpChanToOpDet);
    
    TFile *file;
    TTree *T;
    WireCell::OpflashSelection cosmic_flashes;
    WireCell::OpflashSelection beam_flashes;
    WireCell::OpflashSelection flashes;
    // WireCell::COphitSelection op_hits;

    double gain[32];
    double beam_dt[32];
    TH1F **hraw; // raw
    TH1F **hdecon; // deconvolution and rebin ...
    TH1F **hl1;

    TH1F *h_totPE;
    TH1F *h_mult;
    TH1F *h_l1_mult;
    TH1F *h_l1_totPE;
    
    
  };
}

#endif
