#ifndef WIRECELL_UBOONE_LIGHT_RECO_H
#define WIRECELL_UBOONE_LIGHT_RECO_H


#include "WCPData/Opflash.h"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include <vector>

namespace WCP2dToy{
  class uBooNE_light_reco{
  public:
    uBooNE_light_reco(const char* root_file);
    ~uBooNE_light_reco();

    void load_event_raw(int eve_num); // raw data from swizzler
    void load_event(int eve_num); // from 'saturation' producer w/in uboonecode
    TH1F* get_raw_hist(int ch){return hraw[ch];};
    TH1F* get_decon_hist(int ch){return hdecon[ch];};
    TH1F* get_l1_hist(int ch){return hl1[ch];};

    TClonesArray* get_rawWfm(){return fop_wf;}
    std::vector<short>* get_rawChan(){return fop_femch;}
    std::vector<double>* get_rawTimestamp(){return fop_timestamp;}
    TH1F* get_totPE(){return h_totPE;};
    TH1F* get_mult(){return h_mult;};
    TH1F* get_l1_mult(){return h_l1_mult;};
    TH1F* get_l1_totPE(){return h_l1_totPE;};
    
    WCP::OpflashSelection& get_flashes(){return flashes;};
    WCP::OpflashSelection& get_cosmic_flashes(){return cosmic_flashes;};
    WCP::OpflashSelection& get_beam_flashes(){return beam_flashes;};
    
  protected:
    void Process_beam_wfs();
    void sort_flashes();
    std::pair<double,double> cal_mean_rms(TH1 *hist, int nbin=1500);
    void mergeRawBeam(TClonesArray *hg_wf, std::vector<short> *hg_chan, std::vector<double> *hg_timestamp,
		      TClonesArray *lg_wf, std::vector<short> *lg_chan, std::vector<double> *lg_timestamp,
		      bool beamDisc, std::vector<int> *OpChanToOpDet, int ctr);
    void mergeRawCosmic(TClonesArray *hg_wf, std::vector<short> *hg_chan, std::vector<double> *hg_timestamp,
			TClonesArray *lg_wf, std::vector<short> *lg_chan, std::vector<double> *lg_timestamp,
			bool beamDisc, std::vector<int> *OpChanToOpDet, int ctr); 
    void reorderMerged(TClonesArray *beam_wf, std::vector<short> *beam_femch, std::vector<double> *beam_timestamp,
		       TClonesArray *cosmic_wf, std::vector<short> *cosmic_femch, std::vector<double> *cosmic_timestamp);
    float findScaling(int opdet);
    
    TFile *file;
    TTree *T;
    WCP::OpflashSelection cosmic_flashes;
    WCP::OpflashSelection beam_flashes;
    WCP::OpflashSelection flashes;
    // WCP::COphitSelection op_hits;

    double gain[32];
    double beam_dt[32];
    TH1F **hraw; // raw
    TH1F **hdecon; // deconvolution and rebin ...
    TH1F **hl1;

    TH1F *h_totPE;
    TH1F *h_mult;
    TH1F *h_l1_mult;
    TH1F *h_l1_totPE;
    
    TClonesArray *fop_wf_beam;
    std::vector<short> *fop_femch_beam;
    std::vector<double> *fop_timestamp_beam;
    int ctr_beam;

    TClonesArray *fop_wf_cosmic;
    std::vector<short> *fop_femch_cosmic;
    std::vector<double> *fop_timestamp_cosmic;
    int ctr_cosmic;

    TClonesArray *fop_wf;
    std::vector<short> *fop_femch;
    std::vector<double> *fop_timestamp;
    int ctr;
  };
}

#endif
