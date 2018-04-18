#ifndef WIRECELL_TOYLIGHTRECO_H
#define WIRECELL_TOYLIGHTRECO_H

#include "WireCellData/Opflash.h"

#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TClonesArray.h"
#include <vector>

namespace WireCell2dToy{

  struct pmtDisc{
    short channel;
    double timestamp;
    std::vector<short> wfm;
    bool saturated;
    bool isolated;
    bool highGain;
  };
  typedef std::map<short,pmtDisc> pmtMap;
  struct timeOrder_pmtDisc{
    bool operator()(const pmtDisc& a, const pmtDisc& b) const
    {return a.timestamp < b.timestamp;}
  };
  typedef std::set<pmtDisc,timeOrder_pmtDisc> pmtSet;
  typedef std::map<short,pmtSet> pmtMapSet;
  typedef std::pair<pmtDisc,pmtDisc> pmtPair;
  struct timeOrder_pmtPair{
    bool operator()(const pmtPair& a, const pmtPair& b) const
    {return a.first.timestamp < b.first.timestamp;}
  };
  typedef std::set<pmtPair,timeOrder_pmtPair> pmtSetPair;
  typedef std::map<short,pmtPair> pmtMapPair;
  typedef std::map<short,pmtSetPair> pmtMapSetPair;
  typedef std::vector<std::pair<short,short> > saturationTick;
  
  class ToyLightReco{
  public:
    ToyLightReco(const char* root_file, bool imagingoutput=false);
    ~ToyLightReco();

    void load_event_raw(int eve_num); // raw data from swizzler
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
    
    WireCell::OpflashSelection& get_flashes(){return flashes;};
    WireCell::OpflashSelection& get_cosmic_flashes(){return cosmic_flashes;};
    WireCell::OpflashSelection& get_beam_flashes(){return beam_flashes;};

    void clear_flashes();

  private:
    //    bool delete_status;

  protected:
    void Process_beam_wfs();
    void sort_flashes();
    void update_pmt_map();
    
    pmtMapSet makePmtContainer(bool high, bool beam, TClonesArray *wf, std::vector<short> *chan, std::vector<double> *timestamp);
    pmtMapPair makeBeamPairs(pmtMapSet &high, pmtMapSet &low);
    pmtMapSetPair makeCosmicPairs(pmtMapSet &high, pmtMapSet &low);
    pmtMap mergeBeam(pmtMapPair &beam);
    pmtMapSet mergeCosmic(pmtMapSetPair &cosmic);
    saturationTick findSaturationTick(std::vector<short> &wfm);
    std::vector<short> replaceSaturatedBin(std::vector<short> &high, std::vector<short> &low, saturationTick &st);
    void dumpPmtVec(pmtMap &beam, pmtMapSet &cosmic);

    double findBaselineLg(TH1 *hist, int nbin=1500);
    std::pair<double,double> cal_mean_rms(TH1 *hist, int nbin=1500);

    float findScaling(int opdet);

    TFile *file;
    TTree *T;
    WireCell::OpflashSelection cosmic_flashes;
    WireCell::OpflashSelection beam_flashes;
    WireCell::OpflashSelection flashes;
    WireCell::COphitSelection op_hits;

    TClonesArray* cosmic_hg_wf;
    TClonesArray* cosmic_lg_wf;
    TClonesArray* beam_hg_wf;
    TClonesArray* beam_lg_wf;
    std::vector<short> *cosmic_hg_opch;
    std::vector<short> *cosmic_lg_opch;
    std::vector<short> *beam_hg_opch;
    std::vector<short> *beam_lg_opch;
    std::vector<double> *cosmic_hg_timestamp;
    std::vector<double> *cosmic_lg_timestamp;
    std::vector<double> *beam_hg_timestamp;
    std::vector<double> *beam_lg_timestamp;
    std::vector<float> *op_gain;
    std::vector<float> *op_gainerror;
    

    double triggerTime;
    double gain[32];
    double beam_dt[32];
    TH1F **hraw; // raw
    TH1F **hdecon; // deconvolution and rebin ...
    TH1F **hl1;

    TH1F *h_totPE;
    TH1F *h_mult;
    TH1F *h_l1_mult;
    TH1F *h_l1_totPE;
    
    /* TClonesArray *fop_wf_beam; */
    /* std::vector<short> *fop_femch_beam; */
    /* std::vector<double> *fop_timestamp_beam; */
    /* int ctr_beam; */

    /* TClonesArray *fop_wf_cosmic; */
    /* std::vector<short> *fop_femch_cosmic; */
    /* std::vector<double> *fop_timestamp_cosmic; */
    /* int ctr_cosmic; */

    TClonesArray *fop_wf;
    std::vector<short> *fop_femch;
    std::vector<double> *fop_timestamp;
    int ctr;
  };
}

#endif
