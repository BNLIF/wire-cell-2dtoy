#ifndef WIRECELL2dToy_TOYSIGNALSIMUTRUE_H
#define WIRECELL2dToy_TOYSIGNALSIMUTRUE_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"
#include "TH1F.h"
#include "TGraph.h"

namespace WireCell2dToy {
  class ToySignalSimuTrueFDS : public WireCell::FrameDataSource
  {
  public:
    ToySignalSimuTrueFDS(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, int bins_per_frame1 = 9600, int nframes_total = -1, int flag_smear = 1);
    ~ToySignalSimuTrueFDS();

    virtual int size() const;
    virtual int jump(int frame_number);
    //fixed it ...
    

    void Save();
    
  private:
    WireCell::FrameDataSource* fds;
    const WireCell::GeomDataSource& gds;
    int  max_frames;

    int nwire_u, nwire_v, nwire_w;

    int nbin;
    
    int flag_smear;

    /* TH1F **hu; */
    /* TH1F **hv; */
    /* TH1F **hw; */

    TH1F *hu; 
    TH1F *hv; 
    TH1F *hw; 

    TH1F *hfilter_time_gaus;
    TH1 *hfilter_gaus;
    TF1 *filter_g;
    
  };
}

#endif
