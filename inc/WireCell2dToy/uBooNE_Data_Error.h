#ifndef WIRECELL2dToy_UBOONE_DATA_ERROR_H
#define WIRECELL2dToy_UBOONE_DATA_ERROR_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"
#include "TH2F.h"


namespace WireCell2dToy{
  class uBooNEDataError : public WireCell::FrameDataSource
  {
  public:
    uBooNEDataError(const WireCell::GeomDataSource& gds, TH2F *hu, TH2F *hv, TH2F *hw, int eve_num, int nrebin);
    uBooNEDataError(const WireCell::GeomDataSource& gds, TH2I *hu, TH2I *hv, TH2I *hw, int eve_num, int nrebin);
    ~uBooNEDataError();

    void refresh(TH2F *hu, TH2F *hv, TH2F *hw, int eve_num);
    
    virtual int jump(int frame_number);
    virtual int size() const;

  private:
    const WireCell::GeomDataSource& gds;
    int nwire_u;
    int nwire_v;
    int nwire_w;
    int nrebin;
  };

}
#endif
