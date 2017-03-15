#ifndef WIRECELL2dToy_PD_DATA_FDS_H
#define WIRECELL2dToy_PD_DATA_FDS_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"

#include "TH2F.h"


namespace WireCell2dToy{
  class pdDataFDS : public WireCell::FrameDataSource
  {
  public:
    pdDataFDS(const WireCell::GeomDataSource& gds, TH2I *hu, TH2I *hv, TH2I *hw, int eve_num);
    ~pdDataFDS();

    virtual int jump(int frame_number);
    virtual int size() const;

  private:
    const WireCell::GeomDataSource& gds;
    int nwire_u;
    int nwire_v;
    int nwire_w;
  };

}
#endif
