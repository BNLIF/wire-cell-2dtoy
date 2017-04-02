#ifndef WIRECELL_LOWMEMTILING_H
#define WIRECELL_LOWMEMTILING_H

#include "WireCellTiling/TilingBase.h"
#include "WireCellData/GeomWireCellMap.h"
#include "WireCellNav/SliceDataSource.h"
#include "WireCellSst/GeomDataSource.h"
#include "WireCellNav/DetectorGDS.h"

#include "Rtypes.h"

namespace WireCell2dToy{

  class LowmemTiling : public WireCell::TilingBase {
  public:
    LowmemTiling(const WireCell::Slice& slice,WireCell::GeomDataSource& gds,
		 WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map,
		 std::vector<float>& uplane_rms, std::vector<float>& vplane_rms, std::vector<float>& wplane_rms);
    
    virtual ~LowmemTiling();
  private:
    
    
  };
}



#endif
