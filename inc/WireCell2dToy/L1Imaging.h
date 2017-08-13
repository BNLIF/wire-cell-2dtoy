#ifndef WIRECELL_L1IMAGING_H
#define WIRECELL_L1IMAGING_H

#include "TMatrixD.h"
#include "TVectorD.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"
#include "LowmemTiling.h"

namespace WireCell2dToy{
  class L1Imaging {
  public:
    L1Imaging(const WireCell::GeomDataSource& gds, LowmemTiling& tiling);
    virtual ~L1Imaging();

  protected:
    const WireCell::GeomDataSource& gds;
    const LowmemTiling& tiling;
    
  };

}

#endif
