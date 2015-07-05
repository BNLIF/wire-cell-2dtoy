#ifndef WIRECELL_TRUTHTOYTILING_H
#define WIRECELL_TRUTHTOYTILING_H

#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/ToyDepositor.h"
#include "WireCellNav/GeomDataSource.h"

#include "Rtypes.h"

namespace WireCell2dToy{
  class TruthToyTiling : public WireCell2dToy::ToyTiling {
  public:
    TruthToyTiling(){};
    ~TruthToyTiling(){};
    TruthToyTiling(WireCell2dToy::ToyTiling& tiling, const WireCell::PointValueVector &pvv, int tbin, const WireCell::GeomDataSource& gds, int offset1 = 0);
    
    float charge(const WireCell::GeomCell& cell) const {return cellchargemap.find(&cell)->second;};
    
    WireCell::CellChargeMap ccmap(){return cellchargemap;};

  protected:
    WireCell::CellChargeMap cellchargemap;
    int offset;
    
    ClassDef(TruthToyTiling,1);
  };
}


#endif
