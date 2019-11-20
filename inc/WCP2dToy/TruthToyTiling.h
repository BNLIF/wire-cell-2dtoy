#ifndef WIRECELL_TRUTHTOYTILING_H
#define WIRECELL_TRUTHTOYTILING_H

#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/ToyDepositor.h"
#include "WCPNav/GeomDataSource.h"

#include "Rtypes.h"

namespace WCP2dToy{
  class TruthToyTiling : public WCP2dToy::ToyTiling {
  public:
    TruthToyTiling(){};
    ~TruthToyTiling(){};
    TruthToyTiling(WCP2dToy::ToyTiling& tiling, const WCP::PointValueVector &pvv, int tbin, const WCP::GeomDataSource& gds, int offset1 = 0, float unit_dis = 1.6);
    TruthToyTiling(WCP2dToy::ToyTiling& tiling, const WCP::PointValueVector &pvv, const std::vector<int> &timeoffset, int tbin, const WCP::GeomDataSource& gds, float unit_dis = 1.6);

    TruthToyTiling(WCP2dToy::ToyTiling& tiling, const WCP::PointValueVector &pvv, int tbin, const WCP::DetectorGDS& gds, int offset1 = 0, float unit_dis = 1.6);
    TruthToyTiling(WCP2dToy::ToyTiling& tiling, const WCP::PointValueVector &pvv, const std::vector<int> &timeoffset, int tbin, const WCP::DetectorGDS& gds, float unit_dis = 1.6);

    float charge(const WCP::GeomCell& cell) const {return cellchargemap.find(&cell)->second;};
    
    WCP::CellChargeMap ccmap(){return cellchargemap;};

  protected:
    WCP::CellChargeMap cellchargemap;
    int offset;
    
    ClassDef(TruthToyTiling,1);
  };
}


#endif
