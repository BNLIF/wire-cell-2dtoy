#ifndef WIRECELL_TOYHYPOTHESIS_H
#define WIRECELL_TOYHYPOTHESIS_H

#include "WireCellData/MergeGeomCell.h"

namespace WireCell2dToy{
  class ToyHypothesis {
  public:
    ToyHypothesis();
    ToyHypothesis(WireCell::MergeGeomCell& mcell1, WireCell::MergeGeomCell& mcell2);
    ~ToyHypothesis();
    double CalValue(WireCell::Point p, WireCell::Point p1, WireCell::Point p2);
    bool IsInside(WireCell::Point p);
    bool IsInside(const WireCell::GeomCell& cell);
    
    
  private:
    WireCell::Point p1[2],p2[2],pc[2];
  };
  typedef std::vector<const WireCell2dToy::ToyHypothesis*> HypoSelection;
}


#endif
