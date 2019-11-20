#ifndef WIRECELL_TOYHYPOTHESIS_H
#define WIRECELL_TOYHYPOTHESIS_H

#include "WCPData/MergeGeomCell.h"

namespace WCP2dToy{
  class ToyHypothesis {
  public:
    ToyHypothesis();
    ToyHypothesis(WCP::MergeGeomCell& mcell1, WCP::MergeGeomCell& mcell2);
    ~ToyHypothesis();
    double CalValue(WCP::Point p, WCP::Point p1, WCP::Point p2);
    bool IsInside(WCP::Point p);
    bool IsInside(const WCP::GeomCell& cell);
    
    
  private:
    WCP::Point p1[2],p2[2],pc[2];
  };
  typedef std::vector<const WCP2dToy::ToyHypothesis*> HypoSelection;
}


#endif
