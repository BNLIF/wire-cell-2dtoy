#ifndef TOYMATRIXITERATE_SINGLEWIRE_H
#define TOYMATRIXITERATE_SINGLEWIRE_H

#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/ToyMatrixKalman.h"

namespace WireCell2dToy{
  class ToyMatrixIterate_SingleWire{
  public:
    ToyMatrixIterate_SingleWire(WireCell2dToy::ToyMatrix& toymatrix, WireCell2dToy::MergeToyTiling* mergetiling);
    virtual ~ToyMatrixIterate_SingleWire();
  protected:
    WireCell2dToy::ToyMatrix &toymatrix;
    WireCell2dToy::MergeToyTiling *mergetiling;
    WireCell::WireChargeMap wirechargemap;
    int ncount;
    int nlevel;
    
    void Iterate(WireCell::GeomCellSelection cells, WireCell::GeomCellSelection single_cells, WireCell::GeomCellSelection tried_cell, WireCell::GeomCellMap cellmap, WireCell::GeomWireMap wiremap);

  };
}

#endif
