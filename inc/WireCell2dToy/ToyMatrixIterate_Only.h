#ifndef TOYMATRIXITERATE_ONLY_H
#define TOYMATRIXITERATE_ONLY_H

#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/ToyMatrixKalman.h"

namespace WireCell2dToy{
  class ToyMatrixIterate_Only{
  public:
    ToyMatrixIterate_Only(WireCell2dToy::ToyMatrix& toymatrix, WireCell2dToy::MergeToyTiling* mergetiling);
    ~ToyMatrixIterate_Only();
    

  protected:
     WireCell2dToy::ToyMatrix &toymatrix;
     WireCell2dToy::MergeToyTiling *mergetiling;
     WireCell::WireChargeMap wirechargemap;
     WireCell::GeomCellMap cellmap;
     
     WireCell::GeomCellSelection all_cells;
     
     void Iterate(WireCell::GeomCellSelection remaining_cells, 
		  WireCell::GeomCellSelection good_cells,
		  WireCell::GeomWireSelection used_wires,
		  WireCell::GeomCellSelection bad_cells,
		  WireCell::GeomCellSelection& tried_cells);

     int ncount;
     int nlevel;
     
  };
}

#endif

