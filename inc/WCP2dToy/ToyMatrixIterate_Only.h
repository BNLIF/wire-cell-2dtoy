#ifndef TOYMATRIXITERATE_ONLY_H
#define TOYMATRIXITERATE_ONLY_H

#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/ToyMatrixKalman.h"

namespace WCP2dToy{
  class ToyMatrixIterate_Only{
  public:
    ToyMatrixIterate_Only(WCP2dToy::ToyMatrix& toymatrix, WCP2dToy::MergeToyTiling* mergetiling);
    ~ToyMatrixIterate_Only();
    

  protected:
     WCP2dToy::ToyMatrix &toymatrix;
     WCP2dToy::MergeToyTiling *mergetiling;
     WCP::WireChargeMap wirechargemap;
     WCP::GeomCellMap cellmap;
     
     WCP::GeomCellSelection all_cells;
     
     void Iterate(WCP::GeomCellSelection remaining_cells, 
		  WCP::GeomCellSelection good_cells,
		  WCP::GeomWireSelection used_wires,
		  WCP::GeomCellSelection bad_cells,
		  WCP::GeomCellSelection& tried_cells);

     int ncount;
     int nlevel;
     
  };
}

#endif

