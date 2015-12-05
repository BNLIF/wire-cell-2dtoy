#ifndef TOYMATRIXITERATE_SINGLEWIRE_H
#define TOYMATRIXITERATE_SINGLEWIRE_H

#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/ToyMatrixKalman.h"

namespace WireCell2dToy{
  class ToyMatrixIterate_SingleWire{
  public:
    ToyMatrixIterate_SingleWire(WireCell2dToy::ToyMatrix& toymatrix, WireCell2dToy::MergeToyTiling* mergetiling);
    ToyMatrixIterate_SingleWire(WireCell2dToy::ToyMatrix *toybefore, WireCell2dToy::ToyMatrix *toycur, WireCell2dToy::ToyMatrix *toyafter, WireCell2dToy::MergeToyTiling *mergebefore, WireCell2dToy::MergeToyTiling *mergecur, WireCell2dToy::MergeToyTiling *mergeafter, 
				int recon_t = 2000, float limit = 1e5, double penalty = 5, double penalty_ncpt = 0);  
    
    virtual ~ToyMatrixIterate_SingleWire();
  protected:
    WireCell2dToy::ToyMatrix &toymatrix;
    WireCell2dToy::MergeToyTiling *mergetiling;
    WireCell::WireChargeMap wirechargemap;
    int ncount;
    int nlevel;
    float limit;
    
    double penalty_ncpt;
    std::map<int,double> cell_penal;
    std::map<int,std::vector<int>> cells_ncpt;

    void Iterate(WireCell::GeomCellSelection cells, WireCell::GeomCellSelection single_cells, WireCell::GeomCellSelection tried_cell, WireCell::GeomCellMap cellmap, WireCell::GeomWireMap wiremap);

  };
}

#endif
