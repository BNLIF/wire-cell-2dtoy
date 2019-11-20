#ifndef TOYMATRIXITERATE_SINGLEWIRE_H
#define TOYMATRIXITERATE_SINGLEWIRE_H

#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/ToyMatrixKalman.h"

namespace WCP2dToy{
  class ToyMatrixIterate_SingleWire{
  public:
    ToyMatrixIterate_SingleWire(WCP2dToy::ToyMatrix& toymatrix, WCP2dToy::MergeToyTiling* mergetiling);
    ToyMatrixIterate_SingleWire(WCP2dToy::ToyMatrix *toybefore, WCP2dToy::ToyMatrix *toycur, WCP2dToy::ToyMatrix *toyafter, WCP2dToy::MergeToyTiling *mergebefore, WCP2dToy::MergeToyTiling *mergecur, WCP2dToy::MergeToyTiling *mergeafter, 
				int recon_t = 2000, float limit = 1e5, double penalty = 5, double penalty_ncpt = 0);  
    
    virtual ~ToyMatrixIterate_SingleWire();
  protected:
    WCP2dToy::ToyMatrix &toymatrix;
    WCP2dToy::MergeToyTiling *mergetiling;
    WCP::WireChargeMap wirechargemap;
    int ncount;
    int nlevel;
    float limit;
    
    double penalty_ncpt;
    std::map<int,double> cell_penal;
    std::map<int,std::vector<int>> cells_ncpt;

    void Iterate(WCP::GeomCellSelection cells, WCP::GeomCellSelection single_cells, WCP::GeomCellSelection tried_cell, WCP::GeomCellMap cellmap, WCP::GeomWireMap wiremap);

  };
}

#endif
