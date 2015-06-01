#ifndef WIRECELL_TOYMATRIX_H
#define WIRECELL_TOYMATRIX_H

#include "TMatrixD.h"
#include "TVectorD.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"
#include "WireCell2dToy/MergeToyTiling.h"

namespace WireCell2dToy{
  class ToyMatrix {
  public:
    ToyMatrix(WireCell2dToy::ToyTiling& toytiling, WireCell2dToy::MergeToyTiling& mergetiling);
    virtual ~ToyMatrix();

    int Solve();

    float Get_Chi2(){return chi2;};
    float Get_Solve_Flag(){return solve_flag;};

  protected:
    TMatrixD *MA, *MB, *MAT, *MBT;
    TMatrixD *Vy, *VBy, *Vx, *VBy_inv, *Vx_inv;

    TMatrixD *MC, *MC_inv;

    TVectorD *Wy, *Cx, *dCx;
    float chi2;
    int solve_flag;

    void Buildup_index(WireCell2dToy::MergeToyTiling& mergetiling);

    WireCell::WireIndexMap mwimap, swimap;
    WireCell::CellIndexMap mcimap;

    int mcindex;
    int mwindex;
    int swindex;
  };
}

#endif
