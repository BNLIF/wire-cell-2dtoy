#ifndef WIRECELL_TOYMETRIC_H
#define WIRECELL_TOYMETRIC_H

#include "WireCell2dToy/ToyMatrix.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/GeomCell.h"

namespace WireCell2dToy{
  class ToyMetric {
  public:
    ToyMetric();
    virtual ~ToyMetric();

    void Add(WireCell::GeomCellSelection &allmcell, WireCell2dToy::ToyMatrix& toymatrix, WireCell::CellChargeMap& ccmap);
    void Print();
    
  private:
    int rm_cell_true;  float charge_rm_cell_true;
    int rm_cell_false; float charge_rm_cell_false;
    int el_cell_true;  float charge_el_cell_true;
    int el_cell_false; float charge_el_cell_false;


    int threshold = 2000;
  };
}
#endif
