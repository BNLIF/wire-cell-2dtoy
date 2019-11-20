#ifndef WIRECELL_TOYMETRIC_H
#define WIRECELL_TOYMETRIC_H

#include "WCP2dToy/ToyMatrix.h"
#include "WCPData/MergeGeomCell.h"
#include "WCPData/GeomCell.h"

namespace WCP2dToy{
  class ToyMetric {
  public:
    ToyMetric();
    virtual ~ToyMetric();

    void Add(WCP::GeomCellSelection &allmcell, WCP2dToy::ToyMatrix& toymatrix, WCP::CellChargeMap& ccmap);
    void Print();
    void AddSolve(int cond);
    
  private:
    int rm_cell_true;  float charge_rm_cell_true; float Tcharge_rm_cell_true;
    int rm_cell_false; float charge_rm_cell_false;
    int el_cell_true;  float charge_el_cell_true; float Tcharge_el_cell_true;
    int el_cell_false; float charge_el_cell_false;

    int solve_condition[3];

    int threshold;
  };
}
#endif
