#ifndef WIRECELL_BLOBMETRIC_H
#define WIRECELL_BLOBMETRIC_H

#include "WireCell2dToy/SimpleBlobToyTiling.h"
#include "WireCellData/GeomCell.h"

namespace WireCell2dToy{
  class BlobMetric{
  public:
    BlobMetric();
    virtual ~BlobMetric();
    
    void Add(WireCell2dToy::SimpleBlobToyTiling &blobtiling, WireCell::CellChargeMap& ccmap);
    void Print();
    
  private:
    int rm_cell_true; float charge_rm_cell_true; float Tcharge_rm_cell_true;
    int rm_cell_false; float charge_rm_cell_false;
    int el_cell_true;  float charge_el_cell_true; float Tcharge_el_cell_true;
    int el_cell_false; float charge_el_cell_false;
    
    

  };
}

#endif
