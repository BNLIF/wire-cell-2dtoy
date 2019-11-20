#ifndef WIRECELL_BLOBMETRIC_H
#define WIRECELL_BLOBMETRIC_H

#include "WCP2dToy/SimpleBlobToyTiling.h"
#include "WCPData/GeomCell.h"

namespace WCP2dToy{
  class BlobMetric{
  public:
    BlobMetric();
    virtual ~BlobMetric();
    
    void Add(WCP2dToy::SimpleBlobToyTiling &blobtiling, WCP::CellChargeMap& ccmap);
    void Print();
    
  private:
    int rm_cell_true; float charge_rm_cell_true; float Tcharge_rm_cell_true;
    int rm_cell_false; float charge_rm_cell_false;
    int el_cell_true;  float charge_el_cell_true; float Tcharge_el_cell_true;
    int el_cell_false; float charge_el_cell_false;
    
    

  };
}

#endif
