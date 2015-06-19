#ifndef WIRECELL_SimpleBlobToyTiling_H
#define WIRECELL_SimpleBlobToyTiling_H

#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/ToyMatrix.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"

namespace WireCell2dToy{
  class SimpleBlobToyTiling : public WireCell2dToy::MergeToyTiling {
  public:
    SimpleBlobToyTiling(WireCell2dToy::ToyTiling& toytiling1, WireCell2dToy::MergeToyTiling& mergetiling1, WireCell2dToy::ToyMatrix& toymatrix1, 
			WireCell2dToy::MergeToyTiling& prev_mergetiling, WireCell2dToy::ToyMatrix& prev_toymatrix,
			 WireCell2dToy::MergeToyTiling& next_mergetiling, WireCell2dToy::ToyMatrix& next_toymatrix
			);
    ~SimpleBlobToyTiling();
  private:
    WireCell2dToy::ToyTiling *toytiling;
    WireCell2dToy::MergeToyTiling *mergetiling;
    WireCell2dToy::ToyMatrix *toymatrix;
    
    int nsimple_blob; //smaller than 5 for now
    WireCell::GeomCellSelection corner_mcells[10];
    WireCell::GeomCellSelection corner_smcells[10];
    
  };
}

#endif
