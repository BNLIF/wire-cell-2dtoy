#ifndef WIRECELL_SimpleBlobToyTiling_H
#define WIRECELL_SimpleBlobToyTiling_H

#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/ToyMatrix.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"
#include "WireCell2dToy/ToyHypothesis.h"

namespace WireCell2dToy{
  class SimpleBlobToyTiling : public WireCell2dToy::MergeToyTiling {
  public:
    SimpleBlobToyTiling(WireCell2dToy::ToyTiling& toytiling1, WireCell2dToy::MergeToyTiling& mergetiling1, WireCell2dToy::ToyMatrix& toymatrix1, 
			WireCell2dToy::MergeToyTiling& prev_mergetiling, WireCell2dToy::ToyMatrix& prev_toymatrix,
			 WireCell2dToy::MergeToyTiling& next_mergetiling, WireCell2dToy::ToyMatrix& next_toymatrix
			);
    ~SimpleBlobToyTiling();
  private:
    void Organize(int nsimple_blob);
    void FormHypo();
    void ClearHypo();

    WireCell2dToy::ToyTiling *toytiling;
    WireCell2dToy::MergeToyTiling *mergetiling;
    WireCell2dToy::ToyMatrix *toymatrix;
    
    int nsimple_blob; //smaller than 5 for now
    /* WireCell::GeomCellSelection corner_mcells[10]; */
    /* WireCell::GeomCellSelection corner_smcells[10]; */
        
    //save the first pass results
    WireCell::GeomCellSelectionV corner_mcells;
    WireCell::GeomCellSelectionV corner_smcells;
    WireCell::CellIndexMap cell_rank;

    //save the second pass results to use for hypothesis formation
    WireCell::GeomCellSelectionV first_cell;
    WireCell::GeomCellSelectionV second_cell;
    WireCell::GeomCellSelectionV other_cell;
    std::vector<int> flag_cell;
    
    int ncount;

    std::vector<WireCell2dToy::HypoSelection> cur_hypo;
    
  };
}

#endif
