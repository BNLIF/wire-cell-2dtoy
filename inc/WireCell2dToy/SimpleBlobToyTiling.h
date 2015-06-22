#ifndef WIRECELL_SimpleBlobToyTiling_H
#define WIRECELL_SimpleBlobToyTiling_H

#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/ToyMatrix.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"
#include "WireCell2dToy/ToyHypothesis.h"

#include "Rtypes.h"

namespace WireCell2dToy{
  class SimpleBlobToyTiling : public WireCell2dToy::MergeToyTiling {
  public:
    SimpleBlobToyTiling(WireCell2dToy::ToyTiling& toytiling1, WireCell2dToy::MergeToyTiling& mergetiling1, WireCell2dToy::ToyMatrix& toymatrix1, 
			WireCell2dToy::MergeToyTiling& prev_mergetiling, WireCell2dToy::ToyMatrix& prev_toymatrix,
			 WireCell2dToy::MergeToyTiling& next_mergetiling, WireCell2dToy::ToyMatrix& next_toymatrix
			);
    ~SimpleBlobToyTiling();
    
    
    double Get_Cell_Charge( const WireCell::GeomCell *cell, int flag = 1) ;
    WireCell::GeomCellSelection Get_Cells(){return cell_all_save;};
    WireCell::GeomCellSelection Get_SB_Cells(){return sbcells;};

  private:
    void Organize(int nsimple_blob);
    void FormHypo();
    void ClearHypo();
    void DoTiling();
    double CalChi2();
    void SaveResult();
    void Buildup_index();
    

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
    WireCell::GeomCellSelection sbcells;
    std::vector<int> flag_cell;
    
    int ncount;
    double cur_chi2;
    // cell_all, wire_all, cellmap, wiremap wire_u, wire_v, wire_w 
    std::vector<double> Cx, dCx;
    std::vector<WireCell2dToy::HypoSelection> cur_hypo;
    
    int scindex;
    int swindex;
    WireCell::WireIndexMap swimap;
    WireCell::CellIndexMap scimap;

    //save all the results
    double chi2_save;
    std::vector<double> Cx_save, dCx_save;
    WireCell::GeomWireSelection wire_u_save;
    WireCell::GeomWireSelection wire_v_save;
    WireCell::GeomWireSelection wire_w_save;
    WireCell::GeomWireSelection wire_all_save;
    WireCell::GeomCellSelection cell_all_save;
    WireCell::GeomCellMap cellmap_save;
    WireCell::GeomWireMap wiremap_save;
    std::vector<WireCell2dToy::HypoSelection> hypo_save;
    

    int scindex_save;
    int swindex_save;
    WireCell::WireIndexMap swimap_save;
    WireCell::CellIndexMap scimap_save;
    

      ClassDef(SimpleBlobToyTiling,1);
    
  };
}

#endif
