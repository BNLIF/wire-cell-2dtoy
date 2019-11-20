#ifndef WIRECELL_SimpleBlobToyTiling_H
#define WIRECELL_SimpleBlobToyTiling_H

#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/ToyMatrix.h"
#include "WCPData/GeomCell.h"
#include "WCPData/GeomWire.h"
#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"
#include "WCP2dToy/ToyHypothesis.h"

#include "Rtypes.h"

namespace WCP2dToy{
  class SimpleBlobToyTiling : public WCP2dToy::MergeToyTiling {
  public:
    SimpleBlobToyTiling(WCP2dToy::ToyTiling& toytiling1, WCP2dToy::MergeToyTiling& mergetiling1, WCP2dToy::ToyMatrix& toymatrix1, 
			WCP2dToy::MergeToyTiling& prev_mergetiling, WCP2dToy::ToyMatrix& prev_toymatrix,
			 WCP2dToy::MergeToyTiling& next_mergetiling, WCP2dToy::ToyMatrix& next_toymatrix
			, int recon_t = 2000);
    ~SimpleBlobToyTiling();
    
    
    double Get_Cell_Charge( const WCP::GeomCell *cell, int flag = 1) ;
    WCP::GeomCellSelection Get_Cells(){return cell_all_save;};
    WCP::GeomCellSelection Get_SB_Cells(){return sbcells;};

  private:
    void Organize(int nsimple_blob);
    void FormHypo();
    void ClearHypo();
    void DoTiling();
    double CalChi2();
    void SaveResult();
    void Buildup_index();
    int recon_threshold;
    

    WCP2dToy::ToyTiling *toytiling;
    WCP2dToy::MergeToyTiling *mergetiling;
    WCP2dToy::ToyMatrix *toymatrix;
    
    int nsimple_blob; //smaller than 5 for now
    /* WCP::GeomCellSelection corner_mcells[10]; */
    /* WCP::GeomCellSelection corner_smcells[10]; */
        
    //save the first pass results
    WCP::GeomCellSelectionV corner_mcells;
    WCP::GeomCellSelectionV corner_smcells;
    WCP::CellIndexMap cell_rank;

    //save the second pass results to use for hypothesis formation
    WCP::GeomCellSelectionV first_cell;
    WCP::GeomCellSelectionV second_cell;
    WCP::GeomCellSelectionV other_cell;
    WCP::GeomCellSelection sbcells;
    std::vector<int> flag_cell;
    
    int ncount;
    double cur_chi2;
    // cell_all, wire_all, cellmap, wiremap wire_u, wire_v, wire_w 
    std::vector<double> Cx, dCx;
    std::vector<WCP2dToy::HypoSelection> cur_hypo;
    
    int scindex;
    int swindex;
    WCP::WireIndexMap swimap;
    WCP::CellIndexMap scimap;

    //save all the results
    double chi2_save;
    std::vector<double> Cx_save, dCx_save;
    WCP::GeomWireSelection wire_u_save;
    WCP::GeomWireSelection wire_v_save;
    WCP::GeomWireSelection wire_w_save;
    WCP::GeomWireSelection wire_all_save;
    WCP::GeomCellSelection cell_all_save;
    WCP::GeomCellMap cellmap_save;
    WCP::GeomWireMap wiremap_save;
    std::vector<WCP2dToy::HypoSelection> hypo_save;
    

    int scindex_save;
    int swindex_save;
    WCP::WireIndexMap swimap_save;
    WCP::CellIndexMap scimap_save;
    

      ClassDef(SimpleBlobToyTiling,1);
    
  };
}

#endif
