#ifndef WIRECELL2DTOY_TOYEVENTDISPLAY_H
#define WIRECELL2DTOY_TOYEVENTDISPLAY_H
#include "WireCellData/Point.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/MergeGeomCell.h"


#include "WireCell2dToy/ToyMatrix.h"

#include "WireCellNav/SliceDataSource.h"
#include "WireCellNav/GeomDataSource.h"
#include "TPad.h"
#include "TH2F.h"
#include "TString.h"
#include "TGraph.h"
#include "TPolyLine.h"

namespace WireCell2dToy {

    class ToyEventDisplay {
    private:
    
	TPad &pad;
	const WireCell::GeomDataSource& gds;

	TH2F *h1;
	TGraph *g1;
	TGraph *g2;
	TPolyLine *g3;
    
    public:
	/// Create a ToyEventDisplay drawing into the given TPad using the given GeomDataSource.
	ToyEventDisplay(TPad& pad, const WireCell::GeomDataSource& gds);
	virtual ~ToyEventDisplay();
    
	/// Initialize with extents
	virtual int init(float x_min=4.9, float x_max=6.1, float y_min=-1.1, float y_max=1.1);
    
	/// Draw the cell charges
	virtual int draw_mc(int flag, const WireCell::PointValueVector& cellcharges, TString option);
    
	/// Draw a slice
	virtual int draw_slice(const WireCell::Slice& slice, TString option);

	/// Draw a selection of cells
	virtual int draw_cells(const WireCell::GeomCellSelection& cellall ,TString option,int color=4);
	virtual int draw_mergecells(const WireCell::GeomCellSelection& cellall ,TString option, int flag=0);
	virtual int draw_truthcells(const WireCell::CellChargeMap& ccmap,TString option);
	virtual int draw_truthcells_charge(const WireCell::CellChargeMap& ccmap,TString option, int FI);
	virtual int draw_wires_charge(const WireCell::WireChargeMap& wcmap,TString option, int FI);
	virtual int draw_cells_charge(const WireCell::GeomCellSelection& cellall ,TString option);

	virtual int draw_reconcells(const WireCell::GeomCellSelection& cellall, WireCell2dToy::ToyMatrix *toymatrix ,TString option, int color = 4);

	

	/// Clear visual the event display data.
	void clear();

	Double_t charge_min, charge_max;
    };

}

#endif
