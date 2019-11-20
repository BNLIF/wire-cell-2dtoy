#ifndef WIRECELL2DTOY_TOYEVENTDISPLAY_H
#define WIRECELL2DTOY_TOYEVENTDISPLAY_H
#include "WCPData/Point.h"
#include "WCPData/GeomCell.h"
#include "WCPData/MergeGeomCell.h"


#include "WCP2dToy/ToyMatrix.h"

#include "WCPNav/SliceDataSource.h"
#include "WCPNav/GeomDataSource.h"
#include "WCPNav/DetectorGDS.h"

#include "TPad.h"
#include "TH2F.h"
#include "TString.h"
#include "TGraph.h"
#include "TPolyLine.h"

namespace WCP2dToy {

    class ToyEventDisplay {
    private:
    
	TPad &pad;
	int gds_flag;
	const WCP::GeomDataSource* gds;
	const WCP::DetectorGDS* dgds;

	TH2F *h1;
	TH2F *h2;
	TGraph *g1;
	TGraph *g2;
	TGraph *g2b;
	TPolyLine *g3;
    
    public:
	/// Create a ToyEventDisplay drawing into the given TPad using the given GeomDataSource.
	ToyEventDisplay(TPad& pad, const WCP::GeomDataSource& gds);
	ToyEventDisplay(TPad& pad, const WCP::DetectorGDS& gds);
	virtual ~ToyEventDisplay();
    
	/// Initialize with extents
	virtual int init(float x_min=4.9, float x_max=6.1, float y_min=-1.1, float y_max=1.1);
    
	/// Draw the cell charges
	virtual int draw_mc(int flag, const WCP::PointValueVector& cellcharges, TString option);
    
	/// Draw a slice
	virtual int draw_slice(const WCP::Slice& slice, TString option);
	virtual int draw_wires(WCP::GeomWireSelection& wires, TString option);
	virtual int draw_merged_wires(WCP::GeomWireSelection wires, TString option, int color = 2);

	virtual int draw_points(WCP::PointVector pcells, TString option, int color = 2);

	/// Draw a selection of cells
	virtual int draw_cells(const WCP::GeomCellSelection& cellall ,TString option,int color=4);
	virtual int draw_mergecells(const WCP::GeomCellSelection& cellall ,TString option, int flag=0);
	virtual int draw_truthcells(const WCP::CellChargeMap& ccmap,TString option);
	virtual int draw_truthcells_charge(const WCP::CellChargeMap& ccmap,TString option, int FI);
	virtual int draw_wires_charge(const WCP::WireChargeMap& wcmap,TString option, int FI);
	virtual int draw_cells_charge(const WCP::GeomCellSelection& cellall ,TString option);

	virtual int draw_reconcells(const WCP::GeomCellSelection& cellall, WCP2dToy::ToyMatrix *toymatrix ,TString option, int color = 4);

	void Set_ReconThreshold(int abc){recon_threshold = abc;};
	void Set_TruthThreshold(int abc){truth_threshold = abc;};

	void draw_bad_region(WCP::ChirpMap& chirpmap, int time, int scale, int plane, TString option);
	void draw_bad_cell(WCP::GeomCellSelection& cells);

	/// Clear visual the event display data.
	void clear();

	Double_t charge_min, charge_max;

    private:
	int recon_threshold;
	int truth_threshold;
    };

}

#endif
