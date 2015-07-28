#ifndef WIRECELL2DTOY_CLUSTERDISPLAY_H
#define WIRECELL2DTOY_CLUSTERDISPLAY_H
#include "WireCellData/Point.h"
#include "WireCellData/SpaceCell.h"
#include "WireCellData/MergeSpaceCell.h"

#include "WireCell2dToy/ToyCrawler.h"

#include "TPad.h"
#include "TH2F.h"
#include "TString.h"
#include "TGraph.h"
#include "TPolyLine.h"
#include "TString.h"

namespace WireCell2dToy {

    class ClusterDisplay {
    private:
    
	TPad &pad;
	
    
    public:
	/// Create a ToyEventDisplay drawing into the given TPad using the given GeomDataSource.
	ClusterDisplay(TPad& pad);
	virtual ~ClusterDisplay();
    
	void DrawCluster(WireCell::SpaceCellSelection& cells, TString option = "p0");
	void DrawCluster(WireCell::MergeSpaceCellSelection& mcells, TString option = "p");
	void DrawCrawler(ToyCrawler& toycrawler, TString option = "p");
	void DrawHough(WireCell::SpaceCellSelection& cells, WireCell::Point& p, double dis_near, double dis_far);

    };

}

#endif
