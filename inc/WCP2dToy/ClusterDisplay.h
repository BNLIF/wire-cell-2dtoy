#ifndef WIRECELL2DTOY_CLUSTERDISPLAY_H
#define WIRECELL2DTOY_CLUSTERDISPLAY_H
#include "WCPData/Point.h"
#include "WCPData/SpaceCell.h"
#include "WCPData/MergeSpaceCell.h"

#include "WCP2dToy/ToyCrawler.h"
#include "WCP2dToy/ToyTracking.h"

#include "TPad.h"
#include "TH2F.h"
#include "TString.h"
#include "TGraph.h"
#include "TPolyLine.h"
#include "TString.h"

namespace WCP2dToy {

    class ClusterDisplay {
    private:
    
	TPad &pad;
	
    
    public:
	/// Create a ToyEventDisplay drawing into the given TPad using the given GeomDataSource.
	ClusterDisplay(TPad& pad);
	virtual ~ClusterDisplay();
    
	void DrawCluster(WCP::SpaceCellSelection& cells, TString option = "p0");
	void DrawCluster(WCP::MergeSpaceCellSelection& mcells, TString option = "p");
	void DrawCluster(WCP::MergeSpaceCellSelection& mcells, WCP2dToy::ToyTracking& toytracking);

	void DrawCrawler(ToyCrawler& toycrawler, TString option = "p", int flag = 0);
	void DrawHough(WCP::SpaceCellSelection& cells, WCP::Point& p, double dis_near, double dis_far);

	void DrawShower(WCP::WCShower* shower, TString option = "psame", int color = 2); 

	void DrawVertex(WCP::WCVertexSelection& vertices, TString option = "p");
	
	void DrawTracks(WCP::WCTrackSelection& tracks, TString option = "", int color = 1);
    };

}

#endif
