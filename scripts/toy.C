/** a ROOT script / command line interface to wire-cell-2dtoy
 *
 * $ root -l 2dtoy/scripts/toy.C
 * root [1] t = make_toy("../wire-cell-event/geometry/ChannelWireGeometry.txt")
 * root [2] t->next_slice()
 */



#include "WCP2dToy/FrameDataSource.h"
#include "WCP2dToy/ToyEventDisplay.h"
#include "WCPNav/SliceDataSource.h"
#include "WCPSst/GeomDataSource.h"
#include "WCP2dToy/ToyTiling.h"

#include "TPad.h"

#include <iostream>
using namespace std;

class WC2DToy
{
public:
    WCPSst::GeomDataSource* gds;
    WCP2dToy::FrameDataSource* fds;
    WCP::SliceDataSource* sds;
    WCP2dToy::ToyEventDisplay* display;

    WC2DToy(TPad& _pad, const char* sstgeometryfile, const char* sstrootfile=0) 
	{
	    gds = new WCPSst::GeomDataSource(sstgeometryfile);
	    fds = new WCP2dToy::FrameDataSource(10, *gds);
	    sds = new WCP::SliceDataSource(*fds);
	    display = new WCP2dToy::ToyEventDisplay(_pad, *gds);

	    fds->jump(0);
	}

    void next_frame() {
	int fnum = fds->next();
	cerr << "Frame: " << fnum << endl;
	update();
    }

    void next_slice() {
	int snum = sds->next();
	cerr << "Slice: " << snum << endl;

	update();
    }

    void update() {

	const WCP::PointValueVector& mctruth = fds->cell_charges();
	cerr << "#true cells: " << mctruth.size() << endl;

	display->init();
	cerr << "Drawing..." << endl;
	display->draw_mc(1,mctruth,"");
	display->draw_mc(2,mctruth,"TEXTsame");

	const WCP::Slice& slice = sds->get();
	cerr << "Slice: tbin=" << slice.tbin() << " with " << slice.group().size() << " charges" << endl;

	display->draw_slice(slice, "same");

	WCP2dToy::ToyTiling tt(slice, *gds);

	display->draw_cells(tt.get_allcell(),"*same");
	display->draw_mc(3,mctruth,"*same");
	
    }
};    

WC2DToy* make_toy(const char* sstgeofile, const char* sstrootfile=0);
WC2DToy* make_toy(const char* sstgeofile, const char* sstrootfile) 
{
    TCanvas *canvas = new TCanvas("ToyMC","ToyMC",800,600); // leak
    canvas->Draw();
    WC2DToy* toy = new WC2DToy(*canvas, sstgeofile, sstrootfile);
    return toy;
}

