#include "WireCell2dToy/ToyTiling.h"
#include "WireCellNav/ExampleGDS.h"
#include "WireCellData/Slice.h"
#include "WireCellData/ChannelCharge.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>
using namespace std;

template <typename OK>
void assert(const OK& ok, std::string msg="FAIL")
{
    if (ok) return;
    cerr << msg << endl;
    exit(1);
}

int main()
{
    WireCell::GeomDataSource* gds = WireCell::make_example_gds();
    WireCell::Channel::Group group;

    // fill in some bogus charge 
    const int index_offset = 10; // just pick something away from the edge
    for (int iplane=0; iplane<3; ++iplane) {
	WireCell::WirePlaneType_t plane = static_cast<WireCell::WirePlaneType_t>(iplane);

	for (int ind=0; ind<10; ++ind) {
	    const WireCell::GeomWire* wire = gds->by_planeindex(plane, index_offset+ind);
	    group.push_back(WireCell::Channel::Charge(wire->channel(),100.0));
	}
    }
    WireCell::Slice slice(0, group);
    WireCell2dToy::ToyTiling* tt1 = new WireCell2dToy::ToyTiling(slice, *gds);

    TFile* tfile = TFile::Open("test_tt.root","RECREATE");
    TTree* ttree = new TTree("toy","ToyData");
    ttree->Branch("tt",&tt1);
    ttree->Fill();
    tfile->Write();
    tfile->Close();
    delete tfile;
    delete tt1;
    tt1 = 0;
    ttree = 0;
    tfile = 0;			// it's dead, Jim!

    tfile = TFile::Open("test_tt.root");
    assert(tfile, "Failed to open TFile for reading");

    ttree = dynamic_cast<TTree*>(tfile->Get("toy"));
    assert(ttree, "Failed to read tree 'toy'");

    ttree->SetBranchAddress("tt",&tt1);
    assert(tt1,"Failed to set branch for toy tree");

    assert(1 == ttree->GetEntries(), "Failed to get entries from tree");
    assert(ttree->GetEntry(0), "Failed to load entry 0");
    assert(tt1->get_wire_u().size() == 10, "Failed to get the right number of u wires");
    assert(tt1->get_wire_all().size() == 30, "Failed to get the right number of all wires");

    
}
