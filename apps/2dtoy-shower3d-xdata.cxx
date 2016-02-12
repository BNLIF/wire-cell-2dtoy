/// Convert "shower3D" file to a Wire Cell Exchange Data ROOT format.

#include "WireCellData/GeomCell.h"

#include "WireCellXdataRoot/Writer.h"
#include "WireCellXdataRoot/Reader.h"
#include "WireCellXdataRoot/CloneHelper.h"
#include "WireCellXdataRoot/WireDB.h"

#include "TFile.h"
#include "TTree.h"

#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <chrono>

using namespace std;
using namespace WireCell;
namespace Xdata = WireCellXdataRoot;


// shower3D tree Trun
struct Trun {
    int detector, eventNo, runNo, subRunNo;
    float unit_dis, toffset_uv, toffset_uw, toffset_u;
    int total_time_bin, recon_threshold, frame_length, max_events, eve_num, nrebin;
    float threshold_u, threshold_v, threshold_w;
    int time_offset;
};

// shower3D tree TC
struct TC {
    int time_slice, ncluster, mcell_id;
    double charge, xx, yy, zz;
    int face, apa_no, cryostat_no;
    int u_index, v_index, w_index;
    double u_charge, v_charge, w_charge;
    double u_charge_err, v_charge_err, w_charge_err;    
};

uint64_t intpair(int first, int second)
{
    uint64_t f = first;
    uint64_t s = second;
    return (f<<32)|s;
}
std::pair<int,int> intunpair(uint64_t fs)
{
    int first = (fs>>32);
    int second = (0xffffffff&fs);
    return make_pair(first,second);
}

// A functional object which should be used to map all GeomWires to
// Xdata:Wires.  Create it with a CloneHelper attached to the
// TClonesArray for the cells
// fixme: this probably deserves to be moved into data or nav packages.
struct WireMapper {
    unordered_map<int, Xdata::Wire*> wc2xd;
    Xdata::CloneHelper<Xdata::Wire>& helper;    

    WireMapper(Xdata::CloneHelper<Xdata::Wire>& helper) : helper(helper) {}
    Xdata::Wire* operator()(const GeomWire* gwire) {
	auto it = wc2xd.find(gwire->ident());
	if (it != wc2xd.end()) {
	    return it->second;
	}
	Xdata::Wire* wire = helper.make();
	wc2xd[gwire->ident()] = wire;

	wire->ident = gwire->ident();
	wire->chanid = gwire->channel();
	wire->offset = gwire->index();
	wire->plane = (Xdata::planeid_t)gwire->iplane();
	wire->face = Xdata::kFrontFace;
	if (gwire->face()) {
	    wire->face = Xdata::kBackFace;
	}
	wire->segment = gwire->segment();
	Point p1 = gwire->point1();
	wire->point1 = Xdata::Point(p1.x, p1.y, p1.z);
	Point p2 = gwire->point2();
	wire->point2 = Xdata::Point(p2.x, p2.y, p2.z);
	wire->apaid = gwire->cryo()*1000+gwire->apa(); // pick a random convention

	return wire;
    }
};

// A functional object which should be used to map all GeomCells to
// Xdata::Cells.  Create it with a CloneHelper attached to the
// TClonesArray for the cells and a WireMapper to handle internal wire
// mapping.
// fixme: this probably deserves to be moved into data or nav packages.
struct CellMapper {
    unordered_map<int, Xdata::Cell*> wc2xd;
    Xdata::CloneHelper<Xdata::Cell>& helper;    
    WireMapper& wiremapper;
    int tilingid;		// the algorithm that made the cells.
    
    CellMapper(Xdata::CloneHelper<Xdata::Cell>& helper,
	       WireMapper& wiremapper,
	       int tileingid = 1)
	: helper(helper)
	, wiremapper(wiremapper)
	, tilingid(tilingid)
	{  }

  

    Xdata::Cell* operator()(const GeomCell* gcell) {
	Xdata::Wire* uwire = wiremapper(gcell->get_uwire());
	Xdata::Wire* vwire = wiremapper(gcell->get_vwire());
	Xdata::Wire* wwire = wiremapper(gcell->get_wwire());
	Xdata::cellid_t cellid =  Xdata::cell_ident_pack(uwire->ident, vwire->ident, wwire->ident,tilingid);
	auto it = wc2xd.find(cellid);
	if (it != wc2xd.end()) { // already have it
	    return it->second;
	}
	Xdata::Cell* cell = helper.make();
	wc2xd[cellid] = cell;

	cell->ident = cellid;
	cell->area = gcell->cross_section();
	Point p = gcell->center();
	cell->center = Xdata::Point(p.x,p.y,p.z);

	return cell;
    }
};

// A functional object which should be used to fill  merged cell info into to
// Xdata::Blobs.  
// fixme: this probably deserves to be moved into data or nav packages.
struct BlobMapper {
    unordered_map<int, Xdata::Blob*> wc2xd; // blob id

    Xdata::CloneHelper<Xdata::Blob>& helper;
    CellMapper& cellmapper;
    BlobMapper(Xdata::CloneHelper<Xdata::Blob>& helper,
	       CellMapper& cellmapper)
	: helper(helper)
	, cellmapper(cellmapper)
	{  }
    
    Xdata::Blob* operator()(int blob_id, const GeomCell* gcell,
			    int time_slice,
			    float charge) {
	Xdata::Cell* cell = cellmapper(gcell);

	Xdata::Blob* blob = 0;
	auto it = wc2xd.find(blob_id);
	if (it == wc2xd.end()) { // first time
	    blob = helper.make();
	    wc2xd[blob_id] = blob;
	    blob->ident = blob_id;
	    blob->slice = time_slice;
	    blob->charge = charge;
	    blob->cellids.push_back(cell->ident);
	    return blob;
	}
	blob = it->second;
	for (auto maybe_index : blob->cellids) {
	    if (maybe_index == cell->ident) { // already seen this cell
		return blob;
	    }
	}
	if (blob->slice != time_slice) {
	    cerr << "Corrupt blob data: blob ID=" << blob_id << " cell ID " << cell->ident
		 << " time slice mismatch: " << blob->slice << " != " << time_slice << endl;
	    throw runtime_error("data corruption");
	}
	// add
	blob->cellids.push_back(cell->ident);
	blob->charge += charge;
	return blob;
    }
};

struct DecoMapper {

    unordered_map<uint64_t, Xdata::Deco*> wc2xd; // chan+slice -> deco
    Xdata::CloneHelper<Xdata::Deco>& helper;
    WireMapper& wiremapper;

    DecoMapper(Xdata::CloneHelper<Xdata::Deco>& helper,
	       WireMapper& wiremapper)
	: helper(helper)
	, wiremapper(wiremapper)
	{  }

    Xdata::Deco* operator()(const GeomWire* gwire, int time_slice,
			    float charge, float uncertainty) {
	Xdata::Wire* wire = wiremapper(gwire);
	uint64_t decoid = intpair(wire->chanid, time_slice);
	auto it = wc2xd.find(decoid);
	if (it != wc2xd.end()) {
	    return it->second;
	}
	Xdata::Deco* deco = helper.make();
	wc2xd[decoid] = deco;

	deco->chanid = gwire->channel();
	deco->slice = time_slice;
	deco->charge = charge;
	deco->uncertainty = uncertainty;
	return deco;
    }
};


int main(int argc, const char* argv[])
{
    // larsoft conventins
    vector<string> detectors{"uboone","dune35t","protodune","dune10kt_workspace"};
    auto infile = TFile::Open(argv[1]);

    auto runtree = dynamic_cast<TTree*>(infile->Get("Trun"));
    auto TMC = dynamic_cast<TTree*>(infile->Get("TMC"));
    //
    Trun runobj;		// choo choo
    runtree->SetBranchAddress("detector", &runobj.detector);
    runtree->SetBranchAddress("eventNo", &runobj.eventNo);
    runtree->SetBranchAddress("runNo", &runobj.runNo);
    runtree->SetBranchAddress("subRunNo", &runobj.subRunNo);
    runtree->SetBranchAddress("unit_dis", &runobj.unit_dis);
    runtree->SetBranchAddress("toffset_uv", &runobj.toffset_uv);
    runtree->SetBranchAddress("toffset_uw", &runobj.toffset_uw);
    runtree->SetBranchAddress("toffset_u", &runobj.toffset_u);
    runtree->SetBranchAddress("total_time_bin", &runobj.total_time_bin);
    runtree->SetBranchAddress("recon_threshold", &runobj.recon_threshold);
    runtree->SetBranchAddress("frame_length", &runobj.frame_length);
    runtree->SetBranchAddress("max_events", &runobj.max_events);
    runtree->SetBranchAddress("eve_num", &runobj.eve_num);
    runtree->SetBranchAddress("nrebin", &runobj.nrebin);
    runtree->SetBranchAddress("threshold_u", &runobj.threshold_u);
    runtree->SetBranchAddress("threshold_v", &runobj.threshold_v);
    runtree->SetBranchAddress("threshold_w", &runobj.threshold_w);
    runtree->SetBranchAddress("time_offset", &runobj.time_offset);
    runtree->GetEntry(0); // get immediately

    auto ctree = dynamic_cast<TTree*>(infile->Get("TC"));
    TC cobj;
    const GeomCell* geomcell = 0;
    ctree->SetBranchAddress("cell",&geomcell);
    ctree->SetBranchAddress("time_slice", &cobj.time_slice);
    ctree->SetBranchAddress("ncluster", &cobj.ncluster);
    ctree->SetBranchAddress("mcell_id", &cobj.mcell_id);
    // charge in blob reduced by ratio of cell cross section / blob cross section
    ctree->SetBranchAddress("charge", &cobj.charge);
    ctree->SetBranchAddress("xx", &cobj.xx);//cm
    ctree->SetBranchAddress("yy", &cobj.yy);//cm
    ctree->SetBranchAddress("zz", &cobj.xx);//cm
    ctree->SetBranchAddress("face", &cobj.face);
    ctree->SetBranchAddress("apa_no", &cobj.apa_no);
    ctree->SetBranchAddress("cryostat_no", &cobj.cryostat_no);
    ctree->SetBranchAddress("u_index", &cobj.u_index);
    ctree->SetBranchAddress("v_index", &cobj.v_index);
    ctree->SetBranchAddress("w_index", &cobj.w_index);
    ctree->SetBranchAddress("u_charge", &cobj.u_charge);
    ctree->SetBranchAddress("v_charge", &cobj.v_charge);
    ctree->SetBranchAddress("w_charge", &cobj.w_charge);
    ctree->SetBranchAddress("u_charge_err", &cobj.u_charge_err);
    ctree->SetBranchAddress("v_charge_err", &cobj.v_charge_err);
    ctree->SetBranchAddress("w_charge_err", &cobj.w_charge_err);

    
    struct Ttrue {double x,y,z,q; } trueobj;
    auto truetree = dynamic_cast<TTree*>(infile->Get("T_true"));
    truetree->SetBranchAddress("x",&trueobj.x);
    truetree->SetBranchAddress("y",&trueobj.y);
    truetree->SetBranchAddress("z",&trueobj.z);
    truetree->SetBranchAddress("q",&trueobj.q);

    struct Trec {double x,y,z; } recobj;
    auto rectree = dynamic_cast<TTree*>(infile->Get("T_rec"));
    rectree->SetBranchAddress("x",&recobj.x);
    rectree->SetBranchAddress("y",&recobj.y);
    rectree->SetBranchAddress("z",&recobj.z);

    struct Trecq {double x,y,z,q,nq,chi2,ndf; } recqobj;
    auto recqtree = dynamic_cast<TTree*>(infile->Get("T_rec_charge"));
    recqtree->SetBranchAddress("x",&recqobj.x);
    recqtree->SetBranchAddress("y",&recqobj.y);
    recqtree->SetBranchAddress("z",&recqobj.z);
    recqtree->SetBranchAddress("q",&recqobj.q);
    recqtree->SetBranchAddress("nq",&recqobj.nq);
    recqtree->SetBranchAddress("chi2",&recqobj.chi2);
    recqtree->SetBranchAddress("ndf",&recqobj.ndf);


    Xdata::Writer xwriter(argv[2]);
    xwriter.set_tree_mc(TMC);

    Xdata::Geom& geom = xwriter.geom.obj();
    geom.ident = runobj.detector;
    geom.description = detectors[runobj.detector].c_str();
    cerr << "Detector: " << geom.description << endl;
    //geom.wire_bounds = ...;
    //geom.pitches = ...;

    Xdata::CloneHelper<Xdata::Wire> wireca(*geom.wires);
    WireMapper wiremap(wireca);
    Xdata::CloneHelper<Xdata::Cell> cellca(*geom.cells);
    CellMapper cellmap(cellca, wiremap);

    Xdata::RunInfo& runinfo = xwriter.runinfo.obj();
    uint64_t run64 = runobj.runNo;
    uint64_t subrun64 = runobj.subRunNo;
    uint64_t pack = (run64<<32)|subrun64;
    runinfo.ident = pack;
    xwriter.runinfo.fill();
    cerr << "Run="  << (pack>>32) << " subrun:" << (0xffffffff&pack) << endl;    


    // fixme: in general many triggers will be created.
    Xdata::Trigger& trigger = xwriter.trigger.obj();
    trigger.ident = runobj.eventNo;
    trigger.type = 0;		//fixme: this is reason for trigger.
    trigger.runid = runinfo.ident;
    trigger.second = 0;		// fixme: nothing real to set here
    trigger.nanosecond = 0;	// fixme: ibid
    xwriter.trigger.fill();

    // fixme: in general many frames will be created.
    Xdata::Frame& frame = xwriter.frame.obj();
    frame.ident = 1;
    frame.trigid = trigger.ident;
    frame.geomid = geom.ident;
    frame.toffset = 0;		// fixme: I think this needs to actually be set to something
    //frame.slicespan = ...;
    Xdata::CloneHelper<Xdata::Deco> frameca(*frame.decos);
    DecoMapper decomap(frameca, wiremap);

    // fixme: in general there may be multiple images per frame.
    Xdata::Image& image = xwriter.image.obj();
    image.ident = runobj.eve_num;
    image.frameid = frame.ident;
    Xdata::CloneHelper<Xdata::Blob> blobca(*image.blobs);
    BlobMapper blobmap(blobca, cellmap);
    
    std::chrono::duration<double> duration1, duration2;
    chrono::time_point<std::chrono::system_clock> clock1, clock2, clock3;

    auto nentries = ctree->GetEntries();
    cerr << "Loading " << nentries <<  " entries of TC" << endl;
    for (auto entry=0; entry<nentries; ++entry) {
	clock1 = chrono::system_clock::now();

	ctree->GetEntry(entry);

	clock2 = chrono::system_clock::now();

	cellmap(geomcell);
	decomap(geomcell->get_uwire(), cobj.time_slice, cobj.u_charge, cobj.u_charge_err);
	decomap(geomcell->get_vwire(), cobj.time_slice, cobj.v_charge, cobj.v_charge_err);
	decomap(geomcell->get_wwire(), cobj.time_slice, cobj.w_charge, cobj.w_charge_err);
	// note: sums charge
	blobmap(cobj.mcell_id, geomcell, cobj.time_slice, cobj.charge);

	clock3 = chrono::system_clock::now();
	duration2 += clock3 - clock2;
	duration1 += clock2 - clock1;
    }
    xwriter.geom.fill();
    xwriter.frame.fill();
    xwriter.image.fill();


    Xdata::Field& field = xwriter.field.obj();
    Xdata::CloneHelper<Xdata::FieldPoint> fieldca(*field.points);

    field.ident=1;
    field.trigid = trigger.ident;
    field.geomid = geom.ident;
    field.name = "true charge deposition";
    fieldca.clear();
    nentries = truetree->GetEntries();
    cerr << "Loading " << nentries <<  " entries of T_true" << endl;
    for (auto entry=0; entry<nentries; ++entry) {
	clock1 = chrono::system_clock::now();
	truetree->GetEntry(entry);
	clock2 = chrono::system_clock::now();

	Xdata::FieldPoint* fp = fieldca.make();
	fp->point = Xdata::Point(trueobj.x, trueobj.y, trueobj.z);
	fp->values.push_back(trueobj.q);

	clock3 = chrono::system_clock::now();
	duration1 += clock2 - clock1;
	duration2 += clock3 - clock2;
    }
    xwriter.field.fill();

    field.ident=2;
    field.trigid = trigger.ident;
    field.name = "no charge image";
    fieldca.clear();
    nentries = rectree->GetEntries();
    cerr << "Loading " << nentries <<  " entries of T_rec" << endl;
    for (auto entry=0; entry<nentries; ++entry) {
	clock1 = chrono::system_clock::now();
	rectree->GetEntry(entry);
	clock2 = chrono::system_clock::now();

	Xdata::FieldPoint* fp = fieldca.make();
	fp->point = Xdata::Point(recobj.x, recobj.y, recobj.z);
	//fp->values.push_back(...);//none...

	clock3 = chrono::system_clock::now();
	duration1 += clock2 - clock1;
	duration2 += clock3 - clock2;
    }
    xwriter.field.fill();

    field.ident=3;
    field.trigid = trigger.ident;
    field.name = "charge image";
    fieldca.clear();
    nentries = recqtree->GetEntries();
    cerr << "Loading " << nentries <<  " entries of T_rec_charge" << endl;
    for (auto entry=0; entry<nentries; ++entry) {
	clock1 = chrono::system_clock::now();
	recqtree->GetEntry(entry);
	clock2 = chrono::system_clock::now();

	Xdata::FieldPoint* fp = fieldca.make();
	fp->point = Xdata::Point(recqobj.x, recqobj.y, recqobj.z);
	fp->values.push_back(recqobj.q);
	fp->values.push_back(recqobj.nq);
	fp->values.push_back(recqobj.chi2);
	fp->values.push_back(recqobj.ndf);

	clock3 = chrono::system_clock::now();
	duration1 += clock2 - clock1;
	duration2 += clock3 - clock2;
    }
    xwriter.field.fill();

    cerr << "shower3d read time: " << duration1.count() << "s\n";
    cerr << "xdata packing time: " << duration2.count() << "s\n";

    cerr << "Writing xdata:\n"
	 << "\t" << wireca.size() <<  " wires\n"
	 << "\t" << cellca.size() <<  " cells\n"
	 << "\t" << frameca.size() <<  " decos\n"
	 << "\t" << blobca.size() <<  " blobs\n"
	 << "\t" << fieldca.size() <<  " fields\n";

    clock1 = chrono::system_clock::now();
    xwriter.close();
    clock2 = chrono::system_clock::now();
    duration2 = clock2 - clock1;
    cerr << "Write xdata file in " << duration2.count() << "s\n";

    

    {				// do a readback test....
    	cerr << "Test readback:" << endl;
    	clock1 = chrono::system_clock::now();
	Xdata::Reader xreader(argv[2]);
	while (xreader.runinfo_reader.Next()) { }
	Xdata::WireDB wdb;
	while (xreader.geom_reader.Next()) {
	    Xdata::Geom& g = *xreader.geom;
	    Xdata::CloneHelper<Xdata::Wire> wireca(*g.wires);
	    for (int ind=0; ind < wireca.size(); ++ind) {
		wdb(wireca.get(ind));
	    }
	}
	while (xreader.trigger_reader.Next()) {}
	while (xreader.frame_reader.Next()) { }
	while (xreader.image_reader.Next()) { }
	while (xreader.field_reader.Next()) { }

    	clock2 = chrono::system_clock::now();
    	duration1 = clock2 - clock1;
    	cerr << "... in " << duration1.count() << "s\n";
    }    
    return 0;

}
    
