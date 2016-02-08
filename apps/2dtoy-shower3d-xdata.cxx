/// Convert "shower3D" file to a Wire Cell Exchange Data ROOT format.

#include "WireCellData/GeomCell.h"

#include "WireCellXdataRoot/XdataFile.h"

#include "TFile.h"
#include "TTree.h"

#include <map>
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
   
// grungy reverse mapper from ident to object index in a vector.
struct IdentMapper {
    map<uint64_t,int> wire_id2ind;
    map<uint64_t,int> cell_id2ind; // id is packed wires indices + "context"
    map<uint64_t,int> blob_id2ind;
    map<uint64_t,int> deco_id2ind; // id here is a wire index + timeslice
    
    Xdata::XdataFile& xdata;
    uint64_t context;
    IdentMapper(Xdata::XdataFile& df, uint16_t detector)
	: xdata(df), context(detector) {}


    int fill_wire(const GeomWire& gwire) {
	int ident = gwire.ident();
	auto it = wire_id2ind.find(ident);
	if (it != wire_id2ind.end()) {
	    return it->second;
	}

	Xdata::Wire* wire = new Xdata::Wire;
	wire->wireid = ident;
	wire->chanid = gwire.channel();
	wire->volid = gwire.cryo()*1000+gwire.apa(); // pick a random convention
	wire->plane = gwire.iplane();
	wire->face = gwire.face();
	wire->offset = gwire.index();
	wire->segment = gwire.segment();
	Point p1 = gwire.point1();
	wire->point1 = Xdata::Point(p1.x, p1.y, p1.z);
	Point p2 = gwire.point2();
	wire->point2 = Xdata::Point(p2.x, p2.y, p2.z);

	int ind = xdata.geom().wires.size();
	wire_id2ind[ident] = ind;
	xdata.geom().wires.push_back(wire);
	return ind;
    }

    // Some sugar to return the packed cell ID
    uint64_t cell_ident(const GeomCell& gcell) {
	return Xdata::Cell::ident_pack(fill_wire(*gcell.get_uwire()),
				       fill_wire(*gcell.get_vwire()),
				       fill_wire(*gcell.get_wwire()),
				       context);
    }

    int fill_cell(const GeomCell& gcell) {
	//int ident = gcell.ident();
	uint64_t ident = cell_ident(gcell);
	auto it = cell_id2ind.find(ident);
	if (it != cell_id2ind.end()) {
	    return it->second;
	}
	int ind=-1;
	Xdata::Cell* cell = xdata.image().new_cell(ind);
	cell_id2ind[ident] = ind;
	cell->ident = ident;
	cell->area = gcell.cross_section();
	Point p = gcell.center();
	cell->center = Xdata::Point(p.x,p.y,p.z);
	return ind;
    }
	
    int fill_deco(const GeomWire& gwire, int time_slice,
		  const vector<float>& values) {
	int wind = fill_wire(gwire);
	uint64_t decoid = intpair(wind, time_slice);
	auto it = deco_id2ind.find(decoid);
	if (it != deco_id2ind.end()) {
	    return it->second;
	}
	int ind=-1;
	Xdata::Deco* deco = xdata.image().new_deco(ind);
	deco_id2ind[decoid] = ind;
	deco->wireind = wind;
	deco->slice = time_slice;
	deco->values = values;
	return ind;
    }


    // call this after adding cell
    int fill_blob(int blob_ident, int cell_ident, int time_slice,
		  const vector<float>& values) {
	int cell_index = cell_id2ind[cell_ident];
	auto it = blob_id2ind.find(blob_ident);

	// no such blob yet
	if (it == blob_id2ind.end()) { 
	    int ind = -1;
	    Xdata::Blob* blob = xdata.image().new_blob(ind);
	    blob_id2ind[blob_ident] = ind;
	    blob->ident = blob_ident;
	    blob->slice = time_slice;
	    blob->cellind.push_back(cell_index);
	    blob->values = values;
	    return ind;
	}
	int ind = it->second;
	Xdata::Blob* blob = xdata.image().get_blob(ind);
	for (auto maybe_index : blob->cellind) {
	    if (maybe_index == cell_index) { // already seen this cell
		return ind;
	    }
	}
	blob->cellind.push_back(cell_index);
	for (int ival=0; ival<values.size(); ++ival) {
	    blob->values[ival] += values[ival];
	}
	return ind;	
    }
};

int main(int argc, const char* argv[])
{
    // larsoft conventins
    vector<string> detectors{"uboone","dune35t","protodune","dune10kt_workspace"};
    auto infile = TFile::Open(argv[1]);

    auto runtree = dynamic_cast<TTree*>(infile->Get("Trun"));
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
    runtree->GetEntry(0);

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


    Xdata::XdataFile xdatafile;

    Xdata::RunInfo& ri = xdatafile.runinfo();
    uint64_t run64 = runobj.runNo;
    uint64_t subrun64 = runobj.subRunNo;
    uint64_t pack = (run64<<32)|subrun64;

    ri.ident = pack;
    cerr << "Run="  << (pack>>32) << " subrun:" << (0xffffffff&pack) << endl;    
    ri.detector = detectors[runobj.detector].c_str();
    cerr << "Detector: " << ri.detector << endl;


    IdentMapper idmap(xdatafile, runobj.detector);

    Xdata::Image& img = xdatafile.image();
    img.ident = runobj.eve_num;
    img.slicespan = runobj.nrebin*0.5;

    
    std::chrono::duration<double> in_duration;
    std::chrono::duration<double> out_duration;

    chrono::time_point<std::chrono::system_clock> clock1, clock2, clock3;

    auto nentries = ctree->GetEntries();
    cerr << "Loading " << nentries <<  " entries of TC" << endl;
    for (auto entry=0; entry<nentries; ++entry) {
	clock1 = chrono::system_clock::now();
	ctree->GetEntry(entry);
	clock2 = chrono::system_clock::now();

	idmap.fill_cell(*geomcell);

	// fixme: is narrowing to float okay?
	idmap.fill_deco(*geomcell->get_uwire(), cobj.time_slice,
			vector<float>{(float)cobj.u_charge, (float)cobj.u_charge_err});
	idmap.fill_deco(*geomcell->get_vwire(), cobj.time_slice,
			vector<float>{(float)cobj.v_charge, (float)cobj.v_charge_err});
	idmap.fill_deco(*geomcell->get_wwire(), cobj.time_slice,
			vector<float>{(float)cobj.w_charge, (float)cobj.w_charge_err});

	// reconstruct blob charge density due to this cell
	float cell_area = geomcell->cross_section();
	float charge_density = cobj.charge / cell_area;
	vector<float> values{cell_area, charge_density};

	idmap.fill_blob(cobj.mcell_id, geomcell->ident(), cobj.time_slice, values);

	clock3 = chrono::system_clock::now();

	in_duration += clock2 - clock1;
	out_duration += clock3 - clock2;
    }

    int dummy;			// don't care

    // T_true
    Xdata::Field* field_true = xdatafile.image().new_field(dummy);
    field_true->name = "true charge deposition";
    nentries = truetree->GetEntries();
    cerr << "Loading " << nentries <<  " entries of T_true" << endl;
    for (auto entry=0; entry<nentries; ++entry) {
	clock1 = chrono::system_clock::now();
	truetree->GetEntry(entry);
	clock2 = chrono::system_clock::now();

	Xdata::FieldPoint fp(Xdata::Point(trueobj.x, trueobj.y, trueobj.z),
			     vector<float>{(float)trueobj.q});
	field_true->points.push_back(fp);

	clock3 = chrono::system_clock::now();
	in_duration += clock2 - clock1;
	out_duration += clock3 - clock2;
    }


    Xdata::Field* field_rec = xdatafile.image().new_field(dummy);
    field_rec->name = "no charge image";
    nentries = rectree->GetEntries();
    cerr << "Loading " << nentries <<  " entries of T_rec" << endl;
    for (auto entry=0; entry<nentries; ++entry) {
	clock1 = chrono::system_clock::now();
	rectree->GetEntry(entry);
	clock2 = chrono::system_clock::now();

	Xdata::FieldPoint fp(Xdata::Point(recobj.x, recobj.y, recobj.z),
			     vector<float>{});
	field_rec->points.push_back(fp);

	clock3 = chrono::system_clock::now();
	in_duration += clock2 - clock1;
	out_duration += clock3 - clock2;
    }


    Xdata::Field* field_recq = xdatafile.image().new_field(dummy);
    field_recq->name = "charge image";
    nentries = recqtree->GetEntries();
    cerr << "Loading " << nentries <<  " entries of T_rec_charge" << endl;
    for (auto entry=0; entry<nentries; ++entry) {
	clock1 = chrono::system_clock::now();
	recqtree->GetEntry(entry);
	clock2 = chrono::system_clock::now();

	vector<float> values{
	    (float)recqobj.q,(float)recqobj.nq,(float)recqobj.chi2,(float)recqobj.ndf
		};
	Xdata::FieldPoint fp(Xdata::Point(recqobj.x, recqobj.y, recqobj.z), values);
	field_recq->points.push_back(fp);

	clock3 = chrono::system_clock::now();
	in_duration += clock2 - clock1;
	out_duration += clock3 - clock2;
    }

    
    cerr << "shower3d read time: " << in_duration.count() << "s\n";
    cerr << "xdata packing time: " << out_duration.count() << "s\n";

    cerr << "Writing xdata:\n"
	 << "\t" << xdatafile.geom().wires.size() <<  " wires\n"
	 << "\t" << xdatafile.image().num_decos() <<  " decos\n"
	 << "\t" << xdatafile.image().num_cells() <<  " cells\n"
	 << "\t" << xdatafile.image().num_blobs() <<  " blobs\n"
	 << "\t" << xdatafile.image().num_fields() <<  " fields\n";

    clock1 = chrono::system_clock::now();
    auto siz = xdatafile.write(argv[2]);
    clock2 = chrono::system_clock::now();
    out_duration = clock2 - clock1;
    cerr << "Write "<<siz<<" bytes xdata in " << out_duration.count() << "s\n";

    
    {				// do a readback test....
	cerr << "Test readback:" << endl;
	clock1 = chrono::system_clock::now();
	Xdata::XdataFile readback;
	readback.read(argv[2]);
	clock2 = chrono::system_clock::now();
	in_duration = clock2 - clock1;
	cerr << "... in " << in_duration.count() << "s\n";

	Xdata::RunInfo& ri = readback.runinfo();
	uint64_t rs = ri.ident;
	uint64_t run64 = (rs>>32);
	uint64_t subrun64 = (rs&0xFFFFFFFF);
	cerr << "Run="  << run64 << " subrun:" << subrun64 << endl;
	cerr << "Detector: " << ri.detector << endl;
	cerr << "Read xdata:\n"
	     << "\t" << readback.geom().wires.size() <<  " wires\n"
	     << "\t" << readback.image().num_decos() <<  " decos\n"
	     << "\t" << readback.image().num_cells() <<  " cells\n"
	     << "\t" << readback.image().num_blobs() <<  " blobs\n"
	     << "\t" << readback.image().num_fields() <<  " fields\n";

	if (ri.detector == "") {
	    cerr << "ERROR: failed to read back detector name" << endl;
	}
    }    
    return 0;

}
    
