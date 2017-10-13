#include "WireCellSst/GeomDataSource.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"
#include "WireCell2dToy/ExecMon.h"


#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace WireCell;
using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root " << endl;
    return 1;
  }

  ExecMon em("starting");
  cerr << em("load geometry") << endl;

  WireCellSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;

  cout << "Pitch: " << gds.pitch(WirePlaneType_t(0)) 
       << " " << gds.pitch(WirePlaneType_t(1)) 
       << " " << gds.pitch(WirePlaneType_t(2))
       << endl;
  cout << "Angle: " << gds.angle(WirePlaneType_t(0)) 
       << " " << gds.angle(WirePlaneType_t(1)) 
       << " " << gds.angle(WirePlaneType_t(2))
       << endl;
  TString filename = argv[2];
  TFile *file = new TFile(filename);
  TTree *Trun = (TTree*)file->Get("Trun");

  int run_no, subrun_no, event_no;
  int time_offset;
  int nrebin;
  int frame_length;
  int eve_num;
  float unit_dis;
  
  Trun->SetBranchAddress("eventNo",&event_no);
  Trun->SetBranchAddress("runNo",&run_no);
  Trun->SetBranchAddress("subRunNo",&subrun_no);
  Trun->SetBranchAddress("unit_dis",&unit_dis);
  Trun->SetBranchAddress("frame_length",&frame_length);
  Trun->SetBranchAddress("eve_num",&eve_num);
  Trun->SetBranchAddress("nrebin",&nrebin);
  Trun->SetBranchAddress("time_offset",&time_offset);
  Trun->GetEntry(0);

  // define singleton ... 
  TPCParams& mp = Singleton<TPCParams>::Instance();
  
  double pitch_u = gds.pitch(WirePlaneType_t(0));
  double pitch_v = gds.pitch(WirePlaneType_t(1));
  double pitch_w = gds.pitch(WirePlaneType_t(2));
  double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;

  mp.set_pitch_u(pitch_u);
  mp.set_pitch_v(pitch_v);
  mp.set_pitch_w(pitch_w);
  mp.set_ts_width(time_slice_width);

  // load mcell
  
  TTree *TC = (TTree*)file->Get("TC");
  Int_t cluster_id;
  Int_t time_slice;
  Double_t q, uq, vq, wq, udq, vdq, wdq;
  Int_t nwire_u=0, flag_u; //number of wires, dead?
  Int_t nwire_v=0, flag_v;
  Int_t nwire_w=0, flag_w;
  Int_t wire_index_u[2400];
  Int_t wire_index_v[2400];
  Int_t wire_index_w[3256];
  Double_t wire_charge_u[2400];
  Double_t wire_charge_v[2400];
  Double_t wire_charge_w[2400];
  Double_t wire_charge_err_u[2400];
  Double_t wire_charge_err_v[2400];
  Double_t wire_charge_err_w[2400];

  TC->SetBranchAddress("cluster_id",&cluster_id);
  TC->SetBranchAddress("time_slice",&time_slice);
  TC->SetBranchAddress("q",&q);
  TC->SetBranchAddress("uq",&uq);
  TC->SetBranchAddress("vq",&vq);
  TC->SetBranchAddress("wq",&wq);
  TC->SetBranchAddress("udq",&udq);
  TC->SetBranchAddress("vdq",&vdq);
  TC->SetBranchAddress("wdq",&wdq);
  TC->SetBranchAddress("nwire_u",&nwire_u);
  TC->SetBranchAddress("nwire_v",&nwire_v);
  TC->SetBranchAddress("nwire_w",&nwire_w);
  TC->SetBranchAddress("flag_u",&flag_u);
  TC->SetBranchAddress("flag_v",&flag_v);
  TC->SetBranchAddress("flag_w",&flag_w);
  TC->SetBranchAddress("wire_index_u",wire_index_u);
  TC->SetBranchAddress("wire_index_v",wire_index_v);
  TC->SetBranchAddress("wire_index_w",wire_index_w);
  TC->SetBranchAddress("wire_charge_u",wire_charge_u);
  TC->SetBranchAddress("wire_charge_v",wire_charge_v);
  TC->SetBranchAddress("wire_charge_w",wire_charge_w);
  TC->SetBranchAddress("wire_charge_err_u",wire_charge_err_u);
  TC->SetBranchAddress("wire_charge_err_v",wire_charge_err_v);
  TC->SetBranchAddress("wire_charge_err_w",wire_charge_err_w);


  //load mcell
  
  TTree *TDC = (TTree*)file->Get("TDC");
  int ntime_slice = 0,time_slices[2400];
  TDC->SetBranchAddress("cluster_id",&cluster_id);
  TDC->SetBranchAddress("ntime_slice",&ntime_slice);
  TDC->SetBranchAddress("time_slice",time_slices);
  
  TDC->SetBranchAddress("nwire_u",&nwire_u);
  TDC->SetBranchAddress("nwire_v",&nwire_v);
  TDC->SetBranchAddress("nwire_w",&nwire_w);
  TDC->SetBranchAddress("flag_u",&flag_u);
  TDC->SetBranchAddress("flag_v",&flag_v);
  TDC->SetBranchAddress("flag_w",&flag_w);
  TDC->SetBranchAddress("wire_index_u",wire_index_u);
  TDC->SetBranchAddress("wire_index_v",wire_index_v);
  TDC->SetBranchAddress("wire_index_w",wire_index_w);
  
}
