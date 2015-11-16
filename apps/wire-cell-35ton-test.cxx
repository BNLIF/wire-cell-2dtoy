#include "WireCellNav/DetectorGDS.h"
#include "WireCellNav/DetGenerativeFDS.h"
#include "WireCellNav/FrameDataSource.h"
#include "WireCell2dToy/ToyDepositor.h"
#include "WireCellSst/Util.h"
#include "WireCell2dToy/ToySignalSimu.h"
#include "WireCell2dToy/ToySignalSimuTrue.h"
#include "WireCell2dToy/ToySignalGaus.h"
#include "WireCell2dToy/ToySignalWien.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCellData/MergeGeomCell.h"

#include "WireCell2dToy/ToyEventDisplay.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TBenchmark.h"
#include "TH1F.h"
#include "TLine.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TStyle.h"

#include <iostream>
#include <string>
#include <vector>


using namespace WireCell;
using namespace std;

int main(int argc, char* argv[])
{

  if (argc < 3){
    cerr << "usage: wire-cell-35ton-test /path/to/celltree.root eve_num" << endl;
    return 1;
  }
  
  //build GDS ... 
  DetectorGDS gds;
  gds.set_ncryos(1);
  gds.set_napas(0,4);
  Vector center0(-4.0122, 15.3431, 24.6852);
  Vector halves0(3.26512, 99.7439, 26.7233);
  gds.set_apa(0, 0, 45.705, 44.274, 0.4880488, 0.4880488, 0.4880488, center0, halves0);
  Vector center1(-4.0122, -42.2348, 77.3702);
  Vector halves1(3.26512, 42.2504, 25.9617);
  gds.set_apa(0, 1, 45.705, 44.274, 0.4880488, 0.4880488, 0.4880488, center1, halves1);
  Vector center2(-4.0122, 57.5435, 77.3702);
  Vector halves2(3.26512, 57.5435, 25.9617);
  gds.set_apa(0, 2, 45.705, 44.274, 0.4880488, 0.4880488, 0.4880488, center2, halves2);
  Vector center3(-4.0122, 15.3431, 130.055);
  Vector halves3(3.26512, 99.7439, 26.7235);
  gds.set_apa(0, 3, 45.705, 44.274, 0.4880488, 0.4880488, 0.4880488, center3, halves3);
  gds.buildGDS();

  
  const char* root_file = argv[1];
  const char* tpath = "/Event/Sim";

  int recon_threshold = 2000;
  int max_events = 100;
  int eve_num = atoi(argv[2]);
  float unit_dis = 1.6;  
  float toffset_1=1.647;
  float toffset_2=1.539+1.647;
  float toffset_3=0;

  
  int total_time_bin=9600;
  //  int frame_length = 3200;
  int frame_length = 800;  // hack for now
  int nrebin = 4;

  TFile *tfile = TFile::Open(root_file);
  TTree* sst = dynamic_cast<TTree*>(tfile->Get(tpath));

  int run_no, subrun_no, event_no;
  sst->SetBranchAddress("eventNo",&event_no);
  sst->SetBranchAddress("runNo",&run_no);
  sst->SetBranchAddress("subRunNo",&subrun_no);
  sst->GetEntry(eve_num);

  cout << "Run No: " << run_no << " " << subrun_no << " " << eve_num << endl;
  
  WireCell::FrameDataSource* fds = 0;
  fds = WireCellSst::make_fds(*tfile);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  
  TH1::AddDirectory(kFALSE);
  
  WireCell::ToyDepositor toydep(fds,0,unit_dis,frame_length);
  const PointValueVector& pvv = toydep.depositions(eve_num);

  std::cout << "Points deposited: " << pvv.size() << std::endl;

  DetGenerativeFDS gfds(toydep,gds, 2400,max_events,2.0*1.6*units::millimeter);
  gfds.jump(eve_num);

   WireCellSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons

  // DetGenerativeFDS gfds(toydep, gds,total_time_bin,max_events,0.5*unit_dis*units::millimeter);
  // //gfds.jump(eve_num);
  
  // cout << "Put in Truth " << endl; 
  // WireCell2dToy::ToySignalSimuTrueFDS st_fds(gfds,gds,total_time_bin/nrebin,max_events,0); //truth
  // st_fds.jump(eve_num);
  
  // cout << "Simulate Raw WaveForm " << endl; 
  // WireCell2dToy::ToySignalSimuFDS simu_fds(gfds,gds,total_time_bin,max_events,toffset_1,toffset_2,1); // time offset among different planes for the time electrons travel among different planes
  // simu_fds.jump(eve_num);
  // //simu_fds.Save();
  
  // cout << "Deconvolution with Gaussian filter" << endl;
  // WireCell2dToy::ToySignalGausFDS gaus_fds(simu_fds,gds,total_time_bin/nrebin,max_events,toffset_1,toffset_2); // gaussian smearing for charge estimation
  // gaus_fds.jump(eve_num);
  // //gaus_fds.Save();

  // cout << "Deconvolution with Wiener filter" << endl;
  // WireCell2dToy::ToySignalWienFDS wien_fds(simu_fds,gds,total_time_bin/nrebin,max_events,toffset_1,toffset_2); // weiner smearing for hit identification
  // wien_fds.jump(eve_num);
  // // //wien_fds.Save();


   WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
   WireCell2dToy::MergeToyTiling **mergetiling = new WireCell2dToy::MergeToyTiling*[2400];
   WireCell2dToy::TruthToyTiling **truthtiling = new WireCell2dToy::TruthToyTiling*[2400];
  
   //for (int i=0;i!=2400;i++)
   int i = 1137;
   {
     sds.jump(i);
     WireCell::Slice slice = sds.get();
     if ( slice.group().size() >0){
       cout << i << " " << slice.group().size() << endl;
       toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds);
       //GeomCellSelection allcell = toytiling[i]->get_allcell();
     }
     
     TApplication theApp("theApp",&argc,argv);
    theApp.SetReturnFromRun(true);
    
    TCanvas c1("ToyMC","ToyMC",1200,600);
    c1.Divide(2,1);
    c1.Draw();
    
    float charge_min = 0;
    float charge_max = 1e5;


    WireCell2dToy::ToyEventDisplay display(c1, gds);
    display.charge_min = charge_min;
    display.charge_max = charge_max;


    gStyle->SetOptStat(0);
    
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Int_t MyPalette[NCont];
    Double_t stops[NRGBs] = {0.0, 0.34, 0.61, 0.84, 1.0};
    Double_t red[NRGBs] = {0.0, 0.0, 0.87 ,1.0, 0.51};
    Double_t green[NRGBs] = {0.0, 0.81, 1.0, 0.2 ,0.0};
    Double_t blue[NRGBs] = {0.51, 1.0, 0.12, 0.0, 0.0};
    Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    for (int kk=0;kk!=NCont;kk++) MyPalette[kk] = FI+kk;
    gStyle->SetPalette(NCont,MyPalette);
    
    display.init(-0.03,1.568,-0.845,1.151);
    display.draw_mc(1,WireCell::PointValueVector(),"colz");
    display.draw_slice(slice,"");
  
    theApp.Run();
   }

    
  // TCanvas *c = new TCanvas();
  // c->Range(-5*units::cm, -90*units::cm, 160*units::cm, 120*units::cm);    
  
  // TLine *l = new TLine();
  // l->SetLineWidth(0);
  // c->cd();
  
  
  
  // int colors[] = {2,4,1};
  // for (short cryo = 0; cryo < gds.ncryos(); cryo++) {
  //   for (short apa = 0; apa < gds.napa(cryo); apa++) {
  // 	std::cout << cryo << " " << apa << std::endl;
  
  // 	    const WrappedGDS *apa_gds = gds.get_apaGDS(cryo, apa);
  // 	    for (int iplane=0; iplane<3; ++iplane) {
  // 	        WirePlaneType_t plane = (WirePlaneType_t)iplane;
  // 		GeomWireSelection wip = apa_gds->wires_in_plane(plane);
  // 	        // std::cout<<"\n[CRYO] "<<cryo<<" [APA] "<<apa<<" [PLANE] "<<iplane
  // 		// 	 <<" has "<< wip.size()<<" wires, wire angle is "<<apa_gds->angle(plane)*180/TMath::Pi()<<std::endl;
  // 		//for (auto wit = wip.begin(); wit != wip.end(); ++wit) {
  // 		//  const GeomWire& wire = **wit;
  // 		for (int index=0; index<(int)wip.size(); ++index) {
  // 		    const GeomWire* wire = apa_gds->by_planeindex(plane, index);
  // 		    //if (wire.face() == 0) continue;
  // 		    const Vector& p1 = wire->point1();
  // 		    const Vector& p2 = wire->point2();
  // 		    //  std::cout<<*wire<<" ("<<p1.x<<","<<p1.y<<","<<p1.z<<") ("<<p2.x<<","<<p2.y<<","<<p2.z<<")\n";
  // 		    l->SetLineColor(colors[iplane]);
  // 		    l->DrawLine(p1.z, p1.y, p2.z, p2.y);
  // 		}
  // 	    }	    
  // 	}
  // }
  
  // c->SaveAs("./test_detectorgds_35t.pdf");
}
