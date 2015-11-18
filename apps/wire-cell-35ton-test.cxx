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
  Vector center0(-4.0122*units::cm, 15.3431*units::cm, 24.6852*units::cm);
  Vector halves0(3.26512*units::cm, 99.7439*units::cm, 26.7233*units::cm);
  gds.set_apa(0, 0, 45.705*units::deg, 44.274*units::deg, 0.4880488*units::cm, 0.4880488*units::cm, 0.4880488*units::cm, center0, halves0);
  Vector center1(-4.0122*units::cm, -42.2348*units::cm, 77.3702*units::cm);
  Vector halves1(3.26512*units::cm, 42.2504*units::cm, 25.9617*units::cm);
  gds.set_apa(0, 1, 45.705*units::deg, 44.274*units::deg, 0.4880488*units::cm, 0.4880488*units::cm, 0.4880488*units::cm, center1, halves1);
  Vector center2(-4.0122*units::cm, 57.5435*units::cm, 77.3702*units::cm);
  Vector halves2(3.26512*units::cm, 57.5435*units::cm, 25.9617*units::cm);
  gds.set_apa(0, 2, 45.705*units::deg, 44.274*units::deg, 0.4880488*units::cm, 0.4880488*units::cm, 0.4880488*units::cm, center2, halves2);
  Vector center3(-4.0122*units::cm, 15.3431*units::cm, 130.055*units::cm);
  Vector halves3(3.26512*units::cm, 99.7439*units::cm, 26.7235*units::cm);
  gds.set_apa(0, 3, 45.705*units::deg, 44.274*units::deg, 0.4880488*units::cm, 0.4880488*units::cm, 0.4880488*units::cm, center3, halves3);
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

   float threshold_u = 5.87819e+02 * 4.0;
  float threshold_v = 8.36644e+02 * 4.0;
  float threshold_w = 5.67974e+02 * 4.0;

  float threshold_ug = 755.96;
  float threshold_vg = 822.81;
  float threshold_wg = 510.84;
  
  int time_offset = 0;



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

   int ncount = 0;
  int ncount1 = 0;  
  int ncount2 = 0;

  int ncount_t = 0;


   WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
   WireCell2dToy::MergeToyTiling **mergetiling = new WireCell2dToy::MergeToyTiling*[2400];
   WireCell2dToy::TruthToyTiling **truthtiling = new WireCell2dToy::TruthToyTiling*[2400];


   int start_num = 0 ;
  int end_num = sds.size()-1;
  
   //for (int i=0;i!=2400;i++)
  //int i = 317+800;{
   //int i = 292+800;
  for (int i=start_num;i!=end_num+1;i++){
    sds.jump(i);
     WireCell::Slice slice = sds.get();
     
     // if ( slice.group().size() >0){
     cout << i << " " << slice.group().size() << endl;
     toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds);
     //allcell = toytiling[i]->get_allcell();
     
     truthtiling[i] = new WireCell2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds,frame_length);
     
       //}
     
    //  TApplication theApp("theApp",&argc,argv);
    // theApp.SetReturnFromRun(true);
    
    // TCanvas c1("ToyMC","ToyMC",1200,600);
    // c1.Divide(2,1);
    // c1.Draw();
    
    // float charge_min = 0;
    // float charge_max = 1e5;


    // WireCell2dToy::ToyEventDisplay display(c1, gds);
    // display.charge_min = charge_min;
    // display.charge_max = charge_max;


    // gStyle->SetOptStat(0);
    
    // const Int_t NRGBs = 5;
    // const Int_t NCont = 255;
    // Int_t MyPalette[NCont];
    // Double_t stops[NRGBs] = {0.0, 0.34, 0.61, 0.84, 1.0};
    // Double_t red[NRGBs] = {0.0, 0.0, 0.87 ,1.0, 0.51};
    // Double_t green[NRGBs] = {0.0, 0.81, 1.0, 0.2 ,0.0};
    // Double_t blue[NRGBs] = {0.51, 1.0, 0.12, 0.0, 0.0};
    // Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    // gStyle->SetNumberContours(NCont);
    // for (int kk=0;kk!=NCont;kk++) MyPalette[kk] = FI+kk;
    // gStyle->SetPalette(NCont,MyPalette);
    
    // display.init(-0.03,1.568,-0.845,1.151);
    // display.draw_mc(1,WireCell::PointValueVector(),"colz");
    // display.draw_slice(slice,"");
    // display.draw_cells(toytiling[i]->get_allcell(),"*same");

    // CellChargeMap ccmap = truthtiling[i]->ccmap();
    
    // //std::cout << ccmap.size() << std::endl;
    // display.draw_truthcells(ccmap,"*same");

    // theApp.Run();
   }

    
  TFile *file = new TFile(Form("shower3D_signal_%d.root",eve_num),"RECREATE");
  TTree *t_true = new TTree("T_true","T_true");
  TTree *t_rec = new TTree("T_rec","T_rec");
  TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
  TTree *t_rec_charge_blob = new TTree("T_rec_charge_blob","T_rec_charge_blob");

  Double_t x_save, y_save, z_save;
  Double_t charge_save;
  Double_t ncharge_save;
  Double_t chi2_save;
  Double_t ndf_save;

  t_true->SetDirectory(file);
  t_true->Branch("x",&x_save,"x/D");
  t_true->Branch("y",&y_save,"y/D");
  t_true->Branch("z",&z_save,"z/D");
  t_true->Branch("q",&charge_save,"q/D");
  
  t_rec->SetDirectory(file);
  t_rec->Branch("x",&x_save,"x/D");
  t_rec->Branch("y",&y_save,"y/D");
  t_rec->Branch("z",&z_save,"z/D");
  
  t_rec_charge->SetDirectory(file);
  t_rec_charge->Branch("x",&x_save,"x/D");
  t_rec_charge->Branch("y",&y_save,"y/D");
  t_rec_charge->Branch("z",&z_save,"z/D");
  t_rec_charge->Branch("q",&charge_save,"q/D");
  t_rec_charge->Branch("nq",&ncharge_save,"nq/D");
  t_rec_charge->Branch("chi2",&chi2_save,"chi2/D");
  t_rec_charge->Branch("ndf",&ndf_save,"ndf/D");

  //blob stuff
  t_rec_charge_blob->SetDirectory(file);
  t_rec_charge_blob->Branch("x",&x_save,"x/D");
  t_rec_charge_blob->Branch("y",&y_save,"y/D");
  t_rec_charge_blob->Branch("z",&z_save,"z/D");
  t_rec_charge_blob->Branch("q",&charge_save,"q/D");
  t_rec_charge_blob->Branch("nq",&ncharge_save,"nq/D");
  
  TGraph2D *g = new TGraph2D();
  TGraph2D *gt = new TGraph2D();
  TGraph2D *g_rec = new TGraph2D();
  TGraph2D *g_rec_blob = new TGraph2D();

   //save results 
  for (int i=start_num;i!=end_num+1;i++){
    //truth
    CellChargeMap ccmap = truthtiling[i]->ccmap();
    for (auto it = ccmap.begin();it!=ccmap.end(); it++){
      Point p = it->first->center();
      x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.*4; // *4 is temporary
      y_save = p.y/units::cm;
      z_save = p.z/units::cm;
      charge_save = it->second;
      
      gt->SetPoint(ncount_t,x_save,y_save,z_save);
      t_true->Fill();
            
      ncount_t ++;
    }

     //recon 1
    GeomCellSelection allcell = toytiling[i]->get_allcell();
    for (int j=0;j!=allcell.size();j++){
      Point p = allcell[j]->center();
      x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.*4;
      y_save = p.y/units::cm;
      z_save = p.z/units::cm;
      

      g->SetPoint(ncount,x_save,y_save,z_save);
      t_rec->Fill();

      ncount ++;
    }
  }

  g->Write("shower3D");
  gt->Write("shower3D_truth");
  g_rec->Write("shower3D_charge");
  g_rec_blob->Write("shower3D_charge_blob");

  
  


  TTree *Trun = new TTree("Trun","Trun");
  Trun->SetDirectory(file);

  int detector = 1; // 35 ton
  Trun->Branch("detector",&detector,"detector/I");

  Trun->Branch("eventNo",&event_no,"eventNo/I");
  Trun->Branch("runNo",&run_no,"runNo/I");
  Trun->Branch("subRunNo",&subrun_no,"runRunNo/I");

  Trun->Branch("unit_dis",&unit_dis,"unit_dis/F");
  Trun->Branch("toffset_uv",&toffset_1,"toffset_uv/F");
  Trun->Branch("toffset_uw",&toffset_2,"toffset_uw/F");
  Trun->Branch("toffset_u",&toffset_3,"toffset_u/F");
  Trun->Branch("total_time_bin",&total_time_bin,"total_time_bin/I");
  Trun->Branch("recon_threshold",&recon_threshold,"recon_threshold/I");
  Trun->Branch("frame_length",&frame_length,"frame_length/I");
  Trun->Branch("max_events",&max_events,"max_events/I");
  Trun->Branch("eve_num",&eve_num,"eve_num/I");
  Trun->Branch("nrebin",&nrebin,"nrebin/I");
  Trun->Branch("threshold_u",&threshold_u,"threshold_u/F");
  Trun->Branch("threshold_v",&threshold_v,"threshold_v/F");
  Trun->Branch("threshold_w",&threshold_w,"threshold_w/F");
  Trun->Branch("time_offset",&time_offset,"time_offset/I");



  Trun->Fill();


  file->Write();
  file->Close();

  return 0;

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
