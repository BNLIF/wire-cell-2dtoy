#include "WCPNav/DetectorGDS.h"
#include "WCPNav/DetGenerativeFDS.h"
#include "WCPNav/FrameDataSource.h"
#include "WCP2dToy/ToyDepositor.h"
#include "WCPSst/Util.h"
#include "WCP2dToy/ToySignalSimu.h"
#include "WCP2dToy/ToySignalSimuTrue.h"
#include "WCP2dToy/ToySignalGaus.h"
#include "WCP2dToy/ToySignalWien.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/TruthToyTiling.h"
#include "WCPData/MergeGeomCell.h"

#include "WCP2dToy/ToyEventDisplay.h"

#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/ToyMatrixExclusive.h"
#include "WCP2dToy/ToyMatrixKalman.h"
#include "WCP2dToy/ToyMatrixIterate.h"
#include "WCP2dToy/ToyMatrixMarkov.h"
#include "WCP2dToy/ToyMetric.h"
#include "WCPData/TPCParams.h"
#include "WCPData/Singleton.h"

#include "WCPData/GeomCluster.h"

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


using namespace WCP;
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
  int frame_length = 3200;
  //int frame_length = 800;  // hack for now
  int nrebin = 4;

  TPCParams& mp = Singleton<TPCParams>::Instance();
  
  double pitch_u = gds.get_pitch(0,WirePlaneType_t(0));
  double pitch_v = gds.get_pitch(0,WirePlaneType_t(1));
  double pitch_w = gds.get_pitch(0,WirePlaneType_t(2));
  double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;

  mp.set_pitch_u(pitch_u);
  mp.set_pitch_v(pitch_v);
  mp.set_pitch_w(pitch_w);
  mp.set_ts_width(time_slice_width);
  
  std::cout << "Singleton: " << mp.get_pitch_u() << " " << mp.get_pitch_v() << " " << mp.get_pitch_w() << " " << mp.get_ts_width() << std::endl;
  
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
  
  WCP::FrameDataSource* fds = 0;
  fds = WCPSst::make_fds(*tfile);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  
  TH1::AddDirectory(kFALSE);
  
  WCP::ToyDepositor toydep(fds,0,unit_dis,frame_length);
  const PointValueVector& pvv = toydep.depositions(eve_num);

  std::cout << "Points deposited: " << pvv.size() << std::endl;

  // DetGenerativeFDS gfds(toydep,gds, 2400,max_events,2.0*1.6*units::millimeter);
  DetGenerativeFDS gfds(toydep,gds, total_time_bin,max_events,0.5*1.6*units::millimeter);
  gfds.jump(eve_num);

  tfile->Close("R");
  delete tfile;

  //WCPSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons

  
  
  cout << "Put in Truth " << endl; 
  WCP2dToy::ToySignalSimuTrueFDS st_fds(gfds,gds,total_time_bin/nrebin,max_events,0); //truth
  st_fds.jump(eve_num);
  
  cout << "Simulate Raw WaveForm " << endl; 
  WCP2dToy::ToySignalSimuFDS simu_fds(gfds,gds,total_time_bin,max_events,toffset_1,toffset_2,1); // time offset among different planes for the time electrons travel among different planes
  simu_fds.jump(eve_num);
  //simu_fds.Save();
  
  cout << "Deconvolution with Gaussian filter" << endl;
  WCP2dToy::ToySignalGausFDS gaus_fds(simu_fds,gds,total_time_bin/nrebin,max_events,toffset_1,toffset_2); // gaussian smearing for charge estimation
  gaus_fds.jump(eve_num);
  //gaus_fds.Save();

  cout << "Deconvolution with Wiener filter" << endl;
  WCP2dToy::ToySignalWienFDS wien_fds(simu_fds,gds,total_time_bin/nrebin,max_events,toffset_1,toffset_2); // weiner smearing for hit identification
  wien_fds.jump(eve_num);
  //wien_fds.Save();
  
  
  

  int nwire_u = gds.get_total_nwires(WirePlaneType_t(0));
  int nwire_v = gds.get_total_nwires(WirePlaneType_t(1));
  int nwire_w = gds.get_total_nwires(WirePlaneType_t(2));
  
 
  
  // cin >> abc;

  WCPSst::ToyuBooNESliceDataSource sds(gds,wien_fds,gaus_fds,threshold_u, 
  					    threshold_v, threshold_w, 
  					    threshold_ug, 
  					    threshold_vg, threshold_wg, 
  					    nwire_u, 
  					    nwire_v, nwire_w); 

  WCPSst::ToyuBooNESliceDataSource sds_th(gds,st_fds,st_fds,500, 
  					    500, 500, 
  					    threshold_ug, 
  					    threshold_vg, threshold_wg, 
  					    nwire_u, 
  					    nwire_v, nwire_w); 
  

   int ncount = 0;
  int ncount1 = 0;  
  int ncount2 = 0;

  int ncount_t = 0;


   WCP2dToy::ToyTiling **toytiling = new WCP2dToy::ToyTiling*[2400];
   WCP2dToy::MergeToyTiling **mergetiling = new WCP2dToy::MergeToyTiling*[2400];
   WCP2dToy::TruthToyTiling **truthtiling = new WCP2dToy::TruthToyTiling*[2400];
   WCP2dToy::ToyMatrix **toymatrix = new WCP2dToy::ToyMatrix*[2400];
   
    //save truth ...
   WCP2dToy::ToyTiling **toytiling_th = new WCP2dToy::ToyTiling*[2400];
   WCP2dToy::TruthToyTiling **truthtiling_th = new WCP2dToy::TruthToyTiling*[2400];

   WCP2dToy::ToyMetric toymetric;
   

   //add in cluster
   GeomClusterSet cluster_set, cluster_delset;
   
   int ncount_mcell = 0;
   
   delete fds;

   int start_num = 0 ;
   int end_num = sds.size()-1;
   
   
  //int i = 317+800;{
   //int i = 292+800;
   for (int i=start_num;i!=end_num+1;i++){
     sds.jump(i);
     sds_th.jump(i);
     
     WCP::Slice slice = sds.get();
     WCP::Slice slice_th = sds_th.get();
     cout << i << " " << slice.group().size() << " " << slice_th.group().size() << endl;

     // if ( slice.group().size() >0){
     // cout << i << " " << slice.group().size() << endl;
     toytiling[i] = new WCP2dToy::ToyTiling(slice,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
     //allcell = toytiling[i]->get_allcell();
     GeomCellSelection allcell = toytiling[i]->get_allcell();
     GeomWireSelection allwire = toytiling[i]->get_allwire();
     cout << "Single Cell: " << i << " "  << allcell.size() << " " << allwire.size() << endl;

     mergetiling[i] = new WCP2dToy::MergeToyTiling(gds,*toytiling[i],i); 
     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
     GeomWireSelection allmwire = mergetiling[i]->get_allwire();
     
     cout <<"Blob: " << i << " " << allmcell.size() << " " << allmwire.size() << endl;
     
     
     truthtiling[i] = new WCP2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds,frame_length/nrebin,unit_dis);
     
     
     toymatrix[i] = new WCP2dToy::ToyMatrix(gds,*toytiling[i],*mergetiling[i]);
     if (toymatrix[i]->Get_Solve_Flag()==0){
       WCP2dToy::ToyMatrixIterate toymatrix_it(*toymatrix[i],2000,1e5);
     }
     
     cout << i << " chi2: " << toymatrix[i]->Get_Chi2() <<
       " NDF: " << toymatrix[i]->Get_ndf() << endl;
     
     toytiling_th[i] = new WCP2dToy::ToyTiling(slice_th,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
    truthtiling_th[i] = new WCP2dToy::TruthToyTiling(*toytiling_th[i],pvv,i,gds,frame_length/nrebin,unit_dis);
    
     CellChargeMap ccmap = truthtiling[i]->ccmap();
     if (toymatrix[i]->Get_Solve_Flag()!=0)
       toymetric.Add(allmcell,*toymatrix[i],ccmap);
     
     toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());

     
       //}
     
    //  TApplication theApp("theApp",&argc,argv);
    // theApp.SetReturnFromRun(true);
    
    // TCanvas c1("ToyMC","ToyMC",1200,600);
    // c1.Divide(2,1);
    // c1.Draw();
    
    // float charge_min = 0;
    // float charge_max = 1e5;


    // WCP2dToy::ToyEventDisplay display(c1, gds);
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
    // display.draw_mc(1,WCP::PointValueVector(),"colz");
    // display.draw_slice(slice,"");
    // display.draw_cells(toytiling[i]->get_allcell(),"*same");

    
    // //std::cout << ccmap.size() << std::endl;
    // display.draw_truthcells(ccmap,"*same");
    // display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell

    
    // theApp.Run();
   }
   toymetric.Print();
   


   // figure out how to add chi^2 in the solving ... 
   




   std::cout << "Starting MCMC" << std::endl;
   
   //without  time information
   // for (int i=start_num;i!=end_num+1;i++){
   //   if (toymatrix[i]->Get_Solve_Flag()==0){
   //     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   //     WCP2dToy::ToyMatrixMarkov toymatrix_markov(*toymatrix[i],*mergetiling[i],&allmcell,Good_MCells.at(i-start_num));
   //     //WCP2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i],&allmcell);
   //     CellChargeMap ccmap = truthtiling[i]->ccmap();
   //     if (toymatrix[i]->Get_Solve_Flag()!=0)
   // 	toymetric.Add(allmcell,*toymatrix[i],ccmap);
   //     toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
   //     cout << " chi2: " << i << " " << toymatrix[i]->Get_Chi2() 
   // 	   << " NDF: " << toymatrix[i]->Get_ndf() << endl;
   //   }
   // }
   
   
   //with time information
   if (start_num != end_num){
     int first_solve=-1;
     for (int i=start_num; i!=end_num+1;i++){
       if (toymatrix[i]->Get_Solve_Flag()!=0){
	 first_solve = i;
	 break;
       }
     }
     if (first_solve <0){
       for (int i=start_num;i!=end_num+1;i++){
	 if (toymatrix[i]->Get_Solve_Flag()==0){
	   GeomCellSelection allmcell = mergetiling[i]->get_allcell();
	   WCP2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i],&allmcell);
	   CellChargeMap ccmap = truthtiling[i]->ccmap();
	   if (toymatrix[i]->Get_Solve_Flag()!=0)
	     toymetric.Add(allmcell,*toymatrix[i],ccmap);
	   toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	   cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
	   cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
	 }
       }
     }else{
       for (int i=first_solve+1;i<=end_num-1;i++){
	 if (toymatrix[i]->Get_Solve_Flag()==0){
	   GeomCellSelection allmcell = mergetiling[i]->get_allcell();
	   WCP2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell);
	   CellChargeMap ccmap = truthtiling[i]->ccmap();
	   if (toymatrix[i]->Get_Solve_Flag()!=0)
	     toymetric.Add(allmcell,*toymatrix[i],ccmap);
	   toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	   
	   cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
	   cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
	   
	   //toymetric.Print();
	 }
       }
       
       if (toymatrix[end_num]->Get_Solve_Flag()==0){
	 GeomCellSelection allmcell = mergetiling[end_num]->get_allcell();
	 WCP2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[end_num-1],toymatrix[end_num],0,mergetiling[end_num-1],mergetiling[end_num],0,&allmcell);
	 
	 
	 CellChargeMap ccmap = truthtiling[end_num]->ccmap();
	 if (toymatrix[end_num]->Get_Solve_Flag()!=0)
	   toymetric.Add(allmcell,*toymatrix[end_num],ccmap);
	 toymetric.AddSolve(toymatrix[end_num]->Get_Solve_Flag());
	 
  	cout << "chi2: " << end_num << " " << toymatrix[end_num]->Get_Chi2() << endl;
  	cout << "NDF: " << toymatrix[end_num]->Get_ndf() << endl;
      }
      
      // go to early ones 
      for (int i=first_solve-1;i>=start_num+1;i--){
  	if (toymatrix[i]->Get_Solve_Flag()==0){
  	  GeomCellSelection allmcell = mergetiling[i]->get_allcell();
  	  WCP2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell);
	  
  	  CellChargeMap ccmap = truthtiling[i]->ccmap();
  	  if (toymatrix[i]->Get_Solve_Flag()!=0)
  	    toymetric.Add(allmcell,*toymatrix[i],ccmap);
  	  toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	  
  	  cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
  	  cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
  	}
      }
      
      if (toymatrix[start_num]->Get_Solve_Flag()==0){
  	GeomCellSelection allmcell = mergetiling[start_num]->get_allcell();
  	WCP2dToy::ToyMatrixMarkov toymatrix_markov(0,toymatrix[start_num],toymatrix[start_num+1],0,mergetiling[start_num],mergetiling[start_num+1],&allmcell);
	
	
  	CellChargeMap ccmap = truthtiling[start_num]->ccmap();
  	if (toymatrix[start_num]->Get_Solve_Flag()!=0)
  	  toymetric.Add(allmcell,*toymatrix[start_num],ccmap);
  	toymetric.AddSolve(toymatrix[start_num]->Get_Solve_Flag());
	
  	cout << "chi2: " << start_num << " " << toymatrix[start_num]->Get_Chi2() << endl;
  	cout << "NDF: " << toymatrix[start_num]->Get_ndf() << endl;
      }
    }
  }

     // form a map to illustrate connectivities 
   GeomCellCellsMap cell_prev_map;
   GeomCellCellsMap cell_next_map;
   
   for (int i=start_num;i!=end_num;i++){
     GeomCellSelection curr_mcell = mergetiling[i]->get_allcell();
     GeomCellSelection next_mcell = mergetiling[i+1]->get_allcell();
     for (int j=0;j!=curr_mcell.size();j++){
       const MergeGeomCell *curr_cell = (MergeGeomCell*)curr_mcell.at(j);
       for (int k=0;k!=next_mcell.size();k++){
	 const MergeGeomCell *next_cell = (MergeGeomCell*)next_mcell.at(k);
	 if (curr_cell->Overlap(*next_cell) && curr_cell->get_allcell().at(0)->get_face() == next_cell->get_allcell().at(0)->get_face()){
	   if (cell_next_map.find(curr_cell)==cell_next_map.end()){
	     GeomCellSelection cells;
	     cell_next_map[curr_cell] = cells;
	   }
	   cell_next_map[curr_cell].push_back(next_cell);
	   
	   if (cell_prev_map.find(next_cell)==cell_prev_map.end()){
	     GeomCellSelection cells;
	     cell_prev_map[next_cell] = cells;
	   }
	   cell_prev_map[next_cell].push_back(curr_cell);
	   
	 }
       }
     }
   }
   // 

   // save good cluster_cells;
   GeomCellSelection good_cluster_cells;
   for (int i=start_num;i!=end_num+1;i++){
     GeomCellSelection curr_mcell = mergetiling[i]->get_allcell();
     for (int j=0;j!=curr_mcell.size();j++){
       const MergeGeomCell *curr_cell = (MergeGeomCell*)curr_mcell.at(j);
       int flag_good_cluster_cell = 0;
     
       if (cell_next_map.find(curr_cell) != cell_next_map.end()){
	 for (int k=0;k!=cell_next_map[curr_cell].size();k++){
	   if (toymatrix[i+1]->Get_Cell_Charge(cell_next_map[curr_cell].at(k))> recon_threshold || toymatrix[i+1]->Get_Solve_Flag()==0 ){
	     flag_good_cluster_cell = 1;
	     break;
	   }
	 }
       }
       if (flag_good_cluster_cell == 0){
	 if (cell_prev_map.find(curr_cell) != cell_prev_map.end()){
	   for (int k=0;k!=cell_prev_map[curr_cell].size();k++){
	     if (toymatrix[i-1]->Get_Cell_Charge(cell_prev_map[curr_cell].at(k))> recon_threshold || toymatrix[i-1]->Get_Solve_Flag()==0 ){
	       flag_good_cluster_cell = 1;
	       break;
	     }
	   }
	 }
       }

       if (flag_good_cluster_cell == 1)
	 good_cluster_cells.push_back(curr_cell);
     }
   }
   


   //do clustering ... 
   for (int i=start_num;i!=end_num+1;i++){
     GeomCellSelection pallmcell = mergetiling[i]->get_allcell();
     GeomCellSelection allmcell;
     for (int j=0;j!=pallmcell.size();j++){
       const GeomCell* mcell = pallmcell[j];
       
       // if (toymatrix[i]->Get_Cell_Charge(mcell)> recon_threshold || toymatrix[i]->Get_Solve_Flag()==0){
       // 	 allmcell.push_back(mcell);
       // }
       
        int flag_save_cell = 0;

      if (toymatrix[i]->Get_Solve_Flag()==0){
      	flag_save_cell = 1;
      }else{
      	if (toymatrix[i]->Get_Cell_Charge(mcell)> recon_threshold){
      	  flag_save_cell = 1;
      	}else{
      	  if (i == start_num || i == end_num + 1) continue;
	  
      	  // there are good cells from the prev and next
      	  flag_save_cell = 0;
      	  //	  std::cout << "Xin: " << cell_next_map[mcell].size() <<  " " << cell_prev_map[mcell].size() << std::endl;

      	  for (int k=0;k!=cell_next_map[mcell].size();k++){
	    if (find(good_cluster_cells.begin(),good_cluster_cells.end(),cell_next_map[mcell].at(k))!=good_cluster_cells.end()){
	      //if (toymatrix[i+1]->Get_Cell_Charge(cell_next_map[mcell].at(k))> recon_threshold){
      	      flag_save_cell = 1;
      	      break;
      	    }
      	  }
      	  if (flag_save_cell==1){
      	    flag_save_cell = 0;
      	    for (int k=0;k!=cell_prev_map[mcell].size();k++){
	      if (find(good_cluster_cells.begin(),good_cluster_cells.end(),cell_prev_map[mcell].at(k))!=good_cluster_cells.end()){
		//if (toymatrix[i-1]->Get_Cell_Charge(cell_prev_map[mcell].at(k))> recon_threshold){
      		flag_save_cell = 1;
      		break;
      	      }
      	    }
      	  }
      	}
      }
      
      
      if (flag_save_cell == 1)
      	allmcell.push_back(mcell);
     }
     
     
     
     if (cluster_set.empty()){
       // if cluster is empty, just insert all the mcell, each as a cluster
       for (int j=0;j!=allmcell.size();j++){
  	 GeomCluster *cluster = new GeomCluster(*((MergeGeomCell*)allmcell[j]));
  	 cluster_set.insert(cluster);
       }
     }else{
       for (int j=0;j!=allmcell.size();j++){
  	 int flag = 0;
  	 int flag_save = 0;
  	 GeomCluster *cluster_save = 0;
	 
  	 cluster_delset.clear();
	 
  	 // loop through merged cell
  	 for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  	   //loop through clusters
	   
  	   flag += (*it)->AddCell(*((MergeGeomCell*)allmcell[j]));
  	   if (flag==1 && flag != flag_save){
  	     cluster_save = *it;
  	   }else if (flag>1 && flag != flag_save){
  	     cluster_save->MergeCluster(*(*it));
  	     cluster_delset.insert(*it);
  	   }
  	   flag_save = flag;
  	   
  	 }
	 
  	 for (auto it = cluster_delset.begin();it!=cluster_delset.end();it++){
  	   cluster_set.erase(*it);
  	   delete (*it);
  	 }
	 
  	 if (flag==0){
  	   GeomCluster *cluster = new GeomCluster(*((MergeGeomCell*)allmcell[j]));
  	   cluster_set.insert(cluster);
  	 }
	 
       }
     }

     int ncount_mcell_cluster = 0;
      for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  	ncount_mcell_cluster += (*it)->get_allcell().size();
      }
      ncount_mcell += allmcell.size();
      cout << i << " " << allmcell.size()  << " " << cluster_set.size()  << endl;
   }

   int ncount_mcell_cluster = 0;
   for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
     ncount_mcell_cluster += (*it)->get_allcell().size();
   }
   cout << "Summary: " << ncount << " " << ncount_mcell << " " << ncount_mcell_cluster << endl;

 
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
      
      int cryo = it->first->get_cryo();
      int apa = it->first->get_apa();
      int face = it->first->get_face();
      const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa);
      std::pair<double, double> xmm = apa_gds->minmax(0); 
      
      if (face == 1){
	x_save = (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
      }else if (face == 0){
	x_save = xmm.first/units::cm - (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
      }
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

      int cryo = allcell[j]->get_cryo();
      int apa = allcell[j]->get_apa();
      int face = allcell[j]->get_face();
      const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa);
      std::pair<double, double> xmm = apa_gds->minmax(0); 
      
      if (face == 1){
	x_save = (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
      }else if (face == 0){
	x_save = xmm.first/units::cm - (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
      }

      //      x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.*4;
      y_save = p.y/units::cm;
      z_save = p.z/units::cm;
      

      g->SetPoint(ncount,x_save,y_save,z_save);
      t_rec->Fill();

      ncount ++;
    }

      //recon 2 with charge
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
      if (charge> recon_threshold || toymatrix[i]->Get_Solve_Flag()==0){

    	if (toymatrix[i]->Get_Solve_Flag()==0)
	  charge = toytiling[i]->get_ave_charge();

	//truth
    	for (int k=0;k!=mcell->get_allcell().size();k++){
    	  Point p = mcell->get_allcell().at(k)->center();


	  int cryo = mcell->get_allcell().at(k)->get_cryo();
	  int apa = mcell->get_allcell().at(k)->get_apa();
	  int face = mcell->get_allcell().at(k)->get_face();
	  const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa);
	  std::pair<double, double> xmm = apa_gds->minmax(0); 
	  
	  if (face == 1){
	    x_save = (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
	  }else if (face == 0){
	    x_save = xmm.first/units::cm - (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
	  }
	  
    	  y_save = p.y/units::cm;
    	  z_save = p.z/units::cm;
    	  charge_save = charge/mcell->get_allcell().size();
    	  ncharge_save = mcell->get_allcell().size();
    	  chi2_save = toymatrix[i]->Get_Chi2();
    	  ndf_save = toymatrix[i]->Get_ndf();

    	  g_rec->SetPoint(ncount1,x_save,y_save,z_save);
    	  t_rec_charge->Fill();
	  
    	  ncount1 ++;
    	}
      }
    }
    
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
      if (charge> recon_threshold || toymatrix[i]->Get_Solve_Flag()==0){
	if (toymatrix[i]->Get_Solve_Flag()==0)
	  charge = toytiling[i]->get_ave_charge();

    	for (int k=0;k!=mcell->get_allcell().size();k++){
    	  Point p = mcell->get_allcell().at(k)->center();

	  int cryo = mcell->get_allcell().at(k)->get_cryo();
	  int apa = mcell->get_allcell().at(k)->get_apa();
	  int face = mcell->get_allcell().at(k)->get_face();
	  const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa);
	  std::pair<double, double> xmm = apa_gds->minmax(0); 
	  
	  if (face == 1){
	    x_save = (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
	  }else if (face == 0){
	    x_save = xmm.first/units::cm - (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
	  }

    	  
    	  y_save = p.y/units::cm;
    	  z_save = p.z/units::cm;
    	  charge_save = charge/mcell->get_allcell().size();
    	  ncharge_save = mcell->get_allcell().size();
	  
    	  g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
    	  t_rec_charge_blob->Fill();
	  
    	  ncount2 ++;
    	}
      }
    }

  }

  g->Write("shower3D");
  gt->Write("shower3D_truth");
  g_rec->Write("shower3D_charge");
  g_rec_blob->Write("shower3D_charge_blob");

  
   TTree *ttree1 = new TTree("TC","TC");
  // To save cluster, we need to save
  // 1. time slice
  // 2. single cell
  // 3. charge
  // 4. cluster number
  const GeomCell* cell_save = 0;
  Int_t cluster_num = -1;
  Int_t mcell_id = -1;
  Int_t time_slice;
  
  ttree1->Branch("time_slice",&time_slice,"time_slice/I"); // done
  ttree1->Branch("cell",&cell_save);
  ttree1->Branch("ncluster",&cluster_num,"cluster_num/I"); //done
  ttree1->Branch("mcell_id",&mcell_id,"mcell_id/I");
  ttree1->Branch("charge",&charge_save,"charge/D"); 
  
  Double_t xx,yy,zz;

  ttree1->Branch("xx",&xx,"xx/D");    //don
  ttree1->Branch("yy",&yy,"yy/D");    //don
  ttree1->Branch("zz",&zz,"zz/D");    //don
  // ttree1->Branch("x",&x,"x/D");    //done
  // ttree1->Branch("y",&y,"y/D");
  // ttree1->Branch("z",&z,"z/D");

  // save information to reconstruct the toytiling
  Int_t u_index, v_index, w_index;
  Double_t u_charge, v_charge, w_charge;
  Double_t u_charge_err, v_charge_err, w_charge_err;
  
  Int_t apa_no=0, cryostat_no=0;
  Int_t face;
  ttree1->Branch("apa_no",&apa_no,"apa_no/I");
  ttree1->Branch("cryostat_no",&cryostat_no,"cryostat_no/I");
  ttree1->Branch("face",&face,"face/I");
  
  ttree1->Branch("u_index",&u_index,"u_index/I");
  ttree1->Branch("v_index",&v_index,"v_index/I");
  ttree1->Branch("w_index",&w_index,"w_index/I");
  
  ttree1->Branch("u_charge",&u_charge,"u_charge/D");
  ttree1->Branch("v_charge",&v_charge,"v_charge/D");
  ttree1->Branch("w_charge",&w_charge,"w_charge/D");

  ttree1->Branch("u_charge_err",&u_charge_err,"u_charge_err/D");
  ttree1->Branch("v_charge_err",&v_charge_err,"v_charge_err/D");
  ttree1->Branch("w_charge_err",&w_charge_err,"w_charge_err/D");

  //end save 

  ttree1->SetDirectory(file);
  
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    cluster_num ++;
    //loop merged cell
    for (int i=0; i!=(*it)->get_allcell().size();i++){
      const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
      mcell_id ++;
      time_slice = mcell->GetTimeSlice();

      //      
      //loop single cell
      for (int j=0; j!=mcell->get_allcell().size();j++){
	cell_save = mcell->get_allcell().at(j);

	//fill the information needed for toytiling
	GeomWireSelection wires = toytiling[time_slice]->wires(*cell_save);
	
	//	if (i==0 && j==0) cout << "abc: " << time_slice << " " << toytiling[time_slice]->get_allcell().size() << " " << wires.size() << endl;

	for (int k=0;k!=wires.size();k++){
	  const GeomWire *wire = wires.at(k);
	  WirePlaneType_t plane = wire->plane();
	  if (plane==0){
	    u_index = wire->index();
	    u_charge = toytiling[time_slice]->wcmap()[wire];
	    u_charge_err = toytiling[time_slice]->wcemap()[wire];
	  }else if (plane==1){
	    v_index = wires.at(k)->index();
	    v_charge = toytiling[time_slice]->wcmap()[wire];
	    v_charge_err = toytiling[time_slice]->wcemap()[wire];
	  }else if (plane==2){
	    w_index = wire->index();
	    w_charge = toytiling[time_slice]->wcmap()[wire];
	    w_charge_err = toytiling[time_slice]->wcemap()[wire];
	  }
	}

	//end fill

	Point p = mcell->get_allcell().at(j)->center();

	int cryo = mcell->get_allcell().at(j)->get_cryo();
	apa_no = mcell->get_allcell().at(j)->get_apa();
	face = mcell->get_allcell().at(j)->get_face();
	const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa_no);
	std::pair<double, double> xmm = apa_gds->minmax(0); 
	
	//std::cout << "xin: " << xmm.first/units::cm << " " << xmm.second/units::cm << " " << nrebin << " " << unit_dis << " " << frame_length << std::endl;

	if (face == 1){
	  xx = (time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
	}else if (face == 0){
	  xx = xmm.first/units::cm - (time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
	}
	//xx = time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
	
	charge_save = toymatrix[time_slice]->Get_Cell_Charge(mcell,1)/mcell->cross_section() * cell_save->cross_section();
	yy = p.y/units::cm;
  	zz = p.z/units::cm;
	ttree1->Fill();
	
      }
    }
  }
  ttree1->Write();


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
