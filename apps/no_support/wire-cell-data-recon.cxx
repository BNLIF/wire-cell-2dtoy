#include "WCPSst/GeomDataSource.h"
#include "WCPSst/ToyuBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/TruthToyTiling.h"
#include "WCP2dToy/SimpleBlobToyTiling.h"

#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/ToyMatrixExclusive.h"
#include "WCP2dToy/ToyMatrixKalman.h"
#include "WCP2dToy/ToyMatrixIterate.h"
#include "WCP2dToy/ToyMatrixMarkov.h"
#include "WCP2dToy/ToyMetric.h"
#include "WCP2dToy/BlobMetric.h"

#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"

#include "WCPData/GeomCluster.h"
//#include "WCPNav/SliceDataSource.h"


#include "WCPNav/FrameDataSource.h"
#include "WCPNav/SimDataSource.h"
#include "WCPNav/SliceDataSource.h"
#include "WCPSst/Util.h"
#include "WCPData/SimTruth.h"
#include "WCP2dToy/ToyDepositor.h"
#include "WCPNav/GenerativeFDS.h"
#include "WCP2dToy/ToySignalSimu.h"
#include "WCP2dToy/ToySignalSimuTrue.h"
#include "WCP2dToy/ToySignalGaus.h"
#include "WCP2dToy/ToySignalWien.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include <iostream>

using namespace WCP;
using namespace std;




int main(int argc, char* argv[])
{
  if (argc < 4) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num" << endl;
    return 1;
  }
  
  
  WCPSst::GeomDataSource gds(argv[1]);
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
  
  
  const char* root_file = argv[2];
  const char* tpath = "/Event/Sim";
  
  WCP::FrameDataSource* fds = 0;
  fds = WCPSst::make_fds(root_file);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  
  float unit_dis = 1.605723;

  int max_events = 100;
  int eve_num  = atoi(argv[3]);

  TFile tfile(root_file,"read");
  TTree* sst = dynamic_cast<TTree*>(tfile.Get(tpath));
  WCPSst::ToyuBooNEFrameDataSource data_fds(*sst,gds);
  data_fds.jump(eve_num);
  //data_fds.Save();

  int recon_threshold = 2000;
  
  WCP::ToyDepositor toydep(fds,0,unit_dis);
  const PointValueVector pvv = toydep.depositions(eve_num);
  

  
  // 
  
  //  WCP::GenerativeFDS gfds(toydep,gds,2400,max_events,2.0*1.6*units::millimeter);
  // gfds.jump(eve_num);

  // WCPSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons

  // WCP::ToyDepositor toydep(fds);
  // const PointValueVector pvv = toydep.depositions(eve_num);
  //WCP::GenerativeFDS gfds(toydep,gds,9600,1,0.5*1.605723*units::millimeter); // 87 K at 0.5 kV/cm
  WCP::GenerativeFDS gfds(toydep,gds,9600,max_events,0.5*unit_dis*units::millimeter,unit_dis); // 87 K at 0.5 kV/cm
  //WCP2dToy::ToySignalSimuFDS simu_fds(gfds,gds,9600,1,1.647,1.539+1.647,1); // time offset among different planes for the time electrons travel among different planes
  //simu_fds.jump(eve_num);
  //simu_fds.Save();

 

  cout << "Put in Truth " << endl; 
  WCP2dToy::ToySignalSimuTrueFDS st_fds(gfds,gds,9600/4,max_events,0); //truth
  //WCP::GenerativeFDS st_fds(toydep,gds,9600/4,max_events,2.0*1.60*units::millimeter); // 87 K at 0.5 kV/cm
  st_fds.jump(eve_num);
  //st_fds.Save();
  
  
  int time_offset = -52;
  // //test
  // cout << "Simulate Raw WaveForm " << endl; 
  // //WCP2dToy::ToySignalSimuFDS simu_fds(gfds,gds,9600,max_events,1.647,1.539+1.647,1,-55); // time offset among different planes for the time electrons travel among different planes
  // //WCP2dToy::ToySignalSimuFDS simu_fds(gfds,gds,9600,max_events,0,0,1,-0.5,time_offset); // time offset among different planes for the time electrons travel among different planes
  // WCP2dToy::ToySignalSimuFDS simu_fds(gfds,gds,9600,max_events,0,0,1,0,0); // time offset among different planes for the time electrons travel among different planes
  // simu_fds.jump(eve_num);
  // //simu_fds.Save();
  
  
 

  cout << "Deconvolution with Gaussian filter" << endl;
  //WCP2dToy::ToySignalGausFDS gaus_fds(data_fds,gds,9600/4,max_events,1.647,1.539+1.647); // gaussian smearing for charge estimation
  WCP2dToy::ToySignalGausFDS gaus_fds(data_fds,gds,9600/4,max_events,0,0,-0.5); // gaussian smearing for charge estimation
  gaus_fds.jump(eve_num);
  //gaus_fds.Save();

   cout << "Deconvolution with Wiener filter" << endl; 
   //WCP2dToy::ToySignalWienFDS wien_fds(data_fds,gds,9600/4,max_events,1.647,1.539+1.647); // weiner smearing for hit identification
   WCP2dToy::ToySignalWienFDS wien_fds(data_fds,gds,9600/4,max_events,0,0,-0.5); // weiner smearing for hit identification
   wien_fds.jump(eve_num);
  //wien_fds.Save();
  
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  int nwire_u = wires_u.size();
  int nwire_v = wires_v.size();
  int nwire_w = wires_w.size();
  
  float threshold_u = 5.87819e+02 * 4.0;
  float threshold_v = 8.36644e+02 * 4.0;
  float threshold_w = 5.67974e+02 * 4.0;

  float threshold_ug = 755.96;
  float threshold_vg = 822.81;
  float threshold_wg = 510.84;
  
  // float threshold_u = 1000;
  // float threshold_v = 1000;
  // float threshold_w = 1000;
  

  WCPSst::ToyuBooNESliceDataSource sds(wien_fds,gaus_fds,threshold_u, 
					    threshold_v, threshold_w, 
					    threshold_ug, 
					    threshold_vg, threshold_wg, 
					    nwire_u, 
					    nwire_v, nwire_w); 
  
   WCPSst::ToyuBooNESliceDataSource sds_th(st_fds,st_fds,1, 
					    1, 1, 
					    threshold_ug, 
					    threshold_vg, threshold_wg, 
					    nwire_u, 
					    nwire_v, nwire_w); 
  // const int N = 100000;
  // Double_t x[N],y[N],z[N];
  // Double_t x1[N],y1[N],z1[N], charge_r1[N];
  // Double_t xt[N],yt[N],zt[N], charge_t[N];
  int ncount = 0;
  int ncount1 = 0;
  int ncount2 = 0;

  int ncount_t = 0;
  

  WCP2dToy::ToyTiling **toytiling = new WCP2dToy::ToyTiling*[2400];
  
  WCP2dToy::MergeToyTiling **mergetiling = new WCP2dToy::MergeToyTiling*[2400];
  WCP2dToy::TruthToyTiling **truthtiling = new WCP2dToy::TruthToyTiling*[2400];
  
  WCP2dToy::SimpleBlobToyTiling **blobtiling = new WCP2dToy::SimpleBlobToyTiling*[2400];

  WCP2dToy::ToyMatrix **toymatrix = new WCP2dToy::ToyMatrix*[2400];
  // WCP2dToy::ToyMatrixIterate **toymatrix_it = new WCP2dToy::ToyMatrixIterate*[2400];
  //WCP2dToy::ToyMatrixMarkov **toymatrix_markov = new WCP2dToy::ToyMatrixMarkov*[2400];
  

  //save truth ...
  WCP2dToy::ToyTiling **toytiling_th = new WCP2dToy::ToyTiling*[2400];
  WCP2dToy::TruthToyTiling **truthtiling_th = new WCP2dToy::TruthToyTiling*[2400];

  
  WCP2dToy::ToyMetric toymetric;
  WCP2dToy::BlobMetric blobmetric;

  //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  
  int ncount_mcell = 0;

  delete fds;

  int start_num = 0 ;
  int end_num = sds.size()-1;

  // int start_num = 1170 ;
  // int end_num = 1170 ;

  // int start_num =1117;
  // int end_num = 1119;
  // int end_num = sds.size()-1;

  // int start_num = 400;
  // int end_num = 462;

  //  int i=454;{ // 46, 26
  // int i=329;{  // 18, 6,
  // int i=344;{ // 29 14
  //int i=459;{ // 23, 8,   5e5 
  //int i = 351;{
  //for (int i=0;i!=sds.size();i++){
  for (int i=start_num;i!=end_num+1;i++){
 
    sds.jump(i);
    sds_th.jump(i);
    WCP::Slice slice = sds.get();
    WCP::Slice slice_th = sds_th.get();
    //if ( slice.group().size() >0){
      
    toytiling[i] = new WCP2dToy::ToyTiling(slice,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
    mergetiling[i] = new WCP2dToy::MergeToyTiling(*toytiling[i],i,3,1);
    
    GeomCellSelection allcell = toytiling[i]->get_allcell();
    GeomWireSelection allwire = toytiling[i]->get_allwire();
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    GeomWireSelection allmwire = mergetiling[i]->get_allwire();
    
    cout << i << " " << allmcell.size() << " " << allmwire.size() << endl;
    
    truthtiling[i] = new WCP2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds,800,unit_dis);
    toymatrix[i] = new WCP2dToy::ToyMatrix(*toytiling[i],*mergetiling[i]);
    if (toymatrix[i]->Get_Solve_Flag()==0){
      WCP2dToy::ToyMatrixIterate toymatrix_it(*toymatrix[i]);
    }
    
    cout << "chi2: " << toymatrix[i]->Get_Chi2() << endl;
    cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
    

    toytiling_th[i] = new WCP2dToy::ToyTiling(slice_th,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
    truthtiling_th[i] = new WCP2dToy::TruthToyTiling(*toytiling_th[i],pvv,i,gds,800,unit_dis);

    // cout << slice_th.group().size() << " " << toytiling_th[i] ->get_allcell().size() << " " 
    // 	 << toytiling_th[i] ->get_allwire().size() << " " << truthtiling_th[i]->ccmap().size() << endl;

      // GeomCellSelection calmcell;
      // for (int j=0;j!=allmcell.size();j++){
      // 	MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      // 	double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
      // 	double charge_err = toymatrix[i]->Get_Cell_Charge(mcell,2);
	
      // 	//	cout << "Recon: " << j << " " << charge << " " << charge_err << endl;

      // 	if (charge > 2000) calmcell.push_back(mcell);
      // }
      
      //



      CellChargeMap ccmap = truthtiling[i]->ccmap();
      if (toymatrix[i]->Get_Solve_Flag()!=0)
	toymetric.Add(allmcell,*toymatrix[i],ccmap);

      toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());

      Double_t charge_min = 10000;
      Double_t charge_max = 0;

     
      
      
      // //loop through merged cell and compare with truth cells
      // for (int j=0;j!=allmcell.size();j++){
      // 	MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      // 	if (mcell->CheckContainTruthCell(ccmap)) 
      // 	  cout << "True: " << toymatrix[i]->Get_mcindex(mcell) << " " << mcell->GetTruthCharge()<< endl;
      // 	//cout << mergetiling.wires(*allmcell[j]).size() << endl;
      // }
      
      // WireChargeMap wcmap = toytiling[i]->wcmap();




      
    // TApplication theApp("theApp",&argc,argv);
    // theApp.SetReturnFromRun(true);
    
    // TCanvas c1("ToyMC","ToyMC",800,600);
    // c1.Draw();
    
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

    

    // display.init(0,10.3698,-2.33/2.,2.33/2.);
    // //display.init(1.1,1.8,0.7,1.0);
    
    // display.draw_mc(1,WCP::PointValueVector(),"colz");
    
    

    // // display.draw_slice(slice,""); // draw wire 
    // display.draw_cells(toytiling[i]->get_allcell(),"*same");
    // //display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",1); //0 is normal, 1 is only draw the ones containt the truth cell
    // display.draw_mergecells(calmcell,"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell
    // display.draw_truthcells(ccmap,"*same");
    
    // // display.draw_wires_charge(wcmap,"Fsame",FI);
    // // display.draw_cells_charge(toytiling.get_allcell(),"Fsame");
    // // display.draw_truthcells_charge(ccmap,"lFsame",FI);
    
    
    // theApp.Run();
      // }
  }

  toymetric.Print();
  std::cout << "Starting MCMC" << std::endl;

  // //without  time information
  // for (int i=start_num;i!=end_num+1;i++){
  //   if (toymatrix[i]->Get_Solve_Flag()==0){
  //     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
  //     toymatrix_markov[i] = new WCP2dToy::ToyMatrixMarkov(toymatrix[i],&allmcell);
  //     CellChargeMap ccmap = truthtiling[i]->ccmap();
  //     if (toymatrix[i]->Get_Solve_Flag()!=0)
  // 	toymetric.Add(allmcell,*toymatrix[i],ccmap);
  //     toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
  //     cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
  //     cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
  //   }
  // }
    



  //with time information
  if (start_num != end_num){
    int first_solve = -1;
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


  


  // //do blob thing ... 
  // //use time information
  // std::cout << "Reduce Blob" << std::endl; 
  // for (int i=start_num;i!=end_num+1;i++){
  //   std::cout << "Check Blob " << i << std::endl;
  //   //std::cout << toymatrix[i]->GetSimpleBlobReduction() << std::endl;
  //   if (!mergetiling[i]->GetRemerged()){
  //     toymatrix[i]->JudgeSimpleBlob(*toytiling[i],*mergetiling[i]);
  //   }
  //   //std::cout << toymatrix[i]->GetSimpleBlobReduction() << std::endl;
  //   if (toymatrix[i]->GetSimpleBlobReduction()){
  //     if (i==start_num){
  // 	blobtiling[i] = new WCP2dToy::SimpleBlobToyTiling(*toytiling[i],*mergetiling[i],*toymatrix[i],*mergetiling[i+1],*toymatrix[i+1],*mergetiling[i+1],*toymatrix[i+1]);
  //     }else if (i==end_num){
  // 	blobtiling[i] = new WCP2dToy::SimpleBlobToyTiling(*toytiling[i],*mergetiling[i],*toymatrix[i],*mergetiling[i-1],*toymatrix[i-1],*mergetiling[i-1],*toymatrix[i-1]);
  //     }else{
  // 	blobtiling[i] = new WCP2dToy::SimpleBlobToyTiling(*toytiling[i],*mergetiling[i],*toymatrix[i],*mergetiling[i-1],*toymatrix[i-1],*mergetiling[i+1],*toymatrix[i+1]);
  //     }
  //   }
  //   if (toymatrix[i]->GetSimpleBlobReduction()){
  //     //save stuff
  //     CellChargeMap ccmap = truthtiling[i]->ccmap();
  //     blobmetric.Add(*blobtiling[i],ccmap);
      
  //     std::cout << "Check Blob " << i << std::endl;
  //     WCP2dToy::BlobMetric tempblob;
  //     tempblob.Add(*blobtiling[i],ccmap);
  //     tempblob.Print();
  //   }
  // }
  
  //do clustering ... 
   for (int i=start_num;i!=end_num+1;i++){
     GeomCellSelection pallmcell = mergetiling[i]->get_allcell();
     GeomCellSelection allmcell;
     for (int j=0;j!=pallmcell.size();j++){
       const GeomCell* mcell = pallmcell[j];
       if (toymatrix[i]->Get_Cell_Charge(mcell)> recon_threshold){
	 allmcell.push_back(mcell);
       }
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


   TFile *file = new TFile(Form("shower3D_data_%d.root",eve_num),"RECREATE");
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
    CellChargeMap ccmap = truthtiling_th[i]->ccmap();
    for (auto it = ccmap.begin();it!=ccmap.end(); it++){
      Point p = it->first->center();
      x_save = i*unit_dis/10.*2. - 256;
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
      x_save = i*unit_dis/10.*2- 256 -time_offset*unit_dis/10.;
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
      if (charge> recon_threshold){
	//truth
	for (int k=0;k!=mcell->get_allcell().size();k++){
	  Point p = mcell->get_allcell().at(k)->center();
	  x_save = i*unit_dis/10.*2- 256 -time_offset*unit_dis/10.;
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

    
    //recon 3 with charge and deblob
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
      //if (charge> recon_threshold && !(mcell->IsSimpleBlob() && mcell->IsBlob())){
      if (charge> recon_threshold ){
  	for (int k=0;k!=mcell->get_allcell().size();k++){
  	  Point p = mcell->get_allcell().at(k)->center();
  	  x_save = i*unit_dis/10.*2-256 -time_offset*unit_dis/10.;
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
     // if (toymatrix[i]->GetSimpleBlobReduction()){
     //   for (int j=0;j!=blobtiling[i]->Get_Cells().size();j++){
     // 	 const GeomCell *cell = blobtiling[i]->Get_Cells().at(j);
     // 	 Point p = cell->center();
     // 	 x_save = i*0.32-256-time_shift*0.16;
     // 	 y_save = p.y/units::cm;
     // 	 z_save = p.z/units::cm;
     // 	 charge_save = blobtiling[i]->Get_Cell_Charge(cell,1);
     // 	 ncharge_save = 1;
	 
     // 	 g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
     // 	 t_rec_charge_blob->Fill();
	  
     // 	 ncount2 ++;
     //   }
     // }
    
    //save all results
    // file->Write(Form("toytiling_%d",i),toytiling[i]);
    // file->Write(Form("mergetiling_%d",i),mergetiling[i]);
    // file->Write(Form("truthtiling_%d",i),truthtiling[i]);
    // file->Write(Form("toymatrix_%d",i),toymatrix[i]);

  }
 
  

  g->Write("shower3D");
  gt->Write("shower3D_truth");
  g_rec->Write("shower3D_charge");
  g_rec_blob->Write("shower3D_charge_blob");
  
  // const int N = 100000;
  // Double_t x[N],y[N],z[N];
  Double_t x,y,z;
  //save cluster
  int ncluster = 0;
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    ncount = 0;
    TGraph2D *g1 = new TGraph2D();
    for (int i=0; i!=(*it)->get_allcell().size();i++){
      const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
      for (int j=0; j!=mcell->get_allcell().size();j++){
  	Point p = mcell->get_allcell().at(j)->center();
  	x = mcell->GetTimeSlice()*0.32- 256;
  	y = p.y/units::cm;
  	z = p.z/units::cm;
  	g1->SetPoint(ncount,x,y,z);
  	ncount ++;
      }
    }
    // cout << ncount << endl;
    g1->Write(Form("cluster_%d",ncluster));
    ncluster ++;
  }

  // WCP2dToy::ToyTiling* tt1;
  // WCP2dToy::MergeToyTiling* tt2;
  // WCP2dToy::TruthToyTiling* tt3;

  // WCP2dToy::ToyMatrix* tt4;
  // WCP2dToy::ToyMatrixIterate* tt5;
  // WCP2dToy::ToyMatrixMarkov* tt6;
  // WCP2dToy::SimpleBlobToyTiling *tt7;
  // int time_slice;
  
  // TTree* ttree = new TTree("T","T");
  // ttree->Branch("time_slice",&time_slice,"time_slice/I");
  // ttree->Branch("toytiling",&tt1);
  // //ttree->Branch("truthtiling",&tt3);
  // //ttree->Branch("mergetiling",&tt2);

  
  // //ttree->Branch("toymatrix",&tt4);
  // // ttree->Branch("toymatrix_ite",&tt5);
  // // ttree->Branch("toymatrix_markov",&tt6);
  // //ttree->Branch("blobtiling",&tt7);
  

  // ttree->SetDirectory(file);
  // for (int i=start_num;i!=end_num+1;i++){
  //   tt1 = toytiling[i];
  //   //tt3 = truthtiling[i];

  //   //tt2 = mergetiling[i];
  

  //   //tt4 = toymatrix[i];
  //   // tt5 = toymatrix_it[i];
  //   // tt6 = toymatrix_markov[i];
  //   //tt7 = blobtiling[i];
    
  //   time_slice = i;
  //   ttree->Fill();
  // }

  // // TTree *ttree1 = new TTree("TC","TC");
  // // GeomCluster *cluster;
  // // ttree1->Branch("cluster",&cluster);
  // // ttree1->SetDirectory(file);
  // // for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  // //   cluster = *it;
  // //   ttree1->Fill();
  // // }
  
 


  file->Write();
  file->Close();

  toymetric.Print();
  blobmetric.Print();

  return 0;
  
} // main()
