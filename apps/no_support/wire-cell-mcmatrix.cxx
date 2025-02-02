#include "WCPSst/GeomDataSource.h"
//#include "WCPSst/ToyuBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/TruthToyTiling.h"
#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/ToyMatrixExclusive.h"
#include "WCP2dToy/ToyMatrixKalman.h"
#include "WCP2dToy/ToyMatrixIterate.h"
#include "WCP2dToy/ToyMatrixMarkov.h"
#include "WCP2dToy/ToyMetric.h"

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
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root" << endl;
    return 1;
  }
  
  
  WCPSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;
  
  
  const char* root_file = argv[2];
  const char* tpath = "/Event/Sim";
  
  WCP::FrameDataSource* fds = 0;
  fds = WCPSst::make_fds(root_file);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  

  
 
  int recon_threshold = 2000;
 


  WCP::ToyDepositor toydep(fds);
  const PointValueVector pvv = toydep.depositions(1);
  
  


  WCP::GenerativeFDS gfds(toydep,gds,2400,5,2.0*1.6*units::millimeter);
  gfds.jump(1);

  WCPSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons

  

  // const int N = 100000;
  // Double_t x[N],y[N],z[N];
  // Double_t x1[N],y1[N],z1[N], charge_r1[N];
  // Double_t xt[N],yt[N],zt[N], charge_t[N];
  int ncount = 0;
  int ncount1 = 0;
  int ncount_t = 0;
  

  WCP2dToy::ToyTiling **toytiling = new WCP2dToy::ToyTiling*[2400];
  WCP2dToy::MergeToyTiling **mergetiling = new WCP2dToy::MergeToyTiling*[2400];
  WCP2dToy::TruthToyTiling **truthtiling = new WCP2dToy::TruthToyTiling*[2400];
  WCP2dToy::ToyMatrix **toymatrix = new WCP2dToy::ToyMatrix*[2400];
  WCP2dToy::ToyMatrixIterate **toymatrix_it = new WCP2dToy::ToyMatrixIterate*[2400];
  WCP2dToy::ToyMatrixMarkov **toymatrix_markov = new WCP2dToy::ToyMatrixMarkov*[2400];
  
  WCP2dToy::ToyMetric toymetric;


  //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  
  int ncount_mcell = 0;

  int start_num = 0 ;
  int end_num = sds.size()-1;

  // int start_num =185;
  // int end_num = 185;
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
    WCP::Slice slice = sds.get();
    //if ( slice.group().size() >0){
      
      toytiling[i] = new WCP2dToy::ToyTiling(slice,gds);
      mergetiling[i] = new WCP2dToy::MergeToyTiling(*toytiling[i],i);

      GeomCellSelection allcell = toytiling[i]->get_allcell();
      GeomWireSelection allwire = toytiling[i]->get_allwire();
      GeomCellSelection allmcell = mergetiling[i]->get_allcell();
      GeomWireSelection allmwire = mergetiling[i]->get_allwire();
      
      cout << i << " " << allmcell.size() << " " << allmwire.size() << endl;

      truthtiling[i] = new WCP2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds);
      toymatrix[i] = new WCP2dToy::ToyMatrix(*toytiling[i],*mergetiling[i]);
      if (toymatrix[i]->Get_Solve_Flag()==0)
      	toymatrix_it[i] = new WCP2dToy::ToyMatrixIterate(*toymatrix[i]);
      
      cout << "chi2: " << toymatrix[i]->Get_Chi2() << endl;
      cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
      
      


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
  std::cout << "Starting Time" << std::endl;

  if (start_num != end_num){
    int first_solve;
    for (int i=start_num; i!=end_num+1;i++){
      if (toymatrix[i]->Get_Solve_Flag()!=0){
	first_solve = i;
  	break;
      }
    }
  

    for (int i=first_solve+1;i<=end_num-1;i++){
      if (toymatrix[i]->Get_Solve_Flag()==0){
  	GeomCellSelection allmcell = mergetiling[i]->get_allcell();
	toymatrix_markov[i] = new WCP2dToy::ToyMatrixMarkov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell);
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
      toymatrix_markov[end_num] = new WCP2dToy::ToyMatrixMarkov(toymatrix[end_num-1],toymatrix[end_num],0,mergetiling[end_num-1],mergetiling[end_num],0,&allmcell);

      
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
	toymatrix_markov[i] = new WCP2dToy::ToyMatrixMarkov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell);
  	
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
      toymatrix_markov[start_num] = new WCP2dToy::ToyMatrixMarkov(0,toymatrix[start_num],toymatrix[start_num+1],0,mergetiling[start_num],mergetiling[start_num+1],&allmcell);

      
      CellChargeMap ccmap = truthtiling[start_num]->ccmap();
      if (toymatrix[start_num]->Get_Solve_Flag()!=0)
  	toymetric.Add(allmcell,*toymatrix[start_num],ccmap);
      toymetric.AddSolve(toymatrix[start_num]->Get_Solve_Flag());

      cout << "chi2: " << start_num << " " << toymatrix[start_num]->Get_Chi2() << endl;
      cout << "NDF: " << toymatrix[start_num]->Get_ndf() << endl;
    }
    
  }


  //do clustering ... 
   for (int i=start_num;i!=end_num+1;i++){
     //GeomCellSelection allcell = toytiling[i]->get_allcell();
     GeomCellSelection pallmcell = mergetiling[i]->get_allcell();
     GeomCellSelection allmcell;
     for (int j=0;j!=pallmcell.size();j++){
       const GeomCell* mcell = pallmcell[j];
       if (toymatrix[i]->Get_Cell_Charge(mcell)>recon_threshold){
	 allmcell.push_back(mcell);
       }
     }
     //GeomWireSelection allmwire = mergetiling[i]->get_allwire();
          

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


  TFile *file = new TFile("shower3D.root","RECREATE");
  TTree *t_true = new TTree("T_true","T_true");
  TTree *t_rec = new TTree("T_rec","T_rec");
  TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
  
  Double_t x_save, y_save, z_save;
  Double_t charge_save;
  Double_t ncharge_save;
  Double_t chi2_save;
  Double_t ndf_save;
  Double_t ncell_save;

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
  t_rec_charge->Branch("ncell",&ncell_save,"ncell/D");
  t_rec_charge->Branch("chi2",&chi2_save,"chi2/D");
  t_rec_charge->Branch("ndf",&ndf_save,"ndf/D");

  TGraph2D *g = new TGraph2D();
  TGraph2D *gt = new TGraph2D();
  TGraph2D *g_rec = new TGraph2D();


  //save results 
  for (int i=start_num;i!=end_num+1;i++){
    //truth
    CellChargeMap ccmap = truthtiling[i]->ccmap();
    for (auto it = ccmap.begin();it!=ccmap.end(); it++){
      Point p = it->first->center();
      x_save = i*0.32;
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
      x_save = i*0.32;
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
      if (charge>recon_threshold){
	//truth
	for (int k=0;k!=mcell->get_allcell().size();k++){
	  Point p = mcell->get_allcell().at(k)->center();
	  x_save = i*0.32;
	  y_save = p.y/units::cm;
	  z_save = p.z/units::cm;
	  charge_save = charge/mcell->get_allcell().size();
	  ncharge_save = mcell->get_allcell().size();
	  ncell_save = mcell->get_allcell().size();
	  chi2_save = toymatrix[i]->Get_Chi2();
	  ndf_save = toymatrix[i]->Get_ndf();

	  g_rec->SetPoint(ncount1,x_save,y_save,z_save);
	  t_rec_charge->Fill();
	  
	  ncount1 ++;
	}
      }
    }

    //save all results
    // file->Write(Form("toytiling_%d",i),toytiling[i]);
    // file->Write(Form("mergetiling_%d",i),mergetiling[i]);
    // file->Write(Form("truthtiling_%d",i),truthtiling[i]);
    // file->Write(Form("toymatrix_%d",i),toymatrix[i]);

  }
 
  

  g->Write("shower3D");
  gt->Write("shower3D_truth");
  g_rec->Write("shower3D_charge");

  const int N = 100000;
  Double_t x[N],y[N],z[N];
  //save cluster
  int ncluster = 0;
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    ncount = 0;
    for (int i=0; i!=(*it)->get_allcell().size();i++){
      const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
      for (int j=0; j!=mcell->get_allcell().size();j++){
  	Point p = mcell->get_allcell().at(j)->center();
  	x[ncount] = mcell->GetTimeSlice()*0.32;
  	y[ncount] = p.y/units::cm;
  	z[ncount] = p.z/units::cm;
  	ncount ++;
      }
    }
    // cout << ncount << endl;
    TGraph2D *g1 = new TGraph2D(ncount,x,y,z);
    g1->Write(Form("cluster_%d",ncluster));
    ncluster ++;
  }


  file->Write();
  file->Close();

  toymetric.Print();

  return 0;
  
} // main()
