#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TMath.h"

#include <iostream>
#include <map>

#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"

#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"

#include "TString.h"
#include "TGraph2D.h"
#include "TPolyLine3D.h"
using namespace std;


/*
  short cluster: satisfy "distance"



  Flag: multiple cluters to one track

*/


void set_plot_style()
{
  TStyle* myStyle = new TStyle("myStyle","My ROOT plot style");
  // plot style
  //myStyle->SetPalette(kInvertedDarkBodyRadiator); 
  //myStyle->SetPalette(kRainBow);
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  myStyle->SetNumberContours(NCont);
  myStyle->SetOptStat(0);
    
  myStyle->SetLabelFont(62,"xyz");
  myStyle->SetLabelSize(0.05,"xyz");
  myStyle->SetTitleFont(62, "xyz");
  myStyle->SetTitleSize(0.06,"xyz");
  myStyle->SetTitleOffset(1.3, "y");
  myStyle->SetTitleOffset(0.8, "x");
  myStyle->SetTitleOffset(1.5, "z");
    
  // only 5 in x to avoid label overlaps
  myStyle->SetNdivisions(505, "x");

  // set the margin sizes
  myStyle->SetPadTopMargin(0.05);
  myStyle->SetPadRightMargin(0.2); // increase for colz plots
  myStyle->SetPadBottomMargin(0.1);
  myStyle->SetPadLeftMargin(0.15);
  myStyle->SetFrameBorderMode(0);
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetCanvasColor(0);
  myStyle->SetPadColor(0);

  gROOT->SetStyle("myStyle");
  gROOT->ForceStyle();
}

int main(int argc, char* argv[])
{
  if(argc < 3){
    cout<<"Usage: wire-cell-track-eval input.root outputname.root [option: cluster_id]"<<endl;
    exit(1);
  }
  const char* inputroot = argv[1];
  const char* outputroot = argv[2];
  int cluster_check=-1; // specify cluster_id to check 
  if(argc==4) cluster_check=atoi(argv[3]);

  //cout<<" ---> test track-imaging A"<<endl;
  //cout<<" ---> test track-imaging B"<<endl;

  TString roostr = "";

  // rotation matrix to convert collection plane (default detector) coordinate to inductions plane
  TMatrixD Muplane(3,3);
  TMatrixD Mvplane(3,3);
  double angle = 60./180*TMath::Pi();
  double cosine = TMath::Cos(angle);
  double sine = TMath::Sin(angle);
  Muplane(0,0) = 1.;
  Muplane(0,1) = 0;
  Muplane(0,2) = 0;
  Muplane(1,0) = 0;
  Muplane(1,1) = cosine;
  Muplane(1,2) = sine;
  Muplane(2,0) = 0;
  Muplane(2,1) = -sine;
  Muplane(2,2) = cosine;

  Mvplane(0,0) = 1.;
  Mvplane(0,1) = 0;
  Mvplane(0,2) = 0;
  Mvplane(1,0) = 0;
  Mvplane(1,1) = cosine;
  Mvplane(1,2) = -sine;
  Mvplane(2,0) = 0;
  Mvplane(2,1) = sine;
  Mvplane(2,2) = cosine;

  // read in the point cloud for each cluster_id
  TFile* f = new TFile(inputroot, "READ");
  TTree* clusters = (TTree*)f->Get("T_cluster");
  Double_t x=0;
  Double_t y=0;
  Double_t z=0;
  Int_t cluster_id=0;
  clusters->SetBranchAddress("x",&x);
  clusters->SetBranchAddress("y",&y);
  clusters->SetBranchAddress("z",&z);
  clusters->SetBranchAddress("cluster_id",&cluster_id);

  std::map<int, std::vector<double>> xpt;
  std::map<int, std::vector<double>> ypt;
  std::map<int, std::vector<double>> zpt;

  
  TTree *truth = (TTree*)f->Get("T_true");
  truth->SetBranchAddress("x",&x);
  truth->SetBranchAddress("y",&y);
  truth->SetBranchAddress("z",&z);

  TGraph2D *gh_truth = new TGraph2D();
  gh_truth->SetName("gh_truth");
  
  for(long ientry=0; ientry<truth->GetEntries(); ientry++) {
    truth->GetEntry(ientry);
    gh_truth->SetPoint( ientry, x,y,z );
  }
  
  
  ///////
  std::map<int, TVector3> vc_cluster_dir;
  std::map<int, TVector3> vc_cluster_start;
  std::map<int, TVector3> vc_cluster_end;
  std::map<int, double> vc_cluster_length;
  
  TGraph2D *gh_tracks = new TGraph2D();

  TGraph2D *gh_ghost[100];
  for(int idx=0; idx<100; idx++) {
    gh_ghost[idx] = new TGraph2D();
    roostr = TString::Format("gh_ghost_%03d", idx);
  }

  
  cout<<endl;
  cout<<" -------> Points in cluster, entries: "<<clusters->GetEntries()<<endl;
  cout<<endl;
  
  for(long i=0; i<clusters->GetEntries();i++) {
    clusters->GetEntry(i);
    auto it = xpt.find(cluster_id);
    if( it == xpt.end() ){
      std::vector<double> xvec;
      xvec.push_back(x);
      std::vector<double> yvec;
      yvec.push_back(y);
      std::vector<double> zvec;
      zvec.push_back(z);
	
      xpt[cluster_id] = xvec;
      ypt[cluster_id] = yvec;
      zpt[cluster_id] = zvec;
    }
    else{
      xpt[cluster_id].push_back(x);
      ypt[cluster_id].push_back(y);
      zpt[cluster_id].push_back(z);
    }
    
    gh_tracks->SetPoint(i, x,y,z);
  }

  // PCA (Priciple Component Analysis): Matrix Decomposition --> Coordinate Rotation
  // n-dimentional points in nxn position covariance matrix;
  // diagonalization (rotation);
  // the biggest eigen value corresponds to the main axis --> primary direction of the point cloud
    
  // Step 1: calculate average position
  // Step 2: construct covariance matrix
  // merge the two steps into a single loop of the point cloud
  
  for( auto it=xpt.begin(); it!=xpt.end(); it++ ) { // each cluster
    int i = it->first;
  
    double npts = xpt[i].size();
    double xx=0;
    double yy=0;
    double zz=0;
    double xy=0;
    double yz=0;
    double zx=0;
    double mx=0;
    double my=0;
    double mz=0;
    for(int j=0; j<npts; j++) { // each point      
      double x = xpt[i].at(j);
      double y = ypt[i].at(j);
      double z = zpt[i].at(j);

      xx += x*x/npts;
      yy += y*y/npts;
      zz += z*z/npts;
      xy += x*y/npts;
      yz += y*z/npts;
      zx += z*x/npts;
      mx += x/npts;
      my += y/npts;
      mz += z/npts;

      //polyline[i]->SetPoint(j, x,y,z);
    } //each point   

    TMatrixDSym m(3);
    m(0,0) = xx - mx*mx; 
    m(0,1) = xy - mx*my;
    m(0,2) = zx - mz*mx;
    m(1,0) = xy - mx*my;
    m(1,1) = yy - my*my;
    m(1,2) = yz - my*mz;
    m(2,0) = zx - mz*mx;
    m(2,1) = yz - my*mz;
    m(2,2) = zz - mz*mz;

    TMatrixDSymEigen me(m);
    auto& eigenval = me.GetEigenValues();
    auto& eigenvec = me.GetEigenVectors();

    double maxeval =  eigenval(0);
    int maxevalaxis = 0;

    for(int k=1; k<3; k++) {
      if(eigenval(k)>maxeval)
	{
	  maxevalaxis = k;
	  maxeval = eigenval(k);
	}
    }
        
    // PCA main direction
    double pcadir_x = eigenvec(0, maxevalaxis);
    double pcadir_y = eigenvec(1, maxevalaxis);
    double pcadir_z = eigenvec(2, maxevalaxis);

    // track direction
    TVector3 pca_dir(pcadir_x, pcadir_y, pcadir_z);
    TVector3 track_dir;
    track_dir = (1./pca_dir.Mag())*pca_dir;
      
    // theta_y: 0-90 degree, upgoing
    if(track_dir.Y()<0) track_dir *= -1.0;
      
    // Loop over point cloud to find the edges (start, end) along PCA main direction        
    TVector3 start(xpt[i].at(0), ypt[i].at(0), zpt[i].at(0));
    TVector3 end(start);
    double start_proj = track_dir.Dot(start);
    double end_proj = track_dir.Dot(end);

    for(int j=1; j<npts; j++) {//each point
      TVector3 point(xpt[i].at(j), ypt[i].at(j), zpt[i].at(j));
      double point_proj = track_dir.Dot(point);
      if(point_proj < start_proj)
	{
	  start=point;
	  start_proj = point_proj;
	}
      if(point_proj > end_proj)
	{
	  end=point;
	  end_proj = point_proj;
	} 
    }//each point

    // the start, end points are more precise than the PCA direction in case of a big blob (isochronous track)

    ///////
    TVector3 cluster_dir = end - start;// direction is alway "Up" in Y-axis
    vc_cluster_dir[i] = cluster_dir;
    vc_cluster_dir[i] = 1./vc_cluster_dir[i].Mag() * vc_cluster_dir[i];
    if( vc_cluster_dir[i].Y()<0 ) vc_cluster_dir[i] *= -1.;
    
    if( start.Y()<=end.Y() ) {
      vc_cluster_start[i] = start;
      vc_cluster_end[i] = end;
    }
    else {
      vc_cluster_start[i] = end;
      vc_cluster_end[i] = start;
    }
    vc_cluster_length[i] = cluster_dir.Mag();
      
    ///////
    TVector3 cluster_dir_y = cluster_dir;
    if(cluster_dir_y.Y()<0) cluster_dir_y *= -1.0;
    
    TVector3 cluster_dir_u = Muplane*cluster_dir;
    if(cluster_dir_u.Y()<0) cluster_dir_u *= -1.0;

    TVector3 cluster_dir_v = Mvplane*cluster_dir;
    if(cluster_dir_v.Y()<0) cluster_dir_v *= -1.0;    
  }

  ///////////////////////////////////////////////////////////////
  ///////////// clusters in one track, and ghost ////////////////
  ///////////////////////////////////////////////////////////////

  cout<<endl<<" -------> processing clusters in one track, and ghost"<<endl<<endl;
  
  // std::map<int, TVector3> vc_cluster_dir;
  // std::map<int, TVector3> vc_cluster_start;
  // std::map<int, TVector3> vc_cluster_end;
  // std::map<int, double> vc_cluster_length;


  TTree* track_true = (TTree*)f->Get("T_track");
  double x0,y0,z0;
  double x1,y1,z1;
  int track_id;
  track_true->SetBranchAddress("x0",&x0);
  track_true->SetBranchAddress("y0",&y0);
  track_true->SetBranchAddress("z0",&z0);
  track_true->SetBranchAddress("x1",&x1);
  track_true->SetBranchAddress("y1",&y1);
  track_true->SetBranchAddress("z1",&z1);  
  track_true->SetBranchAddress("track_id",&track_id);
  std::map<int, TVector3> vc_track_dir;
  std::map<int, TVector3> vc_track_start;
  std::map<int, TVector3> vc_track_end;
  std::map<int, double> vc_track_length;
  int num_track_clusters = track_true->GetEntries();
  for(int ientry=0; ientry<num_track_clusters; ientry++) {
    track_true->GetEntry( ientry );
    
    TVector3 v3_start(x0, y0, z0);
    TVector3 v3_end(x1, y1, z1);

    TVector3 dir = v3_end - v3_start;
    vc_track_length[track_id] = dir.Mag();
    dir = 1./dir.Mag() * dir;
    if( dir.Y()<0 ) dir *= -1.;
    vc_track_dir[track_id] = dir;
    
    if( v3_start.Y() < v3_end.Y() ) {
      vc_track_start[track_id] = v3_start;
      vc_track_end[track_id] = v3_end;
    }
    else {
      vc_track_start[track_id] = v3_end;
      vc_track_end[track_id] = v3_start;
    }
  }
  
  int num_clusters = vc_cluster_dir.size();
  cout<<" -------> Cluster numbers: "<<num_clusters<<endl;
  cout<<" -------> Track numbers:   "<<num_track_clusters<<endl;
  cout<<endl;

  /*
    cell_unit = 3 mm
    ---> +/- uncertainty = 6 mm
    
    arc_length = radius * angle(in PI unit)

    delta_arc_length = sqrt( 6*6 + 6*6 ) = 8.5 mm ---> 0.85 cm

    tolerance = delta_arc_length/radius

    sin(0.05) = 0.04998
    sin(0.10) = 0.09983
    sin(0.20) = 0.19867
    sin(0.50) = 0.47943
  */

  const double delta_arc_length = 0.85;
  //double tolerance_parallel = delta_arc_length/vc_track_length[track_id];      
  double tolerance_distance2track = 2; // cm
  double tolerance_contain = 2;// cm
  double short_cluster = 5;
  
  std::map<int, int> vc_ghost_id;
  std::map<int, int> vc_ghost_angle;
  std::map<int, int> vc_ghost_distance;
  std::map<int, int> vc_ghost_contain;
  
  for( auto it_cluster=vc_cluster_dir.begin(); it_cluster!=vc_cluster_dir.end(); it_cluster++ ) {// By default, all the cluster are ghosts
    vc_ghost_id[it_cluster->first] = it_cluster->first;
  }

  for( auto it_track=vc_track_dir.begin(); it_track!=vc_track_dir.end(); it_track++ ) {// each track cluster
    for( auto it_cluster=vc_cluster_dir.begin(); it_cluster!=vc_cluster_dir.end(); it_cluster++ ) {// each cluster cluster

      int track_id = it_track->first;
      int cluster_id = it_cluster->first;
      
      double tolerance_parallel = delta_arc_length/vc_track_length[track_id];      
      tolerance_parallel = 5*tolerance_parallel;

      TVector3 v3_trackXcluster = vc_track_dir[track_id].Cross( vc_cluster_dir[cluster_id] );
      
      if( v3_trackXcluster.Mag()<tolerance_parallel || vc_track_length[cluster_id]<short_cluster ) {// bool: direction || bool:short_cluster
	// cout<<TString::Format(
	// 		      " Tolerance %6.4f, Parallel of t/r  %3d %3d, angle between t/r %6.4f",
	// 		      tolerance_parallel, track_id, cluster_id, v3_trackXcluster.Mag()
	// 		      )<<endl;
	// vc_track_dir[track_id].Print();
	// vc_cluster_dir[cluster_id].Print();

  
	// std::map<int, TVector3> vc_cluster_dir;
	// std::map<int, TVector3> vc_cluster_start;
	// std::map<int, TVector3> vc_cluster_end;
	// std::map<int, double> vc_cluster_length;
 
	// std::map<int, TVector3> vc_track_dir;
	// std::map<int, TVector3> vc_track_start;
	// std::map<int, TVector3> vc_track_end;
	// std::map<int, double> vc_track_length;
	
	TVector3 v3_cluster_start_2_track_end = vc_track_end[track_id] - vc_cluster_start[cluster_id];
	TVector3 v3_startXend = v3_cluster_start_2_track_end.Cross( vc_track_dir[track_id] );
	double distance_startXend = v3_startXend.Mag();

	TVector3 v3_cluster_end_2_track_start = vc_track_start[track_id] - vc_cluster_end[cluster_id];
	TVector3 v3_endXstart = v3_cluster_end_2_track_start.Cross( vc_track_dir[track_id] );
	double distance_endXstart = v3_endXstart.Mag();

	if( distance_startXend<=tolerance_distance2track && distance_endXstart<=tolerance_distance2track ) {// bool: distance

	  double value_track_start = vc_track_start[track_id].Dot( vc_track_dir[track_id] );//A
	  double value_track_end = vc_track_end[track_id].Dot( vc_track_dir[track_id] );//B

	  double value_cluster_start = vc_cluster_start[cluster_id].Dot( vc_track_dir[track_id] );//C
	  double value_cluster_end = vc_cluster_end[cluster_id].Dot( vc_track_dir[track_id] );//D

	  double min = (value_track_start>value_track_end)?(value_track_end-tolerance_contain):(value_track_start-tolerance_contain);
	  double max = (value_track_start>value_track_end)?(value_track_start+tolerance_contain):(value_track_end+tolerance_contain);
	  
	  if( value_cluster_start>min && value_cluster_start<max && value_cluster_end>min && value_cluster_end<max ) {// bool: contain
	    // not ghost
	    vc_ghost_id.erase( cluster_id );
	  }
	  
	}// bool: distance
      }// bool: direction
       
    }// each cluster cluster
  }// each track cluster


  cout<<endl;

  TCanvas *canv = new TCanvas("canv", "canv", 1000, 800);
  TH3D *h3 = new TH3D("h3", "h3", 100, 0, 300, 100, -150, 150, 100, 0, 1500);
  h3->SetStats(0);
  h3->Draw();
  h3->SetXTitle("X axis");
  h3->SetYTitle("Y axis");
  h3->SetZTitle("Z axis");

  int size_ghost = vc_ghost_id.size();
  cout<<" ---> number of ghost "<<size_ghost<<endl;
  for( auto it=vc_ghost_id.begin(); it!=vc_ghost_id.end(); it++ ) {
    int id = it->first;
    cout<<TString::Format( " ---> ghost %3d, length %7.4f", id, vc_cluster_length[id] )<<endl;

    int npts = xpt[id].size();
    for(int j=0; j<npts; j++) { // each point      
      double x = xpt[id].at(j);
      double y = ypt[id].at(j);
      double z = zpt[id].at(j);
      gh_ghost[id]->SetPoint(j, x,y,z);
    }

    gh_ghost[id]->Draw("same p");
  }
  
  gh_truth->SetMarkerColor(kRed);
  gh_truth->Draw("same p");
    
  canv->SaveAs("canv.root");

  
 
  // numbers of cycle = C_n^2 =  num_clusters*(num_clusters-1)/2 

  // int count_it_a = 0;
  // for( auto it_a=vc_cluster_dir.begin(); it_a!=std::prev( vc_cluster_dir.end(), 1); it_a++ ) {
  //   int cluster_id_a = it_a->first;
  //   cout<<TString::Format(" ---> cluster and length: %3d, %6.2f", cluster_id_a, vc_cluster_length[cluster_id_a] )<<endl;

  //   ///////
  //   count_it_a++;
  //   for( auto it_b=std::prev( vc_cluster_dir.end(), vc_cluster_dir.size()-count_it_a ); it_b!=vc_cluster_dir.end(); it_b++ ) {
  //     int cluster_id_b = it_b->first;
  //     cout<<TString::Format("   -------> %3d, %6.2f", cluster_id_b, vc_cluster_length[cluster_id_b])<<endl;
  //   }
  // }

  // int count_it_a = 0;
  // for( auto it_a=vc_cluster_dir.begin(); it_a!=vc_cluster_dir.end(); it_a++ ) {
  //   int cluster_id_a = it_a->first;
  //   //cout<<TString::Format(" ---> cluster and length: %3d, %6.2f", cluster_id_a, vc_cluster_length[cluster_id_a] )<<endl;

  //   TVector3 vector_unit = 1./vc_cluster_dir[cluster_id_a].Mag() * vc_cluster_dir[cluster_id_a];
  //   vector_unit.Print();
    
  //   for( auto it_b=vc_cluster_dir.begin(); it_b!=vc_cluster_dir.end(); it_b++ ) {
  //     int cluster_id_b = it_b->first;
  //     if( cluster_id_a==cluster_id_b ) continue;
  //     //cout<<TString::Format("   -------> %3d, %6.2f", cluster_id_b, vc_cluster_length[cluster_id_b])<<endl;
      
  //     TVector3 vector_a_unit = 1./vc_cluster_dir[cluster_id_a].Mag() * vc_cluster_dir[cluster_id_a];
  //     TVector3 vector_b_unit = 1./vc_cluster_dir[cluster_id_b].Mag() * vc_cluster_dir[cluster_id_b];

  //     TVector3 vector_axb = vector_a_unit.Cross( vector_b_unit );
  //     if( vector_axb.Mag()==0 ) cout<<TString::Format(" ---> Parallel of cluster %3d and %3d", cluster_id_a, cluster_id_b )<<endl;
      
  //   }
  // }


  
  TFile* output = new TFile(outputroot, "RECREATE");
  gh_tracks->Write();
  // for(int idx=0; idx<vc_cluster_dir.size(); idx++) {
  //   polyline[idx]->Write();
  // }
  output->Close();

  
  return 0;

}
