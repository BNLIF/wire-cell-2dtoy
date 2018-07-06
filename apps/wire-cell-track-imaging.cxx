#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TMath.h"
using namespace std;

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
#include "TGraph.h"


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

  //////////////////////////////////////////////////////////////////////////////////////////////// variable histgram start

  roostr = "h1_ghost_track_length";
  TH1D *h1_ghost_track_length = new TH1D(roostr, roostr, 100, 0, 100);
  
  roostr = "graph_ghost_yz";
  TGraph *graph_ghost_yz = new TGraph();
  graph_ghost_yz->SetName(roostr);

  std::map< int, std::vector<int> > map_non_ghost_track_clusters;
  std::map<int,int> map_flag_cluster_used;

  TGraph2D *graph_track_each[21];
  for(int idx=0; idx<21; idx++) {
    roostr = TString::Format("graph_track_each_%2d", idx);
    graph_track_each[idx] = new TGraph2D();
    graph_track_each[idx]->SetName(roostr);
  }
  int count_graph_track_each[21] = {0};

  TGraph2D *graph_cluster_all = new TGraph2D();
  graph_cluster_all->SetName("graph_cluster_all");
  
  //////////////////////////////////////////////////////////////////////////////////////////////// variable histgram end
  
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

  int track_id_debug = 0;
  TTree *truth = (TTree*)f->Get("T_true");
  truth->SetBranchAddress("x",&x);
  truth->SetBranchAddress("y",&y);
  truth->SetBranchAddress("z",&z);
  truth->SetBranchAddress("track_id",&track_id_debug);
  
  TGraph2D *gh_truth = new TGraph2D();
  gh_truth->SetName("gh_truth");
  
  for(long ientry=0; ientry<truth->GetEntries(); ientry++) {
    truth->GetEntry(ientry);
    gh_truth->SetPoint( ientry, x,y,z );

    //////////
    count_graph_track_each[track_id_debug]++;
    graph_track_each[track_id_debug]->SetPoint( count_graph_track_each[track_id_debug]-1, x,y,z );    
  }
  
  
  ///////
  std::map<int, TVector3> vc_cluster_dir;
  std::map<int, TVector3> vc_cluster_start;
  std::map<int, TVector3> vc_cluster_end;
  std::map<int, double> vc_cluster_length;
  
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

    graph_cluster_all->SetPoint(i, x,y,z);
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
  double tolerance_parallel = delta_arc_length/vc_track_length[track_id];      
  tolerance_parallel = 0.06;
  double tolerance_distance2track = 2; // cm
  double tolerance_contain = 2;// cm
  double short_cluster = 5;
  
  std::map<int, int> vc_ghost_id;
  std::map<int, int> vc_ghost_angle;
  std::map<int, int> vc_ghost_distance;
  std::map<int, int> vc_ghost_contain;

  for( auto it_track=vc_track_dir.begin(); it_track!=vc_track_dir.end(); it_track++ ) {
    map_non_ghost_track_clusters[it_track->first].push_back(-1);//
  }
  
  for( auto it_cluster=vc_cluster_dir.begin(); it_cluster!=vc_cluster_dir.end(); it_cluster++ ) {// By default, all the cluster are ghosts
    vc_ghost_id[it_cluster->first] = it_cluster->first;
  }

  for( auto it_cluster=vc_cluster_dir.begin(); it_cluster!=vc_cluster_dir.end(); it_cluster++ ) {// each cluster cluster
    for( auto it_track=vc_track_dir.begin(); it_track!=vc_track_dir.end(); it_track++ ) {// each track cluster
      
      int track_id = it_track->first;
      int cluster_id = it_cluster->first;
      
      // if( it_cluster->first!=9 && it_cluster->first!=26 ) continue;
      // int id = cluster_id;
      // if( id!=13 && id!=16 && id!=18 && id!=19 ) continue;
      
      /// angle
      TVector3 v3_trackXcluster = vc_track_dir[track_id].Cross( vc_cluster_dir[cluster_id] );

      /// distance
      TVector3 v3_cluster_start_2_track_end = vc_track_end[track_id] - vc_cluster_start[cluster_id];
      TVector3 v3_startXend = v3_cluster_start_2_track_end.Cross( vc_track_dir[track_id] );
      double distance_startXend = v3_startXend.Mag();
      TVector3 v3_cluster_end_2_track_start = vc_track_start[track_id] - vc_cluster_end[cluster_id];
      TVector3 v3_endXstart = v3_cluster_end_2_track_start.Cross( vc_track_dir[track_id] );
      double distance_endXstart = v3_endXstart.Mag();
      double distance_center = 0.5*(distance_startXend + distance_endXstart);

      /// contain
      double value_track_start = vc_track_start[track_id].Dot( vc_track_dir[track_id] );//A
      double value_track_end = vc_track_end[track_id].Dot( vc_track_dir[track_id] );//B
      double value_cluster_start = vc_cluster_start[cluster_id].Dot( vc_track_dir[track_id] );//C
      double value_cluster_end = vc_cluster_end[cluster_id].Dot( vc_track_dir[track_id] );//D
      double min = (value_track_start>value_track_end)?(value_track_end-tolerance_contain):(value_track_start-tolerance_contain);
      double max = (value_track_start>value_track_end)?(value_track_start+tolerance_contain):(value_track_end+tolerance_contain);

      // cout<<TString::Format(" cluster %2d , track %2d, angleT %5.3f, angle %5.3f, distance start/end %8.2f %8.2f, contain min/max %8.2f %8.2f, start/end, %8.2f %8.2f",
      // 			    cluster_id, track_id,
      // 			    tolerance_parallel, v3_trackXcluster.Mag(),
      // 			    distance_startXend, distance_endXstart,
      // 			    min, max, value_cluster_start, value_cluster_end )<<endl;

      
      if( v3_trackXcluster.Mag()<tolerance_parallel || vc_cluster_length[cluster_id]<short_cluster ) {// bool: direction || bool:short_cluster

	// std::map<int, TVector3> vc_cluster_dir;
	// std::map<int, TVector3> vc_cluster_start;
	// std::map<int, TVector3> vc_cluster_end;
	// std::map<int, double> vc_cluster_length;
 
	// std::map<int, TVector3> vc_track_dir;
	// std::map<int, TVector3> vc_track_start;
	// std::map<int, TVector3> vc_track_end;
	// std::map<int, double> vc_track_length;
	
	if( distance_center<=tolerance_distance2track ) {// bool: distance

	  if( value_cluster_start>min && value_cluster_start<max && value_cluster_end>min && value_cluster_end<max ) {// bool: contain
	    // not ghost
	    vc_ghost_id.erase( cluster_id );

	    
	    if( map_flag_cluster_used.find( cluster_id )!=map_flag_cluster_used.end() ) {
	      map_flag_cluster_used[cluster_id] += 1;
	      continue;
	    }
	    
	    map_flag_cluster_used[cluster_id] = 1;
	    map_non_ghost_track_clusters[track_id].push_back( cluster_id );
	  }
	  
	}// bool: distance
      }// bool: direction       
    }// each cluster cluster
  }// each track cluster


  map<int, int>flag_track2recon;// 0: inefficiency, 1: broken, 2: good
  for( auto it=map_non_ghost_track_clusters.begin(); it!=map_non_ghost_track_clusters.end(); it++ ) {
    int user_size = it->second.size() - 1;

    cout<<it->first<<"\t"<<user_size<<endl;
      
    if( user_size==0 ) flag_track2recon[it->first] = 0;
    if( user_size>1  ) flag_track2recon[it->first] = 1;
    
    if( user_size==1 ) {
      flag_track2recon[it->first] = 2;

      /// get cluster_id
      int cluster_id = it->second.at(1);
      double cluster_length = vc_cluster_length[cluster_id];

      if( cluster_length<5 ) {// ghost ---> inefficinecy
	flag_track2recon[it->first] = 0;
	vc_ghost_id[cluster_id] = cluster_id;// pay attention
      }
      else if( cluster_length<95 ) {// broken
	flag_track2recon[it->first] = 1;
      }
      
    }
    
  }

  
  for( auto it=map_flag_cluster_used.begin(); it!=map_flag_cluster_used.end(); it++ ) {
    if(it->second > 1) {
      cout<<" -----------> WARNING : cluster "<<it->first<<" used "<<it->second<< " times." <<endl;
    }
  }
  
  cout<<endl;

  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  
  TCanvas *canv = new TCanvas("canv", "canv", 1000, 800);
  TH3D *h3 = new TH3D("h3", "h3", 100, 0, 300, 100, -150, 150, 100, 0, 1500);
  h3->SetStats(0);
  h3->Draw();
  h3->SetXTitle("X axis");
  h3->SetYTitle("Y axis");
  h3->SetZTitle("Z axis");

  int size_ghost = vc_ghost_id.size();
  cout<<" ---> number of ghost "<<size_ghost<<endl;

  long count_graph_ghost_yz = 0;
  for( auto it=vc_ghost_id.begin(); it!=vc_ghost_id.end(); it++ ) {
    int id = it->first;
    cout<<TString::Format( " ---> ghost %3d, length %7.4f", id, vc_cluster_length[id] )<<endl;

    int npts = xpt[id].size();
    for(int j=0; j<npts; j++) { // each point      
      double x = xpt[id].at(j);
      double y = ypt[id].at(j);
      double z = zpt[id].at(j);
      gh_ghost[id]->SetPoint(j, x,y,z);

      count_graph_ghost_yz++;
      graph_ghost_yz->SetPoint(count_graph_ghost_yz-1, z,y);      
    }


    // if( id!=9 && id!=26 ) continue;
    // if( id!=13 && id!=16 && id!=18 && id!=19 ) continue;
    gh_ghost[id]->Draw("same p");

    /////// ghost histgram
    /////// ghost histgram
    h1_ghost_track_length->Fill( vc_cluster_length[id] );
  }
  
  gh_truth->SetMarkerColor(kRed);
  gh_truth->Draw("same p");
    
  canv->SaveAs("canv.root");

  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  roostr = "canv_visuliaztion_tracks_cluster";
  TCanvas *canv_visuliaztion_tracks_cluster = new TCanvas(roostr, roostr, 1000, 800);

  roostr = "h3_visuliaztion_tracks_cluster";
  TH3D *h3_visuliaztion_tracks_cluster = new TH3D(roostr, roostr, 100, 0, 300, 100, -150, 150, 100, 0, 1500);
  h3_visuliaztion_tracks_cluster->SetStats(0);
  h3_visuliaztion_tracks_cluster->Draw();
  h3_visuliaztion_tracks_cluster->SetXTitle("X axis");
  h3_visuliaztion_tracks_cluster->SetYTitle("Y axis");
  h3_visuliaztion_tracks_cluster->SetZTitle("Z axis");

  int colors_vis[3] = {kRed, kOrange-3, kGreen};// inefficiency, broken, good
  
  for( auto it=vc_track_dir.begin(); it!=vc_track_dir.end(); it++ ) {
    int id = it->first;
    
    graph_track_each[id]->Draw("same p");
    graph_track_each[id]->SetMarkerColor( colors_vis[flag_track2recon[id]] );
    graph_track_each[id]->SetMarkerStyle(4);
    graph_track_each[id]->SetMarkerSize(0.5);
  }

  graph_cluster_all->Draw("same p");
  graph_cluster_all->SetMarkerColor(kGray+2);
  
  for( auto it=vc_ghost_id.begin(); it!=vc_ghost_id.end(); it++ ) {
    gh_ghost[it->first]->Draw("same p");
    gh_ghost[it->first]->SetMarkerColor(kBlue);
  }
  
  canv_visuliaztion_tracks_cluster->SaveAs("canv_visuliaztion_tracks_cluster.root");
  
  /////////////////////////////////////////////////////////////////////////////////// WRITE FILE
  /////////////////////////////////////////////////////////////////////////////////// WRITE FILE
  
  TFile* output = new TFile(outputroot, "RECREATE");
  
  h1_ghost_track_length->Write();
  graph_ghost_yz->Write();

  
  output->Close();

  
  return 0;

 
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
  
}
