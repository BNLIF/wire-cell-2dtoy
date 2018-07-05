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

using namespace std;

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

  cout<<" ---> test track-imaging"<<endl;

  return 0;
}
