#include "WireCell2dToy/ToyTiling.h"
#include "WireCellData/MergeGeomCell.h"

#include "TApplication.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include <iostream>

using namespace WireCell;
using namespace std;

int main(int argc, char* argv[])
{
   if (argc < 2) {
    cerr << "usage: wire-cell-cluster /path/to/shower_3D.root" << endl;
    return 1;
  }
   
  TString filename = argv[1];
  TFile *file = new TFile(filename);
  TTree *T = (TTree*)file->Get("T");
  TTree *TC = (TTree*)file->Get("TC");
  
  cout << T->GetEntries() << " " << TC->GetEntries() << endl;
  
}
