#include "WireCellNav/GeomDataSource.h"
#include "WireCellSst/GeomWireReader.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCellData/MergeGeomCell.h"

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/SimDataSource.h"
#include "WireCellNav/SliceDataSource.h"
#include "WireCellSst/Util.h"
#include "WireCellData/SimTruth.h"
#include "WireCell2dToy/ToyDepositor.h"
#include "WireCellNav/GenerativeFDS.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TSystem.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include <iostream>

using namespace WireCell;
using namespace std;

int main(int argc, char* argv[])
{
  // Setup environment
  gErrorIgnoreLevel = kError;

  // Application instructions
  if (argc < 3) {
    cerr << "usage:  wire-cell-uboone-truth-MM /path/to/ChannelWireGeometry.txt /path/to/celltree.root" << endl;
    return 1;
  }

  // Set additional user-defined parameters
  const Int_t runMode = 0;
  const Int_t firstSlice = 460;
  const Int_t lastSlice = 560;
  const Int_t sliceStep = 10;
  //const Double_t lowZ = 0.0;
  //const Double_t highZ = 10.3698;
  //const Double_t lowY = -1.165;
  //const Double_t highY = 1.165;
  const Double_t lowZ = 1.5;
  const Double_t highZ = 2.5;
  const Double_t lowY = 0.5;
  const Double_t highY = 1.0;
  const Int_t binsPerFrame = 2400;
  const Int_t totalFrames = 5;
  const Int_t elecThreshold = 2000;

  // Define diagnostic histograms
  TH2F *eigenValHist = new TH2F("eigenValHist","",15,0,30,15,0,30);
  eigenValHist->GetXaxis()->SetTitle("Number of Cells");
  eigenValHist->GetYaxis()->SetTitle("Number of Non-zero Eigenvalues");
  TH1F *passChargeRecoRes = new TH1F("passChargeRecoRes","",40,-2.0,2.0);
  passChargeRecoRes->GetXaxis()->SetTitle("(recoCharge-trueCharge)/trueCharge");
  passChargeRecoRes->GetYaxis()->SetTitle("# of Events");
  TH1F *failChargeRecoRes = new TH1F("failChargeRecoRes","",40,-2.0,2.0);
  failChargeRecoRes->GetXaxis()->SetTitle("(recoCharge-trueCharge)/trueCharge");
  failChargeRecoRes->GetYaxis()->SetTitle("# of Events");
  TH1F *totalChargeRecoRes = new TH1F("totalChargeRecoRes","",40,-2.0,2.0);
  totalChargeRecoRes->GetXaxis()->SetTitle("(recoCharge-trueCharge)/trueCharge");
  totalChargeRecoRes->GetYaxis()->SetTitle("# of Events");

  // Setup display options
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
  
  // Get GDS from geometry text file
  WireCellSst::GeomWireReader reader(argv[1]);
  WireCell::GeomDataSource gds;
  gds.use_wires(reader);
  
  // Get FDS from ROOT file (Chao's data format for now)
  const char* root_file = argv[2];
  const char* tpath = "/Event/Sim";
  WireCell::FrameDataSource* fds = 0;
  fds = WireCellSst::make_fds(root_file);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  
  // Make GenerativeFDS to represent reco info using SimChannels for now (eventually replace with actual CalWire output)
  WireCell::ToyDepositor toydep(fds);
  const PointValueVector pvv = toydep.depositions(1);
  WireCell::GenerativeFDS gfds(toydep,gds,binsPerFrame,totalFrames,2.0*1.6*units::millimeter);
  gfds.jump(1); // NB: look at first frame only for now

  // Create MicroBooNE SDS using GenerativeFDS (used for per-slice tiling)
  WireCellSst::ToyuBooNESliceDataSource sds(gfds,elecThreshold);

  // Start interactive application used for plotting during program execution
  TApplication theApp("theApp",0,0); 
  theApp.SetReturnFromRun(true);

  // Create toy display for plotting
  TCanvas c1("ToyMC","ToyMC",1200,600);
  WireCell2dToy::ToyEventDisplay display(c1,gds);

  // Main loop over slices
  const int N = 100000;
  Double_t x[N],y[N],z[N];
  Double_t xt[N],yt[N],zt[N];
  int ncount = 0;
  int ncount_t = 0;
  int numSolved1 = 0;
  int numSolved2 = 0;
  int numTotal = 0;
  for (int i = 0; i != sds.size(); i++){
    // Get next slice
    sds.jump(i);
    WireCell::Slice slice = sds.get();

    // Look at slice if there is above-threshold activity
    if (slice.group().size() > 0){
      // Do tiling for all possible hits (fake and real hits)
      WireCell2dToy::ToyTiling toytiling(slice,gds);
      // Do merged tiling (Xin's blob idea) for all possible hits (fake and real hits)
      WireCell2dToy::MergeToyTiling mergetiling(toytiling);
      // Do tiling using truth info (real hits only) - associates 2D points to closest {U,V,Y} wires
      WireCell2dToy::TruthToyTiling truthtiling(toytiling,pvv,i,gds);
     
      // Collect cell and wire information
      GeomCellSelection allcell = toytiling.get_allcell();
      GeomWireSelection allwire = toytiling.get_allwire();
      GeomCellSelection allmcell = mergetiling.get_allcell();
      GeomWireSelection allmwire = mergetiling.get_allwire();

      // Get wire-to-cells mapping
      GeomWireMap wmap = toytiling.wmap();

      // Get cell-to-wires mapping
      GeomCellMap cmap = toytiling.cmap();

      // Get wire charge map (using only SimChannel information for now)
      WireChargeMap wcmap = toytiling.wcmap();

      // Get true cell charge map (real hits) using truth info tiling
      CellChargeMap ccmap = truthtiling.ccmap();

      // Block for solving matrix equation
      if((allcell.size() > 0) && (allwire.size() > 0)){
        // Initialize matrices for Wire-Cell inversion
        TMatrixD wireMat;
        wireMat.ResizeTo(allwire.size(),1);
        TMatrixD cellMat;
        cellMat.ResizeTo(allcell.size(),1);
        TMatrixD geoMat;
        geoMat.ResizeTo(allwire.size(),allcell.size());
        TMatrixD transpGeoMat;
        transpGeoMat.ResizeTo(allcell.size(),allwire.size());
        TMatrixD tempGeoMat1;
        tempGeoMat1.ResizeTo(allcell.size(),allcell.size());
        TMatrixD tempInvGeoMat1;
        tempInvGeoMat1.ResizeTo(allcell.size(),allcell.size());
        TMatrixD tempGeoMat2;
        tempGeoMat2.ResizeTo(allwire.size(),allwire.size());
        TMatrixD tempInvGeoMat2;
        tempInvGeoMat2.ResizeTo(allwire.size(),allwire.size());
        TMatrixD inverseGeoMat;
        inverseGeoMat.ResizeTo(allcell.size(),allwire.size());
        TMatrixD crosscheckMat;
        crosscheckMat.ResizeTo(allwire.size(),1);
        TMatrixD resultMat;
        resultMat.ResizeTo(allcell.size(),1);
        
        // Fill matrices
        for (int k = 0; k < allwire.size(); k++){
          if (k < allwire.size()){
            wireMat[k][0] = wcmap[allwire[k]];
          }
          else{
            wireMat[k][0] = 0.0;
          }
        
          for (int h = 0; h < allcell.size(); h++){
            if (k == 0){
              if (h < allcell.size()){
                cellMat[h][0] = ccmap[allcell[h]];
              }
              else{
                cellMat[h][0] = 0.0;
              }
            }
        
            if ((k < allwire.size()) && (h < allcell.size()) && ((allwire[k] == cmap[allcell[h]].at(0)) || (allwire[k] == cmap[allcell[h]].at(1)) || (allwire[k] == cmap[allcell[h]].at(2)))){
              geoMat[k][h] = 1.0;
            }
            else{
              geoMat[k][h] = 0.0;
            }
          }
        }
        
        // Do matrix calculations and check eigenvalues
        transpGeoMat.Transpose(geoMat);
        tempGeoMat1 = transpGeoMat*geoMat;
        tempGeoMat2 = geoMat*transpGeoMat;
        TMatrixDEigen eigenValMat1(tempGeoMat1);
        TVectorD eigenValsRe1 = eigenValMat1.GetEigenValuesRe();
        TVectorD eigenValsIm1 = eigenValMat1.GetEigenValuesIm();
        Int_t numNonzeroEigenVals1 = 0;
        for (int k = 0; k < allcell.size(); k++){
          if (fabs(eigenValsRe1[k]) > 0.0001){
            numNonzeroEigenVals1++;
          }
        }
        TMatrixDEigen eigenValMat2(tempGeoMat2);
        TVectorD eigenValsRe2 = eigenValMat2.GetEigenValuesRe();
        TVectorD eigenValsIm2 = eigenValMat2.GetEigenValuesIm();
        Int_t numNonzeroEigenVals2 = 0;
        for (int k = 0; k < allwire.size(); k++){
          if (fabs(eigenValsRe2[k]) > 0.0001){
            numNonzeroEigenVals2++;
          }
        }

        // Do matrix inversion
        //if(numNonzeroEigenVals1 >= allcell.size()){ // Over-constrained method (left inverse)
          tempInvGeoMat1 = tempGeoMat1;
          tempInvGeoMat1.Invert();
          inverseGeoMat = tempInvGeoMat1*transpGeoMat;
          resultMat = inverseGeoMat*wireMat;
          if(numNonzeroEigenVals1 >= allcell.size()){
            numSolved1++;
          }
        //}
        //else{// if(numNonzeroEigenVals2 >= allwire.size()){ // Under-constrained Method (right inverse)
        //  tempInvGeoMat2 = tempGeoMat2;
        //  tempInvGeoMat2.Invert();
        //  inverseGeoMat = transpGeoMat*tempInvGeoMat2;
        //  resultMat = inverseGeoMat*wireMat;
        //  numSolved2++;
        //}
        crosscheckMat = geoMat*cellMat;
        numTotal++;
        
        // Do validation of matrix multiplication/inversion
        cout << "NUMWIRES:  " << allwire.size() << endl;
        cout << "NUMCELLS:  " << allcell.size() << endl;
        cout << "NONZEROEIGENVALNUM1:  " << numNonzeroEigenVals1 << " out of " << allcell.size() << endl;
        cout << "NONZEROEIGENVALNUM2:  " << numNonzeroEigenVals2 << " out of " << allwire.size() << endl;
        for (int k = 0; k < allwire.size(); k++){
          cout << wireMat[k][0] << " " << crosscheckMat[k][0] << endl;
        }
        for (int k = 0; k < allwire.size(); k++){
          for (int h = 0; h < allcell.size(); h++){
            cout << geoMat[k][h] << " ";
          }
          cout << endl;
        }
        for (int h = 0; h < allcell.size(); h++){
          cout << cellMat[h][0] << " " << resultMat[h][0] << endl;
        }

        // Fill diagnostic plots
        eigenValHist->Fill(allcell.size(),numNonzeroEigenVals1);
        for (int h = 0; h < allcell.size(); h++){
          if (cellMat[h][0] > 0.0){
            totalChargeRecoRes->Fill((resultMat[h][0]-cellMat[h][0])/cellMat[h][0]);
	  }
        }
        if(numNonzeroEigenVals1 >= allcell.size()){
          for (int h = 0; h < allcell.size(); h++){
            if (cellMat[h][0] > 0.0){
              passChargeRecoRes->Fill((resultMat[h][0]-cellMat[h][0])/cellMat[h][0]);
	    }
          }
	}
	else{
          for (int h = 0; h < allcell.size(); h++){
            if (cellMat[h][0] > 0.0){
              failChargeRecoRes->Fill((resultMat[h][0]-cellMat[h][0])/cellMat[h][0]);
	    }
          }
	}
      }

      // Grab {x,y,z} points in TPC for all possible hits (fake and real hits)
      for (int j = 0; j != allcell.size(); j++){
        // Get center coordinates of cell
	Point p = allcell[j]->center();
	x[ncount] = i*0.32*units::cm; // NB: drift distance for 2 us spacing
	y[ncount] = p.y/units::cm;
	z[ncount] = p.z/units::cm;
	ncount++;
      }

      // Grab {x,y,z} points in TPC using truth info (real hits only)
      Double_t charge_min = 10000000.0;
      Double_t charge_max = -10000000.0;
      for (auto it = ccmap.begin(); it != ccmap.end(); it++){
        // Get center coordinates of cell
	Point p = it->first->center();
      	xt[ncount_t] = i*0.32*units::cm; // NB: drift distance for 2 us spacing
      	yt[ncount_t] = p.y/units::cm;
      	zt[ncount_t] = p.z/units::cm;
      	ncount_t++;

        // Check for minimal and maximal charge values for plotting purposes
	float charge = it->second;
	if (charge > charge_max) charge_max = charge;
	if (charge < charge_min) charge_min = charge;
      }
      
      // Do plotting interactively
      if((((i >= firstSlice) && (i <= lastSlice) && ((i-firstSlice) % sliceStep == 0)) && (runMode == 2)) || ((i == firstSlice) && (runMode == 1))){
        // Update charge bounds
        display.charge_min = charge_min;
        display.charge_max = charge_max;
        display.init(lowZ,highZ,lowY,highY);

        // Draw cells and wires
        c1.Draw();
	display.draw_mc(1,WireCell::PointValueVector(),"colz");           
        display.draw_slice(slice,"");
        display.draw_cells(toytiling.get_allcell(),"*same");
        display.draw_mergecells(mergetiling.get_allcell(),"*same");
        display.draw_truthcells(ccmap,"*same");

        // Update canvas and wait for input
        c1.Update();
        if(runMode == 1)
	  theApp.Run();
        else if(runMode == 2)
          c1.WaitPrimitive();
      }
    }
  }

  // Save output for 3D display
  TGraph2D *g = new TGraph2D(ncount,x,y,z);
  g->SetName("g");
  TGraph2D *gt = new TGraph2D(ncount_t,xt,yt,zt);
  gt->SetName("gt");
  TFile *file = new TFile("eventInfo.root","RECREATE");
  g->Write("recoHits");
  gt->Write("truthHits");
  eigenValHist->Write();
  passChargeRecoRes->Write();
  failChargeRecoRes->Write();
  totalChargeRecoRes->Write();
  file->Write();
  file->Close();

  // End of application
  cout << "///////////////////////////////////////////////////////" << endl;
  cout << "//// NUM SOLVED1:  " << numSolved1 << endl;
  cout << "//// NUM SOLVED2:  " << numSolved2 << endl;
  cout << "//// NUM TOTAL:    " << numTotal << endl;
  cout << "///////////////////////////////////////////////////////" << endl;
  return 0;
}
