#include "WireCellSst/GeomDataSource.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"
//#include "WireCell2dToy/ToyLightReco.h"

#include "WireCell2dToy/ToyLibPe.h"

#include "TH1F.h"
#include "TH2F.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace WireCell;
using namespace std;

int main(int argc, char* argv[]){

  float shiftX = 147, shiftY = -1000, shiftZ = -5000;
  
  if(argc < 3){
    cerr << "usage: wire-cell-photonlib-pred /path/to/ChannelWireGeometry.txt sim.root ";
    return 1;
  }

  // global switch disabling the reference to histo (i.e. not in list of in-memory objects)
  TH1::AddDirectory(kFALSE);

  int setX =  128, minX = 0, maxX = 256, allX=1;
  
  for(int i=1; i!=argc; i++){
    switch(argv[i][1]){
    case 'a':
      setX = atoi(&argv[i][2]);
      break;
    case 'b':
      minX = atoi(&argv[i][2]);
      break;
    case 'c':
      maxX = atoi(&argv[i][2]);
      break;
    case 'd':
      allX = atoi(&argv[i][2]);
      break;
    }
  }

  
  //int eve_num = atoi(argv[3]);
  int eve_num = 50;

  const char* root_file = argv[2];
  TFile *file = new TFile(root_file);
  TTree *tout = (TTree*)file->Get("tout");
  vector<int> *trackId = new vector<int>;
  vector<float> *energy = new vector<float>;
  vector<float> *numElectrons = new vector<float>;
  vector<float> *x = new vector<float>;
  vector<float> *y = new vector<float>;
  vector<float> *z = new vector<float>;
  tout->SetBranchAddress("trackId",&trackId);
  tout->SetBranchAddress("energy",&energy);
  tout->SetBranchAddress("numElectrons",&numElectrons);
  tout->SetBranchAddress("x",&x);
  tout->SetBranchAddress("y",&y);
  tout->SetBranchAddress("z",&z);
  tout->GetEntry(eve_num);

  WireCell2dToy::ToyLibPe simCharge(root_file);

  std::vector<double> shifted_v = simCharge.getShiftedXYZQ(x,y,z,numElectrons,shiftX,shiftY,shiftZ);
  std::cout<<"shifted x="<<shifted_v[0]<<" shifted y="<<shifted_v[1]<<" shifted z="<<shifted_v[2]<<" q="<<shifted_v[3]<<std::endl;
    
  std::pair<float,float> shiftedMinMax = simCharge.getShiftedXminXmax(x,shiftX);
  std::cout<<"minimum shifted x: "<<shiftedMinMax.first<<" maximum shifted x:"<<shiftedMinMax.second<<std::endl;
  
  bool fidVol = simCharge.xInFidVol(shiftedMinMax);
  bool atBoundary = simCharge.xAtBoundary(shiftedMinMax);
  std::cout<<"Is cluster inside fiducial volume? \n"<< fidVol 
	   <<"\nIs cluster at boundary? \n"<< atBoundary<<std::endl;

  std::vector<double> pe_pred_v(32,0);
  pe_pred_v= simCharge.fromClusterFindPeDist(x,y,z,numElectrons,shiftX,shiftY,shiftZ);
  
 
  
  TFile *fout = new TFile("libexp.root","RECREATE");
  TTree *T_exp = new TTree("T_exp","T_exp");

  Double_t totalEnergy = 0;
  Double_t pe_pred[32] = {0.};
  Double_t xyzq[4] = {0.};

  T_exp->Branch("totalEnergy",&totalEnergy,"totalEnergy/D");
  T_exp->Branch("pe_pred",pe_pred,"pe_pred[32]/D");
  T_exp->Branch("xyzq",xyzq,"xyzq[4]/D");

  for(int a=0; a<(int)shifted_v.size(); a++){
    xyzq[a] = shifted_v.at(a);
  }
  for(int b=0; b<(int)pe_pred_v.size(); b++){
    pe_pred[b] = pe_pred_v.at(b);
  }

  int nE = (int)energy->size();
  for(int c=0; c<nE; c++){ totalEnergy += energy->at(c); }
  
  T_exp->Fill();
  fout->Write();
  fout->Close();

} // main
