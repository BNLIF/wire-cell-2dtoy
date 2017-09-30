#include "WireCell2dToy/uBooNE_light_reco.h"

//using namespace WireCell;
#include "TH1S.h"
#include <iostream>



WireCell2dToy::uBooNE_light_reco::uBooNE_light_reco(const char* root_file){
  file = new TFile(root_file);
  T = (TTree*)file->Get("/Event/Sim");
  //  T->AddFile(root_file);
  // std::cout << T->GetEntries() << std::endl;

 
}

WireCell2dToy::uBooNE_light_reco::~uBooNE_light_reco(){
  delete T;
  delete file;
}

void WireCell2dToy::uBooNE_light_reco::load_event(int eve_num){

  TClonesArray* op_wf = new TClonesArray;
  std::vector<int> *op_femch = new std::vector<int>;
  std::vector<double> *op_timestamp = new std::vector<double>;
  std::vector<double> *op_gain = new std::vector<double>;
  std::vector<double> *op_gainerror = new std::vector<double>;
  double triggerTime;

  T->SetBranchAddress("op_femch",&op_femch);
  T->SetBranchAddress("op_gain",&op_gain);
  T->SetBranchAddress("op_gainerror",&op_gainerror);
  T->SetBranchAddress("op_timestamp",&op_timestamp);
  T->SetBranchAddress("op_wf",&op_wf);
  T->SetBranchAddress("triggerTime",&triggerTime);

  T->SetBranchStatus("*",0);
  T->SetBranchStatus("op_femch",1);
  T->SetBranchStatus("op_gain",1);
  T->SetBranchStatus("op_gainerror",1);
  T->SetBranchStatus("op_timestamp",1);
  T->SetBranchStatus("op_wf",1);
  T->SetBranchStatus("triggerTime",1);

  //  std::cout << T->GetEntries() << std::endl;
  T->GetEntry(eve_num);
  //  std::cout << op_femch->size() << " " << op_gain->size() << std::endl;
  
  
}
