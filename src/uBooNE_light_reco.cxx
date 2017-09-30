#include "WireCell2dToy/uBooNE_light_reco.h"


#include "TH1S.h"
#include <iostream>

using namespace WireCell;


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
  
  for (Int_t i=32;i!=op_femch->size();i++){
    COphit *op_hit = new COphit(op_femch->at(i), (TH1S*)op_wf->At(i), op_timestamp->at(i) - triggerTime, op_gain->at(op_femch->at(i)), op_gainerror->at(op_femch->at(i)));
    std::cout << op_hit->get_PE() << " " << op_hit->get_PE_err() << " " << op_hit->get_gain() << std::endl;
    // std::cout << op_hit->get_ch_no() << " " << op_hit->get_peak() << " " << op_hit->get_baseline() << " " << op_hit->get_time() << " " << op_hit->get_integral() << " " << op_hit->get_gain() << std::endl;
    op_hits.push_back(op_hit);
  }
  
}
