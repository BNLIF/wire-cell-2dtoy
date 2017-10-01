#include "WireCell2dToy/uBooNE_light_reco.h"


#include "TH1S.h"
#include <iostream>

using namespace WireCell;


WireCell2dToy::uBooNE_light_reco::uBooNE_light_reco(const char* root_file){
  file = new TFile(root_file);
  T = (TTree*)file->Get("/Event/Sim");
  //  T->AddFile(root_file);
  // std::cout << T->GetEntries() << std::endl;
  hraw = new TH1F*[32];
  hdecon = new TH1F*[32];
  for (int i=0;i!=32;i++){
    hraw[i] = new TH1F(Form("hraw_%d",i),Form("hraw_%d",i),1500,0,1500);
    hdecon[i] = new TH1F(Form("hdecon_%d",i),Form("hdecon_%d",i),250,0,250);
  }
}

WireCell2dToy::uBooNE_light_reco::~uBooNE_light_reco(){
  for (int i=0;i!=32;i++){
    delete hraw[i];
    delete hdecon[i];
  }
  delete hraw;
  delete hdecon;
  
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

  std::vector<COphitSelection> ophits_group;
  COphitSelection left_ophits;
  
  for (int i=32;i!=op_femch->size();i++){
    COphit *op_hit = new COphit(op_femch->at(i), (TH1S*)op_wf->At(i), op_timestamp->at(i) - triggerTime, op_gain->at(op_femch->at(i)), op_gainerror->at(op_femch->at(i)));
    if (op_hit->get_type()){
      bool flag_used = false;
      if (ophits_group.size()==0){
	COphitSelection ophits;
	ophits.push_back(op_hit);
	ophits_group.push_back(ophits);
	flag_used = true;
      }else{
	for (size_t j=0; j!=ophits_group.size();j++){
	  for (size_t k=0; k!= ophits_group.at(j).size(); k++){
	    if (fabs(op_hit->get_time() - ophits_group.at(j).at(k)->get_time()) < 0.1 ){
	      ophits_group.at(j).push_back(op_hit);
	      flag_used = true;
	      break;
	    }
	  }
	  if (flag_used)
	    break;
	}
      }
      if (!flag_used){
	COphitSelection ophits;
	ophits.push_back(op_hit);
	ophits_group.push_back(ophits);
      }
    }else{
      left_ophits.push_back(op_hit);
    }
    //std::cout << op_hit->get_PE() << " " << op_hit->get_PE_err() << " " << op_hit->get_gain() << " " << op_hit->get_time() << std::endl;
    // std::cout << op_hit->get_ch_no() << " " << op_hit->get_peak() << " " << op_hit->get_baseline() << " " << op_hit->get_time() << " " << op_hit->get_integral() << " " << op_hit->get_gain() << std::endl;
    // op_hits.push_back(op_hit);
  }

  //int count = 0;
  for (size_t i=0;i!=left_ophits.size();i++){
    bool flag_used = false;
    for (size_t j=0; j!=ophits_group.size();j++){
      for (size_t k=0; k!= ophits_group.at(j).size(); k++){
	//std::cout << left_ophits.at(i)->get_time() - ophits_group.at(j).at(k)->get_time() << std::endl;
	if (fabs(left_ophits.at(i)->get_time() - ophits_group.at(j).at(k)->get_time())<0.1){
	  
	  ophits_group.at(j).push_back(left_ophits.at(i));
	  flag_used = true;
	  //	  count ++;
	  break;
	}
      }
      if (flag_used)
	break;
    }
  }
  //  std::cout << count << std::endl;
  // std::cout << ophits_group.size() << " " << left_ophits.size() << std::endl;
  //for (size_t i=0; i!=ophits_group.size();i++){
    //std::cout << i << " " << ophits_group.at(i).size() << std::endl;
    //  if (ophits_group.at(i).size()==1)
    //std::cout << ophits_group.at(i).at(0)->get_time() << " " << ophits_group.at(i).at(0)->get_PE() << std::endl;
  //}

  
  // form flash ....
  for (size_t j=0; j!=ophits_group.size();j++){
    Opflash *flash = new Opflash(ophits_group.at(j));
    //std::cout << flash->get_num_fired() << " " << flash->get_time() << " " << flash->get_total_PE() << std::endl;
    flashes.push_back(flash);
  }

  // fill in the rawl ... 
  
  for (int i=0;i!=32;i++){
    TH1S *hsignal = (TH1S*)op_wf->At(i);
    for (int j=0;j!=1500;j++){
      hraw[i]->SetBinContent(j+1,hsignal->GetBinContent(j+1)-2050);
    }
  }
  
  Process_beam_wfs();
  
}

void WireCell2dToy::uBooNE_light_reco::Process_beam_wfs(){
  // correct the baseline ...
  TH1F h1("h1","h1",200,-100,100);
  for (int i=0;i!=32;i++){
    h1.Reset();
    for (int j=0;j!=1500;j++){
      h1.Fill(hraw[i]->GetBinContent(j+1));
    }
    // double xq = 0.5;
    // double baseline;
    // h1.GetQuantiles(1,&baseline,&xq);
    double baseline = h1.GetMaximumBin()-100;
    //std::cout << h1.GetMaximum() << " " << baseline << std::endl;
    for (int j=0;j!=1500;j++){
      hraw[i]->SetBinContent(j+1,hraw[i]->GetBinContent(j+1)-baseline);
    }
  }

  // deconvolution ...

  // ROI finding

  // rebin ... 
  
}
