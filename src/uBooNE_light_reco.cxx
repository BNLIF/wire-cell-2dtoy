#include "WireCell2dToy/uBooNE_light_reco.h"
#include "WireCellRess/LassoModel.h"
#include "WireCellRess/ElasticNetModel.h"
#include <Eigen/Dense>

#include "TH1S.h"
#include "TF1.h"
#include "TVirtualFFT.h"
#include <iostream>

using namespace WireCell;
using namespace Eigen;


WireCell2dToy::uBooNE_light_reco::uBooNE_light_reco(const char* root_file){
  file = new TFile(root_file);
  T = (TTree*)file->Get("/Event/Sim");
  //  T->AddFile(root_file);
  // std::cout << T->GetEntries() << std::endl;
  hraw = new TH1F*[32];
  hdecon = new TH1F*[32];
  hl1 = new TH1F*[32];
  for (int i=0;i!=32;i++){
    hraw[i] = new TH1F(Form("hraw_%d",i),Form("hraw_%d",i),1500,0,1500);
    hdecon[i] = new TH1F(Form("hdecon_%d",i),Form("hdecon_%d",i),250,0,250);
    hl1[i] = new TH1F(Form("hl1_%d",i),Form("hl1_%d",i),250,0,250);
  }
  h_totPE = new TH1F("h_totPE","h_totPE",250,0,250);
  h_mult = new TH1F("h_mult","h_mult",250,0,250);
  h_l1_mult = new TH1F("h_l1_mult","h_l1_mult",250,0,250);
  h_l1_totPE = new TH1F("h_l1_totPE","h_l1_totPE",250,0,250);
}

WireCell2dToy::uBooNE_light_reco::~uBooNE_light_reco(){
  for (int i=0;i!=32;i++){
    delete hraw[i];
    delete hdecon[i];
    delete hl1[i];
  }
  delete hraw;
  delete hdecon;
  delete hl1;

  delete h_totPE;
  delete h_mult;
  delete h_l1_mult;
  delete h_l1_totPE;
  
  delete T;
  delete file;
}

void WireCell2dToy::uBooNE_light_reco::load_event_raw(int eve_num){

  TClonesArray* cosmic_hg_wf = new TClonesArray;
  TClonesArray* cosmic_lg_wf = new TClonesArray;
  TClonesArray* beam_hg_wf = new TClonesArray;
  TClonesArray* beam_lg_wf = new TClonesArray;
  std::vector<short> *cosmic_hg_opch = new std::vector<short>;
  std::vector<short> *cosmic_lg_opch = new std::vector<short>;
  std::vector<short> *beam_hg_opch = new std::vector<short>;
  std::vector<short> *beam_lg_opch = new std::vector<short>;
  std::vector<double> *cosmic_hg_timestamp = new std::vector<double>;
  std::vector<double> *cosmic_lg_timestamp = new std::vector<double>;
  std::vector<double> *beam_hg_timestamp = new std::vector<double>;
  std::vector<double> *beam_lg_timestamp = new std::vector<double>;
  std::vector<float> *op_gain = new std::vector<float>;
  std::vector<float> *op_gainerror = new std::vector<float>;
  std::vector<int> *opch_to_opdet = new std::vector<int>;
  double triggerTime;

  T->SetBranchAddress("cosmic_hg_wf",&cosmic_hg_wf);
  T->SetBranchAddress("cosmic_lg_wf",&cosmic_lg_wf);
  T->SetBranchAddress("beam_hg_wf",&beam_hg_wf);
  T->SetBranchAddress("beam_lg_wf",&beam_lg_wf);
  T->SetBranchAddress("cosmic_hg_opch",&cosmic_hg_opch);
  T->SetBranchAddress("cosmic_lg_opch",&cosmic_lg_opch);
  T->SetBranchAddress("beam_hg_opch",&beam_hg_opch);
  T->SetBranchAddress("beam_lg_opch",&beam_lg_opch);
  T->SetBranchAddress("cosmic_hg_timestamp",&cosmic_hg_timestamp);
  T->SetBranchAddress("cosmic_lg_timestamp",&cosmic_lg_timestamp);
  T->SetBranchAddress("beam_hg_timestamp",&beam_hg_timestamp);
  T->SetBranchAddress("beam_lg_timestamp",&beam_lg_timestamp);
  T->SetBranchAddress("op_gain",&op_gain);
  T->SetBranchAddress("op_gainerror",&op_gainerror);
  T->SetBranchAddress("opch_to_opdet",&opch_to_opdet);
  T->SetBranchAddress("triggerTime",&triggerTime);

  T->GetEntry(eve_num);
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
    cosmic_flashes.push_back(flash);
  }
 

  // fill in the raw ... 
  for (int i=0;i!=32;i++){
    TH1S *hsignal = (TH1S*)op_wf->At(i);
    for (int j=0;j!=1500;j++){
      hraw[i]->SetBinContent(j+1,hsignal->GetBinContent(j+1)-2050);
    }
    gain[i] = op_gain->at(i);
    beam_dt[i] = op_timestamp->at(i) - triggerTime;
  }
  
  
  Process_beam_wfs();

  // sort the flashes ... 
  sort_flashes();
}
void WireCell2dToy::uBooNE_light_reco::sort_flashes(){
  OpFlashSet cosmic_set;
  for (auto it= cosmic_flashes.begin(); it!= cosmic_flashes.end(); it++){
    cosmic_set.insert(*it);
  }
  cosmic_flashes.clear();
  std::copy(cosmic_set.begin(), cosmic_set.end(), std::back_inserter(cosmic_flashes));

  OpFlashSet beam_set;
  for (auto it=beam_flashes.begin(); it!=beam_flashes.end(); it++){
    beam_set.insert(*it);
  }
  beam_flashes.clear();
  std::copy(beam_set.begin(), beam_set.end(), std::back_inserter(beam_flashes));

  // std::cout << flashes.size() << std::endl;
  OpFlashSet all_set;
  for (auto it=flashes.begin(); it!=flashes.end(); it++){
    all_set.insert(*it);
  }
  flashes.clear();
  std::copy(all_set.begin(), all_set.end(), std::back_inserter(flashes));
  // std::cout << flashes.size() << std::endl;
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
  double rc_tau[32];
  for (int i=0;i!=32;i++){
    rc_tau[i] = 800;
  }
  rc_tau[28] = 28.6;
  
  TF1 f1("f1","1./19.6348*pow(x/[0],4)*exp(-x/[0])",0,1500); //RC-CR4 response 
  double par[4];
  par[0] = 8.18450e-01;
  f1.SetParameters(par);

  TF1 f2("f2","exp(-pow(x/[0],[1]))",0,1); // high frequency filter
  // par[0] = 0.72;
  // par[1] = 5.086;
  par[0] = 0.45;
  par[1] = 3.07;
  f2.SetParameters(par);

  TF1 f3("f3","(1-exp(-pow(x/[0],2)))*exp(-pow(x/[1],[2]))",0,1);
  par[0] = 0.05;
  par[1] = 0.45;
  par[2] = 3.07;
  f3.SetParameters(par);
  
  
  TH1F *hrc = new TH1F("hrc","hrc",1500,0,1500);
  TH1F *hspe = new TH1F("hspe","hspe",1500,0,1500);
  for (int j=0;j!=32;j++){
    //std::cout << gain[j] << std::endl;
    hrc->Reset();
    hspe->Reset();
    for (int i=0;i!=1500;i++){
      double x = hrc->GetBinCenter(i+1);
      hspe->SetBinContent(i+1,f1.Eval(x)*gain[j]);
      double content = -1./rc_tau[j] * exp(-x/rc_tau[j]);
      if (i==0) content += 1;
      hrc->SetBinContent(i+1,content);
      // hraw[j]->SetBinContent(i+1,content);
    }
    // deconvolution ...
    TH1 *hm = hraw[j]->FFT(0,"MAG");
    TH1 *hp = hraw[j]->FFT(0,"PH");

    TH1 *hm_rc = hrc->FFT(0,"MAG");
    TH1 *hp_rc = hrc->FFT(0,"PH");

    TH1 *hm_spe = hspe->FFT(0,"MAG");
    TH1 *hp_spe = hspe->FFT(0,"PH");

    double value_re[1500], value_im[1500];
    double value_re1[1500], value_im1[1500];
    for (int i=0;i!=1500;i++){
      double freq;
      if (i<=750){
	freq = i/1500.*2;
      }else{
	freq = (1500-i)/1500.*2;
      }
      double rho = hm->GetBinContent(i+1)/ hm_rc->GetBinContent(i+1) / hm_spe->GetBinContent(i+1)  ;
      double phi = hp->GetBinContent(i+1) - hp_rc->GetBinContent(i+1) - hp_spe->GetBinContent(i+1);
      if (i==0) rho = 0;
      value_re[i] = rho * f3.Eval(freq)* cos(phi)/1500.;
      value_im[i] = rho * f3.Eval(freq)* sin(phi)/1500.;

      value_re1[i] = rho * cos(phi)/1500.* f2.Eval(freq);
      value_im1[i] = rho * sin(phi)/1500.* f2.Eval(freq);
    }

    // ROI finding
    int n = 1500;
    TVirtualFFT *ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    TH1 *fb = TH1::TransformHisto(ifft,0,"Re");
    // calcumate rms and mean ... 
    std::pair<double,double> results = cal_mean_rms(fb);
    //std::cout << results.first << " " << results.second << std::endl;
    TH1F *hflag = new TH1F("hflag","hflag",1500,0,1500);
    for (int i=0;i!=1500;i++){
      double content = fb->GetBinContent(i+1);
      if (fabs(content-results.first)>5*results.second){
	for (int j=-20;j!=20;j++){
	  hflag->SetBinContent(i+j+1,1);
	}
      }
    }
    delete fb;

    // solve for baseline 
    ifft->SetPointsComplex(value_re1,value_im1);
    ifft->Transform();
    fb = TH1::TransformHisto(ifft,0,"Re");
    double A11 = 0, A12 = 0, A21=0, A22=0;
    double B1 = 0, B2 = 0;
    double a=0, b=0;
    for (int i=0;i!=1500;i++){
      if (hflag->GetBinContent(i+1)==0){
	B2 += fb->GetBinContent(i+1);
	B1 += fb->GetBinContent(i+1) * fb->GetBinCenter(i+1);
	A11 += pow(fb->GetBinCenter(i+1),2);
	A12 += fb->GetBinCenter(i+1);
	A21 += fb->GetBinCenter(i+1);
	A22 += 1;
      }
    }
    
    if (A22>0){
      a = (B1*A22-B2*A12)/(A11*A22-A21*A12);
      b = (B1*A21-B2*A11)/(A22*A11-A12*A21);
    }
    // std::cout << a << " " << b << std::endl;
    for (int i=0;i!=1500;i++){
      fb->SetBinContent(i+1,fb->GetBinContent(i+1) - a * fb->GetBinCenter(i+1) -b);
    }
    results = cal_mean_rms(fb);
    for (int i=0;i!=1500;i++){
      fb->SetBinContent(i+1,fb->GetBinContent(i+1)-results.first+0.01);
    }

    // prepare L1 fit ... 
    TH1F *hrebin = new TH1F("hrebin","hrebin",250,0,250);
    
    for (int i=0;i!=250;i++){
      hrebin->SetBinContent(i+1,
			       fb->GetBinContent(6*i+1) +
			       fb->GetBinContent(6*i+2) +
			       fb->GetBinContent(6*i+3) +
			       fb->GetBinContent(6*i+4) +
			       fb->GetBinContent(6*i+5) +
			       fb->GetBinContent(6*i+6) );
    }
    
    // work on the L1 ... 
    std::vector<double> vals_y;
    std::vector<double> vals_x;
    std::vector<int> vals_bin;
    
    for (int i=0;i!=250;i++){
      double content = hrebin->GetBinContent(i+1);
      if (content>0.3){
	vals_y.push_back(content);
	vals_x.push_back(hrebin->GetBinCenter(i+1));
	vals_bin.push_back(i);
      }
      // if (content <0) content =0;
      //W(i) = content;
    }
    int nbin_fit = vals_x.size();
    VectorXd W = VectorXd::Zero(nbin_fit);
    MatrixXd G = MatrixXd::Zero(nbin_fit,nbin_fit);
    
    for (int i=0;i!=nbin_fit;i++){
      W(i) = vals_y.at(i) / sqrt(vals_y.at(i));
      double t1 = vals_x.at(i); // measured time
      for (int k=0;k!=nbin_fit;k++){
	double t2 = vals_x.at(k); // real time 
	if (t1>t2) {
	  G(i,k) = (0.75 * (exp(-((t1-t2)*6*15.625/1000.-3*15.625/1000.)/1.5)-exp(-((t1-t2)*6*15.625/1000.+3*15.625/1000.)/1.5))) / sqrt(vals_y.at(i));
	}else if (t1==t2){
	  G(i,k) = (0.25 + 0.75 *(1-exp(-3*15.625/1000./1.5))) / sqrt(vals_y.at(i));
	}else{
	  continue;
	}
	//std::cout << i << " " << k << " " << G(i,k) << std::endl;
      }
    }
    
    double lambda = 5;//1/2.;
    WireCell::LassoModel m2(lambda, 100000, 0.05);
    m2.SetData(G, W);
    m2.Fit();
    VectorXd beta = m2.Getbeta();
    
    // double sum1 = 0, sum2 = 0;
    // for (int i=0;i!=nbin_fit;i++){
    //   sum1 += beta(i);
    // }
    // for (int i=0;i!=250;i++){
    //   sum2 += hrebin->GetBinContent(i+1);
    // }
    // std::cout << j << " " << sum1 << " " << sum2 << std::endl;
    
    
    for (int i=0;i!=250;i++){
      // 
      hdecon[j]->SetBinContent(i+1,hrebin->GetBinContent(i+1));
    }
    for (int i=0;i!=nbin_fit;i++){
      hl1[j]->SetBinContent(vals_bin.at(i)+1,beta(i));
    }

    delete hrebin;
    delete hflag;
    
    // rebin ...

    delete fb;
    delete ifft;
    
    delete hm;
    delete hp;
    delete hm_rc;
    delete hp_rc;
    delete hm_spe;
    delete hp_spe;
  }
  
  // Now need to define the flash, 8 us, 85 bins out of 250 bins
  // PE, multiplicity threshold  
  // L1 PE, multiplicity ...
  //TH1F *h_totPE = new TH1F("h_totPE","h_totPE",250,0,250);
  for (int i=0;i!=32;i++){
    for (int j=0;j!=250;j++){
      double content = hdecon[i]->GetBinContent(j+1);
      if (content >0.2)
	h_totPE->SetBinContent(j+1,h_totPE->GetBinContent(j+1) + content);
      if (content > 1.5) // ~2 PE threshold ... 
	h_mult->SetBinContent(j+1,h_mult->GetBinContent(j+1)+1);
      
      content = hl1[i]->GetBinContent(j+1);
      h_l1_totPE->SetBinContent(j+1,h_l1_totPE->GetBinContent(j+1)+content);
      if (content > 1) // 1 PE threshold
	h_l1_mult->SetBinContent(j+1,h_l1_mult->GetBinContent(j+1)+1);
    }
  }

  // Now working on the flashes ... 
  // [-2,78)

  std::vector<int> flash_time;
  std::vector<double> flash_pe;
  
  for (int i=0;i!=250;i++){
    double pe = h_totPE->GetBinContent(i+1);
    double mult = h_mult->GetBinContent(i+1);
    // careteria: multiplicity needs to be higher than 3, PE needs to be higher than 6
    //std::cout << pe << " " << mult << std::endl;
    if (pe >= 6 && mult >= 3){
      if (flash_time.size()==0){
	flash_time.push_back(i);
	flash_pe.push_back(pe);
      }else{
	if (i - flash_time.back() >= 78){
	  flash_time.push_back(i);
	  flash_pe.push_back(pe);
	  // start one, and open a window of 8 us, if anything above it, not the next 2 bin
	  // if find one is bigger than this, save a flash ... start a new flash?
	}else if (i-flash_time.back() > 4 && pe > flash_pe.back()){
	  flash_time.push_back(i);
	  flash_pe.push_back(pe);
	}
      }
    }
  }
  
  //  std::cout << flash_time.size() << " " << flash_pe.size() << std::endl;
  //  for a flash, examine the L1 one to decide if add in more time ...?
  for (size_t i=0; i!=flash_time.size(); i++){
    //std::cout << flash_time.at(i) << " " << flash_pe.at(i) << std::endl;
    int start_bin = flash_time.at(i)-2;
    if (start_bin <0) start_bin = 0;

    int end_bin = start_bin + 78; // default
    if (end_bin > 250) end_bin = 250;
    if (i+1<flash_time.size()){
      if (end_bin > flash_time.at(i+1)-2)
	end_bin = flash_time.at(i+1)-2;
    }
    //  std::cout << start_bin << " " << end_bin << std::endl;
    //check with the next bin content ...
    Opflash *flash = new Opflash(hdecon, beam_dt[0], start_bin, end_bin);
    flash->Add_l1info(h_l1_totPE, h_l1_mult, beam_dt[0], start_bin, end_bin);
    //std::cout << flash->get_time() << " " <<flash->get_total_PE() << " " << flash->get_num_fired() << std::endl;
    beam_flashes.push_back(flash);
  }

  for (size_t i=0; i!=cosmic_flashes.size();i++){
    Opflash *cflash = cosmic_flashes.at(i);
    bool save = true;
    for (size_t j=0; j!=beam_flashes.size();j++){
      Opflash *bflash = beam_flashes.at(j);
      if (cflash->get_time() >= bflash->get_low_time() &&
	  cflash->get_time() <= bflash->get_high_time()){
	save = false;
	break;
      }
    }
    if (save)
      flashes.push_back(cflash);
  }
  for (size_t j=0; j!=beam_flashes.size();j++){
    Opflash *bflash = beam_flashes.at(j);
    flashes.push_back(bflash);
  }
  // std::cout << cosmic_flashes.size() << " " << beam_flashes.size() << " " << flashes.size() << std::endl;
  
  
  delete hrc;
  delete hspe;
}

std::pair<double,double> WireCell2dToy::uBooNE_light_reco::cal_mean_rms(TH1 *hist, int nbin){
  TH1F *h4 = new TH1F("h4","h4",2000,-10,10);
  double mean, rms;
  for (int i=0;i!=nbin;i++){
    double content = hist->GetBinContent(i+1);
    if (fabs(content)<10)
      h4->Fill(content);
  }
  mean = h4->GetBinCenter(h4->GetMaximumBin()+1);

  double xq = 0.5;
  double par[3];
  h4->GetQuantiles(1,&par[1],&xq);
  xq = 0.5 - 0.34;
  h4->GetQuantiles(1,&par[0],&xq);
  xq = 0.5 + 0.34;
  h4->GetQuantiles(1,&par[2],&xq);

  rms = sqrt((pow(par[0]-par[1],2)+pow(par[2]-par[1],2))/2.);
  
  delete h4;
  return std::make_pair(mean,rms);
}
