#include "WireCell2dToy/ToyLightReco.h"
#include "WireCellRess/LassoModel.h"
#include "WireCellRess/ElasticNetModel.h"
#include <Eigen/Dense>

#include "TH1S.h"
#include "TF1.h"
#include "TVirtualFFT.h"
#include <iostream>

using namespace WireCell;
using namespace Eigen;


WireCell2dToy::ToyLightReco::ToyLightReco(const char* root_file, bool imagingoutput){
  file = new TFile(root_file);
  if(imagingoutput){
    T = (TTree*)file->Get("Trun");
  }
  else{
    T = (TTree*)file->Get("/Event/Sim");
  }

  cosmic_hg_wf = new TClonesArray;
  cosmic_lg_wf = new TClonesArray;
  beam_hg_wf = new TClonesArray;
  beam_lg_wf = new TClonesArray;
  cosmic_hg_opch = new std::vector<short>;
  cosmic_lg_opch = new std::vector<short>;
  beam_hg_opch = new std::vector<short>;
  beam_lg_opch = new std::vector<short>;
  cosmic_hg_timestamp = new std::vector<double>;
  cosmic_lg_timestamp = new std::vector<double>;
  beam_hg_timestamp = new std::vector<double>;
  beam_lg_timestamp = new std::vector<double>;
  op_gain = new std::vector<float>;
  op_gainerror = new std::vector<float>;
  

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
  T->SetBranchAddress("triggerTime",&triggerTime);


  

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

  // fop_wf_beam = new TClonesArray("TH1S");
  // fop_femch_beam = new std::vector<short>;
  // fop_timestamp_beam = new std::vector<double>;
  // ctr_beam = 0;

  // fop_wf_cosmic = new TClonesArray("TH1S");
  // fop_femch_cosmic = new std::vector<short>;
  // fop_timestamp_cosmic = new std::vector<double>;
  // ctr_cosmic = 0;

  fop_wf = new TClonesArray("TH1S");
  fop_femch = new std::vector<short>;
  fop_timestamp = new std::vector<double>;
  ctr = 0;

  delete_status = true;
}

WireCell2dToy::ToyLightReco::~ToyLightReco(){
  if(delete_status){
  for (int i=0;i!=32;i++){
    delete hraw[i];
    hraw[i] = nullptr;
    delete hdecon[i];
    hdecon[i] = nullptr;
    delete hl1[i];
    hl1[i] = nullptr;
  }
  delete hraw;
  delete hdecon;
  delete hl1;

  delete h_totPE;
  delete h_mult;
  delete h_l1_mult;
  delete h_l1_totPE;

  // delete fop_wf_beam;
  // delete fop_femch_beam;
  // delete fop_timestamp_beam;

  // delete fop_wf_cosmic;
  // delete fop_femch_cosmic;
  // delete fop_timestamp_cosmic;

  delete fop_wf;
  delete fop_femch;
  delete fop_timestamp;
  
  delete T;
  delete file;

  hraw = nullptr;
  hdecon = nullptr;
  hl1 = nullptr;
  
  h_totPE = nullptr;
  h_mult = nullptr;
  h_l1_mult = nullptr;
  h_l1_totPE = nullptr;

  // fop_wf_beam = nullptr;
  // fop_femch_beam = nullptr;
  // fop_timestamp_beam = nullptr;

  // fop_wf_cosmic = nullptr;
  // fop_femch_cosmic = nullptr;
  // fop_timestamp_cosmic = nullptr;

  fop_wf = nullptr;
  fop_femch = nullptr;
  fop_timestamp = nullptr;

  T = nullptr;
  file = nullptr;
  }
  delete_status = false;
}


void WireCell2dToy::ToyLightReco::clear_flashes(){
  

  
  // clear flashes and actually delete them ... 
  beam_flashes.clear();
  cosmic_flashes.clear();
  flashes.clear(); 
}

void WireCell2dToy::ToyLightReco::load_event_raw(int eve_num){
  

  T->GetEntry(eve_num);

  WireCell2dToy::pmtMapSet beamHG = makePmtContainer(true, true, beam_hg_wf, beam_hg_opch, beam_hg_timestamp);
  WireCell2dToy::pmtMapSet beamLG = makePmtContainer(false, true, beam_lg_wf, beam_lg_opch, beam_lg_timestamp);
  WireCell2dToy::pmtMapSet cosmicHG = makePmtContainer(true, false, cosmic_hg_wf, cosmic_hg_opch, cosmic_hg_timestamp);
  WireCell2dToy::pmtMapSet cosmicLG = makePmtContainer(false, false, cosmic_lg_wf, cosmic_lg_opch, cosmic_lg_timestamp);

  WireCell2dToy::pmtMapPair beamPair = makeBeamPairs(beamHG, beamLG);
  WireCell2dToy::pmtMapSetPair cosmicPair = makeCosmicPairs(cosmicHG, cosmicLG);

  WireCell2dToy::pmtMap beamMerge = mergeBeam(beamPair);
  WireCell2dToy::pmtMapSet cosmicMerge = mergeCosmic(cosmicPair);

  fop_wf->Clear();
  fop_timestamp->clear();
  fop_femch->clear();
  
  dumpPmtVec(beamMerge, cosmicMerge);

  std::vector<COphitSelection> ophits_group;
  COphitSelection left_ophits;
  for (int i=32;i!=fop_femch->size();i++){
    COphit *op_hit = new COphit(fop_femch->at(i), (TH1S*)fop_wf->At(i), fop_timestamp->at(i) - triggerTime, op_gain->at(fop_femch->at(i)), op_gainerror->at(fop_femch->at(i)));
    if (op_hit->get_type()){ // what type  good baseline ???
      bool flag_used = false;
      if (ophits_group.size()==0){
	COphitSelection ophits;
	ophits.push_back(op_hit);
	ophits_group.push_back(ophits);
	flag_used = true;
      }else{
	for (size_t j=0; j!=ophits_group.size();j++){
	  for (size_t k=0; k!= ophits_group.at(j).size(); k++){
	    if (fabs(op_hit->get_time() - ophits_group.at(j).at(k)->get_time()) < 0.1 ){  // time unit??? 0.1 us?
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
  }
  
  for (size_t i=0;i!=left_ophits.size();i++){
    bool flag_used = false;
    for (size_t j=0; j!=ophits_group.size();j++){
      for (size_t k=0; k!= ophits_group.at(j).size(); k++){
	if (fabs(left_ophits.at(i)->get_time() - ophits_group.at(j).at(k)->get_time())<0.1){ // time unit??? 0.1 us?
	  ophits_group.at(j).push_back(left_ophits.at(i));
	  flag_used = true;
	  break;
	}
      }
      if (flag_used)
	break;
    }
  }

  for (size_t j=0; j!=ophits_group.size();j++){
    Opflash *flash = new Opflash(ophits_group.at(j));
    cosmic_flashes.push_back(flash);
  }
 
  for (int i=0;i!=32;i++){
    TH1S *hsignal = (TH1S*)fop_wf->At(i);
    for (int j=0;j!=1500;j++){
      hraw[i]->SetBinContent(j+1,hsignal->GetBinContent(j+1)-2050);
    }
    gain[i] = op_gain->at(i);
    beam_dt[i] = fop_timestamp->at(i) - triggerTime;
  }
  
  Process_beam_wfs();
  sort_flashes();

  // update map
  update_pmt_map();
  
}

void WireCell2dToy::ToyLightReco::update_pmt_map(){
  // std::cout << "Update map!" << std::endl;
  for (auto it=flashes.begin(); it!=flashes.end(); it++){
    Opflash *flash = *it;
    flash->swap_channels(); // hard coded for now ... 
  }
}


    
void WireCell2dToy::ToyLightReco::sort_flashes(){
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

void WireCell2dToy::ToyLightReco::Process_beam_wfs(){
  // correct the baseline ...
  TH1F h1("h1","h1",200,-100,100);
  for (int i=0;i!=32;i++){
    h1.Reset();
    for (int j=0;j!=20;j++){
      h1.Fill(hraw[i]->GetBinContent(j+1));
    }
    // double xq = 0.5;
    // double baseline;
    // h1.GetQuantiles(1,&baseline,&xq);
    double baseline = h1.GetMaximumBin()-100;
    if (fabs(baseline)>=8) baseline = 0;
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

  
  // prepare for global L1 fit ...
  // std::vector<double> global_vals_y;
  // std::vector<double> global_vals_x;
  // std::vector<int> global_vals_bin;
  // std::vector<int> global_vals_pmtid;
  // std::vector<double> global_sols;
  
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
      if (i<1500-4){
	fb->SetBinContent(i+1,fb->GetBinContent(i+1)-results.first+0.01);
      }else{
	fb->SetBinContent(i+1,0);
      }
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
    for (int i=0;i!=250;i++){
      hdecon[j]->SetBinContent(i+1,hrebin->GetBinContent(i+1));
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

	// global_vals_y.push_back(content);
	// global_vals_x.push_back(hrebin->GetBinCenter(i+1));
	// global_vals_bin.push_back(i);
	// global_vals_pmtid.push_back(j);
      }
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
      }
    }
    
    double lambda = 5;//1/2.;
    WireCell::LassoModel m2(lambda, 100000, 0.05);
    m2.SetData(G, W);
    m2.Fit();
    VectorXd beta = m2.Getbeta();
    for (int i=0;i!=nbin_fit;i++){
      //global_sols.push_back(beta(i));
      hl1[j]->SetBinContent(vals_bin.at(i)+1,beta(i));
    }
    

    delete hrebin;
    delete hflag;
    delete fb;
    delete ifft;
    delete hm;
    delete hp;
    delete hm_rc;
    delete hp_rc;
    delete hm_spe;
    delete hp_spe;
  }

  // std::set<int> global_fired_bin_set;
  // for (size_t i=0;i!=global_vals_bin.size();i++){
  //   if (global_fired_bin_set.find(global_vals_bin.at(i))==global_fired_bin_set.end()){
  //     global_fired_bin_set.insert(global_vals_bin.at(i));
  //   }
  // }
  // std::vector<int> global_fired_bin;
  // std::copy(global_fired_bin_set.begin(), global_fired_bin_set.end(), std::back_inserter(global_fired_bin));
  // std::vector<double> global_sum_sols;
  // global_sum_sols.resize(global_fired_bin.size(),0);
  // for (size_t i=0;i!=global_fired_bin.size();i++){
  //   for (size_t j=0;j!=global_vals_x.size();j++){
  //     if (global_vals_bin.at(j)==global_fired_bin.at(i))
  // 	global_sum_sols.at(i)+=global_sols.at(j);
  //   }
  // }
  // global_sols.insert(global_sols.end(),global_sum_sols.begin(),global_sum_sols.end());

  
  // //std::cout << global_vals_x.size() << std::endl;
  // VectorXd W_tot = VectorXd::Zero(global_vals_x.size());
  // MatrixXd G_tot = MatrixXd::Zero(global_vals_x.size(),global_vals_x.size() + global_fired_bin.size()); // unknown to measurement
  // //MatrixXd G_tot = MatrixXd::Zero(global_vals_x.size(),global_vals_x.size()); // unknown to measurement
  // MatrixXd F_tot = MatrixXd::Zero(global_fired_bin.size(),global_vals_x.size()+global_fired_bin.size()); //unknow to constraints 

  // for (size_t i=0;i!=global_vals_x.size();i++){
  //   W_tot(i) = global_vals_y.at(i) / sqrt(global_vals_y.at(i));
  //   double t1 = global_vals_x.at(i); // measured time
  //   for (size_t k=0;k!=global_vals_x.size();k++){
  //     if (global_vals_pmtid.at(i)!=global_vals_pmtid.at(k)) continue;
  //     double t2 = global_vals_x.at(k); // real time 
  //     if (t1>t2) {
  // 	G_tot(i,k) = (0.75 * (exp(-((t1-t2)*6*15.625/1000.-3*15.625/1000.)/1.5)-exp(-((t1-t2)*6*15.625/1000.+3*15.625/1000.)/1.5))) / sqrt(global_vals_y.at(i));
  //     }else if (t1==t2){
  // 	G_tot(i,k) = (0.25 + 0.75 *(1-exp(-3*15.625/1000./1.5))) / sqrt(global_vals_y.at(i));
  //     }else{
  // 	continue;
  //     }
  //   }
  // }
  // for (size_t i=0;i!=global_fired_bin.size();i++){
  //   for (size_t j=0;j!=global_vals_x.size();j++){
  //     if (global_vals_bin.at(j)==global_fired_bin.at(i)){
  // 	//std::cout << i << " " << j << " " << global_vals_pmtid.at(j) << " " << global_fired_bin.at(i) << std::endl;
  // 	F_tot(i,j)=1./0.001;
  //     }
  //   }
  //   F_tot(i,global_vals_x.size()+i)=-1./0.001;
  // }
  
  // // now the global fit ... 
  // MatrixXd GT_tot = G_tot.transpose();
  // MatrixXd FT_tot = F_tot.transpose();
  // VectorXd Wp = GT_tot * W_tot;
  // MatrixXd Gp = GT_tot * G_tot + FT_tot * F_tot;

  // double lambda = 0.1;//1/2.;
  // WireCell::LassoModel m3(lambda, 100000, 0.05);
  // m3.SetData(Gp, Wp);
  // for (size_t i=0;i!=global_vals_x.size();i++){
  //   m3.SetLambdaWeight(i,0.0);
  // }
  // m3.Set_init_values(global_sols);
  // m3.Fit();
  // VectorXd beta = m3.Getbeta();
  // for (size_t i=0;i!=global_vals_x.size();i++){
  //   //std::cout << global_vals_pmtid.at(i) << " " << global_vals_bin.at(i) << " " << beta(i) << std::endl;
  //   hl1[global_vals_pmtid.at(i)]->SetBinContent(global_vals_bin.at(i)+1,beta(i));
  // }
  // // finish the global fit ...
  
  
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

  TH1F *prev_hpe = new TH1F("prev_hpe","prev_hpe",32,0,32);
  TH1F *curr_hpe = new TH1F("curr_hpe","curr_hpe",32,0,32);
  
  for (int i=0;i!=250;i++){
    double pe = h_totPE->GetBinContent(i+1);
    double mult = h_mult->GetBinContent(i+1);
    // careteria: multiplicity needs to be higher than 3, PE needs to be higher than 6
    //std::cout << pe << " " << mult << std::endl;
    if (pe >= 6 && mult >= 3){
      bool flag_save = false;
      if (flash_time.size()==0){
	flag_save = true;
	for (int j=0;j!=32;j++){
	  prev_hpe->SetBinContent(j+1,hdecon[j]->GetBinContent(i+1));
	}
      }else{
	for (int j=0;j!=32;j++){
	  curr_hpe->SetBinContent(j+1,hdecon[j]->GetBinContent(i+1));
	}
	if (i - flash_time.back() >= 78){
	  flag_save = true;
	  // start one, and open a window of 8 us, if anything above it, not the next 2 bin
	  // if find one is bigger than this, save a flash ... start a new flash?
	}else if (i-flash_time.back() > 4 && pe > flash_pe.back()){
	  if (i-flash_time.back()>15){
	    flag_save = true;
	  }else{
	    if (curr_hpe->KolmogorovTest(prev_hpe,"M")>0.1){
	      flag_save = true;
	    }
	  }
	}
	for (int j=0;j!=32;j++){
	  prev_hpe->SetBinContent(j+1,hdecon[j]->GetBinContent(i+1));
	}
      }

      if (flag_save){
	flash_time.push_back(i);
	flash_pe.push_back(pe);
      }else{
	if (i - flash_time.back()<=6 && pe > flash_pe.back())
	  flash_pe.back()=pe;
      }
    }
  }

  delete prev_hpe;
  delete curr_hpe;
  
  
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

  Opflash *prev_cflash = 0;
  
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
    if (save){
      flashes.push_back(cflash);
      if (prev_cflash==0){
	prev_cflash = cflash;
      }else{
	if (cflash->get_time() < 0 && cflash->get_time() > prev_cflash->get_time())
	  prev_cflash = cflash;
      }
    }
  }
  for (size_t j=0; j!=beam_flashes.size();j++){
    Opflash *bflash = beam_flashes.at(j);
    if (prev_cflash!=0){
      if (bflash->get_time() - prev_cflash->get_time() < 2.4 && // veto for 3 us
	  bflash->get_total_PE() < 0.7 * prev_cflash->get_total_PE())
	continue;
      // std::cout << bflash->get_time() << " " << prev_cflash->get_time() << " " << bflash->get_total_PE() << " " << prev_cflash->get_total_PE() << std::endl;
    }
    flashes.push_back(bflash);
  }
  // std::cout << cosmic_flashes.size() << " " << beam_flashes.size() << " " << flashes.size() << std::endl;
  
  
  delete hrc;
  delete hspe;
}

WireCell2dToy::pmtMapSet WireCell2dToy::ToyLightReco::makePmtContainer(bool high, bool beam, TClonesArray *wf, std::vector<short> *chan, std::vector<double> *timestamp){
  short discSize=1500, baseline=2050, saturationThreshold=4080, refSize=1000;
  WireCell2dToy::pmtMapSet result;
  int num = wf->GetEntriesFast();
  if(beam == false){ discSize = 40; }
  std::array<float,32> scalePMT; scalePMT.fill(10.0);

  if(high == true){
    for(int i=0; i<num; i++){
      WireCell2dToy::pmtDisc disc;
      if(beam == true){
	if(chan->at(i)>=32) continue;
	disc.channel = chan->at(i);
      }
      else if(beam == false){
	disc.channel = chan->at(i)-200;
      }
      TH1S *h = (TH1S*)wf->At(i);
      disc.wfm.resize(discSize);
      int sat = 0;
      for(int j=0; j<discSize; j++){
	short content = h->GetBinContent(j+1);
	disc.wfm.at(j) = content;
	if(content>saturationThreshold){ sat++; }
      }
      if(sat >= 3){ disc.saturated = true; }
      else{ disc.saturated = false; }
      disc.timestamp = timestamp->at(i);
      disc.isolated = true;
      disc.highGain = true;
      result[disc.channel].insert(disc);
      h->Delete();
    }
  }
  else if(high == false){
    for(int i=0; i<num; i++){
      WireCell2dToy::pmtDisc disc;
      if(chan->at(i)>=132) continue;
      TH1S *h = (TH1S*)wf->At(i);
      if( (beam == true && h->GetNbinsX()<refSize) ||
	  (beam == false && h->GetNbinsX()>refSize) ){
	h->Delete();
	continue;
      }
      disc.channel = chan->at(i)-100;
      disc.saturated = false;
      disc.timestamp = timestamp->at(i);
      disc.isolated = true;
      disc.highGain = false;
      disc.wfm.resize(discSize);
      if(beam == true){
	baseline = 2050;
	double temp_baseline = findBaselineLg(h);
	//std::cout << temp_baseline << " " << baseline << std::endl;
	if (fabs(temp_baseline-baseline)<=8)
	  baseline = temp_baseline;
      }
      if(beam == false){ baseline = h->GetBinContent(1); }
      for(int j=0; j<discSize; j++){
	// is 2050 a good approximation??? 
	//disc.wfm.at(j) = (h->GetBinContent(j+1)-baseline)*scalePMT[disc.channel]+baseline;
	//if (j==0) std::cout << baseline << std::endl;
	disc.wfm.at(j) = (h->GetBinContent(j+1)-baseline)*findScaling(disc.channel)+baseline;
      }
      //      std::cout <<  disc.channel << " A " << disc.wfm.at(0) << std::endl;
      result[disc.channel].insert(disc);
      h->Delete();
    }
  }
  return result;
}

WireCell2dToy::pmtMapPair WireCell2dToy::ToyLightReco::makeBeamPairs(WireCell2dToy::pmtMapSet &high, WireCell2dToy::pmtMapSet &low){
  WireCell2dToy::pmtMapPair result;

  for(auto h=high.begin(); h!=high.end(); h++){
    WireCell2dToy::pmtPair tempPair;
    tempPair.first = *h->second.begin();
    tempPair.second = *low[h->first].begin();
    result[h->first] = tempPair;
  }
  return result;
}

WireCell2dToy::pmtMapSetPair WireCell2dToy::ToyLightReco::makeCosmicPairs(WireCell2dToy::pmtMapSet &high, WireCell2dToy::pmtMapSet &low){
  WireCell2dToy::pmtMapSetPair result;
  int countIsoLg = 0, tickWindow=20;
  float tick = 0.015625;

  // can use the 1D kd tree to search??? 
  for(auto h=high.begin(); h!=high.end(); h++){
    WireCell2dToy::pmtSetPair tempSetPair;

    for(auto hs : h->second){
      WireCell2dToy::pmtPair tempPair;
      bool found = false;

      for(auto ls : low[h->first]){
	if(std::fabs(hs.timestamp-ls.timestamp) > tickWindow*tick) continue;
	found = true;
	hs.isolated = false;
	ls.isolated = false;
	tempPair.first = hs;
	tempPair.second = ls;
	low[h->first].erase(ls);
	break;
      }
      if(found == false){
	tempPair.first = hs;
	tempPair.second = hs;
      }
      tempSetPair.insert(tempPair);
    }

    for(auto ls : low[h->first]){
      WireCell2dToy::pmtPair tempPair;
      if(ls.isolated == false) continue;
      tempPair.first = ls;
      tempPair.second = ls;
      tempSetPair.insert(tempPair);
    }
    result[h->first] = tempSetPair;
  }
  return result;
}

WireCell2dToy::pmtMap WireCell2dToy::ToyLightReco::mergeBeam(WireCell2dToy::pmtMapPair &beam){
  WireCell2dToy::pmtMap result;
  for(auto b=beam.begin(); b!=beam.end(); b++){
    if(b->second.first.saturated == true){
      //  std::cout << "Xin: " << b->second.first.channel <<  " " << b->second.first.wfm.at(0) << " " << b->second.second.wfm.at(0) << std::endl;
      WireCell2dToy::saturationTick tickVec = findSaturationTick(b->second.first.wfm);
      b->second.first.wfm = replaceSaturatedBin(b->second.first.wfm,b->second.second.wfm,tickVec);
      //b->second.first.wfm = b->second.second.wfm;
    }
    result[b->first] = b->second.first;
  }
  return result;
}

WireCell2dToy::pmtMapSet WireCell2dToy::ToyLightReco::mergeCosmic(WireCell2dToy::pmtMapSetPair &cosmic){
  WireCell2dToy::pmtMapSet result;
  for(auto c=cosmic.begin(); c!=cosmic.end(); c++){
    for(auto p : c->second){

      // std::cout << p.first.isolated << " " << p.first.saturated << " " << p.first.highGain << " " << p.first.timestamp-triggerTime << std::endl;
      
      if(p.first.isolated == true){ // isolated cos. disc.
	if(p.first.highGain == true){ // isolated HG
	  result[c->first].insert(p.first);
	}
	else{ // isolated LG
	  result[c->first].insert(p.first);
	}
      }
      else{ // paired cos. disc.
	if(p.first.saturated == true){ // saturated HG
	  //WireCell2dToy::saturationTick tickVec = findSaturationTick(p.first.wfm);
	  //p.first.wfm = replaceSaturatedBin(p.first.wfm,p.second.wfm,tickVec);
	  p.first.wfm = p.second.wfm;
	  result[c->first].insert(p.first);
	}
	else{ // non-saturated HG
	  result[c->first].insert(p.first);
	}
      }
    }
  }
  return result;
}

WireCell2dToy::saturationTick WireCell2dToy::ToyLightReco::findSaturationTick(std::vector<short> &wfm){
  WireCell2dToy::saturationTick result;
  bool saturatedStatus = false;
  std::pair<short,short> tempPair;
  short saturationThreshold=4080;
  
  for(int i=0; i<(int)wfm.size(); i++){
    if(wfm.at(i)>saturationThreshold){
      if(saturatedStatus == false){ tempPair.first = i; }
      saturatedStatus = true;
    }
    if(wfm.at(i)<saturationThreshold && saturatedStatus == true){
      saturatedStatus = false;
      tempPair.second = i;
      result.push_back(tempPair);
    }
  }
  return result;
}

std::vector<short> WireCell2dToy::ToyLightReco::replaceSaturatedBin(std::vector<short> &high, std::vector<short> &low, saturationTick &st){

  for(int i=0; i<(int)st.size(); i++){
    for(short j=st.at(i).first; j<st.at(i).second; j++){
      high.at(j) = low.at(j);
    }
  }
  return high;
}

void WireCell2dToy::ToyLightReco::dumpPmtVec(WireCell2dToy::pmtMap &beam, WireCell2dToy::pmtMapSet &cosmic){
  int z=0;
  for(auto b=beam.begin(); b!=beam.end(); b++){
    fop_femch->push_back( b->second.channel );
    fop_timestamp->push_back( b->second.timestamp );
    int nbins = (int)b->second.wfm.size();
    TH1S *h = new((*fop_wf)[z]) TH1S("","",nbins,0,nbins);
    for(int i=0; i<nbins; i++){
      h->SetBinContent(i+1,b->second.wfm.at(i));
    }
    z++;
  }
  for(auto c=cosmic.begin(); c!=cosmic.end(); c++){
    for(auto cs : c->second){
      fop_femch->push_back( cs.channel );
      fop_timestamp->push_back( cs.timestamp );
      int nbins = (int)cs.wfm.size();
      TH1S *h = new((*fop_wf)[z]) TH1S("","",nbins,0,nbins);
      for(int i=0; i<nbins; i++){
	h->SetBinContent(i+1,cs.wfm.at(i));
      }
      z++;
    }
  }
  return;
}

double WireCell2dToy::ToyLightReco::findBaselineLg(TH1 *hist, int nbin){
  TH1F *h = new TH1F("h","",1000,1500-0.5,2500-0.5);
  double baseline=0;
  for(int i=0; i!=20; i++){
    double content = hist->GetBinContent(i+1);
    //    baseline += content;
    if(content>1500 && content<2500){ h->Fill(content); }
  }
  //  baseline /= 6.;
  baseline = h->GetBinCenter(h->GetMaximumBin()+1);
  delete h;
  return baseline;
}

std::pair<double,double> WireCell2dToy::ToyLightReco::cal_mean_rms(TH1 *hist, int nbin){
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

float WireCell2dToy::ToyLightReco::findScaling(int opdet){
  if(opdet == 0){ return 10.13; }
  if(opdet == 1){ return 10.20; }
  if(opdet == 2){ return 10.13; }
  if(opdet == 3){ return 10.05; }
  if(opdet == 4){ return 9.96; }
  if(opdet == 5){ return 9.95; }
  if(opdet == 6){ return 10.04; }
  if(opdet == 7){ return 9.58; }
  if(opdet == 8){ return 9.42; }
  if(opdet == 9){ return 9.81; }
  if(opdet == 10){ return 9.25; }
  if(opdet == 11){ return 9.61; }
  if(opdet == 12){ return 9.56; }
  if(opdet == 13){ return 9.35; }
  if(opdet == 14){ return 9.99; }
  if(opdet == 15){ return 9.66; }
  if(opdet == 16){ return 10.26; }
  if(opdet == 17){ return 9.82; }
  if(opdet == 18){ return 10.32; }
  if(opdet == 19){ return 10.08; }
  if(opdet == 20){ return 9.77; }
  if(opdet == 21){ return 9.64; }
  if(opdet == 22){ return 10.14; }
  if(opdet == 23){ return 9.74; }
  if(opdet == 24){ return 9.76; }
  if(opdet == 25){ return 10.10; }
  if(opdet == 26){ return 10.48; }
  if(opdet == 27){ return 9.81; }
  if(opdet == 28){ return 9.99; }
  if(opdet == 29){ return 9.79; }
  if(opdet == 30){ return 9.93; }
  if(opdet == 31){ return 10.01; }
  return 0;
}
