#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"

#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char* argv[]){

  if(argc < 2){
    cerr << "USAGE: test-light-mismatch /path/to/match.root" << endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);

  const int nChan = 32;
  Double_t ks_m = 0.;
  Double_t chisq = 0.;
  Double_t ndf = 0.;
  
  TString fInName = argv[1];
  TFile *fIn = new TFile(fInName,"update");

  TTree *info = (TTree*)fIn->Get("Trun");
  unsigned int bit;
  int eventNo;
  int runNo;
  int subRunNo;
  info->SetBranchAddress("eventNo",&eventNo);
  info->SetBranchAddress("runNo",&runNo);
  info->SetBranchAddress("subRunNo",&subRunNo);
  info->SetBranchAddress("triggerBits",&bit);
  info->GetEntry(0);
  double low=0., high=0.;
  if(bit == 2048){ low = 3.0; high = 5.0; }
  if(bit == 512){ low = 3.45; high = 5.45; }
  
  TTree *flashT = (TTree*)fIn->Get("T_flash");
  int type;
  int f_flashId;
  double time;
  double totPE;
  flashT->SetBranchAddress("type",&type);
  flashT->SetBranchAddress("flash_id",&f_flashId);
  flashT->SetBranchAddress("time",&time);
  flashT->SetBranchAddress("total_PE",&totPE);
  for(int nf = 0; nf<flashT->GetEntries(); nf++){
    flashT->GetEntry(nf);
    if(time > low && time < high && type == 2) break;
  }

  TH1D *h_m = new TH1D("h_m","; Channel; PE",32,0,32);
  TH1D *h_p = new TH1D("h_p","",32,0,32);
  
  TTree *matchT = (TTree*)fIn->Get("T_match");
  int m_flashId;
  int m_clusterId;
  vector<int> v_clusterId;
  double pe_m[nChan];
  double pe_m_err[nChan];
  double pe_p[nChan];
  matchT->SetBranchAddress("flash_id",&m_flashId);
  matchT->SetBranchAddress("tpc_cluster_id",&m_clusterId);
  matchT->SetBranchAddress("pe_meas",pe_m);
  matchT->SetBranchAddress("pe_meas_err",pe_m_err);
  matchT->SetBranchAddress("pe_pred",pe_p);

  int beam_cluster_id = 10000;
  for(int nm = 0; nm<matchT->GetEntries(); nm++){
    matchT->GetEntry(nm);
    if(f_flashId == m_flashId){
      v_clusterId.push_back(m_clusterId);
      if(m_clusterId < beam_cluster_id){
	beam_cluster_id = m_clusterId;
      }
    }
  }

  TH1D *meas = new TH1D("meas","; Channel; PE",32,0,32);
  TH1D *pred = new TH1D("pred","",32,0,32);
  
  int keepId = 0;
  double ks_temp = 2.0;
  double ks_val = 2.0;
  for(int nm=0; nm<matchT->GetEntries(); nm++){
    matchT->GetEntry(nm);
    for(int i=0; i<(int)v_clusterId.size(); i++){
      if(f_flashId = m_flashId && v_clusterId.at(i) == m_clusterId){
	for(int chan=0; chan < nChan; chan++){
	  h_m->SetBinContent(chan+1, pe_m[chan]);
	  h_m->SetBinError(chan+1, pe_m_err[chan]);
	  h_p->SetBinContent(chan+1, pe_p[chan]);
	}
      }
      
      if(h_m->Integral()==0) continue;
      
      
      ks_temp = h_m->KolmogorovTest(h_p,"M");
      if(ks_temp < ks_val){
	ks_val = ks_temp;
	keepId = v_clusterId.at(i);
	meas = (TH1D*)h_m->Clone("meas");
	pred = (TH1D*)h_p->Clone("pred");
      }
    }
  }

  float m5 = 0., m10 = 0., m20 = 0., m40 = 0.;
  float p5 = 0., p10 = 0., p20 = 0., p40 = 0.;
  double predPE = 0., measPE = 0.;
  for(int i=0; i<meas->GetNbinsX(); i++){
    double p = pred->GetBinContent(i+1);
    if(p > 5){ p5++;
      if(p > 10){ p10++;
	if(p > 20){ p20++;
	  if(p > 40){ p40++; }
	}
      }
    }
    predPE += p;

    double m = meas->GetBinContent(i+1);
    if(m > 5){ m5++;
      if(m > 10){ m10++;
	if(m > 20){ m20++;
	  if(m > 40){ m40++; }
	}
      }
    }
    measPE += p;    
  }
  /*  
  for(int nm = 0; nm<matchT->GetEntries(); nm++){
    matchT->GetEntry(nm);
    if(f_flashId == m_flashId && beam_cluster_id == m_clusterId){
      for(int chan = 0; chan < nChan; chan++){
	h_m->SetBinContent(chan+1,pe_m[chan]);
	h_m->SetBinError(chan+1,pe_m_err[chan]);
	h_p->SetBinContent(chan+1,pe_p[chan]);

	chisq += pow(pe_p[chan]-pe_m[chan],2)/pow(pe_m_err[chan],2);
	predPE += pe_p[chan];
	
	if(pe_m[chan] > 5){ m5++;
	  if(pe_m[chan] > 10){ m10++;
	    if(pe_m[chan] > 20){ m20++;
	      if(pe_m[chan] > 40){ m40++;
	      }
	    }
	  }
	}
	if(pe_p[chan] > 5){ p5++;
	  if(pe_p[chan] > 10){ p10++;
	    if(pe_p[chan] > 20){ p20++;
	      if(pe_p[chan] > 40){ p40++;
	      }
	    }
	  }
	}
	
	if(pe_m[chan]==0 && pe_p[chan]==0){
	}
	else{
	  ndf++;
	}	  
      }
      break;
    }
  }
  */
  if(m5 == 0){ m5 = 1; }
  if(m10 == 0){ m10 = 1; }
  if(m20 == 0){ m20 = 1; }
  if(m40 == 0){ m40 = 1; }
  if(p5 == 0){ p5 = 1; }
  if(p10 == 0){ p10 = 1; }
  if(p20 == 0){ p20 = 1; }
  if(p40 == 0){ p40 = 1; }

  bool pass = true;
  
  ks_m = meas->KolmogorovTest(pred,"M");
  meas->SetLineColor(kBlack);
  pred->SetLineColor(kRed);

  double r_pe10 = m10/p10;
  double d_pe20 = (m20/totPE)/(p20/predPE);
  
  if(ks_m>0.5){ pass = false; }
  if(m10/p10>3){ pass = false; }
  if(totPE>1000 && (m20/totPE)/(p20/predPE)<0.2){ pass = false; }
  if(ks_m<0.1){ pass = true; }
  if(meas->Integral() == 0){ pass = false; }
  
  TString output;
  //output.Form("nuMuCandidate_%d_%d_%d",runNo,subRunNo,eventNo);
  output.Form("categorizedWrongLight_%d_%d_%d",runNo,subRunNo,eventNo);
  double cn = totPE/chisq/ndf;
  TString title;
  if(pass == true){ title.Form("PASS: KS_{#Deltadistance}: %f  PE_{total}: %f  R_{PE>10}: %f  D_{PE>20}: %f", ks_m, totPE, r_pe10, d_pe20); }
  if(pass == false){ title.Form("FAIL: KS_{#Deltadistance}: %f  PE_{total}: %f  R_{PE>10}: %f  D_{PE>20}: %f", ks_m, totPE, r_pe10, d_pe20); }
  TCanvas *c = new TCanvas("c","",500,500);
  meas->SetTitle(title);
  meas->SetStats(0);
  meas->Draw();
  pred->Draw("same");
  c->SaveAs(output+".png");

  std::cout<<"Run "<<runNo<<"  Subrun "<<subRunNo<<"  Event "<<eventNo<<"  total PE "<<totPE<<std::endl;
  std::cout<<"KS max diff: "<<ks_m<<std::endl;
  std::cout<<"N5: " <<m5/p5<<"  N10: "<<m10/p10<<"  N20: "<<m20/p20<<"  N40: "<<m40/p40<<std::endl;
  std::cout<<"L5: " <<(m5/totPE)/(p5/predPE)<<"  L10: "<<(m10/totPE)/(p10/predPE)<<"  L20: "<<(m20/totPE)/(p20/predPE)<<"  L40: "<<(m40/totPE)/(p40/predPE)<<std::endl;
  std::cout<<"---------------------------------------"<<std::endl;
  
  TFile *fout = new TFile("output_lightMismatch.root","RECREATE");
  /*  TTree *tout  = new TTree("tout","tout");
  tout->SetDirectory(fout);
  tout->Branch("ks_m",&ks_m,"ks_m/D");
  tout->Branch("chisq",&chisq,"chisq/D");
  tout->Branch("ndf",&ndf,"ndf/D");

  tout->Fill();
  */
  //h_m->SetDirectory(fout);
  //h_p->SetDirectory(fout);
  meas->SetDirectory(fout);
  pred->SetDirectory(fout);
  fout->Write();
  fout->Close();
  return 0;
}
