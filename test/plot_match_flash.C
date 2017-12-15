#include <set>

void plot_match_flash(int flash_obj = 0){
  TFile *file = new TFile("match_5326_17_900.root");
  TTree *T = (TTree*)file->Get("T_match");
  TTree *T1 = (TTree*)file->Get("T_cluster");
  
  Double_t pe_pred[32];
  Double_t pe_meas[32];
  Double_t pe_meas_err[32];
  Int_t ncluster=0;
  Int_t flash_id;
  Double_t strength;
  TGraph *g1 = new TGraph();
  TGraph *g2 = new TGraph();
  TGraph *g3 = new TGraph();
  TGraph *g4 = new TGraph();
  TGraph *gpos = new TGraph();
  Int_t nps = 0;
  T->SetBranchAddress("tpc_cluster_id",&ncluster);
  T->SetBranchAddress("flash_id",&flash_id);
  T->SetBranchAddress("strength",&strength);
  T->SetBranchAddress("pe_pred",pe_pred);
  T->SetBranchAddress("pe_meas",pe_meas);
  T->SetBranchAddress("pe_meas_err",pe_meas_err);

  Double_t x, y, z;
  Int_t cluster_id;
  T1->SetBranchAddress("x",&x);
  T1->SetBranchAddress("y",&y);
  T1->SetBranchAddress("z",&z);
  T1->SetBranchAddress("cluster_id",&cluster_id);
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  
  // Int_t det_ch_map[32];
  // det_ch_map[0] = 3; det_ch_map[8] = 11; det_ch_map[16] = 18; det_ch_map[24] = 23;
  // det_ch_map[1] = 5; det_ch_map[9] = 7; det_ch_map[17] = 15; det_ch_map[25] = 26;
  // det_ch_map[2] = 1; det_ch_map[10] = 12; det_ch_map[18] = 16; det_ch_map[26] = 29;
  // det_ch_map[3] = 6; det_ch_map[11] = 8; det_ch_map[19] = 21; det_ch_map[27] = 30;
  // det_ch_map[4] = 0; det_ch_map[12] = 10; det_ch_map[20] = 22; det_ch_map[28] = 25;
  // det_ch_map[5] = 2; det_ch_map[13] = 14; det_ch_map[21] = 19; det_ch_map[29] = 31; 
  // det_ch_map[6] = 4; det_ch_map[14] = 17; det_ch_map[22] = 24; det_ch_map[30] = 27;
  // det_ch_map[7] = 9; det_ch_map[15] = 13; det_ch_map[23] = 20; det_ch_map[31] = 28;


  Double_t pe_measure[32],pe_predict[32];
  for (int i=0;i!=32;i++){
    pe_measure[i] = 0;
    pe_predict[i] = 0;
  }
  
  std::set<int> clusters;
  for (int i=0;i!=T->GetEntries();i++){
    T->GetEntry(i);
    if (flash_obj==flash_id){
      for (int j=0;j!=32;j++){
	pe_measure[j] = pe_meas[j];
	pe_predict[j] += pe_pred[j];
      }
      std::cout << "flash id: " << flash_id << " tpc id: " << ncluster << std::endl;
      clusters.insert(ncluster);
    }
  }

  c1->cd(1);
  for (int j=0;j!=32;j++){
    g1->SetPoint(j,j+1,pe_measure[j]);
    g2->SetPoint(j,j+1,pe_predict[j]);
  }
  g1->Draw("AL*");
  g2->Draw("L*same");
  g2->SetMarkerColor(2);
  g1->SetMarkerStyle(20);
  g2->SetMarkerStyle(21);
  g1->GetXaxis()->SetTitle("OpCh + 1");
  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->AddEntry(g1,"Measurement","lp");
  le1->AddEntry(g2,"Prediction","lp");
  le1->Draw();

  c1->cd(2);
  for (int i=0;i!=T1->GetEntries();i++){
    T1->GetEntry(i);
    if (clusters.find(cluster_id)!=clusters.end()){
      gpos->SetPoint(nps,z,y);
      nps++;
    }
  }
  gpos->Draw("AP");
  
  //   T->Scan("tpc_cluster_id",Form("flash_id==%d",flash_id));
  
}
