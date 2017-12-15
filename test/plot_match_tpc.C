void plot_match_tpc(int tpc_obj = 0){
  TFile *file = new TFile("match_5326_17_900.root");
  TTree *T = (TTree*)file->Get("T_match");
  TTree *T1 = (TTree*)file->Get("T_cluster");
  //  TChain *T = new TChain("T_match","T_match");
  //TChain *T1 = new TChain("T_match","T_match");
  //T->AddFile("");
  // T->AddFile("match_9009_199_9952.root");
  // T->AddFile("match_5946_11_593.root");
  // T->AddFile("match_5946_11_568.root");
  // T->AddFile("match_5122_20_1034.root");
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
  T->SetBranchAddress("tpc_cluster_id",&ncluster);
  T->SetBranchAddress("flash_id",&flash_id);
  T->SetBranchAddress("strength",&strength);
  T->SetBranchAddress("pe_pred",pe_pred);
  T->SetBranchAddress("pe_meas",pe_meas);
  T->SetBranchAddress("pe_meas_err",pe_meas_err);

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  //TCanvas *c2 = new TCanvas("c2","c2",800,600);
  c1->Divide(2,1);
  
  Int_t det_ch_map[32];
  det_ch_map[0] = 3; det_ch_map[8] = 11; det_ch_map[16] = 18; det_ch_map[24] = 23;
  det_ch_map[1] = 5; det_ch_map[9] = 7; det_ch_map[17] = 15; det_ch_map[25] = 26;
  det_ch_map[2] = 1; det_ch_map[10] = 12; det_ch_map[18] = 16; det_ch_map[26] = 29;
  det_ch_map[3] = 6; det_ch_map[11] = 8; det_ch_map[19] = 21; det_ch_map[27] = 30;
  det_ch_map[4] = 0; det_ch_map[12] = 10; det_ch_map[20] = 22; det_ch_map[28] = 25;
  det_ch_map[5] = 2; det_ch_map[13] = 14; det_ch_map[21] = 19; det_ch_map[29] = 31; 
  det_ch_map[6] = 4; det_ch_map[14] = 17; det_ch_map[22] = 24; det_ch_map[30] = 27;
  det_ch_map[7] = 9; det_ch_map[15] = 13; det_ch_map[23] = 20; det_ch_map[31] = 28;
  
  for (int i=tpc_obj;i!=tpc_obj+1;i++){
    T->GetEntry(i);
    // flash_id
    
    c1->cd(1);
    for (int j=0;j!=32;j++){
      g1->SetPoint(j,j+1,pe_meas[j]);
      g2->SetPoint(j,j+1,pe_pred[j]);
      g3->SetPoint(det_ch_map[j],det_ch_map[j],pe_meas[j]);
      g4->SetPoint(det_ch_map[j],det_ch_map[j],pe_pred[j]);
      //g4->SetPoint(det_ch_map[j],det_ch_map[j],pe_pred[j]);
    }
    g1->Draw("AL*");
    g2->Draw("L*same");
    g2->SetMarkerColor(2);
    g1->SetMarkerStyle(20);
    g2->SetMarkerStyle(21);
    TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
    le1->AddEntry(g1,"Measurement","lp");
    le1->AddEntry(g2,"Prediction","lp");
    le1->Draw();
    g1->SetTitle(Form("Flash %d, TPC %d, Strength, %3.2f",flash_id,ncluster,strength));
    g1->GetXaxis()->SetTitle("OpCh+1");
    //c1->Update();
    c1->cd(2);
    T1->Draw("y:z",Form("cluster_id==%d",tpc_obj));
    // g3->Draw("AL*");
    // g3->SetMarkerColor(1);
    // g3->SetMarkerStyle(20);
    // g3->GetXaxis()->SetTitle("OpDet");
    // g4->Draw("L*same");
    // g4->SetMarkerColor(2);
    // g4->SetMarkerStyle(21);
    T->Scan("tpc_cluster_id",Form("flash_id==%d",flash_id));
  }
}
