void plot_charge(){
  TFile *file = new TFile("temp_l1sp.root");
  //TFile *file = new TFile("nsp_2D_display_3493_821_41075.root");
  TH2F* hv_decon = (TH2F*)file->Get("hv_decon");
  Int_t nrebin = 4;
  // Int_t wire_min = 3550-2400;
  // Double_t wmin_tmin = 1250;
  // Double_t wmin_tmax = 1500;
  
  // Int_t wire_max = 4020-2400;
  // Double_t wmax_tmin = 500;
  // Double_t wmax_tmax = 700;

   Int_t wire_min = 3850-2400;
  Double_t wmin_tmin = 0;
  Double_t wmin_tmax = 100;
  
  Int_t wire_max = 4190-2400;
  Double_t wmax_tmin = 500;
  Double_t wmax_tmax = 630;

  TGraph *g1 = new TGraph();
  
  for (Int_t i=wire_min; i!=wire_max; i++){
    Double_t sum = 0;
    Int_t tbin_min = (wmin_tmin + (i-wire_min*1.0)/(wire_max-wire_min) * (wmax_tmin - wmin_tmin))/nrebin;
    Int_t tbin_max = (wmin_tmax + (i-wire_min*1.0)/(wire_max-wire_min) * (wmax_tmax - wmin_tmax))/nrebin;

    

    for (Int_t j=tbin_min; j!=tbin_max; j++){
      sum += hv_decon->GetBinContent(i+1,j+1);
    }
    std::cout << i+2400 << " " << tbin_min << " " << tbin_max << " " << sum << std::endl;
    
    g1->SetPoint(i-wire_min,i+2400,sum);
    
  }
  g1->Draw("A*");
  
}
