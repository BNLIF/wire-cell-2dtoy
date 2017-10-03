void plot_beam(int run_no=11110, int subrun_no=135, int event_no=6780){
  TFile *file = new TFile(Form("flash_%d_%d_%d.root",run_no, subrun_no, event_no));
  
  TTree *T_data = (TTree*)file->Get("T_data");
  TClonesArray* op_wf = new TClonesArray;
  std::vector<int> *op_femch = new std::vector<int>;
  std::vector<double> *op_timestamp = new std::vector<double>;
  std::vector<double> *op_gain = new std::vector<double>;
  std::vector<double> *op_gainerror = new std::vector<double>;
  double triggerTime;
  
  T_data->SetBranchAddress("op_femch",&op_femch);
  T_data->SetBranchAddress("op_gain",&op_gain);
  T_data->SetBranchAddress("op_gainerror",&op_gainerror);
  T_data->SetBranchAddress("op_timestamp",&op_timestamp);
  T_data->SetBranchAddress("op_wf",&op_wf);
  T_data->SetBranchAddress("triggerTime",&triggerTime);

  T_data->GetEntry(0);

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->SetLogz(1);
  TH2F *hraw = (TH2F*)file->Get("hraw");
  TH2F *hdecon = (TH2F*)file->Get("hdecon");
  hraw->Draw("COLZ");

   TTree *T_flash = (TTree*)file->Get("T_flash");
  int type;
  double low_time, high_time;
  double time;
  double total_PE;
  double PE[32],PE_err[32];
  std::vector<int> *fired_channels = new std::vector<int>;
  std::vector<double> *l1_fired_time = new std::vector<double>;
  std::vector<double> *l1_fired_pe = new std::vector<double>;
  T_flash->SetBranchAddress("type",&type);
  T_flash->SetBranchAddress("low_time",&low_time);
  T_flash->SetBranchAddress("high_time",&high_time);
  T_flash->SetBranchAddress("time",&time);
  T_flash->SetBranchAddress("total_PE",&total_PE);
  T_flash->SetBranchAddress("PE",PE);
  T_flash->SetBranchAddress("PE_err",PE_err);
  T_flash->SetBranchAddress("fired_channels",&fired_channels);
  T_flash->SetBranchAddress("l1_fired_time",&l1_fired_time);
  T_flash->SetBranchAddress("l1_fired_pe",&l1_fired_pe);
  
  for (int i=0;i!=T_flash->GetEntries();i++){
    T_flash->GetEntry(i);
    if (type==2){
      int bin = (time - op_timestamp->at(0)+triggerTime)/(15.625/1000.);
      int low_bin = (low_time - op_timestamp->at(0)+triggerTime)/(15.625/1000.);
      int high_bin = (high_time - op_timestamp->at(0)+triggerTime)/(15.625/1000.);
      //std::cout << bin << std::endl;
      TLine *l1 = new TLine(bin,0,bin,32);
      l1->Draw("same");
      l1->SetLineColor(1);
      TLine *l2 = new TLine(low_bin,0,low_bin,32);
      l2->Draw("same");
      l2->SetLineColor(2);
      TLine *l3 = new TLine(high_bin,0,high_bin,32);
      l3->Draw("same");
      l3->SetLineColor(2);
      for (size_t j=0;j!=fired_channels->size();j++){
	TLine *l4 = new TLine(low_bin,fired_channels->at(j)+0.5,high_bin,fired_channels->at(j)+0.5);
	l4->Draw("same");
      }
      
      //std::cout << l1_fired_time->size() << std::endl;
      for (size_t j=0;j!=l1_fired_time->size();j++){
	int bin1 = (l1_fired_time->at(j) - op_timestamp->at(0)+triggerTime)/(15.625/1000.);
	TLine *l4 = new TLine(bin1,0,bin1,32);
	l4->Draw("same");
	l4->SetLineColor(6);
      }
    }
  }
  

  TH1F *h1 = new TH1F("h1","h1",1500,0,1500);
  TH1F *h2 = new TH1F("h2","h2",250,0,250);
  for (int i=0;i!=1500;i++){
    h1->SetBinContent(i+1,hraw->GetBinContent(i+1,4));
  }
  for (int i=0;i!=250;i++){
    h2->SetBinContent(i+1,hdecon->GetBinContent(i+1,4));
  }
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  c2->cd(1);
  //h1->Rebin(6);
  h1->Draw();
  c2->cd(2);
  h2->Draw();
  
 

  //std::cout << h1->Integral(1440,1480)/op_gain->at(6)+80/op_gain->at(3) << " " << h2->Integral(239,246) << std::endl;
}
