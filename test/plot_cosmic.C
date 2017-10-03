void plot_cosmic(int run_no=11110, int subrun_no=135, int event_no=6780){
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
  
  //std::cout << T_data->GetEntries() << std::endl;

  TTree *T_flash = (TTree*)file->Get("T_flash");
  int type;
  double low_time, high_time;
  double time;
  double total_PE;
  double PE[32],PE_err[32];
  T_flash->SetBranchAddress("type",&type);
  T_flash->SetBranchAddress("low_time",&low_time);
  T_flash->SetBranchAddress("high_time",&high_time);
  T_flash->SetBranchAddress("time",&time);
  T_flash->SetBranchAddress("total_PE",&total_PE);
  T_flash->SetBranchAddress("PE",PE);
  T_flash->SetBranchAddress("PE_err",PE_err);

  //TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  // c1->Divide(2,1);
  TH2F *h1 = new TH2F("h1","h1",10000,-3200,4800,32,0,32);
  TH2F *h2 = new TH2F("h2","h2",10000,-3200,4800,32,0,32);
  T_data->GetEntry(0);
  for (size_t i=0; i!=op_femch->size();i++){
    h1->Fill(op_timestamp->at(i)-triggerTime,op_femch->at(i),100);
  }
  // c1->cd(1);
  h1->Draw("COLZ");

  for (int i=0;i!=T_flash->GetEntries();i++){
    T_flash->GetEntry(i);
    // TLine *l1 = new TLine(time,0,time,32);
    // l1->Draw("same");
    //std::cout << type << std::endl;
    // if (i==0){
    for (int j=0;j!=32;j++){
      if (PE[j] >0){
	h2->Fill(time,j,50);
      }
      // 	std::cout << time << " " << PE[j] << std::endl;
    //   }
    }
  }
  h2->GetZaxis()->SetRangeUser(0,100);
  //h2->Fill(0.,0.,100);
  // c1->cd(2);
  h2->Draw("COLZsame");

  // int nhit = 42;
  // std::cout << std::endl;
  
  // TH1S *h2 = (TH1S*)op_wf->At(nhit);
  // std::cout << op_femch->at(nhit) << " " <<  op_timestamp->at(nhit) - triggerTime << std::endl;
  // double sum = 0;
  // for (int i=0;i!=40;i++){
  //   sum += h2->GetBinContent(i+1)-h2->GetBinContent(1);
  // }
  // std::cout << (sum/op_gain->at(op_femch->at(nhit))-0.15)*2 << std::endl;
  
  
}
