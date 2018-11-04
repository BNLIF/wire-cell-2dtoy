void plot_charge1(){
  TFile *file = new TFile("temp_l1sp.root");
  TH2F* hv_raw = (TH2F*)file->Get("hv_raw");
  Int_t ch = 3998-2400;
  TH1F *h1 = new TH1F("h1","h1",300,500,800);
  for (Int_t i=500;i!=800;i++){
    h1->SetBinContent(i+1-500,hv_raw->GetBinContent(ch+1,i+1));
  }
  h1->Draw();
  std::cout << h1->Integral(100,116) << std::endl;

#include "./2dtoy/src/data_70_2D_11.txt"
  TGraph *gw = new TGraph(5000,w_2D_g_0_x,w_2D_g_0_y);
  int scaling =   14 * 1.2 * 4096/2000.;

  TH1F *h2 = new TH1F("h2","h2",40,-10,10);
  for (Int_t i=0;i!=40;i++){
    Double_t x = h2->GetBinCenter(i+1);
    h2->SetBinContent(i+1,gw->Eval(x)*scaling);
  }
  h2->Draw();

  std::cout << h2->GetSum() << std::endl;
  
//   //missing a factor of 2 ... 
//   std::cout << h2->GetSum() << " " << h1->Integral(100,116) / h2->GetSum() << std::endl;

  // int ch = 6474-4800;
  // TH2F *hw_raw = (TH2F*)file->Get("hw_raw");
  // TH2F *hw_decon = (TH2F*)file->Get("hw_decon");
  // TH1F *h1 = new TH1F("h1","h1",400,4600,5000);
  // TH1F *h3 = new TH1F("h3","h3",100,4600,5000);
  // for (Int_t i=0;i!=400;i++){
  //   h1->SetBinContent(i+1,hw_raw->GetBinContent(ch+1,4600+i+1));
  // }
  // h1->Draw();
  // for (Int_t i=0;i!=100;i++){
  //   h3->SetBinContent(i+1,hw_decon->GetBinContent(ch+1,1150+i+1));
  // }
  //h3->Draw();
  // std::cout << h1->Integral(200,300) << " " << h1->GetSum() << " " << h3->GetSum() << std::endl;
}
