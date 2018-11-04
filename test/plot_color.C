void plot_color(){
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);
  TFile *file = new TFile("match_5122_6_324.root");
  TTree *T = (TTree*)file->Get("T_proj");
  TH2F *hw = new TH2F("hw","hw",300,6300,6600,400,1500,1900);
  TH2F *hu = new TH2F("hu","hu",400,900,1300,400,1500,1900);
  TH2F *hv = new TH2F("hv","hv",600,3300,3900,400,1500,1900);
  T->Project("hu","time_slice:channel","charge*(cluster_id==19)");
  T->Project("hv","time_slice:channel","charge*(cluster_id==19)");
  T->Project("hw","time_slice:channel","charge*(cluster_id==19)");
  TCanvas *c1 = new TCanvas("c1","c1",600,1200);
  c1->Divide(1,3);
  c1->cd(1);
  hu->Draw("COLZ");
  c1->cd(2);
  hv->Draw("COLZ");
  c1->cd(3);
  hw->Draw("COLZ");
}
