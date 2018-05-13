void plot_dead(){
  TString temp;
  Double_t y, z;
  ifstream infile("1.log");
  TGraph *g1 = new TGraph();
  for (int i=0;i!=74249;i++){
    infile >> temp >> y >> z;
    g1->SetPoint(i,z,y);
  }
  g1->Draw("AP");
}
