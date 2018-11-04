void plot_test(){
  ifstream infile("temp");
  TGraph *g1 = new TGraph();
  Double_t x,y;
  TString temp;
  double temp1;
  for(int i=0;i!=1530;i++){
    infile >> temp >> y >> x >> temp1 >> temp1;
    g1->SetPoint(i,x,y);
  }
  g1->Draw("AL");
  // ifstream infile("2.dat");
  // Double_t t,u,v,w,c;
  // TString temp;
  // infile >> t;
  // TGraph *gu = new TGraph();
  // TGraph *gv = new TGraph();
  // TGraph *gw = new TGraph();

  // for (int i=0;i!=2;i++){
  //   infile >> temp >> u >> c;
  //   gu->SetPoint(i,t,u);
  // }
  // for (int i=0;i!=2;i++){
  //   infile >> temp >> v >> c;
  //   gv->SetPoint(i,t,v);
  // }
  // for (int i=0;i!=3;i++){
  //   infile >> temp >> w >> c;
  //   gw->SetPoint(i,t,w);
  // }
  
  // TGraph *gu1 = new TGraph();
  // TGraph *gv1 = new TGraph();
  // TGraph *gw1 = new TGraph();
  // for (int i=0;i!=10;i++){
  //   infile >> temp >> u >> t >> c;
  //   gu1->SetPoint(i,t,u);
  // }
  // for (int i=0;i!=19;i++){
  //   infile >> temp >> v >> t >> c;
  //   gv1->SetPoint(i,t,v);
  // }
  // for (int i=0;i!=21;i++){
  //   infile >> temp >> w >> t >> c;
  //   gw1->SetPoint(i,t,w);
  // }
  
  // TCanvas *c1 = new TCanvas("c1","c1",1200,400);
  // c1->Divide(3,1);
  // c1->cd(1);
  // // TH2F *h1 = new TH2F("h1","h1",100,1394,1401,100,1410,1430);
  // // h1->Draw();
  // gu1->Draw("A*");
  // gu1->SetMarkerColor(1);
  // gu->Draw("*same");
  // gu->SetMarkerColor(2);
  // gu1->GetYaxis()->SetRangeUser(1050,1100);
  // c1->cd(2);
  
  // gv1->Draw("A*");
  // gv1->SetMarkerColor(1);
  // gv->Draw("*same");
  // gv->SetMarkerColor(2);
  // gv1->GetYaxis()->SetRangeUser(4000,4050);
  // c1->cd(3);
  
  // gw1->Draw("A*");
  // gw1->SetMarkerColor(1);
  // gw->Draw("*same");
  // gw->SetMarkerColor(2);
  // gw1->GetYaxis()->SetRangeUser(6800,6830);
  
  // // gu1->GetXaxis()->SetRangeUser(1394,1401);
  // // gu1->GetYaxis()->SetRangeUser(1410,1425);
}
