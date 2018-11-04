#include <vector>

int pnpoly(std::vector<double>& vertx, std::vector<double>& verty, double testx, double testy)
{
  int i, j, c = 0;
  for (i = 0, j = int(vertx.size())-1; i < int(vertx.size()); j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
      c = !c;
  }
  return c;
}

void plot_boundary(){
  // TH2F *h1 = new TH2F("h1","h1",210,50,260,20,-116,-96);
  // h1->Draw();
  // TLine *l1 = new TLine(77,-113,253,-96);
  // l1->Draw();

  std::vector<double> vertx = {3,77,253,253,97,3};
  std::vector<double> verty = {-113,-113,-96,100,115,114};

  TH2F *h1 = new TH2F("h1","h1",256,0,256,240,-120,120);

  for (int i=0;i!=256;i++){
    double x = h1->GetXaxis()->GetBinCenter(i+1);
    for (int j=0;j!=240;j++){
      double y = h1->GetYaxis()->GetBinCenter(j+1);
      h1->SetBinContent(i+1,j+1,pnpoly(vertx,verty,x,y));
    }
  }
  h1->Draw("COLZ");
  
  //  double x = 166.627;
  // double y = -104.599;

  //std::cout << pnpoly(vertx,verty,x,y) << std::endl;
  //  TGraph *g1 = new TGraph();
  // for (int i=0;i!=240;i++){
  //   double y = -120+i;
  //   double c = pnpoly(vertx,verty,x,y);
  //   g1->SetPoint(i,y,c);
  // }
  // g1->Draw("A*");
  
  // Double_t x[7]={3,77,253,253,97,3,3};
  // Double_t y[7]={-113,-113,-96,100,115,114,-113};
  //TGraph *g1 = new TGraph(7,x,y);
  //g1->Draw("AL");
}


