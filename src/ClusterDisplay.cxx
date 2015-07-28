#include "WireCell2dToy/ClusterDisplay.h"
#include "TGraph2D.h"
#include "TVector3.h"
#include "TRandom.h"

using namespace WireCell;

WireCell2dToy::ClusterDisplay::ClusterDisplay(TPad& pad)
  :pad(pad)
{
}

WireCell2dToy::ClusterDisplay::~ClusterDisplay(){
}

void WireCell2dToy::ClusterDisplay::DrawHough(SpaceCellSelection& cells, Point& p, double dis_near, double dis_far){
  TH2F *h1 = new TH2F("h1","h1",180,0,3.1415926,360,-3.1415926,3.1415926);
  Double_t x,y,z;
  Double_t x0 = p.x;
  Double_t y0 = p.y;
  Double_t z0 = p.z;

  for (int i=0;i!=cells.size();i++){
    SpaceCell *cell = cells.at(i);
    x = cell->x();
    y = cell->y();
    z = cell->z();
    
    TVector3 vec(x-x0,y-y0,z-z0);
    if (vec.Mag() < dis_far && vec.Mag() >= dis_near){
      h1->Fill(vec.Theta(),vec.Phi());
    }
    
  }
  
  h1->Draw("COLZ");
}

void WireCell2dToy::ClusterDisplay::DrawCrawler(WireCell2dToy::ToyCrawler& toycrawler, TString option){
  
  std::cout << "Draw Crawler " << " " << toycrawler.Get_allCT().size() << std::endl;

  int color[7]={3,4,5,6,7,8,9};
  //int style[7]={24,25,32,30,27,28,31};
  int style[7]={20,21,23,29,33,34,22};
  int num = 0;
  
  for (int i=0;i!=toycrawler.Get_allCT().size();i++){
    ClusterTrack* clustertrack = toycrawler.Get_allCT().at(i);
    Double_t x, y, z;
    TGraph2D *g1 = new TGraph2D();
    int n= 0;
    for (int j=0;j!=clustertrack->Get_allmcells().size();j++){
      MergeSpaceCell* mcell = clustertrack->Get_allmcells().at(j);
      Point center = mcell->Get_Center();
      x = center.x/units::cm;// + gRandom->Uniform(-0.15,0.15);
      y = center.y/units::cm;// + gRandom->Uniform(-0.15,0.15);
      z = center.z/units::cm;// + gRandom->Uniform(-0.15,0.15);
      
      g1->SetPoint(n,x,y,z);
      n++;
    }
    std::cout << i << " " << n << std::endl;

    //std::cout << i << " " << n << " " << x << " " << y << " " << z << std::endl;
    
    
    g1->Draw(option);
    if (num == 7 ) num = 0;
    g1->SetMarkerColor(color[num]);
    g1->SetMarkerStyle(style[num]);
    num++;
  }

}
					     

void WireCell2dToy::ClusterDisplay::DrawCluster(SpaceCellSelection& mcells, TString option){
  Double_t x, y, z;
  TGraph2D *g1 = new TGraph2D();
  
  for (int i=0;i!=mcells.size();i++){
    
    SpaceCell *cell = mcells.at(i);
      x = cell->x()/units::cm;
      y = cell->y()/units::cm;
      z = cell->z()/units::cm;
      g1->SetPoint(i,x,y,z);

    // if (i==0) 
      //std::cout << x << " " << y << " " << z << std::endl;
   
  }
  // std::cout << mcells.size() << std::endl;
  g1->Draw(option);
 
}


void WireCell2dToy::ClusterDisplay::DrawCluster(MergeSpaceCellSelection& mcells,TString option){
  Double_t x, y, z;
  TGraph2D *g1 = new TGraph2D();
  
  int n=0;
  for (int i=0;i!=mcells.size();i++){
    
    MergeSpaceCell *mcell = mcells.at(i);
    
    for (int j=0;j!=mcell->Get_all_spacecell().size();j++){
      SpaceCell *cell = mcell->Get_all_spacecell().at(j);
    
      x = cell->x()/units::cm;
      y = cell->y()/units::cm;
      z = cell->z()/units::cm;
      g1->SetPoint(n,x,y,z);
      n++;
    }
    // if (i==0) 
    //std::cout << x << " " << y << " " << z << std::endl;
    
  }
  // std::cout << n << std::endl;
  g1->Draw(option);

  //test 
  // g1->GetXaxis()->SetRangeUser(95,105);
  // g1->GetYaxis()->SetRangeUser(195,210);
  // g1->GetZaxis()->SetRangeUser(335,345);

  std::cout << mcells.size() << std::endl;

  TGraph2D *g2 = new TGraph2D();
  for (int i=0;i!=mcells.size();i++){
    MergeSpaceCell *mcell = mcells.at(i);
    Point center = mcell->Get_Center();
    x = center.x/units::cm;
    y = center.y/units::cm;
    z = center.z/units::cm;
    g2->SetPoint(i,x,y,z);
  }
  g2->Draw("Psame");
  g2->SetMarkerStyle(22);
  g2->SetMarkerColor(2);
  
}
