#include "WireCell2dToy/ClusterDisplay.h"
#include "TGraph2D.h"
#include "TVector3.h"

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
					     

void WireCell2dToy::ClusterDisplay::DrawCluster(SpaceCellSelection& mcells){
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
  std::cout << mcells.size() << std::endl;
  g1->Draw("p0");
 
}


void WireCell2dToy::ClusterDisplay::DrawCluster(MergeSpaceCellSelection& mcells){
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
  std::cout << n << std::endl;
  g1->Draw("p0");
  
}
