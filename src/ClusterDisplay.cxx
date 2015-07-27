#include "WireCell2dToy/ClusterDisplay.h"
#include "TGraph2D.h"

using namespace WireCell;

WireCell2dToy::ClusterDisplay::ClusterDisplay(TPad& pad)
  :pad(pad)
{
}

WireCell2dToy::ClusterDisplay::~ClusterDisplay(){
}

void WireCell2dToy::ClusterDisplay::DrawCluster(SpaceCellSelection& cells){
  Double_t x, y, z;
  TGraph2D *g1 = new TGraph2D();
  
  for (int i=0;i!=cells.size();i++){
    SpaceCell *cell = cells.at(i);
    x = cell->x()/units::cm;
    y = cell->y()/units::cm;
    z = cell->z()/units::cm;
    g1->SetPoint(i,x,y,z);
    
  }
  g1->Draw("p");
  
}
