#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCellData/Units.h"
#include "WireCellData/Point.h"
#include "TGraph.h"
#include "TLine.h"
#include <iostream>

using namespace WireCell2dToy;

ToyEventDisplay::ToyEventDisplay(TPad& pad, const WireCell::GeomDataSource& gds)
    : pad(pad)
    , gds(gds)
    , h1(0)
    , g1(0)
    , g2(0)
{
}

ToyEventDisplay::~ToyEventDisplay()
{
    this->clear();
}
void ToyEventDisplay::clear()
{
    if (h1) {
      delete h1;
    }
    if (g1) {
      delete g1;
    }
    if (g2) { 
      delete g2;
    }
}

int ToyEventDisplay::init(float x_min, float x_max, float y_min, float y_max)
{
  this->clear();

  h1 = new TH2F("h1","h1",1000,x_min,x_max,1000,y_min,y_max);
  h1->SetTitle("Wires and True Hits");
  h1->GetYaxis()->SetNdivisions(506);
  h1->GetXaxis()->SetNdivisions(506);
  h1->SetXTitle("Z (m)");
  h1->SetYTitle("Y (m)");
  h1->SetBinContent(0,0,1);
  h1->GetZaxis()->SetRangeUser(charge_min,charge_max);

  return 0;
}

int ToyEventDisplay::draw_mc(int flag, const WireCell::PointValueVector& mctruth, TString option)
{
  pad.cd();

  if (flag==1){
    h1->Draw(option);
    
  }else if (flag==2){
    h1->Reset();
    for (int i=0;i!=mctruth.size();i++){
      h1->Fill(mctruth[i].first.z/units::m,mctruth[i].first.y/units::m-0.02,int(mctruth[i].second*10)/10.);
    }
    h1->Draw(option);
  }else if (flag==3){
      //if (g1) delete g1;
    g1 = new TGraph();
    for (int i=0;i!=mctruth.size();i++){
      g1->SetPoint(i,mctruth[i].first.z/units::m,mctruth[i].first.y/units::m);
    }
    g1->SetMarkerColor(2);
    g1->SetMarkerSize(0.6);
    g1->Draw(option);
    g1->SetMarkerStyle(20);
  }
  
  return 0;
}

int ToyEventDisplay::draw_slice(const WireCell::Slice& slice, TString option)
{
  pad.cd();

  WireCell::Channel::Group group = slice.group();
  //std::cout << group.size() << std::endl;
  
  for (int i=0;i!=group.size();i++){
    //std::cout << group.at(i).first << std::endl;
    const WireCell::GeomWire *wire = gds.by_channel_segment(group.at(i).first,0);
    // std::cout << wire->channel() << std::endl;
    // if ( wire->channel() ==1429 || wire->channel() ==4461){
    
      float pitch = gds.pitch(wire->plane());
      float angle = gds.angle(wire->plane());
      
      //std::cout << pitch << " " << angle << std::endl;
      
      // std::cout << wire->point1().y << " " << wire->point1().z << std::endl;
      TLine *l1 = new TLine(wire->point1().z/units::m - pitch/2./units::m/std::cos(angle/units::radian) ,wire->point1().y/units::m,
			    wire->point2().z/units::m - pitch/2./units::m/std::cos(angle/units::radian) ,wire->point2().y/units::m);
      l1->Draw(option);
      
      TLine *l2 = new TLine(wire->point1().z/units::m + pitch/2./units::m/std::cos(angle/units::radian) ,wire->point1().y/units::m,
			    wire->point2().z/units::m + pitch/2./units::m/std::cos(angle/units::radian) ,wire->point2().y/units::m);
      l2->Draw(option);
      
      TLine *l3 = new TLine(wire->point1().z/units::m  ,wire->point1().y/units::m,
			    wire->point2().z/units::m  ,wire->point2().y/units::m);
      l3->SetLineColor(2);
      l3->Draw(option);
    }
  // }
  return 0;
}


int ToyEventDisplay::draw_cells(const WireCell::GeomCellSelection& cellall, TString option)
{
  pad.cd();

  //if (g2) delete g2;
  g2 = new TGraph();
  for (int i=0;i!=cellall.size();i++){
    g2->SetPoint(i,cellall[i]->center().z/units::m,cellall[i]->center().y/units::m);
  }
  g2->SetMarkerColor(4);
  g2->SetMarkerSize(0.8);
  g2->Draw(option);
  g2->SetMarkerStyle(21);
  
  return 0;
}


int ToyEventDisplay::draw_mergecells(const WireCell::GeomCellSelection& cellall, TString option, int flag)
{
  pad.cd();
  
  int np = 0;
  g2 = new TGraph();
  for (int i=0;i!=cellall.size();i++){
    if (flag==0){
      g2->SetPoint(np,cellall[i]->center().z/units::m,cellall[i]->center().y/units::m);
      np++;
    }else if (flag==1){
      WireCell::MergeGeomCell *mcell = (WireCell::MergeGeomCell*)cellall[i];
      if (mcell->GetContainTruthCell()){
	g2->SetPoint(np,cellall[i]->center().z/units::m,cellall[i]->center().y/units::m);
	np++;
      }
    }
  }
  g2->SetMarkerColor(2);
  g2->SetMarkerSize(0.8);
  g2->Draw(option);
  g2->SetMarkerStyle(24);
  
  return 0;
}



int ToyEventDisplay::draw_truthcells(const WireCell::CellChargeMap& ccmap, TString option)
{
  pad.cd();
  
  g2 = new TGraph();
  int i=0;
  for (auto it = ccmap.begin();it!=ccmap.end(); it++){
    WireCell::Point p = it->first->center();
    if (it->second > 2000){
      g2->SetPoint(i,p.z/units::m,p.y/units::m);
      i++;
    }
  }
  g2->SetMarkerColor(8);
  g2->SetMarkerSize(0.8);
  g2->Draw(option);
  g2->SetMarkerStyle(26);
  
  return 0;
}

int ToyEventDisplay::draw_truthcells_charge(const WireCell::CellChargeMap& ccmap, TString option, int FI)
{
  pad.cd();

  //  g2 = new TGraph();
  //int i=0;

  Double_t x[100],y[100];

  for (auto it = ccmap.begin();it!=ccmap.end(); it++){
    //WireCell::Point p = it->first->center();
    
    WireCell::PointVector boundary = it->first->boundary();
    int n = 0;
    for (int i=0;i!=boundary.size();i++){
      x[n] = boundary[i].z/units::m;
      y[n] = boundary[i].y/units::m;
      n++;
    }
    x[n] = x[0];
    y[n] = y[0];
    n++;
    g3 = new TPolyLine(n,x,y);

    float charge = it->second;
    int color = (charge - charge_min)/(charge_max-charge_min)*254;
    if (color <0) color = 0;
    if (color>254) color =254;
    

    
    color += FI;
    g3->Draw(option);
    g3->SetFillColor(color);
    g3->SetLineColor(2);
    g3->SetLineWidth(1.5);
    g3->Draw("fsame");
    
  }
  
  return 0;
}


int ToyEventDisplay::draw_wires_charge(const WireCell::WireChargeMap& wcmap, TString option, int FI)
{
  pad.cd();
  Double_t x[100],y[100];

  for (auto it = wcmap.begin();it!=wcmap.end(); it++){

    const WireCell::GeomWire *wire = it->first;
    float charge = it->second;
    
    float pitch = gds.pitch(wire->plane());
    float angle = gds.angle(wire->plane());
    
    
    x[0] = wire->point1().z/units::m - pitch/2./units::m/std::cos(angle/units::radian); 
    y[0] = wire->point1().y/units::m;
    x[1] = wire->point2().z/units::m - pitch/2./units::m/std::cos(angle/units::radian);
    y[1] = wire->point2().y/units::m;
    x[3] = wire->point1().z/units::m + pitch/2./units::m/std::cos(angle/units::radian);
    y[3] = wire->point1().y/units::m;
    x[2] = wire->point2().z/units::m + pitch/2./units::m/std::cos(angle/units::radian);
    y[2] = wire->point2().y/units::m;
    x[4] = x[0];
    y[4] = y[0];
    g3 = new TPolyLine(5,x,y);

    int color = (charge - charge_min)/(charge_max-charge_min)*254;
    if (color <0) color = 0;
    if (color>254) color =254;
    color += FI;
    g3->Draw(option);
    g3->SetFillColor(color);
  
  }
  return 0;
}



int ToyEventDisplay::draw_cells_charge(const WireCell::GeomCellSelection& cellall, TString option)
{
  pad.cd();

 
  Double_t x[100],y[100];
  for (int i=0;i!=cellall.size();i++){
    
    WireCell::PointVector boundary = cellall[i]->boundary();
    int n = 0;
    for (int i=0;i!=boundary.size();i++){
      x[n] = boundary[i].z/units::m;
      y[n] = boundary[i].y/units::m;
      n++;
    }
    x[n] = x[0];
    y[n] = y[0];
    n++;
    g3 = new TPolyLine(n,x,y);
    g3->Draw(option);
    g3->SetFillColor(10);
  }
 
  
  return 0;
}
