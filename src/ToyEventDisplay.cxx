#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCellData/Units.h"
#include "WireCellData/Point.h"
#include "WireCellData/MergeGeomWire.h"
//#include "WireCellData/MergeGeomCell.h"
#include "TGraph.h"
#include "TLine.h"
#include <iostream>

using namespace WireCell2dToy;

ToyEventDisplay::ToyEventDisplay(TPad& pad, const WireCell::GeomDataSource& gds)
    : pad(pad)
    , gds_flag(0)
    , gds(&gds)
    , dgds(0)
    , h1(0)
    , h2(0)
    , g1(0)
    , g2(0)
{
  recon_threshold = 1500;
  truth_threshold = 2000;
}

ToyEventDisplay::ToyEventDisplay(TPad& pad, const WireCell::DetectorGDS& gds)
    : pad(pad)
    , gds_flag(1)
    , gds(0)
    , dgds(&gds)
    , h1(0)
    , h2(0)
    , g1(0)
    , g2(0)
{
  recon_threshold = 1500;
  truth_threshold = 2000;
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
    if (h2) {
      delete h2;
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
  if (gds_flag == 0){
    h1->SetTitle("Wires and True Hits ");
  }else{
    h1->SetTitle("Wires and True Hits (B face +x)");
  }
  h1->GetYaxis()->SetNdivisions(506);
  h1->GetXaxis()->SetNdivisions(506);
  h1->SetXTitle("Z (m)");
  h1->SetYTitle("Y (m)");
  h1->SetBinContent(0,0,1);
  h1->GetZaxis()->SetRangeUser(charge_min,charge_max);

  if (gds_flag == 1){
    h2 = new TH2F("h2","h2",1000,x_min,x_max,1000,y_min,y_max);
    h2->SetTitle("Wires and True Hits (A face -x)");
    h2->GetYaxis()->SetNdivisions(506);
    h2->GetXaxis()->SetNdivisions(506);
    h2->SetXTitle("Z (m)");
    h2->SetYTitle("Y (m)");
    h2->SetBinContent(0,0,1);
    h2->GetZaxis()->SetRangeUser(charge_min,charge_max);
  }

  return 0;
}

int ToyEventDisplay::draw_mc(int flag, const WireCell::PointValueVector& mctruth, TString option)
{
  if (gds_flag == 0 ){
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
  }else{
    pad.cd(1);
    h1->Draw(option);
    
    pad.cd(2);
    h2->Draw(option);
  }
  
  return 0;
}


void ToyEventDisplay::draw_bad_cell(WireCell::GeomCellSelection& cells){
  pad.cd();
  for (int i=0;i!=cells.size();i++){
    double x[100];
    double y[100];
    
    for (int j=0;j!=cells.at(i)->boundary().size();j++){
      x[j] = cells.at(i)->boundary().at(j).z/units::m;
      y[j] = cells.at(i)->boundary().at(j).y/units::m;
      //std::cout << x[j] << " " << y[j] << std::endl;
    }
    x[cells.at(i)->boundary().size()] = x[0];
    y[cells.at(i)->boundary().size()] = y[0];
    // std::cout << x[4] << " " << y[4] << std::endl;
    // std::cout << std::endl;

    TPolyLine *pline = new TPolyLine(cells.at(i)->boundary().size()+1,x,y);
    pline->SetFillColor(5);
    pline->SetLineColor(5);
    //pline->SetLineWidth(1);
    pline->Draw("fsame");
    pline->Draw("same");
  }
}

void ToyEventDisplay::draw_bad_region(WireCell::ChirpMap& chirpmap, int time, int scale, int plane, TString option){
  pad.cd();

  // int nwire_u = gds->wires_in_plane(WirePlaneType_t(0)).size();
  // int nwire_v = gds->wires_in_plane(WirePlaneType_t(1)).size();
  // int nwire_w = gds->wires_in_plane(WirePlaneType_t(2)).size();
  
  //  std::cout << chirpmap.size() << std::endl;

  for (auto it = chirpmap.begin(); it!=chirpmap.end(); it++){
    if (time >= it->second.first/scale && time <= it->second.second/scale){
      const WireCell::GeomWire *wire = gds->by_planeindex(WireCell::WirePlaneType_t(plane),it->first);
      float pitch = gds->pitch(wire->plane());
      float angle = gds->angle(wire->plane());
      TLine *l3 = new TLine(wire->point1().z/units::m  ,wire->point1().y/units::m,
			    wire->point2().z/units::m  ,wire->point2().y/units::m);
      l3->SetLineColor(5);
      l3->Draw(option);
    }
  }

}


int ToyEventDisplay::draw_merged_wires(WireCell::GeomWireSelection wires, TString option, int color){
   for (int j=0;j!=wires.size();j++){
      
     WireCell::MergeGeomWire *wire = (WireCell::MergeGeomWire*)wires.at(j);
     const WireCell::GeomWire *wire1 = wire->get_allwire().front();
     const WireCell::GeomWire *wire2 = wire->get_allwire().back();
          
      TLine *l3 = new TLine(wire1->point1().z/units::m  ,wire1->point1().y/units::m,
      			    wire1->point2().z/units::m  ,wire1->point2().y/units::m);
      l3->SetLineColor(color);
      l3->Draw(option);
      l3->SetLineStyle(1);

      TLine *l1 = new TLine(wire2->point1().z/units::m  ,wire2->point1().y/units::m,
      			    wire2->point2().z/units::m  ,wire2->point2().y/units::m);
      l1->SetLineColor(color);
      l1->Draw(option);
      l1->SetLineStyle(1);
    }
    
    // break;
    



  
  return 0;
}

int ToyEventDisplay::draw_wires(WireCell::GeomWireSelection& wires, TString option){
  
 
    for (int j=0;j!=wires.size();j++){
      
      const WireCell::GeomWire *wire = wires.at(j);
      
      
      
       float pitch = gds->pitch(wire->plane());
      float angle = gds->angle(wire->plane());
      //std::cout << pitch/units::mm << " " << angle << std::endl;
      
      TLine *l1 = new TLine(wire->point1().z/units::m - pitch/2./units::m/std::cos(angle/units::radian) ,wire->point1().y/units::m,
			    wire->point2().z/units::m - pitch/2./units::m/std::cos(angle/units::radian) ,wire->point2().y/units::m);
      l1->Draw(option);
      l1->SetLineStyle(1);
      
      TLine *l2 = new TLine(wire->point1().z/units::m + pitch/2./units::m/std::cos(angle/units::radian) ,wire->point1().y/units::m,
			    wire->point2().z/units::m + pitch/2./units::m/std::cos(angle/units::radian) ,wire->point2().y/units::m);
      l2->Draw(option);
      l2->SetLineStyle(1);
      
      TLine *l3 = new TLine(wire->point1().z/units::m  ,wire->point1().y/units::m,
			    wire->point2().z/units::m  ,wire->point2().y/units::m);
      l3->SetLineColor(2);
      l3->Draw(option);
      l3->SetLineStyle(1);
    }
    
    // break;
    



  
  return 0;
}

int ToyEventDisplay::draw_slice(const WireCell::Slice& slice, TString option)
{
  

  WireCell::Channel::Group group = slice.group();
  //std::cout << group.size() << std::endl;
  
  for (int i=0;i!=group.size();i++){

    if (gds_flag == 1 ){
      const WireCell::GeomWireSelection& wires = dgds->by_channel(group.at(i).first);
      //std::cout << wires.at(0)->plane() << " " << group.at(i).first << std::endl;
       for (int j=0;j!=wires.size();j++){
	 
	 const WireCell::GeomWire *wire = wires.at(j);
	 
	 // std::cout << wire->point1().z/units::m  << " " << wire->point1().y/units::m << " " << 
	 //   wire->point2().z/units::m  << " " << wire->point2().y/units::m << std::endl;
	 //std::cout << wire->cryo() << " " << wire->apa() << " " << wire->face() << std::endl;
	 float pitch;
	 float angle;
	 if (wire->face()==1){
	   pad.cd(1);
	   pitch = dgds->get_pitch(wire->cryo(),wire->plane());
	   angle = -dgds->get_angle(wire->cryo(),wire->plane());
	 }else if (wire->face()==0){
	   pad.cd(2);
	   pitch = dgds->get_pitch(wire->cryo(),wire->plane());
	   angle = dgds->get_angle(wire->cryo(),wire->plane());
	 }	 
	 //	 std::cout << pitch/units::mm << " " << angle << std::endl;
	 
	 TLine *l1 = new TLine(wire->point1().z/units::m - pitch/2./units::m/std::cos(angle/units::radian) ,wire->point1().y/units::m,
			       wire->point2().z/units::m - pitch/2./units::m/std::cos(angle/units::radian) ,wire->point2().y/units::m);
	 l1->Draw(option);
	 l1->SetLineStyle(7);

	 TLine *l2 = new TLine(wire->point1().z/units::m + pitch/2./units::m/std::cos(angle/units::radian) ,wire->point1().y/units::m,
			       wire->point2().z/units::m + pitch/2./units::m/std::cos(angle/units::radian) ,wire->point2().y/units::m);
	 l2->Draw(option);
	 l2->SetLineStyle(7);
	 
	 TLine *l3 = new TLine(wire->point1().z/units::m  ,wire->point1().y/units::m,
			       wire->point2().z/units::m  ,wire->point2().y/units::m);
	 l3->SetLineColor(2);
	 l3->Draw(option);
	 l3->SetLineStyle(7);
       }
       
       // break;

    }else if (gds_flag == 0 ){
      pad.cd();
      //std::cout << group.at(i).first << std::endl;
      const WireCell::GeomWire *wire = gds->by_channel_segment(group.at(i).first,0);
      // std::cout << wire->channel() << std::endl;
      // if ( wire->channel() ==1429 || wire->channel() ==4461){
      
      float pitch = gds->pitch(wire->plane());
      float angle = gds->angle(wire->plane());
      
      //std::cout << pitch << " " << angle << std::endl;
      
      // std::cout << wire->point1().y << " " << wire->point1().z << std::endl;
      TLine *l1 = new TLine(wire->point1().z/units::m - pitch/2./units::m/std::cos(angle/units::radian) ,wire->point1().y/units::m,
			    wire->point2().z/units::m - pitch/2./units::m/std::cos(angle/units::radian) ,wire->point2().y/units::m);
      l1->Draw(option);
      l1->SetLineStyle(7);

      TLine *l2 = new TLine(wire->point1().z/units::m + pitch/2./units::m/std::cos(angle/units::radian) ,wire->point1().y/units::m,
			    wire->point2().z/units::m + pitch/2./units::m/std::cos(angle/units::radian) ,wire->point2().y/units::m);
      l2->Draw(option);
      l2->SetLineStyle(7);
      
      TLine *l3 = new TLine(wire->point1().z/units::m  ,wire->point1().y/units::m,
			    wire->point2().z/units::m  ,wire->point2().y/units::m);
      l3->SetLineColor(2);
      l3->Draw(option);
      l3->SetLineStyle(7);
    }
  }
    // }
  return 0;
}

int ToyEventDisplay::draw_points(WireCell::PointVector pcells, TString option, int color){
  TGraph *g1 = new TGraph();
  for (int i=0; i!=pcells.size(); i++){
    g1->SetPoint(i,pcells.at(i).z/units::m, pcells.at(i).y/units::m);
  }
  g1->Draw(option);
  g1->SetMarkerColor(color);
  g1->SetMarkerSize(1.0);
  g1->SetMarkerStyle(20);
}


int ToyEventDisplay::draw_cells(const WireCell::GeomCellSelection& cellall, TString option, int color)
{
  if (gds_flag == 1){
    g2 = new TGraph();
    g2b = new TGraph();
    int nf = 0;
    int nb = 0;
    for (int i=0;i!=cellall.size();i++){
      if (cellall[i]->get_face()==1){
	g2->SetPoint(nf,cellall[i]->center().z/units::m,cellall[i]->center().y/units::m);
	nf ++;
      }else{
	g2b->SetPoint(nb,cellall[i]->center().z/units::m,cellall[i]->center().y/units::m);
	nb ++;
      }
    }
    pad.cd(1);
    g2->SetMarkerColor(color);
    g2->SetMarkerSize(0.8);
    g2->Draw(option);
    g2->SetMarkerStyle(21);

    pad.cd(2);
    g2b->SetMarkerColor(color);
    g2b->SetMarkerSize(0.8);
    if (g2b->GetN()!=0)
      g2b->Draw(option);
    g2b->SetMarkerStyle(21);
  }else{
    pad.cd();
    
    //if (g2) delete g2;
    g2 = new TGraph();
    for (int i=0;i!=cellall.size();i++){
      g2->SetPoint(i,cellall[i]->center().z/units::m,cellall[i]->center().y/units::m);
    }
    g2->SetMarkerColor(color);
    g2->SetMarkerSize(0.8);
    g2->Draw(option);
    g2->SetMarkerStyle(21);
  }
  return 0;
}


int ToyEventDisplay::draw_mergecells(const WireCell::GeomCellSelection& cellall, TString option, int flag)
{
 
  
  int np = 0;
  int npb = 0;
  g2 = new TGraph();
  g2b = new TGraph();
  
  for (int i=0;i!=cellall.size();i++){
    const WireCell::MergeGeomCell* mcell = (const WireCell::MergeGeomCell*)cellall.at(i);
    int face = mcell->get_allcell().at(0)->get_face();

    if (flag==0){
      
      if (face == 1){
	g2->SetPoint(np,cellall[i]->center().z/units::m,cellall[i]->center().y/units::m);
	//std::cout << cellall[i]->center().z/units::m << " " << cellall[i]->center().y/units::m << std::endl;
	np++;
      }else{
	g2b->SetPoint(npb,cellall[i]->center().z/units::m,cellall[i]->center().y/units::m);
	//std::cout << cellall[i]->center().z/units::m << " " << cellall[i]->center().y/units::m << std::endl;
	npb++;
      }
    }else if (flag==1){
      WireCell::MergeGeomCell *mcell = (WireCell::MergeGeomCell*)cellall[i];
      if (mcell->GetContainTruthCell()){
	if (face == 1){
	  g2->SetPoint(np,cellall[i]->center().z/units::m,cellall[i]->center().y/units::m);
	  np++;
	}else{
	  g2b->SetPoint(npb,cellall[i]->center().z/units::m,cellall[i]->center().y/units::m);
	  npb++;
	}
      }
    }
  }

  
  pad.cd(1);
  g2->SetMarkerColor(2);
  g2->SetMarkerSize(0.8);
  if (g2->GetN()!=0)
    g2->Draw(option);
  g2->SetMarkerStyle(24);
  
  pad.cd(2);
  g2b->SetMarkerColor(2);
  g2b->SetMarkerSize(0.8);
  if (g2b->GetN()!=0)
    g2b->Draw(option);
  g2b->SetMarkerStyle(24);

  return 0;
}



int ToyEventDisplay::draw_truthcells(const WireCell::CellChargeMap& ccmap, TString option)
{
  if (gds_flag==1){
    g2 = new TGraph();
    g2b = new TGraph();
    int nf = 0;
    int nb = 0;
    
    int i=0;
    for (auto it = ccmap.begin();it!=ccmap.end(); it++){
      WireCell::Point p = it->first->center();
      if (it->first->get_face()==1){
	if (it->second > truth_threshold){
	  g2->SetPoint(nf,p.z/units::m,p.y/units::m);
	  nf++;
	}
      }else{
	if (it->second > truth_threshold){
	  g2b->SetPoint(nb,p.z/units::m,p.y/units::m);
	  nb++;
	}
      }
      
    }
    pad.cd(1);
    g2->SetMarkerColor(8);
    //g2->SetMarkerColor(1);
    g2->SetMarkerSize(0.8);
    g2->Draw(option);
    g2->SetMarkerStyle(26);

    pad.cd(2);
    g2b->SetMarkerColor(8);
    //g2->SetMarkerColor(1);
    g2b->SetMarkerSize(0.8);
    if (g2b->GetN()!=0)
      g2b->Draw(option);
    g2b->SetMarkerStyle(26);

  }else{
    pad.cd();
    
    g2 = new TGraph();
    int i=0;
    for (auto it = ccmap.begin();it!=ccmap.end(); it++){
      WireCell::Point p = it->first->center();
      if (it->second > truth_threshold){
	g2->SetPoint(i,p.z/units::m,p.y/units::m);
	i++;
      }
    }
    g2->SetMarkerColor(8);
    //g2->SetMarkerColor(1);
    g2->SetMarkerSize(0.8);
    g2->Draw(option);
    g2->SetMarkerStyle(26);
  }
  
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
    
    float pitch = gds->pitch(wire->plane());
    float angle = gds->angle(wire->plane());
    
    
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
    g3->SetLineColor(1);
    g3->SetLineWidth(1);
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

int ToyEventDisplay::draw_reconcells(const WireCell::GeomCellSelection& cellall, WireCell2dToy::ToyMatrix *toymatrix ,TString option, int color){
  
  pad.cd();
  g2 = new TGraph();
  int ncount = 0;
  for (int i=0;i!=cellall.size();i++){
    WireCell::MergeGeomCell *mcell = (WireCell::MergeGeomCell*)cellall[i];
    if (toymatrix->Get_Cell_Charge(mcell)>recon_threshold){
      WireCell::GeomCellSelection acell = mcell->get_allcell();
      for (int j=0;j!=acell.size();j++){
	g2->SetPoint(ncount,acell[j]->center().z/units::m,acell[j]->center().y/units::m);
	ncount ++;
      }
    }
  }
  
  g2->SetMarkerColor(color);
  g2->SetMarkerSize(0.8);
  g2->Draw(option);
  g2->SetMarkerStyle(21);
  
  return 0;

  
}
