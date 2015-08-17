#include "WireCell2dToy/ClusterDisplay.h"
#include "TGraph2D.h"
#include "TVector3.h"
#include "TPolyLine3D.h"
#include "TRandom.h"

using namespace WireCell;

WireCell2dToy::ClusterDisplay::ClusterDisplay(TPad& pad)
  :pad(pad)
{
}

WireCell2dToy::ClusterDisplay::~ClusterDisplay(){
}


void WireCell2dToy::ClusterDisplay::DrawVertex(WCVertexSelection& vertices, TString option){
  TGraph2D *g1 = new TGraph2D();
  int n = 0;
  double x,y,z;
  for (int i=0;i!=vertices.size();i++){
    Point center = vertices.at(i)->Center();
    x = center.x/units::cm;// + gRandom->Uniform(-0.15,0.15);
    y = center.y/units::cm;// + gRandom->Uniform(-0.15,0.15);
    z = center.z/units::cm;// + gRandom->Uniform(-0.15,0.15);
    
    g1->SetPoint(n,x,y,z);
    n++;

    double x1[2],y1[2],z1[2];
    WCTrackSelection& tracks = vertices.at(i)->get_tracks();
    if (tracks.size()>=1){
      for (int j=0;j!=tracks.size();j++){
    	WCTrack *track = tracks.at(j);
    	Point p1 = track->get_end_scells().at(0)->Get_Center();
    	Point p2 = track->get_end_scells().at(1)->Get_Center();
    	x1[0] = p1.x/units::cm;
    	y1[0] = p1.y/units::cm;
    	z1[0] = p1.z/units::cm;
	
    	x1[1] = p2.x/units::cm;
    	y1[1] = p2.y/units::cm;
    	z1[1] = p2.z/units::cm;
    	TPolyLine3D *l1 = new TPolyLine3D(2,x1,y1,z1);
    	l1->Draw("same");
    	l1->SetLineColor(6);
    	// std::cout << x1[0] << " " << y1[0] << " " << z1[0] << " " 
    	// 	  << x1[1] << " " << y1[1] << " " << z1[1] << std::endl;
      }
    }

    // if (tracks.size()>1 ){
    //   for (int j=0;j!=tracks.size();j++){
    // 	double ky = vertices.at(i)->get_ky(j);
    // 	double kz = vertices.at(i)->get_kz(j);
	
    // 	x1[0] = x-2.5;
    // 	y1[0] = ky * (x1[0]-x) + y;
    // 	z1[0] = kz * (x1[0]-x) + z;
	
    // 	x1[1] = x+2.5;
    // 	y1[1] = ky * (x1[1]-x) + y;
    // 	z1[1] = kz * (x1[1]-x) + z;
    // 	TPolyLine3D *l1 = new TPolyLine3D(2,x1,y1,z1);
    // 	l1->Draw("same");
    // 	l1->SetLineColor(4);
    //   }
    // }
    

  }
  
  std::cout <<"abc " << n << " " <<std::endl;

  g1->SetMarkerColor(2);
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(1.0);
  g1->Draw(option);
  
  // Draw tracks associated with it
  
  
  
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

void WireCell2dToy::ClusterDisplay::DrawCrawler(WireCell2dToy::ToyCrawler& toycrawler, TString option, int flag){
  
  std::cout << "Draw Crawler " << " " << toycrawler.Get_allCT().size() << " " << toycrawler.Get_allMCT().size()<< std::endl;

  int color[9]={3,4,5,6,7,8,9,2,1};
  //int style[7]={24,25,32,30,27,28,31};
  int style[7]={20,21,23,29,33,34,22};
  int num = 0;
  

  if (flag == 0 ){
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
      //std::cout << i << " " << n << std::endl;
      // Point center = clustertrack->Get_FirstMSCell()->Get_Center();
      // clustertrack->SC_Hough(center);
      // std::cout << i << " " << clustertrack->Get_Theta() << " " << clustertrack->Get_Phi() << " " << center.x/units::cm << " " << center.y/units::cm << " " << center.z/units::cm << std::endl;
      //std::cout << i << " " << n << " " << x << " " << y << " " << z << std::endl;
      
      //   if (i<12)
      
      
      
      //if (n> 7) 
      g1->Draw(option);
      // if (num == 9 ) num = 0;
      // g1->SetMarkerColor(color[num]);
      // g1->SetMarkerStyle(style[num]);

      int num1 = gRandom->Uniform(0,8.9);
      int num2 = gRandom->Uniform(0,6.9);

      g1->SetMarkerColor(color[num1]);
      g1->SetMarkerStyle(style[num2]);

      num++;
    }
  }else if (flag == 1){
    for (int i=0;i!=toycrawler.Get_allMCT().size();i++){
      MergeClusterTrack* clustertrack = toycrawler.Get_allMCT().at(i);
      
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
      //std::cout << i << " " << n << std::endl;
      
      //if (i<9)
      //if (n>5)
      g1->Draw(option);
      if (num ==7 ) num = 0;
      
      gRandom->SetSeed(0);
      int num1 = gRandom->Uniform(0,7.9);
      int num2 = gRandom->Uniform(0,6.9);

      g1->SetMarkerColor(color[num1]);
      g1->SetMarkerStyle(style[num2]);
      num++;
    }
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
      // std::cout << x << " " << y << " " << z << std::endl;
   
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
    //std::cout << "Xin1: " << x << " " << y << " " << z << std::endl;
    
  }
  // std::cout << n << std::endl;
  g1->Draw(option);
  //g1->Draw("p0");

  //test 
  // g1->GetXaxis()->SetRangeUser(50,70);
  // g1->GetYaxis()->SetRangeUser(-70,-30);
  // g1->GetZaxis()->SetRangeUser(280,294);

  //  std::cout << mcells.size() << std::endl;

  TGraph2D *g2 = new TGraph2D();
  for (int i=0;i!=mcells.size();i++){
    MergeSpaceCell *mcell = mcells.at(i);
    Point center = mcell->Get_Center();
    x = center.x/units::cm;
    y = center.y/units::cm;
    z = center.z/units::cm;
    g2->SetPoint(i,x,y,z);
    //std::cout << "Xin2: " << x << " " << y << " " << z << std::endl;
    
  }
  //g2->Draw("Psame");
  g2->SetMarkerStyle(22);
  g2->SetMarkerColor(2);
  
}
