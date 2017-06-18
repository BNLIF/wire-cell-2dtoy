#include "WireCellSst/GeomDataSource.h"
#include "WireCellRess/LassoModel.h"
#include "WireCellRess/ElasticNetModel.h"


#include "TFile.h"
#include "TH2F.h"
#include "TGraph.h"

#include <Eigen/Dense>
using namespace Eigen;

using namespace WireCell;
using namespace std;




int main(int argc, char* argv[])
{
  if (argc < 4) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/magnify.root ch_id " << endl;
  }
  
  WireCellSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;

  cout << "Pitch: " << gds.pitch(WirePlaneType_t(0)) 
       << " " << gds.pitch(WirePlaneType_t(1)) 
       << " " << gds.pitch(WirePlaneType_t(2))
       << endl;
  cout << "Angle: " << gds.angle(WirePlaneType_t(0)) 
       << " " << gds.angle(WirePlaneType_t(1)) 
       << " " << gds.angle(WirePlaneType_t(2))
       << endl;
  
  const char* root_file = argv[2];
  int chid = atoi(argv[3]);

  TFile *file = new TFile(root_file);
  
  TH2F *hu_raw, *hv_raw, *hw_raw;
  hu_raw = (TH2F*)file->Get("hu_raw");
  hv_raw = (TH2F*)file->Get("hv_raw");
  hw_raw = (TH2F*)file->Get("hw_raw");
  const int nbins = hu_raw->GetNbinsY();
  int nwire_u = hu_raw->GetNbinsX();
  int nwire_v = hv_raw->GetNbinsX();
  int nwire_w = hw_raw->GetNbinsX();
  
  TH2F *htemp;
  if (chid < nwire_u){
    htemp = hu_raw;
  }else if (chid < nwire_v+nwire_u){
    htemp = hv_raw;
    chid -= nwire_u;
  }else{
    htemp = hw_raw;
    chid -= nwire_u + nwire_v;
  }
  
  const int nbin_fit = 70;

  TH1F *hsig = new TH1F("hsig","hsig",nbins,0,nbins);
  TH1F *hsig1 = new TH1F("hsig1","hsig1",nbin_fit,0,nbin_fit);
  TH1F *hsig_w = new TH1F("hsig_w","hsig_w",nbin_fit,0,nbin_fit);
  TH1F *hsig_v = new TH1F("hsig_v","hsig_v",nbin_fit,0,nbin_fit);
  
  
  for (int i=0;i!=nbins;i++){
    hsig->SetBinContent(i+1,htemp->GetBinContent(chid+1,i+1));
  }
  for (int i=0;i!=nbin_fit;i++){
    hsig1->SetBinContent(i+1,hsig->GetBinContent(3450+i+1));
  }
  

  #include "/home/xqian/uboone/matrix_inversion/work/wire-cell/2dtoy/src/data_70_2D_11.txt"

  // read in the response functions ...
  // collection response ...
  TGraph **gw_2D_g = new TGraph*[11];
  gw_2D_g[0] = new TGraph(5000,w_2D_g_0_x,w_2D_g_0_y);
  gw_2D_g[1] = new TGraph(5000,w_2D_g_1_x,w_2D_g_1_y);
  gw_2D_g[2] = new TGraph(5000,w_2D_g_2_x,w_2D_g_2_y);
  gw_2D_g[3] = new TGraph(5000,w_2D_g_3_x,w_2D_g_3_y);
  gw_2D_g[4] = new TGraph(5000,w_2D_g_4_x,w_2D_g_4_y);
  gw_2D_g[5] = new TGraph(5000,w_2D_g_5_x,w_2D_g_5_y);
  gw_2D_g[6] = new TGraph(5000,w_2D_g_6_x,w_2D_g_6_y);
  gw_2D_g[7] = new TGraph(5000,w_2D_g_7_x,w_2D_g_7_y);
  gw_2D_g[8] = new TGraph(5000,w_2D_g_8_x,w_2D_g_8_y);
  gw_2D_g[9] = new TGraph(5000,w_2D_g_9_x,w_2D_g_9_y);
  gw_2D_g[10] = new TGraph(5000,w_2D_g_10_x,w_2D_g_10_y);

  // induction add all 11 thing ...
  TGraph **gv_2D_g = new TGraph*[11];
  gv_2D_g[0] = new TGraph(5000,v_2D_g_0_x,v_2D_g_0_y);
  gv_2D_g[1] = new TGraph(5000,v_2D_g_1_x,v_2D_g_1_y);
  gv_2D_g[2] = new TGraph(5000,v_2D_g_2_x,v_2D_g_2_y);
  gv_2D_g[3] = new TGraph(5000,v_2D_g_3_x,v_2D_g_3_y);
  gv_2D_g[4] = new TGraph(5000,v_2D_g_4_x,v_2D_g_4_y);
  gv_2D_g[5] = new TGraph(5000,v_2D_g_5_x,v_2D_g_5_y);
  gv_2D_g[6] = new TGraph(5000,v_2D_g_6_x,v_2D_g_6_y);
  gv_2D_g[7] = new TGraph(5000,v_2D_g_7_x,v_2D_g_7_y);
  gv_2D_g[8] = new TGraph(5000,v_2D_g_8_x,v_2D_g_8_y);
  gv_2D_g[9] = new TGraph(5000,v_2D_g_9_x,v_2D_g_9_y);
  gv_2D_g[10] = new TGraph(5000,v_2D_g_10_x,v_2D_g_10_y);

  TGraph *gw = new TGraph();
  TGraph *gv = new TGraph();
  
  for (Int_t i=0;i!=5000;i++){
    double x,y;
    gw_2D_g[0]->GetPoint(i,x,y);
    double sum1 = gw_2D_g[0]->Eval(x);
    double sum2 = gv_2D_g[0]->Eval(x);
    for (int j=1;j!=11;j++){
      sum1 += gw_2D_g[j]->Eval(x)*2.;
      sum2 += gv_2D_g[j]->Eval(x)*2.;
    }
    gw->SetPoint(i,x,sum1);
    gv->SetPoint(i,x,sum2);
  }

  // read in the waveform ... 
  VectorXd W = VectorXd::Zero(nbin_fit);
  for(int i=0;i!=nbin_fit;i++){
    W(i) = hsig1->GetBinContent(i+1);
    //std::cout << W(i) << std::endl;
  }

  int scaling = 4096/2000.*14.*1.2*500;  // as a scale of 500 electrons

  // form matrix G
  MatrixXd G = MatrixXd::Zero(nbin_fit,nbin_fit*2);
  for (int i=0;i!=nbin_fit;i++){ 
    // X 
    Double_t t1 = i/2.; // us, measured time  
    for (int j=0;j!=nbin_fit;j++){
      // Y ... 
      Double_t t2 = j/2.; // us, real signal time
      double delta_t = t1 - t2;
      if (delta_t >-84 && delta_t < 15.8){
	G(i,j) = gw->Eval(delta_t+3) * scaling;
	G(i,nbin_fit+j) = gv->Eval(delta_t) * scaling;
      }
    }
  }
  // solve ... 
  double lambda = 0.1;
  WireCell::LassoModel m2(lambda, 100000, 1e-3);
  m2.SetData(G, W);
  m2.Fit();
  
  VectorXd beta = m2.Getbeta();
  int nbeta = beta.size();
  for (int i=0;i!=nbin_fit;i++){
    hsig_w->SetBinContent(i+1,beta(i));
    hsig_v->SetBinContent(i+1,beta(nbin_fit+i));
  }

  // hsig->Draw();

  TFile *file1 = new TFile("L1_sp.root","RECREATE");
  hsig->SetDirectory(file1);
  hsig1->SetDirectory(file1);
  
  hsig_w->SetDirectory(file1);
  hsig_v->SetDirectory(file1);

  gw->Write("gw");
  gv->Write("gv");
  file1->Write();
  file1->Close();


  return 0;
}
