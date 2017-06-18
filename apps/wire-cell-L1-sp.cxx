#include "WireCellSst/GeomDataSource.h"
#include "WireCellRess/LassoModel.h"
#include "WireCellRess/ElasticNetModel.h"


#include "TFile.h"
#include "TH2F.h"

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
  

  TH1F *hsig = new TH1F("hsig","hsig",nbins,0,nbins);
  for (int i=0;i!=nbins;i++){
    hsig->SetBinContent(i+1,htemp->GetBinContent(chid+1,i+1));
  }

  

  // hsig->Draw();

  TFile *file1 = new TFile("L1_sp.root","RECREATE");
  hsig->SetDirectory(file1);
  file1->Write();
  file1->Close();


  return 0;
}
