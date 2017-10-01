#include "WireCellSst/GeomDataSource.h"
#include "WireCell2dToy/uBooNE_light_reco.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>

using namespace WireCell;
using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 4) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num " << endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);

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
  int eve_num = atoi(argv[3]);

  WireCell2dToy::uBooNE_light_reco uboone_flash(root_file);
  uboone_flash.load_event(eve_num);

  TFile *file = new TFile("temp.root","RECREATE");
  TH2F *h1 = new TH2F("h1","h1",1500,0,1500,32,0,32);
  TH2F *h2 = new TH2F("h2","h2",250,0,250,32,0,32);
  h1->SetDirectory(file);
  h2->SetDirectory(file);
  for (int i=0;i!=32;i++){
    TH1F *h3 = uboone_flash.get_raw_hist(i);
    TH1F *h4 = uboone_flash.get_decon_hist(i);
    for (int j=0;j!=1500;j++){
      h1->SetBinContent(j+1,i+1,h3->GetBinContent(j+1));
    }
    for (int j=0;j!=250;j++){
      h2->SetBinContent(j+1,i+1,h4->GetBinContent(j+1));
    }
  }
  file->Write();
  file->Close();
  
  return 1;
}
