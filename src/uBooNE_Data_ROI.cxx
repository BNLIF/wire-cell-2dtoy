#include "WireCell2dToy/uBooNE_Data_ROI.h"

#include "WireCellData/GeomWire.h"


using namespace WireCell;

WireCell2dToy::uBooNEDataROI::uBooNEDataROI(WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds, WireCell::ChirpMap& umap, WireCell::ChirpMap& vmap, WireCell::ChirpMap& wmap)
  : fds(fds)
  , gds(gds)
  , umap(umap)
  , vmap(vmap)
  , wmap(wmap)
{
  find_ROI_by_itself();
}

WireCell2dToy::uBooNEDataROI::~uBooNEDataROI()
{
  
}

void WireCell2dToy::uBooNEDataROI::find_ROI_by_itself(){
  
}



void WireCell2dToy::uBooNEDataROI::restore_baseline(TH1F *htemp){
  //correct baseline 
  double max = htemp->GetMaximum();
  double min = htemp->GetMinimum();
  int nbin_b = max - min;
  if (nbin_b ==0) nbin_b = 1;
  TH1F *h1 = new TH1F("h1","h1",nbin_b,min,max);
  int nbin = htemp->GetNbinsX();
  for (int j=0;j!=nbin;j++){
    h1->Fill(htemp->GetBinContent(j+1));
  }
  float ped = h1->GetMaximumBin()*(max-min)/(nbin_b*1.) + min;
  float ave=0,ncount = 0;
  
  for (int j=0;j!=nbin;j++){
    if (fabs(htemp->GetBinContent(j+1)-ped)<400){
      ave +=htemp->GetBinContent(j+1);
      ncount ++;
      }
    }
    if (ncount==0) ncount=1;
    ave = ave/ncount;
    
    for (int j=0;j!=nbin;j++){
      double content = htemp->GetBinContent(j+1);
      content -= ave;
      htemp->SetBinContent(j+1,content);
    }
    delete h1;
}
