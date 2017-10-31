#include "WireCell2dToy/ToyMatching.h"

#include "TChain.h"

using namespace WireCell;

void WireCell2dToy::tpc_light_match(int time_offset, std::map<WireCell::PR3DCluster*,std::vector<std::pair<WireCell::PR3DCluster*,double>>>& group_clusters, WireCell::OpflashSelection& flashes){
  TChain *T = new TChain("/pmtresponse/PhotonLibraryData","/pmtresponse/PhotonLibraryData");
  T->AddFile("uboone_photon_library.root");
  std::cout << T->GetEntries();
}
