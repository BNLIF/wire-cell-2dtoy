#ifndef WCP2dToy_TOYDATAQUALITY_H
#define WCP2dToy_TOYDATAQUALITY_H

#include "WCPData/GeomWire.h"

#include "TH2F.h"
#include <iostream>

namespace WCP2dToy{
  int Noisy_Event_ID(TH2F *hu_decon, TH2F *hv_decon, TH2F *hw_decon, std::vector<float>& uplane_rms, std::vector<float>& vplane_rms, std::vector<float>& wplane_rms, WCP::ChirpMap& uplane_map, WCP::ChirpMap& vplane_map, WCP::ChirpMap& wplane_map, TH2F *hu_decon_g, TH2F *hv_decon_g, TH2F *hw_decon_g, int nrebin, TH2F *hv_raw, bool flag_corr = false);

  void Organize_Dead_Channels(WCP::ChirpMap& uplane_map, WCP::ChirpMap& vplane_map, WCP::ChirpMap& wplane_map, int nbins, int nrebin, int n_div=5);
}

#endif
