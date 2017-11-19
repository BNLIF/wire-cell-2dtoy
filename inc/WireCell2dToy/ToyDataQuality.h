#ifndef WireCell2dToy_TOYDATAQUALITY_H
#define WireCell2dToy_TOYDATAQUALITY_H

#include "TH2F.h"
#include <iostream>

namespace WireCell2dToy{
  bool Noisy_Event_ID(TH2F *hu_decon, TH2F *hv_decon, TH2F *hw_decon, TH1F *hu_threshold, TH1F *hv_threshold, TH1F *hw_threshold);
}

#endif
