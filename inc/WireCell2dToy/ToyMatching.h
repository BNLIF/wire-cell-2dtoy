#ifndef WireCell2dToy_TOYMATCHING_H
#define WireCell2dToy_TOYMATCHING_H

#include "WireCellData/PR3DCluster.h"
//#include "WireCell2dToy/uBooNE_light_reco.h"
#include "WireCell2dToy/ToyLightReco.h"
#include "WireCellData/FlashTPCBundle.h"
#include "WireCellData/PhotonLibrary.h"

namespace WireCell2dToy{

  /* class Photon_Library { */
  /*   public: */
  /*     std::map<int,int> map_lib_pmt, map_pmt_lib; */
  /*     std::vector<std::list<std::pair<int,float>>> library; */
  /*     double scaling_light_mag, rel_light_yield_err; */

  /*     Photon_Library(Int_t run_no = 0, bool flag_data = true, bool flag_add_light_yield_err = false); */
  /* }; */

  // time_offset in us
  int convert_xyz_voxel_id(WireCell::Point& p);

  void calculate_pred_pe(int run_no, int time_offset, int nrebin, double time_slice_width, WireCell::Photon_Library *pl, WireCell::FlashTPCBundle* bundle, std::vector<double>* pred_pmt_light, 
			std::vector<std::pair<WireCell::PR3DCluster*,double>>* additional_clusters, WireCell::PR3DClusterSelection* other_clusters, WireCell::PR3DClusterSelection* more_clusters, bool &flag_good_bundle, bool flag_data);

  WireCell::FlashTPCBundleSelection tpc_light_match(int time_offset, int nrebin, WireCell::Photon_Library *pl, std::map<WireCell::PR3DCluster*,std::vector<std::pair<WireCell::PR3DCluster*,double>>>& group_clusters, WireCell::OpflashSelection& flashes, Int_t runno = 0, bool flag_data = true, bool flag_add_light_yield_err = false);

  // WireCell::FlashTPCBundleSelection tpc_light_match_ana(int time_offset, int nrebin, std::map<WireCell::PR3DCluster*,std::vector<std::pair<WireCell::PR3DCluster*,double>>>& group_clusters, WireCell::OpflashSelection& flashes);

  void organize_matched_bundles(WireCell::FlashTPCBundleSelection& results_bundles, Double_t *cos_pe_low, Double_t *cos_pe_mid, std::map<std::pair<WireCell::Opflash*,WireCell::PR3DCluster*>,WireCell::FlashTPCBundle*>& fc_bundles_map);
  
  
}

#endif
