#ifndef WCP2dToy_TOYMATCHING_H
#define WCP2dToy_TOYMATCHING_H

#include "WCPData/PR3DCluster.h"
//#include "WCP2dToy/uBooNE_light_reco.h"
#include "WCP2dToy/ToyLightReco.h"
#include "WCPData/FlashTPCBundle.h"
#include "WCPData/PhotonLibrary.h"

namespace WCP2dToy{

  /* class Photon_Library { */
  /*   public: */
  /*     std::map<int,int> map_lib_pmt, map_pmt_lib; */
  /*     std::vector<std::list<std::pair<int,float>>> library; */
  /*     double scaling_light_mag, rel_light_yield_err; */

  /*     Photon_Library(Int_t run_no = 0, bool flag_data = true, bool flag_add_light_yield_err = false); */
  /* }; */

  // time_offset in us
  int convert_xyz_voxel_id(WCP::Point& p);

  void calculate_pred_pe(int run_no, int time_offset, int nrebin, double time_slice_width, WCP::Photon_Library *pl, WCP::FlashTPCBundle* bundle, std::vector<double>* pred_pmt_light, 
			std::vector<std::pair<WCP::PR3DCluster*,double>>* additional_clusters, WCP::PR3DClusterSelection* other_clusters, WCP::PR3DClusterSelection* more_clusters, bool &flag_good_bundle, bool flag_data);

  WCP::FlashTPCBundleSelection tpc_light_match(int time_offset, int nrebin, WCP::Photon_Library *pl, std::map<WCP::PR3DCluster*,std::vector<std::pair<WCP::PR3DCluster*,double>>>& group_clusters, WCP::OpflashSelection& flashes, Int_t runno = 0, bool flag_data = true, bool flag_add_light_yield_err = false);

  // WCP::FlashTPCBundleSelection tpc_light_match_ana(int time_offset, int nrebin, std::map<WCP::PR3DCluster*,std::vector<std::pair<WCP::PR3DCluster*,double>>>& group_clusters, WCP::OpflashSelection& flashes);

  void organize_matched_bundles(WCP::FlashTPCBundleSelection& results_bundles, Double_t *cos_pe_low, Double_t *cos_pe_mid, std::map<std::pair<WCP::Opflash*,WCP::PR3DCluster*>,WCP::FlashTPCBundle*>& fc_bundles_map);
  
  
}

#endif
