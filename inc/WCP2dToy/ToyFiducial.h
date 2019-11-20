#ifndef WIRECELL_UBOONE_TOYFIDUCIAL_H
#define WIRECELL_UBOONE_TOYFIDUCIAL_H

#include "WCP2dToy/ToyMatching.h"

#include "WCPData/Units.h"
#include "WCPData/Point.h"
#include "WCPData/SlimMergeGeomCell.h"
#include "WCPData/FlashTPCBundle.h"
#include "WCPData/ToyCTPointCloud.h"
#include "WCPData/LMBDT.h"

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <vector>
#include <map>


#include "WCPData/Point.h"
#include "WCPData/ToyCTPointCloud.h"
#include "WCP2dToy/ImprovePR3DCluster.h"
#include "WCPData/PR3DCluster.h"

namespace WCP2dToy{
  class ToyFiducial{
  public:
    ToyFiducial(int dead_region_ch_ext = 3, double offset_t=800, double offset_u=0, double offset_v=0, double offset_w=0, double slope_t=1./2*units::mm, double slope_u=1./(3*units::mm), double slope_v=1./(3*units::mm), double slope_w=1./(3*units::mm), double angle_u=-1.047198, double angle=1.047198, double angle_w=0,
		double boundary_dis_cut=2*units::cm, double top=117*units::cm, double bottom=-116*units::cm, double upstream=0*units::cm, double downstream=1037*units::cm, double anode = 0*units::cm, double cathode=256*units::cm, int flag_data=1);
    ~ToyFiducial();

    void set_offset_t(double value){offset_t=value;};
    
    bool check_neutrino_candidate(WCP::PR3DCluster *main_cluster, WCP::WCPointCloud<double>::WCPoint& wcp1 ,WCP::WCPointCloud<double>::WCPoint& wcp2, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud, bool flag_2view_check = true);
    bool inside_fiducial_volume(WCP::Point& p, double offset_x=0, std::vector<double>* tolerance_vec=NULL);
    bool inside_dead_region(WCP::Point& p);
    bool check_dead_volume(WCP::Point& p, TVector3& dir, double step = 1.0*units::cm, double offset_x=0);
    
    bool check_signal_processing(WCP::Point& p, TVector3& dir, WCP::ToyCTPointCloud& ct_point_cloud, double step = 1.0*units::cm, double offset_x=0);

    bool check_tgm(WCP::FlashTPCBundle *bundle, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud,std::map<WCP::PR3DCluster*, WCP::PR3DCluster*>& old_new_cluster_map, int flag = 1);
    
    bool check_fully_contained(WCP::FlashTPCBundle *bundle, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud,std::map<WCP::PR3DCluster*, WCP::PR3DCluster*>& old_new_cluster_map, unsigned int* fail_mode=nullptr, int flag = 1);
    
    int check_LM(WCP::FlashTPCBundle *bundle, double& cluster_length);
    int check_LM_cuts(WCP::FlashTPCBundle *bundle, double& cluster_length);
    int check_LM_bdt(WCP::FlashTPCBundle *bundle, double& cluster_length);
      
    void AddDeadRegion(WCP::SlimMergeGeomCell* mcell, std::vector<int>& time_slices);

    std::vector<double> get_boundary_SCB_xy_x(WCP::Point& p){
      int index_z = floor(p.z/units::m);
      if(index_z<0){index_z=0;} else if(index_z>9){index_z=9;}
      return boundary_SCB_xy_x_array[index_z];
    }
    std::vector<double> get_boundary_SCB_xy_y(WCP::Point& p){
      int index_z = floor(p.z/units::m);
      if(index_z<0){index_z=0;} else if(index_z>9){index_z=9;}
      return boundary_SCB_xy_y_array[index_z];
    }
    std::vector<double> get_boundary_SCB_xz_x(WCP::Point& p){
      int index_y = floor((p.y/units::cm+116)/24);
      if(index_y<0){index_y=0;} else if(index_y>9){index_y=9;}
      return boundary_SCB_xz_x_array[index_y];
    }
    std::vector<double> get_boundary_SCB_xz_z(WCP::Point& p){
      int index_y = floor((p.y/units::cm+116)/24);
      if(index_y<0){index_y=0;} else if(index_y>9){index_y=9;}
      return boundary_SCB_xz_z_array[index_y];
    }
    int check_boundary(std::vector<std::vector<WCP::WCPointCloud<double>::WCPoint>> extreme_points, double offset_x, std::vector<double>* tol_vec);
    void cosmic_tagger(WCP::OpflashSelection& flashes,WCP::FlashTPCBundleSelection *matched_bundles, WCP::FlashTPCBundle* main_bundle, WCP::Photon_Library *pl,
      int time_offset, int nrebin, float unit_dis, WCP::ToyCTPointCloud& ct_point_cloud,
      std::map<WCP::PR3DCluster*, WCP::PR3DCluster*>& old_new_cluster_map, int run_no, int subrun_no, int event_no, bool flag_data, bool debug_tagger=false);
    
  protected:
    // boundary
    double m_top; // top distance
    double m_bottom; // bottom distance
    double m_upstream;
    double m_downstream;
    double m_anode;
    double m_cathode;
    
    // space charge boundary
    double m_sc_bottom_1_x, m_sc_bottom_1_y;
    double m_sc_bottom_2_x, m_sc_bottom_2_y;

    double m_sc_top_1_x, m_sc_top_1_y;
    double m_sc_top_2_x, m_sc_top_2_y;

    double m_sc_upstream_1_x, m_sc_upstream_1_z;
    double m_sc_upstream_2_x, m_sc_upstream_2_z;

    double m_sc_downstream_1_x, m_sc_downstream_1_z;
    double m_sc_downstream_2_x, m_sc_downstream_2_z;

    std::vector<double> boundary_xy_x, boundary_xy_y;
    std::vector<double> boundary_xz_x, boundary_xz_z;
    std::vector<double> boundary_SCB_xy_x, boundary_SCB_xy_y;
    std::vector<double> boundary_SCB_xz_x, boundary_SCB_xz_z;

    std::vector<std::vector<double>> boundary_xy_x_array, boundary_xy_y_array;
    std::vector<std::vector<double>> boundary_xz_x_array, boundary_xz_z_array;
    std::vector<std::vector<double>> boundary_SCB_xy_x_array, boundary_SCB_xy_y_array;
    std::vector<std::vector<double>> boundary_SCB_xz_x_array, boundary_SCB_xz_z_array;
    
    // dead regions ... 
    WCP::SMGCSelection mcells;
    std::map<WCP::SlimMergeGeomCell*, std::pair<int,int>> mcell_time_map;
    std::map<int, std::set<WCP::SlimMergeGeomCell*>> ch_mcell_set_map;
    
    // conversion between positions to the channel and time ???

    // convert time into a position
    // (time_slice - offset_t) / slope_t = position_x 
    double offset_t, slope_t;
    // convert u wire number into a position
    // (u_index -offset_u) / slope_u = position_u
    double offset_u, slope_u;
    // convert v wire number into a position
    double offset_v, slope_v;
    // convert w wire number into a position 
    double offset_w, slope_w;
    double angle_u, angle_v, angle_w;

    int dead_region_ch_ext;
    
  };
}

#endif
