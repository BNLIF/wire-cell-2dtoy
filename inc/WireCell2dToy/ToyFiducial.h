#ifndef WIRECELL_UBOONE_TOYFIDUCIAL_H
#define WIRECELL_UBOONE_TOYFIDUCIAL_H

#include "WireCellData/Units.h"
#include "WireCellData/Point.h"
#include "WireCellData/SlimMergeGeomCell.h"
#include "WireCellData/FlashTPCBundle.h"
#include "WireCellData/ToyCTPointCloud.h"
#include "WireCellData/LMBDT.h"

#include "TVector3.h"

#include <vector>
#include <map>

namespace WireCell2dToy{
  class ToyFiducial{
  public:
    ToyFiducial(int dead_region_ch_ext = 3, double offset_t=800, double offset_u=0, double offset_v=0, double offset_w=0, double slope_t=1./2*units::mm, double slope_u=1./(3*units::mm), double slope_v=1./(3*units::mm), double slope_w=1./(3*units::mm), double angle_u=-1.047198, double angle=1.047198, double angle_w=0,
		double boundary_dis_cut=2*units::cm, double top=117*units::cm, double bottom=-116*units::cm, double upstream=0*units::cm, double downstream=1037*units::cm, double anode = 0*units::cm, double cathode=256*units::cm);
    ~ToyFiducial();

    void set_offset_t(double value){offset_t=value;};
    
    bool check_neutrino_candidate(WireCell::PR3DCluster *main_cluster, WireCell::WCPointCloud<double>::WCPoint& wcp1 ,WireCell::WCPointCloud<double>::WCPoint& wcp2, double offset_x, WireCell::ToyCTPointCloud& ct_point_cloud, bool flag_2view_check = true);
    bool inside_fiducial_volume(WireCell::Point& p, double offset_x=0);
    bool inside_dead_region(WireCell::Point& p);
    bool check_dead_volume(WireCell::Point& p, TVector3& dir, double step = 1.0*units::cm, double offset_x=0);
    
    bool check_signal_processing(WireCell::Point& p, TVector3& dir, WireCell::ToyCTPointCloud& ct_point_cloud, double step = 1.0*units::cm, double offset_x=0);

    bool check_tgm(WireCell::FlashTPCBundle *bundle, double offset_x, WireCell::ToyCTPointCloud& ct_point_cloud,std::map<WireCell::PR3DCluster*, WireCell::PR3DCluster*>& old_new_cluster_map);

    int check_LM(WireCell::FlashTPCBundle *bundle, double& cluster_length);
    int check_LM_cuts(WireCell::FlashTPCBundle *bundle, double& cluster_length);
    int check_LM_bdt(WireCell::FlashTPCBundle *bundle, double& cluster_length);
      
    void AddDeadRegion(WireCell::SlimMergeGeomCell* mcell, std::vector<int>& time_slices);
    
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
    
    // dead regions ... 
    WireCell::SMGCSelection mcells;
    std::map<WireCell::SlimMergeGeomCell*, std::pair<int,int>> mcell_time_map;
    std::map<int, std::set<WireCell::SlimMergeGeomCell*>> ch_mcell_set_map;
    
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
