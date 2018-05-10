#ifndef WIRECELL_UBOONE_TOYFIDUCIAL_H
#define WIRECELL_UBOONE_TOYFIDUCIAL_H

#include "WireCellData/Units.h"
#include "WireCellData/Point.h"
#include "WireCellData/SlimMergeGeomCell.h"

#include <vector>
#include <map>

namespace WireCell2dToy{
  class ToyFiducial{
  public:
    ToyFiducial(int dead_region_ch_ext = 3, double offset_t=800, double offset_u=0, double offset_v=0, double offset_w=0, double slope_t=1./2*units::mm, double slope_u=1./3*units::mm, double slope_v=1./3*units::mm, double slope_w=1./3*units::mm, double angle_u=-1.047198, double angle=1.047198, double angle_w=0,
		double boundary_dis_cut=3*units::cm, double top=117*units::cm, double bottom=-116*units::cm, double upstream=0*units::cm, double downstream=1037*units::cm, double anode = 0*units::cm, double cathode=256*units::cm);
    ~ToyFiducial();

    bool inside_fiducial_volume(WireCell::Point& p, double offset_x=0);
    bool inside_dead_region(WireCell::Point& p);
      
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
