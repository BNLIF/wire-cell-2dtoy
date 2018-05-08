#ifndef WIRECELL_UBOONE_TOYFIDUCIAL_H
#define WIRECELL_UBOONE_TOYFIDUCIAL_H

#include "WireCellData/Units.h"

namespace WireCell2dToy{
  class ToyFiducial{
  public:
    ToyFiducial(double top=117*units::cm, double bottom=-116*units::cm, double upstream=0*units::cm, double downstream=1037*units::cm, double anode = 0*units::cm, double cathode=256*units::cm);
    ~ToyFiducial();

    
    //void AddDeadRegion();
    
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
    
    // dead regions

    
  };
}

#endif
