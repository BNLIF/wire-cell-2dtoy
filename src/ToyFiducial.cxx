#include "WireCell2dToy/ToyFiducial.h"



WireCell2dToy::ToyFiducial::ToyFiducial(double top, double bottom, double upstream, double downstream, double anode, double cathode)
  : m_top(top)
  , m_bottom(bottom)
  , m_upstream(upstream)
  , m_downstream(downstream)
  , m_anode(anode)
  , m_cathode(cathode)
{
  m_sc_bottom_1_y=-116*units::cm;
  m_sc_bottom_1_x=80*units::cm;

  m_sc_bottom_2_y=-99*units::cm;
  m_sc_bottom_2_x=256*units::cm;

  m_sc_top_1_y = 118*units::cm;
  m_sc_top_1_x = 100*units::cm;

  m_sc_top_2_y = 103*units::cm;
  m_sc_top_2_x = 256*units::cm;

  m_sc_upstream_1_z = 0*units::cm;
  m_sc_upstream_1_x = 120*units::cm;

  m_sc_upstream_2_z = 11*units::cm;
  m_sc_upstream_2_x = 256*units::cm;

  m_sc_downstream_1_z=1037*units::cm;
  m_sc_downstream_1_x=120*units::cm;

  m_sc_downstream_2_z=1026*units::cm;
  m_sc_downstream_2_x=256*units::cm;
}

WireCell2dToy::ToyFiducial::~ToyFiducial(){
}
