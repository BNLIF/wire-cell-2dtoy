#include "WireCell2dToy/ToyFiducial.h"

int pnpoly(std::vector<double>& vertx, std::vector<double>& verty, double testx, double testy)
{
  int i, j, c = 0;
  for (i = 0, j = int(vertx.size())-1; i < int(vertx.size()); j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
      c = !c;
  }
  return c;
}


WireCell2dToy::ToyFiducial::ToyFiducial(double boundary_dis_cut, double top, double bottom, double upstream, double downstream, double anode, double cathode)
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

  //
  boundary_xy_x.clear(); boundary_xy_y.clear();
  boundary_xy_x.push_back(m_anode + boundary_dis_cut); boundary_xy_y.push_back(m_bottom + boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_bottom_1_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_bottom_1_y+boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_bottom_2_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_bottom_2_y+boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_top_2_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_top_2_y-boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_top_1_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_top_1_y-boundary_dis_cut);
  boundary_xy_x.push_back(m_anode + boundary_dis_cut); boundary_xy_y.push_back(m_top - boundary_dis_cut);

  
  boundary_xz_x.clear(); boundary_xz_z.clear();
  boundary_xz_x.push_back(m_anode + boundary_dis_cut); boundary_xz_z.push_back(m_upstream + boundary_dis_cut);
  boundary_xz_x.push_back(m_sc_upstream_1_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_upstream_1_z + boundary_dis_cut);
  boundary_xz_x.push_back(m_sc_upstream_2_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_upstream_2_z + boundary_dis_cut);
  boundary_xz_x.push_back(m_sc_downstream_2_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_downstream_2_z - boundary_dis_cut);
  boundary_xz_x.push_back(m_sc_downstream_1_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_downstream_1_z - boundary_dis_cut);
  boundary_xz_x.push_back(m_anode + boundary_dis_cut); boundary_xz_z.push_back(m_downstream - boundary_dis_cut);
  
}

WireCell2dToy::ToyFiducial::~ToyFiducial(){
}


bool WireCell2dToy::ToyFiducial::inside_fiducial_volume(WireCell::Point &p){
  int c1 = pnpoly(boundary_xy_x, boundary_xy_y, p.x, p.y);
  int c2 = pnpoly(boundary_xz_x, boundary_xz_z, p.x, p.z);
  
  if (c1 && c2){
    return true;
  }else{
    return false;
  }
}
