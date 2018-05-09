#include "WireCell2dToy/ToyFiducial.h"

using namespace WireCell;

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


WireCell2dToy::ToyFiducial::ToyFiducial(int dead_region_ch_ext, double offset_t, double offset_u, double offset_v, double offset_w, double slope_t, double slope_u, double slope_v, double slope_w, double angle_u, double angle_v, double angle_w, double boundary_dis_cut, double top, double bottom, double upstream, double downstream, double anode, double cathode)
  : dead_region_ch_ext(dead_region_ch_ext)
  , offset_t(offset_t)
  , offset_u(offset_u)
  , offset_v(offset_v)
  , offset_w(offset_w)
  , slope_t(slope_t)
  , slope_u(slope_u)
  , slope_v(slope_v)
  , slope_w(slope_w)
  , angle_u(angle_u)
  , angle_v(angle_v)
  , angle_w(angle_w)
  , m_top(top)
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


bool WireCell2dToy::ToyFiducial::inside_fiducial_volume(WireCell::Point &p, double offset_x){
  int c1 = pnpoly(boundary_xy_x, boundary_xy_y, p.x-offset_x, p.y);
  int c2 = pnpoly(boundary_xz_x, boundary_xz_z, p.x-offset_x, p.z);
  
  if (c1 && c2){
    return true;
  }else{
    return false;
  }
}


void WireCell2dToy::ToyFiducial::AddDeadRegion(WireCell::SlimMergeGeomCell* mcell, std::vector<int>& time_slices){

  mcells.push_back(mcell);
  int start_time = time_slices.front() - dead_region_ch_ext ;
  int end_time = time_slices.back() + dead_region_ch_ext;
  mcell_time_map[mcell] = std::make_pair(start_time, end_time);

  GeomWireSelection& uwires = mcell->get_uwires();
  GeomWireSelection& vwires = mcell->get_vwires();
  GeomWireSelection& wwires = mcell->get_wwires();
  
  std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();

  if (find(bad_planes.begin(), bad_planes.end(), WirePlaneType_t(0))!=bad_planes.end()){
    int start_ch = uwires.front()->channel() - dead_region_ch_ext;
    if (start_ch <0) start_ch = 0;
    int end_ch = uwires.back()->channel() + dead_region_ch_ext;
    if (end_ch>=2400) end_ch = 2399;
    for (int i = start_ch; i<=end_ch;i++){
      if (ch_mcell_set_map.find(i)==ch_mcell_set_map.end()){
	std::set<SlimMergeGeomCell*> mcells_set;
	mcells_set.insert(mcell);
	ch_mcell_set_map[i] = mcells_set;
      }else{
	ch_mcell_set_map[i].insert(mcell);
      }
    }
  }
  if (find(bad_planes.begin(), bad_planes.end(), WirePlaneType_t(1))!=bad_planes.end()){
    int start_ch = vwires.front()->channel() - dead_region_ch_ext;
    if (start_ch <2400) start_ch = 2400;
    int end_ch = vwires.back()->channel() + dead_region_ch_ext;
    if (end_ch>=4800) end_ch = 4799;
    for (int i = start_ch; i<=end_ch;i++){
      if (ch_mcell_set_map.find(i)==ch_mcell_set_map.end()){
	std::set<SlimMergeGeomCell*> mcells_set;
	mcells_set.insert(mcell);
	ch_mcell_set_map[i] = mcells_set;
      }else{
	ch_mcell_set_map[i].insert(mcell);
      }
    }
  }
  if (find(bad_planes.begin(), bad_planes.end(), WirePlaneType_t(2))!=bad_planes.end()){
    int start_ch = wwires.front()->channel() - dead_region_ch_ext;
    if (start_ch <4800) start_ch = 4800;
    int end_ch = wwires.back()->channel() + dead_region_ch_ext;
    if (end_ch>=7256) end_ch = 7255;
    for (int i = start_ch; i<=end_ch;i++){
      if (ch_mcell_set_map.find(i)==ch_mcell_set_map.end()){
	std::set<SlimMergeGeomCell*> mcells_set;
	mcells_set.insert(mcell);
	ch_mcell_set_map[i] = mcells_set;
      }else{
	ch_mcell_set_map[i].insert(mcell);
      }
    }
  }
  
  
}
