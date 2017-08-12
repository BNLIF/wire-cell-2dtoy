#ifndef WIRECELL_LOWMEMTILING_H
#define WIRECELL_LOWMEMTILING_H

#include "WireCellTiling/TilingBase.h"
#include "WireCellData/GeomWireCellMap.h"
#include "WireCellData/MergeGeomWire.h"
#include "WireCellData/SlimMergeGeomCell.h"
#include "WireCellNav/SliceDataSource.h"
#include "WireCellSst/GeomDataSource.h"
#include "WireCellNav/DetectorGDS.h"
#include "WireCell2dToy/WireCellHolder.h"

#include "Rtypes.h"

namespace WireCell2dToy{

  class LowmemTiling : public WireCell::TilingBase {
  public:
    LowmemTiling(int time_slice, int nrebin,WireCell::GeomDataSource& gds, WireCell2dToy::WireCellHolder& holder);
    
    virtual ~LowmemTiling();
    
    WireCell::GeomWireSelection wires(const WireCell::GeomCell& cell) const;
    WireCell::GeomCellSelection cells(const WireCell::GeomWire& wire) const;
    const WireCell::GeomCell* cell(const WireCell::GeomWireSelection& wires) const;

    int get_time_slice(){return time_slice;};
    void form_bad_merge_wires(WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map);

    void form_fired_merge_wires(const WireCell::Slice& slice, const WireCell::Slice& slice_err);
    
    void form_two_bad_cells();
    void init_bad_cells(WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map);
    void check_bad_cells(LowmemTiling* tiling,WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map);
    void init_good_cells(const WireCell::Slice& slice, const WireCell::Slice& slice_err, std::vector<float>& uplane_rms, std::vector<float>& vplane_rms, std::vector<float>& wplane_rms);
    
    WireCell::GeomCellSelection& get_two_bad_wire_cells(){return two_bad_wire_cells;};
    WireCell::GeomCellSelection& get_three_good_wire_cells(){return three_good_wire_cells;};
    

    WireCell::GeomWireSelection& get_bad_wire_u(){return bad_wire_u;};
    WireCell::GeomWireSelection& get_bad_wire_v(){return bad_wire_v;};
    WireCell::GeomWireSelection& get_bad_wire_w(){return bad_wire_w;};

    WireCell::SlimMergeGeomCell* create_slim_merge_cell(WireCell::MergeGeomWire *uwire, WireCell::MergeGeomWire *vwire, WireCell::MergeGeomWire *wwire);

    bool test_point(WireCell::PointVector& pcell, float dis_u, float dis_v,
		    float bmin_w, float bmax_w, float u_pitch, float v_pitch,
		    float& dis1);
    void test_cross(WireCell::PointVector& pcell, float dis, WireCell::WirePlaneType_t plane, float low_limit, float high_limit, float pitch, float pitch1, float w_pitch, const WireCell::GeomWire* wire1, const WireCell::GeomWire *wire2, int dir, float bmin_w, float bmax_w); // pitch corresponds to plane

    /* void test_crossing(WireCell::PointVector& pcell, float dis_u, float dis_v, float bmin_w, float bmax_w, float u_pitch, float v_pitch, float w_pitch, const WireCell::GeomWire* u_wire1, const WireCell::GeomWire* u_wire2, const WireCell::GeomWire *v_wire1, const WireCell::GeomWire *v_wire2, int dir_u, int dir_v); */
    
    bool check_crossing(const WireCell::GeomWire* wire1, const WireCell::GeomWire* wire2, float pitch1, float pitch2, WireCell::WirePlaneType_t plane, float min, float max, float tolerance);
    
    WireCell::GeomCellSelection create_single_cells();
    WireCell::GeomCellSelection create_single_cells(WireCell::SlimMergeGeomCell * mcell);
    void create_one_good_wire_cells();
    void MergeWires(); // main algorithm to merge things together ... 
    
    void DivideWires(int wire_limit = 2, int min_wire = 2); // do one round of dividing for all planes
    
    WireCell::GeomCellMap get_cell_wires_map(){
      return cell_wires_map;
    }
    //map wire to cell
    WireCell::GeomWireMap get_wire_cells_map(){
      return wire_cells_map;
    }
    WireCell::GeomWireSelection get_all_good_wires();
    WireCell::GeomWireSelection get_all_bad_wires();

    WireCell::PointVector get_all_cell_centers();

  protected:
    WireCell::GeomDataSource& gds;
    WireCell2dToy::WireCellHolder& holder;

    WireCell::PointVector points;

    // original group of wires
    WireCell::GeomWireSelection bad_wire_u; 
    WireCell::GeomWireSelection bad_wire_v;
    WireCell::GeomWireSelection bad_wire_w;

    WireCell::GeomWireSelection fired_wire_u;
    WireCell::GeomWireSelection fired_wire_v;
    WireCell::GeomWireSelection fired_wire_w;

    // current cells ... 
    WireCell::GeomCellSelection three_good_wire_cells;
    WireCell::GeomCellSelection two_good_wire_cells;
    WireCell::GeomCellSelection one_good_wire_cells;

    WireCell::GeomCellSelection not_used_one_good_wire_cells;
    WireCell::GeomCellSelection two_bad_wire_cells;

    bool remove_cell(WireCell::SlimMergeGeomCell *cell);
    bool remove_wire(WireCell::MergeGeomWire *wire);
    bool replace_wire(WireCell::MergeGeomWire *old_wire, WireCell::MergeGeomWire *wire);

    
    int further_mergewire(WireCell::GeomWireSelection &allwire);

   

    
    //map wire --> bad or good
    std::map<const WireCell::GeomWire*,bool> wire_type_map;
    //map cell to wire
    WireCell::GeomCellMap cell_wires_map;
    //map wire to cell
    WireCell::GeomWireMap wire_cells_map;
    //map wire to parent wire
    WireCell::GeomWireWireMap wire_pwire_map;
    //map parent wire to wire
    WireCell::GeomWireWiresMap pwire_wires_map;
    
    

    int nrebin;
    int time_slice;

    int nwire_u;
    int nwire_v;
    int nwire_w;
    
    WireCell::WireChargeMap wirechargemap;
    WireCell::WireChargeMap wirecharge_errmap;
    
  };
}



#endif
