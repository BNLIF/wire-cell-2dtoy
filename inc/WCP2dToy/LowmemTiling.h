#ifndef WIRECELL_LOWMEMTILING_H
#define WIRECELL_LOWMEMTILING_H

#include "WCPTiling/TilingBase.h"
#include "WCPData/GeomWCPMap.h"
#include "WCPData/MergeGeomWire.h"
#include "WCPData/SlimMergeGeomCell.h"
#include "WCPNav/SliceDataSource.h"
#include "WCPSst/GeomDataSource.h"
#include "WCPNav/DetectorGDS.h"
#include "WCP2dToy/WCPHolder.h"

#include "Rtypes.h"

namespace WCP2dToy{

  class LowmemTiling : public WCP::TilingBase {
  public:
    LowmemTiling(int time_slice, WCP::GeomDataSource& gds, WCP2dToy::WCPHolder& holder);
    LowmemTiling(int time_slice, int nrebin, WCP::GeomDataSource& gds, WCP2dToy::WCPHolder& holder);
    
    virtual ~LowmemTiling();
    
    WCP::GeomWireSelection wires(const WCP::GeomCell& cell) const;
    WCP::GeomCellSelection cells(const WCP::GeomWire& wire) const;
    const WCP::GeomCell* cell(const WCP::GeomWireSelection& wires) const;

    void re_establish_maps();
    
    int get_time_slice(){return time_slice;};
    void form_bad_merge_wires(WCP::ChirpMap& uplane_map, WCP::ChirpMap& vplane_map, WCP::ChirpMap& wplane_map);

    void form_fired_merge_wires(const WCP::Slice& slice, const WCP::Slice& slice_err);
    void form_fired_merge_wires(std::map<int,std::set<int>>& u_time_chs, std::map<int,std::set<int>>& v_time_chs, std::map<int,std::set<int>>& w_time_chs);
    void form_fired_merge_wires_with_charge(std::map<int,std::set<int>>& u_time_chs, std::map<int,std::set<int>>& v_time_chs, std::map<int,std::set<int>>& w_time_chs, std::map<std::pair<int,int>,double>& time_ch_charge_map, std::map<std::pair<int,int>,double>& time_ch_charge_err_map);
    
    void form_two_bad_cells();
    void init_bad_cells(WCP::ChirpMap& uplane_map, WCP::ChirpMap& vplane_map, WCP::ChirpMap& wplane_map);
    void check_bad_cells(LowmemTiling* tiling,WCP::ChirpMap& uplane_map, WCP::ChirpMap& vplane_map, WCP::ChirpMap& wplane_map);
    void reset_cells();
    void init_good_cells(const WCP::Slice& slice, const WCP::Slice& slice_err, std::vector<float>& uplane_rms, std::vector<float>& vplane_rms, std::vector<float>& wplane_rms);
    
    void init_good_cells(std::map<int,std::set<int>>& u_time_chs, std::map<int,std::set<int>>& v_time_chs, std::map<int,std::set<int>>& w_time_chs);
    
    void init_good_cells_with_charge(std::map<int,std::set<int>>& u_time_chs, std::map<int,std::set<int>>& v_time_chs, std::map<int,std::set<int>>& w_time_chs, std::map<std::pair<int,int>,double>& time_ch_charge_map, std::map<std::pair<int,int>,double>& time_ch_charge_err_map);
    
    WCP::GeomCellSelection& get_two_bad_wire_cells(){return two_bad_wire_cells;};
    WCP::GeomCellSelection& get_three_good_wire_cells(){return three_good_wire_cells;};
    WCP::GeomCellSelection& get_two_good_wire_cells(){return two_good_wire_cells;};
    WCP::GeomCellSelection& get_one_good_wire_cells(){return one_good_wire_cells;};
    

    WCP::GeomWireSelection& get_bad_wire_u(){return bad_wire_u;};
    WCP::GeomWireSelection& get_bad_wire_v(){return bad_wire_v;};
    WCP::GeomWireSelection& get_bad_wire_w(){return bad_wire_w;};

    WCP::SlimMergeGeomCell* create_slim_merge_cell(WCP::MergeGeomWire *uwire, WCP::MergeGeomWire *vwire, WCP::MergeGeomWire *wwire);

    bool test_point(WCP::PointVector& pcell, float dis_u, float dis_v,
		    float bmin_w, float bmax_w, float u_pitch, float v_pitch,
		    float& dis1);
    void test_cross(WCP::PointVector& pcell, float dis, WCP::WirePlaneType_t plane, float low_limit, float high_limit, float pitch, float pitch1, float w_pitch, const WCP::GeomWire* wire1, const WCP::GeomWire *wire2, int dir, float bmin_w, float bmax_w); // pitch corresponds to plane

    /* void test_crossing(WCP::PointVector& pcell, float dis_u, float dis_v, float bmin_w, float bmax_w, float u_pitch, float v_pitch, float w_pitch, const WCP::GeomWire* u_wire1, const WCP::GeomWire* u_wire2, const WCP::GeomWire *v_wire1, const WCP::GeomWire *v_wire2, int dir_u, int dir_v); */
    
    bool check_crossing(const WCP::GeomWire* wire1, const WCP::GeomWire* wire2, float pitch1, float pitch2, WCP::WirePlaneType_t plane, float min, float max, float tolerance);
    
    WCP::GeomCellSelection create_single_cells();
    WCP::GeomCellSelection create_single_cells(WCP::SlimMergeGeomCell * mcell);
    void create_one_good_wire_cells();
    void MergeWires(); // main algorithm to merge things together ... 
    
    void DivideWires(int wire_limit = 2, int min_wire = 2); // do one round of dividing for all planes
    
    WCP::GeomCellMap& get_cell_wires_map(){
      return cell_wires_map;
    }
    //map wire to cell
    WCP::GeomWireMap& get_wire_cells_map(){
      return wire_cells_map;
    }
    WCP::GeomWireSelection get_all_good_wires();
    WCP::GeomWireSelection get_all_bad_wires();

    WCP::PointVector get_all_cell_centers();

    WCP::WireChargeMap& get_wire_charge_map(){return wirechargemap;};
    WCP::WireChargeMap& get_wire_charge_error_map(){return wirecharge_errmap;};

    void Erase_Cell(WCP::SlimMergeGeomCell *cell);

    void Print_maps();
    WCP::GeomCellSelection local_deghosting(std::set<WCP::SlimMergeGeomCell*>& potential_good_mcells, std::set<WCP::SlimMergeGeomCell*>& good_mcells, bool flag_del = false);

    void local_deghosting1(std::set<WCP::SlimMergeGeomCell*>& good_mcells, std::map<WCP::SlimMergeGeomCell*, double>& map_mcell_charge);

    WCP::GeomWireSelection find_L1SP_wires();

    std::map<const WCP::GeomWire*,bool> get_wire_type_map(){return wire_type_map;};

    WCP::GeomWireWireMap get_wire_pwire_map(){ return wire_pwire_map;};
    //map parent wire to wire
    WCP::GeomWireWiresMap get_pwire_wires_map(){return pwire_wires_map;};

    bool get_regen_two_bad_wire_cells(){return regen_two_bad_wire_cells;};
    
  protected:
    WCP::GeomDataSource& gds;
    WCP2dToy::WCPHolder& holder;

    WCP::PointVector points;

    // original group of wires
    WCP::GeomWireSelection bad_wire_u; 
    WCP::GeomWireSelection bad_wire_v;
    WCP::GeomWireSelection bad_wire_w;

    WCP::GeomWireSelection fired_wire_u;
    WCP::GeomWireSelection fired_wire_v;
    WCP::GeomWireSelection fired_wire_w;

    // current cells ... 
    WCP::GeomCellSelection three_good_wire_cells;
    WCP::GeomCellSelection two_good_wire_cells;
    WCP::GeomCellSelection one_good_wire_cells;

    WCP::GeomCellSelection not_used_one_good_wire_cells;
    WCP::GeomCellSelection two_bad_wire_cells;
    bool regen_two_bad_wire_cells;

    bool remove_cell(WCP::SlimMergeGeomCell *cell);
    bool remove_wire(WCP::MergeGeomWire *wire);
    bool remove_wire_clear(WCP::MergeGeomWire *wire);
    bool replace_wire(WCP::MergeGeomWire *old_wire, WCP::MergeGeomWire *wire);

    
    int further_mergewire(WCP::GeomWireSelection &allwire);

    void calculate_merged_wire_charge();

    
    
    
    
    //map wire --> bad or good
    std::map<const WCP::GeomWire*,bool> wire_type_map;
    //map cell to wire
    WCP::GeomCellMap cell_wires_map;
    //map wire to cell
    WCP::GeomWireMap wire_cells_map;
    //map wire to parent wire
    WCP::GeomWireWireMap wire_pwire_map;
    //map parent wire to wire
    WCP::GeomWireWiresMap pwire_wires_map;
    
    

    int nrebin;
    int time_slice;

    int nwire_u;
    int nwire_v;
    int nwire_w;
    
    WCP::WireChargeMap wirechargemap;
    WCP::WireChargeMap wirecharge_errmap;
    
  };
}



#endif
