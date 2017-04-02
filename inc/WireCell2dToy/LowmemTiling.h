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
    LowmemTiling(int time_slice, int nrebin, const WireCell::Slice& slice,WireCell::GeomDataSource& gds, WireCell2dToy::WireCellHolder& holder,
		 std::vector<float>& uplane_rms, std::vector<float>& vplane_rms, std::vector<float>& wplane_rms);
    
    virtual ~LowmemTiling();
    
    WireCell::GeomWireSelection wires(const WireCell::GeomCell& cell) const;
    WireCell::GeomCellSelection cells(const WireCell::GeomWire& wire) const;
    const WireCell::GeomCell* cell(const WireCell::GeomWireSelection& wires) const;

    int get_time_slice(){return time_slice;};
    void form_bad_merge_wires(WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map);

    void form_fired_merge_wires(const WireCell::Slice& slice);
    
    void form_two_bad_cells();
    void init_bad_cells(WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map);
    void check_bad_cells(LowmemTiling* tiling,WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map);
    
    WireCell::GeomCellSelection& get_two_bad_wire_cells(){return two_bad_wire_cells;};
    WireCell::GeomWireSelection& get_bad_wire_u(){return bad_wire_u;};
    WireCell::GeomWireSelection& get_bad_wire_v(){return bad_wire_v;};
    WireCell::GeomWireSelection& get_bad_wire_w(){return bad_wire_w;};

  protected:
    WireCell::GeomDataSource& gds;
    WireCell2dToy::WireCellHolder& holder;

    WireCell::GeomWireSelection bad_wire_u;
    WireCell::GeomWireSelection bad_wire_v;
    WireCell::GeomWireSelection bad_wire_w;

    WireCell::GeomWireSelection fired_wire_u;
    WireCell::GeomWireSelection fired_wire_v;
    WireCell::GeomWireSelection fired_wire_w;

    WireCell::GeomCellSelection three_good_wire_cells;
    WireCell::GeomCellSelection two_good_wire_cells;
    WireCell::GeomCellSelection two_bad_wire_cells;

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
