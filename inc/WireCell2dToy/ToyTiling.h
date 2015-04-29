#ifndef WIRECELL_TOYTILING_H
#define WIRECELL_TOYTILING_H

#include "WireCellTiling/TilingBase.h"
#include "WireCellData/GeomWireCellMap.h"
#include "WireCellNav/SliceDataSource.h"
#include "WireCellSst/GeomDataSource.h"

namespace WireCell2dToy {
    /**
     *  A bogus tiling class that doesn't do anything.
     */
    class ToyTiling : public WireCell::TilingBase { 
    public:
      ToyTiling();
	// fixme: tiling should not know about slices data.
      ToyTiling(WireCell::Slice slice,WireCellSst::GeomDataSource gds);
      virtual ~ToyTiling();
      
      WireCell::GeomWireSelection wires(const WireCell::GeomCell& cell) const;
      WireCell::GeomCellSelection cells(const WireCell::GeomWire& wire) const;
      const WireCell::GeomCell* cell(const WireCell::GeomWireSelection& wires) const;
      
      WireCell::GeomCellSelection get_allcell(){ return cell_all;}
      WireCell::GeomWireSelection get_allwire(){ return wire_all;}
      int get_ncell(){return ncell;}

      WireCell::GeomWireSelection get_wire_u(){ return wire_u;}
      WireCell::GeomWireSelection get_wire_v(){ return wire_v;}
      WireCell::GeomWireSelection get_wire_w(){return wire_w;}
      WireCell::GeomWireSelection get_wire_all(){ return wire_all;}

    protected:
      WireCell::GeomWireSelection wire_u;
      WireCell::GeomWireSelection wire_v;
      WireCell::GeomWireSelection wire_w;
      WireCell::GeomWireSelection wire_all;
      int ncell;

      WireCell::GeomCellSelection cell_all;
      
      WireCell::GeomCellMap cellmap;
      WireCell::GeomWireMap wiremap;
      
      

      //

      
      
    };
}
#endif
