#ifndef WIRECELL_TOTALTILING_H
#define WIRECELL_TOTALTILING_H
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/GeomWireCellMap.h"

namespace WireCell2dToy {
  class TotalTiling {
    public:
    TotalTiling();
    ~TotalTiling();
    
    void AddCellWire(const WireCell::GeomCell *cell, const WireCell::GeomWire *uwire, const WireCell::GeomWire *vwire, const WireCell::GeomWire *wwire);
    
    WireCell::GeomWireSelection wires(const WireCell::GeomCell& cell) const;
    WireCell::GeomCellSelection cells(const WireCell::GeomWire& wire) const;
    const WireCell::GeomCell* cell(const WireCell::GeomWireSelection& wires);
    void Clear();
  
    WireCell::GeomWireSelection wire_u;
    WireCell::GeomWireSelection wire_v;
    WireCell::GeomWireSelection wire_w;

    WireCell::GeomWireSelection wire_all;
    
    WireCell::GeomCellSelection cell_all;

    WireCell::GeomCellMap cellmap;
    WireCell::GeomWireMap wiremap;
  };
  

}

#endif
