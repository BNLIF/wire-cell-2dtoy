#ifndef WIRECELL_TOTALTILING_H
#define WIRECELL_TOTALTILING_H
#include "WCPData/GeomCell.h"
#include "WCPData/GeomWire.h"
#include "WCPData/GeomWCPMap.h"

namespace WCP2dToy {
  class TotalTiling {
    public:
    TotalTiling();
    ~TotalTiling();
    
    void AddCellWire(const WCP::GeomCell *cell, const WCP::GeomWire *uwire, const WCP::GeomWire *vwire, const WCP::GeomWire *wwire);
    
    WCP::GeomWireSelection wires(const WCP::GeomCell& cell) const;
    WCP::GeomCellSelection cells(const WCP::GeomWire& wire) const;
    const WCP::GeomCell* cell(const WCP::GeomWireSelection& wires);
    void Clear();
  
    WCP::GeomWireSelection wire_u;
    WCP::GeomWireSelection wire_v;
    WCP::GeomWireSelection wire_w;

    WCP::GeomWireSelection wire_all;
    
    WCP::GeomCellSelection cell_all;

    WCP::GeomCellMap cellmap;
    WCP::GeomWireMap wiremap;
  };
  

}

#endif
