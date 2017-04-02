#ifndef WireCell2dToy_WireCellHolder_h
#define WireCell2dToy_WireCellHolder_h

#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWireCellMap.h"

namespace WireCell2dToy{
  class WireCellHolder{
  public:
    WireCellHolder();
    ~WireCellHolder();

    void AddWire(WireCell::GeomWire *wire);
    void AddCell(WireCell::GeomCell *cell);
    
    int get_ncell(){return ncell;};
    int get_nwire(){return nwire;};

  protected:
    int ncell;
    int nwire;

    WireCell::GeomWireSelection wires;
    WireCell::GeomCellSelection cells;
  };
}

#endif
