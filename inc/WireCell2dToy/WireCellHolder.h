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

    void AddWire_No(){cell_no++;};
    void AddCell_No(){wire_no++;};
    
    int get_cell_no(){return cell_no;};
    int get_wire_no(){return wire_no;};

  protected:
    int ncell;
    int nwire;

    int cell_no;
    int wire_no;

    WireCell::GeomWireSelection wires;
    WireCell::GeomCellSelection cells;
  };
}

#endif
