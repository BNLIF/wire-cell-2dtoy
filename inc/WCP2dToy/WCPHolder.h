#ifndef WCP2dToy_WCPHolder_h
#define WCP2dToy_WCPHolder_h

#include "WCPData/GeomCell.h"
#include "WCPData/GeomWCPMap.h"

namespace WCP2dToy{
  class WCPHolder{
  public:
    WCPHolder();
    ~WCPHolder();

    void AddWire(WCP::GeomWire *wire);
    void AddCell(WCP::GeomCell *cell);
    
    int get_ncell(){return ncell;};
    int get_nwire(){return nwire;};

    void AddWire_No(){cell_no++;};
    void AddCell_No(){wire_no++;};
    
    int get_cell_no(){return cell_no;};
    int get_wire_no(){return wire_no;};

    void clear_cell();
    void clear_wire();

    WCP::GeomCellSelection& get_cells(){return cells;};
    
  protected:
    int ncell;
    int nwire;

    int cell_no;
    int wire_no;

    WCP::GeomWireSelection wires;
    WCP::GeomCellSelection cells;
  };
}

#endif
