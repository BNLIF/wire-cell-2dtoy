#ifndef WIRECELL_BADTILING_H
#define WIRECELL_BADTILING_H

#include "WireCellTiling/TilingBase.h"
#include "WireCellData/MergeGeomWire.h"
#include "WireCellData/GeomCell.h"
#include "WireCellSst/GeomDataSource.h"

namespace WireCell2dToy {
  class BadTiling {
  public :
    BadTiling(int time, int scale, WireCell::ChirpMap& uplane_map, 
	      WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map, WireCell::GeomDataSource& gds);

    ~BadTiling();

    WireCell::GeomCellSelection& get_cell_all(){return cell_all;};
    
  protected :
    WireCell::GeomWireSelection wire_u;
    WireCell::GeomWireSelection wire_v;
    WireCell::GeomWireSelection wire_w;
    
    
    WireCell::GeomCellSelection cell_all;
  };
}


#endif
