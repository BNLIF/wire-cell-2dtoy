#ifndef WIRECELL_BADTILING_H
#define WIRECELL_BADTILING_H

#include "WCPTiling/TilingBase.h"
#include "WCPData/MergeGeomWire.h"
#include "WCPData/GeomCell.h"
#include "WCPSst/GeomDataSource.h"
#include "WCPData/GeomWCPMap.h"

namespace WCP2dToy {
  class BadTiling {
  public :
    BadTiling(int time, int scale, WCP::ChirpMap& uplane_map, 
	      WCP::ChirpMap& vplane_map, WCP::ChirpMap& wplane_map, WCP::GeomDataSource& gds, int flag_1plane = 0, int flag_all = 0);
    void BadTiling1(int time, int scale, WCP::ChirpMap& uplane_map, WCP::ChirpMap& vplane_map, WCP::ChirpMap& wplane_map, WCP::GeomDataSource& gds, int flag_all = 0);

    WCP::GeomCellMap& cmap(){return cellmap;};
    ~BadTiling();

    WCP::GeomCellSelection& get_cell_all(){return cell_all;};
    
  protected :
    WCP::GeomWireSelection wire_u;
    WCP::GeomWireSelection wire_v;
    WCP::GeomWireSelection wire_w;
    
    WCP::GeomCellMap cellmap;

    WCP::GeomCellSelection cell_all;
  };
}


#endif
