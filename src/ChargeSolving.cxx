#include "WireCell2dToy/ChargeSolving.h"


using namespace WireCell;

WireCell2dToy::ChargeSolving::ChargeSolving(const WireCell::GeomDataSource& gds, LowmemTiling& tiling)
  : gds(gds)
  , tiling(tiling)
{
  GeomCellSelection& one_wire_cells = tiling.get_one_good_wire_cells();
  GeomCellSelection& two_wire_cells = tiling.get_two_good_wire_cells();
  GeomCellSelection& three_wire_cells = tiling.get_three_good_wire_cells();

  GeomWireSelection good_wires = tiling.get_all_good_wires();
  GeomWireSelection bad_wires = tiling.get_all_bad_wires();

  GeomCellMap cell_wire_map = tiling.get_cell_wires_map();
  GeomWireMap wire_cell_map = tiling.get_wire_cells_map();

  std::cout << one_wire_cells.size() << " " << two_wire_cells.size() << " " << three_wire_cells.size() << " " << good_wires.size() << " " << bad_wires.size() << " " << cell_wire_map.size() << " " << wire_cell_map.size() << std::endl;

  // divide all wires into seperate groups ...

  // start one good wire, find all cells move in, all the good wires move in
}

WireCell2dToy::ChargeSolving::~ChargeSolving(){
}
