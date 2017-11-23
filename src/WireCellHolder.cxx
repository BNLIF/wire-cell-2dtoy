#include "WireCell2dToy/WireCellHolder.h"

using namespace WireCell;

WireCell2dToy::WireCellHolder::WireCellHolder(){
  ncell = 0;
  nwire = 0;

  cell_no = 0;
  wire_no = 0;
}

WireCell2dToy::WireCellHolder::~WireCellHolder(){
  // for (int i=0;i!=wires.size();i++){
  //   delete wires.at(i);
  // }
  // wires.clear();
  // for (int i=0;i!=cells.size();i++){
  //   delete cells.at(i);
  // }
  // cells.clear();
}


void WireCell2dToy::WireCellHolder::AddWire(GeomWire *wire){
  nwire++;
  wires.push_back(wire);
}

void WireCell2dToy::WireCellHolder::AddCell(GeomCell *cell){
  ncell++;
  cells.push_back(cell);
}
