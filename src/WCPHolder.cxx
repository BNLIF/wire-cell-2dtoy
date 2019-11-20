#include "WCP2dToy/WCPHolder.h"

using namespace WCP;

WCP2dToy::WCPHolder::WCPHolder(){
  ncell = 0;
  nwire = 0;

  cell_no = 0;
  wire_no = 0;
}

WCP2dToy::WCPHolder::~WCPHolder(){
  for (int i=0;i!=wires.size();i++){
    delete wires.at(i);
  }
  wires.clear();
  for (int i=0;i!=cells.size();i++){
    delete cells.at(i);
  }
  cells.clear();
}


void WCP2dToy::WCPHolder::clear_cell(){
  for (int i=0;i!=cells.size();i++){
    delete cells.at(i);
  }
  cells.clear();
}

void WCP2dToy::WCPHolder::clear_wire(){
  for (int i=0;i!=wires.size();i++){
    delete wires.at(i);
  }
  wires.clear();
}

void WCP2dToy::WCPHolder::AddWire(GeomWire *wire){
  nwire++;
  wires.push_back(wire);
}

void WCP2dToy::WCPHolder::AddCell(GeomCell *cell){
  ncell++;
  cells.push_back(cell);
}
