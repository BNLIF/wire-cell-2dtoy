#include "WireCell2dToy/TotalTiling.h"

using namespace WireCell;

WireCell2dToy::TotalTiling::TotalTiling(){
}

WireCell2dToy::TotalTiling::~TotalTiling(){
}

void WireCell2dToy::TotalTiling::AddCellWire(const GeomCell *cell, const GeomWire *uwire, const GeomWire *vwire, const GeomWire *wwire){
  cell_all.push_back(cell);
  
  GeomWireSelection wires;
  wires.push_back(uwire);
  wires.push_back(vwire);
  wires.push_back(wwire);
  cellmap[cell] = wires;

  auto it_u = find(wire_u.begin(),wire_u.end(),uwire);
  if (it_u == wire_u.end()){
    wire_u.push_back(uwire);
    wire_all.push_back(uwire);
    GeomCellSelection cells;
    cells.push_back(cell);
    wiremap[uwire] = cells; 
  }else{
    wiremap[uwire].push_back(cell);
  }

  auto it_v = find(wire_v.begin(),wire_v.end(),vwire);
  if (it_v == wire_v.end()){
    wire_v.push_back(vwire);
    wire_all.push_back(vwire);
    GeomCellSelection cells;
    cells.push_back(cell);
    wiremap[vwire] = cells; 
  }else{
    wiremap[vwire].push_back(cell);
  }
  
  auto it_w = find(wire_w.begin(),wire_w.end(),wwire);
  if (it_w == wire_w.end()){
    wire_w.push_back(wwire);
    wire_all.push_back(wwire);
    GeomCellSelection cells;
    cells.push_back(cell);
    wiremap[wwire] = cells; 
  }else{
    wiremap[wwire].push_back(cell);
  }

  
}


void WireCell2dToy::TotalTiling::Clear(){
  wire_u.clear();
  wire_v.clear();
  wire_w.clear();
  wire_all.clear();
  cell_all.clear();
  cellmap.clear();
  wiremap.clear();
}


GeomWireSelection WireCell2dToy::TotalTiling::wires(const GeomCell& cell) const
{
  if (cellmap.find(&cell) == cellmap.end()){
    //not found 
    return GeomWireSelection();
  }else{
    //found
    return cellmap.find(&cell)->second;
  }
    
}
	
GeomCellSelection WireCell2dToy::TotalTiling::cells(const GeomWire& wire) const
{
  if (wiremap.find(&wire) == wiremap.end()){
    return GeomCellSelection();
  }else{
    return wiremap.find(&wire)->second;
  }
}


const GeomCell* WireCell2dToy::TotalTiling::cell(const GeomWireSelection& wires)
{
  if (wires.size()!=3) return 0;
  const GeomWire *wire1 = wires[0];
  const GeomWire *wire2 = wires[1];
  const GeomWire *wire3 = wires[2];

  if (wire1->plane() == wire2->plane() ||
      wire1->plane() == wire3->plane() || 
      wire2->plane() == wire3->plane()) return 0;

  GeomCellSelection cells1 = cells(*wire1);
  // GeomCellSelection cells2 = cells(*wire2);
  // GeomCellSelection cells3 = cells(*wire3);
  if (cells1.size() < cells(*wire2).size())
    cells1 = cells(*wire2);
  if (cells1.size() < cells(*wire3).size())
    cells1 = cells(*wire3);
  
  for (int i = 0; i!=cells1.size(); i++){
    const GeomCell *cell1 = cells1[i];
    // for (int j =0; j!=cells2.size(); j++){
    //   const GeomCell *cell2 = cells2[j];
    //   if (*cell1==*cell2){
    auto it1 = find(cellmap[cell1].begin(),cellmap[cell1].end(),wire1);
    auto it2 = find(cellmap[cell1].begin(),cellmap[cell1].end(),wire2);
    auto it3 = find(cellmap[cell1].begin(),cellmap[cell1].end(),wire3);
    if (it1 != cellmap[cell1].end() && it2 != cellmap[cell1].end() && it3 != cellmap[cell1].end()){
      return cell1;
    }
    // for (int k=0;k!=cells3.size();k++){
    //   const GeomCell *cell3 = cells3[k];
    //   if (*cell1 == *cell3){
    //     //there is a problem here, not sure what to do 
    //     return cell1;
    //     //return 0;
    //   }
    // }
    // }
    // }
  }

  return 0;

}

