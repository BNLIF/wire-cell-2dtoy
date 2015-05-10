#include "WireCell2dToy/ToyTiling.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include <cmath>

using namespace WireCell;

WireCell2dToy::ToyTiling::ToyTiling()
{
}

WireCell2dToy::ToyTiling::ToyTiling(WireCell::Slice slice,WireCellSst::GeomDataSource gds){
  WireCell::Channel::Group group = slice.group();
  for (int i=0;i!=group.size();i++){
    const WireCell::GeomWire *wire = gds.by_channel_segment(group.at(i).first,0);
    float charge = group.at(i).second;
    
    if (wirechargemap.find(wire) == wirechargemap.end()){
      //not found
      wirechargemap[wire] = charge;
    }else{
      wirechargemap[wire] += charge;
    }


    wire_all.push_back(wire);
    if (wire->plane() == kUwire){
      wire_u.push_back(wire);
    }else if (wire->plane() == kVwire){
      wire_v.push_back(wire);
    }else if (wire->plane() == kYwire){
      wire_w.push_back(wire);
    }
  }

  float dis_u[2],dis_v[2],dis_w[2],
    dis_puv[4],dis_puw[4],dis_pwv[4];
  ncell = 1;

  //std::cout << "Wire Counts: " << wire_u.size() << " " << wire_v.size() << " " << wire_w.size() << std::endl;


  for (int i=0;i!=wire_u.size();i++){
    dis_u[0] = gds.wire_dist(*wire_u[i]) - gds.pitch(kUwire)/2.;
    dis_u[1] = gds.wire_dist(*wire_u[i]) + gds.pitch(kUwire)/2.;
    for (int j=0;j!=wire_v.size();j++){
      dis_v[0] = gds.wire_dist(*wire_v[j]) - gds.pitch(kVwire)/2.;
      dis_v[1] = gds.wire_dist(*wire_v[j]) + gds.pitch(kVwire)/2.;
      
      //four vertices around
      PointVector puv(4);
      
      gds.crossing_point(dis_u[0],dis_v[0],kUwire,kVwire, puv[0]);
      gds.crossing_point(dis_u[0],dis_v[1],kUwire,kVwire, puv[1]);
      gds.crossing_point(dis_u[1],dis_v[1],kUwire,kVwire, puv[2]);
      gds.crossing_point(dis_u[1],dis_v[0],kUwire,kVwire, puv[3]);
      
      for (int k=0;k!=4;k++){
	dis_puv[k] = gds.wire_dist(puv[k],kYwire);
      }
      
      for (int k=0;k!=wire_w.size();k++){
	int flag = 0;
	PointVector pcell;
	dis_w[0] = gds.wire_dist(*wire_w[k]) - gds.pitch(kYwire)/2.;
  	dis_w[1] = gds.wire_dist(*wire_w[k]) + gds.pitch(kYwire)/2.;	
	
	for (int m = 0;m!=4;m++){
	  if (dis_puv[m] > dis_w[0] && dis_puv[m] < dis_w[1]){
	    flag = 1;
	    pcell.push_back(puv[m]);
	  }
	}
	
	if (flag==1 ) {
	  PointVector puw(4);
	  gds.crossing_point(dis_u[0],dis_w[0],kUwire,kYwire, puw[0]);
	  gds.crossing_point(dis_u[0],dis_w[1],kUwire,kYwire, puw[1]);
	  gds.crossing_point(dis_u[1],dis_w[1],kUwire,kYwire, puw[2]);
	  gds.crossing_point(dis_u[1],dis_w[0],kUwire,kYwire, puw[3]);

	  for (int k=0;k!=4;k++){
	    dis_puw[k] = gds.wire_dist(puw[k],kVwire);
	    if (dis_puw[k] > dis_v[0] && dis_puw[k] < dis_v[1]){
	      pcell.push_back(puw[k]);
	    }
	  }
	  PointVector pwv(4);
	  gds.crossing_point(dis_v[0],dis_w[0],kVwire,kYwire, pwv[0]);
	  gds.crossing_point(dis_v[0],dis_w[1],kVwire,kYwire, pwv[1]);
	  gds.crossing_point(dis_v[1],dis_w[1],kVwire,kYwire, pwv[2]);
	  gds.crossing_point(dis_v[1],dis_w[0],kVwire,kYwire, pwv[3]);
	  for (int k=0;k!=4;k++){
	    dis_pwv[k] = gds.wire_dist(pwv[k],kUwire);
	    if (dis_pwv[k] > dis_u[0] && dis_pwv[k] < dis_u[1]){
	      pcell.push_back(pwv[k]);
	    }
	  }

	  
	  
	  //order all the points by phi angle
	  GeomCell *cell = new GeomCell(ncell,pcell);
	  
	  //std::cout << i << " " << j << " " << k << " " << pcell.size() << " " << cell->center().z/units::m << " " << cell->center().y/units::m << std::endl;

	  // fill cellmap
	  GeomWireSelection wiresel;
	  wiresel.push_back(wire_u[i]);
	  wiresel.push_back(wire_v[j]);
	  wiresel.push_back(wire_w[k]);	
	  cellmap[cell]=wiresel;
	  
	  //fillwiremap
	  
	  if (wiremap.find(wire_u[i]) == wiremap.end()){
	    //not found
	    GeomCellSelection cellsel;
	    cellsel.push_back(cell);
	    wiremap[wire_u[i]]=cellsel;
	  }else{
	    //found
	    wiremap[wire_u[i]].push_back(cell);
	  }
	  
	  if (wiremap.find(wire_v[j]) == wiremap.end()){
	    //not found
	    GeomCellSelection cellsel;
	    cellsel.push_back(cell);
	    wiremap[wire_v[j]]=cellsel;
	  }else{
	    //found
	    wiremap[wire_v[j]].push_back(cell);
	  }
	  
	  if (wiremap.find(wire_w[k]) == wiremap.end()){
	    //not found
	    GeomCellSelection cellsel;
	    cellsel.push_back(cell);
	    wiremap[wire_w[k]]=cellsel;
	  }else{
	    //found
	    wiremap[wire_w[k]].push_back(cell);
	  }
	  
	  cell_all.push_back(cell);
	  ncell++;
	}
      }
      
      // initialize uw and vw points
      // if (flag==1){
	// check order 
	//	pcell = cell->boundary();
	// std::cout << "Cell Count: " << pcell.size() << " " << cell->cross_section() << std::endl;
	// for (int k=0;k!=pcell.size();k++){
	//   std::cout << pcell[k].y << " " << pcell[k].z << " " << std::atan2(pcell[k].z - cell->center().z, pcell[k].y-cell->center().y) << std::endl;
	// }
      // }
    }
  }
  
  //std::cout << wire_u.size() << " " << wire_v.size() << " " << wire_w.size() << std::endl;
  
}


WireCell2dToy::ToyTiling::~ToyTiling()
{
  //delete all the cells
  for (int i=0;i!=cell_all.size();i++){
    cell_all[i] = 0;
  }
}

GeomWireSelection WireCell2dToy::ToyTiling::wires(const GeomCell& cell) const
{
  if (cellmap.find(&cell) == cellmap.end()){
    //not found 
    return GeomWireSelection();
  }else{
    //found
    return cellmap.find(&cell)->second;
  }
    
}
	
GeomCellSelection WireCell2dToy::ToyTiling::cells(const GeomWire& wire) const
{
  if (wiremap.find(&wire) == wiremap.end()){
    return GeomCellSelection();
  }else{
    return wiremap.find(&wire)->second;
  }
}


const GeomCell* WireCell2dToy::ToyTiling::cell(const GeomWireSelection& wires) const
{
  if (wires.size()!=3) return 0;
  const GeomWire *wire1 = wires[0];
  const GeomWire *wire2 = wires[1];
  const GeomWire *wire3 = wires[2];

  if (wire1->plane() == wire2->plane() ||
      wire1->plane() == wire3->plane() || 
      wire2->plane() == wire3->plane()) return 0;

  GeomCellSelection cells1 = cells(*wire1);
  GeomCellSelection cells2 = cells(*wire2);
  GeomCellSelection cells3 = cells(*wire3);
  
  for (int i = 0; i!=cells1.size(); i++){
    const GeomCell *cell1 = cells1[i];
    for (int j =0; j!=cells2.size(); j++){
      const GeomCell *cell2 = cells2[j];
      if (*cell1==*cell2){
	for (int k=0;k!=cells3.size();k++){
	  const GeomCell *cell3 = cells3[k];
	  if (*cell1 == *cell3){
	    //there is a problem here, not sure what to do 
	    return cell1;
	    //return 0;
	  }
	}
      }
    }
  }

  return 0;

}
