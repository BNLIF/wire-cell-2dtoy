#include "WireCell2dToy/ToyHypothesis.h"
using namespace WireCell;

WireCell2dToy::ToyHypothesis::ToyHypothesis(MergeGeomCell& mcell1, MergeGeomCell& mcell2){
  pc[0] = mcell1.center();
  pc[1] = mcell2.center();
  
  double_t max = 0, min = 0;
  // find p1[0] and p2[0] from mcell1
  GeomCellSelection cells = mcell1.get_allcell();
  for (int i=0;i!=cells.size();i++){
    const GeomCell* cell = cells.at(i);
    PointVector pv = cell->boundary();
    for (int j = 0; j!=pv.size();j++){
      double val = CalValue(pv.at(j),pc[0],pc[1]);
      if (val > max){
	max = val;
	p1[0] = pv.at(j);
      }
      if (val < min){
	min = val;
	p2[0] = pv.at(j);
      }
    }
  }
  
  max = 0; min = 0;
  //find p1[1] and p2[1] from mcell2
  cells = mcell2.get_allcell();
  for (int i=0;i!=cells.size();i++){
    const GeomCell* cell = cells.at(i);
    PointVector pv = cell->boundary();
    for (int j = 0; j!=pv.size(); j++){
      double val = CalValue(pv.at(j),pc[0],pc[1]);
      if (val > max){
	max = val;
	p1[1] = pv.at(j);
      }
      if (val < min){
	min = val;
	p2[1] = pv.at(j);
      }
    }
  }
}

WireCell2dToy::ToyHypothesis::~ToyHypothesis(){
}

double WireCell2dToy::ToyHypothesis::CalValue(Point p, Point p1, Point p2){
  double val;

  double  y = p.y, z = p.z;
  double y1 = p1.y, z1 = p1.z;
  double y2 = p2.y, z2 = p2.z;
  
  val = (z2 - z1) * (y - y1) - (y2 - y1) * (z - z1);

  return val;
}

bool WireCell2dToy::ToyHypothesis::IsInside(Point p){
  bool val = false;
  
  double dis1 = CalValue(p,p1[0],p1[1]);
  double dis2 = CalValue(p,p2[0],p2[1]);
  if (dis1* dis2 < 0 ) val = true;
  
  return val;
}

bool WireCell2dToy::ToyHypothesis::IsInside(GeomCell& cell){
  bool val = false;
  PointVector pv = cell.boundary();
  for (int i=0;i!=pv.size();i++){
    val = IsInside(pv.at(i));
    if (val) break;
  }
  return val;
}
