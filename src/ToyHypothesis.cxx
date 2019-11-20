#include "WCP2dToy/ToyHypothesis.h"
using namespace WCP;
WCP2dToy::ToyHypothesis::ToyHypothesis(){
  pc[0] = Point();
  pc[1] = Point();
  p1[0] = Point();
  p1[1] = Point();
  p2[0] = Point();
  p2[1] = Point();
  
}

WCP2dToy::ToyHypothesis::ToyHypothesis(MergeGeomCell& mcell1, MergeGeomCell& mcell2){
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
      
      // std::cout << i << " " << j << " " << val << std::endl;

      if (val > max){
  	max = val;
  	p1[0] = pv.at(j);
      }
      if (val < min){
  	min = val;
  	p2[0] = pv.at(j);
      }
    }

    // //try center
    // double val = CalValue(cell->center(),pc[0],pc[1]);
    // if (val > max){
    //   max = val;
    //   p1[0] = cell->center();
    // }
    // if (val < min){
    //   min = val;
    //   p2[0] = cell->center();
    // }
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

    //try center
    // double val = CalValue(cell->center(),pc[0],pc[1]);
    // if (val > max){
    //   max = val;
    //   p1[0] = cell->center();
    // }
    // if (val < min){
    //   min = val;
    //   p2[0] = cell->center();
    // }


  }

  // for (int i=0;i!=mcell1.get_allcell().size();i++){
  //   std::cout << mcell1.get_allcell().at(i)->center().y << " " << mcell1.get_allcell().at(i)->center().z << std::endl;
  // }

  // for (int i=0;i!=mcell2.get_allcell().size();i++){
  //   std::cout << mcell2.get_allcell().at(i)->center().y << " " << mcell2.get_allcell().at(i)->center().z << std::endl;
  // }

  // std::cout << pc[0].y << " " << pc[0].z << " " 
  // 	    << pc[1].y << " " << pc[1].z << " " 
  // 	    << p1[0].y << " " << p1[0].z << " " 
  // 	    << p1[1].y << " " << p1[1].z << " " 
  // 	    << p2[0].y << " " << p2[0].z << " " 
  // 	    << p2[1].y << " " << p2[1].z << " "
  // 	    << std::endl;
 
}

WCP2dToy::ToyHypothesis::~ToyHypothesis(){
}

double WCP2dToy::ToyHypothesis::CalValue(Point p, Point p1, Point p2){
  double val;

  double  y = p.y, z = p.z;
  double y1 = p1.y, z1 = p1.z;
  double y2 = p2.y, z2 = p2.z;

  // std::cout << y << " " << z << " " 
  // 	    << y1 << " " << z1 << " " 
  // 	    << y2 << " " << z2 << std::endl;
  
  val = (z2 - z1) * (y - y1) - (y2 - y1) * (z - z1);

  return val;
}

bool WCP2dToy::ToyHypothesis::IsInside(Point p){
  bool val = false;
  
  double dis1 = CalValue(p,p1[0],p1[1]);
  double dis2 = CalValue(p,p2[0],p2[1]);
  if ( dis1* dis2 < 0 ) val = true;
  
  return val;
}

bool WCP2dToy::ToyHypothesis::IsInside(const GeomCell& cell){
  bool val = false;
  PointVector pv = cell.boundary();
  for (int i=0;i!=pv.size();i++){
    val = IsInside(pv.at(i));
    if (val) break;
  }
  
  // //try center 
  // val = IsInside(cell.center());
  

  return val;
}
