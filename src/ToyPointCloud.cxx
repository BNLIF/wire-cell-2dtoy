#include "WireCell2dToy/ToyPointCloud.h"

using namespace WireCell;

WireCell2dToy::ToyPointCloud::ToyPointCloud(){
  index = 0;
}

WireCell2dToy::ToyPointCloud::~ToyPointCloud(){
  if (index!=0)
    delete index;

  cloud.pts.clear();
}

void WireCell2dToy::ToyPointCloud::AddPoint(Point& p, SlimMergeGeomCell *mcell){
  WCPointCloud<float>::WCPoint point;
  point.x = p.x;
  point.y = p.y;
  point.z = p.z;
  point.mcell = mcell;
  cloud.pts.push_back(point);
}

void WireCell2dToy::ToyPointCloud::AddPoints(PointVector& ps, SlimMergeGeomCell *mcell){
  size_t current_size = cloud.pts.size();
  cloud.pts.resize(current_size + ps.size());
  for (size_t i=0;i!=ps.size();i++){
    cloud.pts[current_size+i].x = ps.at(i).x;
    cloud.pts[current_size+i].y = ps.at(i).y;
    cloud.pts[current_size+i].z = ps.at(i).z;
    cloud.pts[current_size+i].mcell = mcell;
  }
}

void WireCell2dToy::ToyPointCloud::build_index(){
  if (index!=(my_kd_tree_t*)0){
    delete index;
  }
  index = new my_kd_tree_t(3 /*dim*/, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index->buildIndex();
}
