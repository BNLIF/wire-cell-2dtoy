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
  WCPointCloud<double>::WCPoint point;
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

void WireCell2dToy::ToyPointCloud::build_kdtree_index(){
  if (index!=(my_kd_tree_t*)0){
    delete index;
  }
  index = new my_kd_tree_t(3 /*dim*/, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index->buildIndex();
}

std::vector<std::pair<size_t,double>> WireCell2dToy::ToyPointCloud::get_closest(WireCell::Point& p, int N){
  std::vector<size_t> ret_index(N);
  std::vector<double> out_dist_sqr(N);
  
  double query_pt[3];
  query_pt[0] = p.x;
  query_pt[1] = p.y;
  query_pt[2] = p.z;
  
  N = index->knnSearch(&query_pt[0], N, &ret_index[0], &out_dist_sqr[0]);
  ret_index.resize(N);
  out_dist_sqr.resize(N);
  
  std::vector<std::pair<size_t,double>> results(N);
  for(size_t i=0; i!=N; i++){
    results.at(i) = std::make_pair(ret_index.at(i),out_dist_sqr.at(i));
  }
  
  return results;
}

std::vector<std::pair<size_t,double>> WireCell2dToy::ToyPointCloud::get_closest(WireCell::Point& p, double search_radius){
  double query_pt[3];
  query_pt[0] = p.x;
  query_pt[1] = p.y;
  query_pt[2] = p.z;
  std::vector<std::pair<size_t,double> >   ret_matches;
  nanoflann::SearchParams params;
  const size_t nMatches = index->radiusSearch(&query_pt[0], search_radius, ret_matches, params);
  
  return ret_matches;
}
