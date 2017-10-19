#ifndef WireCell2dToy_TOYPOINTCLOUD_H
#define WireCell2dToy_TOYPOINTCLOUD_H


#include "WireCellNanoflann/nanoflann.h" 
#include "WireCellData/WCPointCloud.h"
#include "WireCellQuickhull/QuickHull.h"
#include "WireCellQuickhull/MathUtils.h"

#include  <map>

typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, WireCell::WCPointCloud<double> > ,
    WireCell::WCPointCloud<double>,
    3 /* dim */
  > my_kd_tree_t;

namespace WireCell2dToy{
  
  class ToyPointCloud {
  public:
    ToyPointCloud();
    ~ToyPointCloud();

    void AddPoint(WireCell::Point& p, WireCell::SlimMergeGeomCell *mcell);
    void AddPoints(WireCell::PointVector& ps, WireCell::SlimMergeGeomCell *mcell);
    void build_kdtree_index();
    std::vector<std::pair<size_t,double>> get_closest_index(WireCell::Point& p, int N);
    std::vector<std::pair<size_t,double>> get_closest_index(WireCell::Point& p, double radius);

    std::map<WireCell::SlimMergeGeomCell*, WireCell::Point> get_closest_mcell(WireCell::Point& p, int N);

    std::map<WireCell::SlimMergeGeomCell*, WireCell::Point> get_closest_mcell(WireCell::Point& p, double radius);
    
    std::vector<std::pair<WireCell::SlimMergeGeomCell*, WireCell::Point>> get_hull(); 
    
    
    
  protected:
    WireCell::WCPointCloud<double> cloud;
    my_kd_tree_t *index;
  };

   
}

#endif
