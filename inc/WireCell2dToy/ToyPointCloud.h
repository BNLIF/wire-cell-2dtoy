#ifndef WireCell2dToy_TOYPOINTCLOUD_H
#define WireCell2dToy_TOYPOINTCLOUD_H


#include "WireCellNanoflann/nanoflann.h" 
#include "WireCellData/WCPointCloud.h"
#include "WireCellQuickhull/QuickHull.h"
#include "WireCellQuickhull/MathUtils.h"

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
    std::vector<std::pair<size_t,double>> get_closest(WireCell::Point& p, int N);
    std::vector<std::pair<size_t,double>> get_closest(WireCell::Point& p, double radius);
    

    
    
    
  protected:
    WireCell::WCPointCloud<double> cloud;
    my_kd_tree_t *index;
  };

   
}

#endif
