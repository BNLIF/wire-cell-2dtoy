#ifndef WireCell2dToy_TOYPOINTCLOUD_H
#define WireCell2dToy_TOYPOINTCLOUD_H


#include "WireCellNanoflann/nanoflann.h" 
#include "WireCellData/WCPointCloud.h"
#include "WireCellQuickhull/QuickHull.h"
#include "WireCellQuickhull/MathUtils.h"

typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<float, WireCell::WCPointCloud<float> > ,
    WireCell::WCPointCloud<float>,
    3 /* dim */
  > my_kd_tree_t;

namespace WireCell2dToy{
  
  class ToyPointCloud {
  public:
    ToyPointCloud();
    ~ToyPointCloud();

    void AddPoint(WireCell::Point& p, WireCell::SlimMergeGeomCell *mcell);
    void AddPoints(WireCell::PointVector& ps, WireCell::SlimMergeGeomCell *mcell);
    void build_index();
    
    
  protected:
    WireCell::WCPointCloud<float> cloud;
    my_kd_tree_t *index;
  };

   
}

#endif
