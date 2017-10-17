#include "WireCell2dToy/CalcPoints.h"

using namespace WireCell;

void WireCell2dToy::calc_boundary_points(WireCell::GeomDataSource& gds, WireCell::PR3DCluster* cluster){
  SMGCSelection mcells = cluster->get_mcells();
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    calc_boundary_points(gds,*it);
  }
}

void WireCell2dToy::calc_sampling_points(WireCell::GeomDataSource& gds, WireCell::PR3DCluster* cluster){
  SMGCSelection mcells = cluster->get_mcells();
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    calc_sampling_points(gds,*it);
  }
}

void WireCell2dToy::calc_boundary_points(WireCell::GeomDataSource& gds, WireCell::SlimMergeGeomCell* mcell){

}

void WireCell2dToy::calc_sampling_points(WireCell::GeomDataSource& gds, WireCell::SlimMergeGeomCell* mcell){

}
