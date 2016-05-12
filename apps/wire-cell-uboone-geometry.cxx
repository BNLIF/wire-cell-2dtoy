#include "WireCellSst/GeomDataSource.h"
#include "WireCellSst/DatauBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/BadTiling.h"

#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCell2dToy/SimpleBlobToyTiling.h"

#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/ToyMatrixExclusive.h"
#include "WireCell2dToy/ToyMatrixKalman.h"
#include "WireCell2dToy/ToyMatrixIterate.h"
#include "WireCell2dToy/ToyMatrixIterate_SingleWire.h"
#include "WireCell2dToy/ToyMatrixIterate_Only.h"


#include "WireCell2dToy/ToyMatrixMarkov.h"
#include "WireCell2dToy/ToyMetric.h"
#include "WireCell2dToy/BlobMetric.h"

#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"

#include "WireCellData/GeomCluster.h"
//#include "WireCellNav/SliceDataSource.h"


#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/SimDataSource.h"
#include "WireCellNav/SliceDataSource.h"
#include "WireCellSst/Util.h"
#include "WireCellData/SimTruth.h"
#include "WireCell2dToy/ToyDepositor.h"
#include "WireCellNav/GenerativeFDS.h"
#include "WireCell2dToy/ToySignalSimu.h"
#include "WireCell2dToy/ToySignalSimuTrue.h"
#include "WireCell2dToy/DataSignalGaus.h"
#include "WireCell2dToy/DataSignalWien.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include <iostream>

using namespace WireCell;
using namespace std;




int main(int argc, char* argv[])
{
  if (argc < 2) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt -u[u channle #] -v[v channel #] -w[w channel #] -x[position cm] -y[position cm] -z[position cm]" << endl;
    return 1;
  }

  WireCellSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  
  int uwire = -1;
  int vwire = -1;
  int wwire = -1;
  double x_pos = 0;
  double y_pos = -100*units::m;
  double z_pos = -100*units::m;

  //  float unit_dis = 1.14753;  // 70 KV @ 226.5 V/cm
  float unit_dis = 1.119;  // 70 KV @ 273 V/cm
  //float unit_dis = 1.6;  // test

  int num_uvw=0;
  int num_xyz=0;
  
  for (Int_t i = 1; i != argc; i++){
    switch(argv[i][1]){
    case 'u':
      uwire = atoi(&argv[i][2]);
      num_uvw++;
      break;
    case 'v':
      vwire = atoi(&argv[i][2])-2400;
      num_uvw++;
      break;
    case 'w':
      wwire = atoi(&argv[i][2])-2400-2400;
      num_uvw++;
      break;
    case 'x':
      x_pos = atof(&argv[i][2])*units::cm; 
      num_xyz ++;
      break;
    case 'y':
      y_pos = atof(&argv[i][2])*units::cm; 
      num_xyz ++;
      break;
    case 'z':
      z_pos = atof(&argv[i][2])*units::cm; 
      num_xyz ++;
      break;
    }
  }

  x_pos -= 52 * unit_dis/10.*units::cm;

  if (num_uvw  > num_xyz){
    
    // convert wire num into position
    Vector p;
    if (uwire >=0 && vwire >=0){
      const GeomWire *u_wire = gds.by_planeindex(WirePlaneType_t(0),uwire);
      const GeomWire *v_wire = gds.by_planeindex(WirePlaneType_t(1),vwire);
      gds.crossing_point(*u_wire,*v_wire,p);
      
      std::cout << "UV (x,y,z): " << p.x/units::cm << ", " << p.y/units::cm << ", " << p.z/units::cm  << " cm" << std::endl;
    }else if (uwire >=0 && wwire >=0){
      const GeomWire *u_wire = gds.by_planeindex(WirePlaneType_t(0),uwire);
      const GeomWire *w_wire = gds.by_planeindex(WirePlaneType_t(2),wwire);
      gds.crossing_point(*u_wire,*w_wire,p);
      
      std::cout << "UW (x,y,z): " << p.x/units::cm << ", " << p.y/units::cm << ", " << p.z/units::cm  << " cm" << std::endl;
    }else if (vwire >=0 && wwire >=0){
      const GeomWire *w_wire = gds.by_planeindex(WirePlaneType_t(2),wwire);
      const GeomWire *v_wire = gds.by_planeindex(WirePlaneType_t(1),vwire);
      gds.crossing_point(*w_wire,*v_wire,p);
      
      std::cout << "VW (x,y,z): " << p.x/units::cm << ", " << p.y/units::cm << ", " << p.z/units::cm  << " cm" << std::endl;
    }
  }else{

    std::cout << x_pos/units::cm << " " << y_pos/units::cm << " " << z_pos/units::cm << " cm" << std::endl;

    // convert position into wire number
    Point p(x_pos,y_pos,z_pos);
    if (gds.contained_yz(p)){
      const GeomWire *u_wire = gds.closest(p,WirePlaneType_t(0));
      const GeomWire *v_wire = gds.closest(p,WirePlaneType_t(1));
      const GeomWire *w_wire = gds.closest(p,WirePlaneType_t(2));
      std::cout << "U: " << u_wire->index() << std::endl;
      std::cout << "V: " << v_wire->index() + 2400 << std::endl;
      std::cout << "W: " << w_wire->index() + 4800 << std::endl;
      std::cout << "T: " << (x_pos)/(unit_dis*units::mm)*2+3200 << " " << (x_pos)/(unit_dis*units::mm)*2/4.+800 << std::endl;
    }else{
      std::cout << "Point is outside the boundary! " << std::endl;
    }
  }
  

  return 0;
}
