#include "WCPSst/GeomDataSource.h"
#include "WCPSst/DatauBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/BadTiling.h"

#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/TruthToyTiling.h"
#include "WCP2dToy/SimpleBlobToyTiling.h"

#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/ToyMatrixExclusive.h"
#include "WCP2dToy/ToyMatrixKalman.h"
#include "WCP2dToy/ToyMatrixIterate.h"
#include "WCP2dToy/ToyMatrixIterate_SingleWire.h"
#include "WCP2dToy/ToyMatrixIterate_Only.h"


#include "WCP2dToy/ToyMatrixMarkov.h"
#include "WCP2dToy/ToyMetric.h"
#include "WCP2dToy/BlobMetric.h"

#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"

#include "WCPData/GeomCluster.h"
//#include "WCPNav/SliceDataSource.h"


#include "WCPNav/FrameDataSource.h"
#include "WCPNav/SimDataSource.h"
#include "WCPNav/SliceDataSource.h"
#include "WCPSst/Util.h"
#include "WCPData/SimTruth.h"
#include "WCP2dToy/ToyDepositor.h"
#include "WCPNav/GenerativeFDS.h"
#include "WCP2dToy/ToySignalSimu.h"
#include "WCP2dToy/ToySignalSimuTrue.h"
#include "WCP2dToy/DataSignalGaus.h"
#include "WCP2dToy/DataSignalWien.h"

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

using namespace WCP;
using namespace std;




int main(int argc, char* argv[])
{
  if (argc < 2) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt " << endl;
    return 1;
  }

  WCPSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  
  int uwire = -1;
  int vwire = -1;
  int wwire = -1;
  double x_pos = 0;
  double y_pos = -100*units::m;
  double z_pos = -100*units::m;

  float unit_dis = 1.14753;  // 70 KV @ 226.5 V/cm

  int num_uvw=0;
  int num_xyz=0;
  
  for (Int_t i = 1; i != argc; i++){
    switch(argv[i][1]){
    case 'u':
      uwire = atoi(&argv[i][2]);
      num_uvw++;
      break;
    case 'v':
      vwire = atoi(&argv[i][2]);
      num_uvw++;
      break;
    case 'w':
      wwire = atoi(&argv[i][2]);
      num_uvw++;
      break;
    case 'x':
      x_pos = atof(&argv[i][2])*units::m; 
      num_xyz ++;
      break;
    case 'y':
      y_pos = atof(&argv[i][2])*units::m; 
      num_xyz ++;
      break;
    case 'z':
      z_pos = atof(&argv[i][2])*units::m; 
      num_xyz ++;
      break;
    }
  }

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
    // convert position into wire number
    Point p(x_pos,y_pos,z_pos);
    if (gds.contained_yz(p)){
      const GeomWire *u_wire = gds.closest(p,WirePlaneType_t(0));
      const GeomWire *v_wire = gds.closest(p,WirePlaneType_t(1));
      const GeomWire *w_wire = gds.closest(p,WirePlaneType_t(2));
      std::cout << "U: " << u_wire->index() << std::endl;
      std::cout << "V: " << v_wire->index() << std::endl;
      std::cout << "W: " << w_wire->index() << std::endl;
      std::cout << "T: " << (x_pos)/(unit_dis*units::mm)*2+3200<< std::endl;
    }else{
      std::cout << "Point is outside the boundary! " << std::endl;
    }
  }
  

  return 0;
}
