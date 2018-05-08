#include "WireCell2dToy/ToyFiducial.h"

using namespace std;

int main(int argc, char * argv[]){
  WireCell2dToy::ToyFiducial *fid = new WireCell2dToy::ToyFiducial();

  WireCell::Point p(30.0*units::cm,30*units::cm,30*units::cm);
  WireCell::Point p1(110.0*units::cm,0*units::cm,30*units::cm);
  
  std::cout << fid->inside_fiducial_volume(p) << " " << fid->inside_fiducial_volume(p1) << std::endl;;
  
  return 0;
}
