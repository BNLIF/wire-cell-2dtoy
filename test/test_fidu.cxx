#include "WCP2dToy/ToyFiducial.h"

using namespace std;

int main(int argc, char * argv[]){
  WCP2dToy::ToyFiducial *fid = new WCP2dToy::ToyFiducial();

  WCP::Point p(30.0*units::cm,30*units::cm,30*units::cm);
  WCP::Point p1(110.0*units::cm,0*units::cm,30*units::cm);
  
  std::cout << fid->inside_fiducial_volume(p) << " " << fid->inside_fiducial_volume(p1) << std::endl;;
  
  return 0;
}
