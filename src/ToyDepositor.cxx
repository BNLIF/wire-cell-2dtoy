#include "WCP2dToy/ToyDepositor.h"
#include "WCPData/Units.h"
#include "TMath.h"

using namespace::WCP;

ToyDepositor::ToyDepositor(WCP::FrameDataSource* fds1, int flag, float unit_dis, int toffset, float x_center, float y_center, float z_center, float rotate_angle)
  : fds(fds1)
  , flag(flag)
  , unit_dis(unit_dis)
  , toffset(toffset)
  , x_center(x_center)
  , y_center(y_center)
  , z_center(z_center)
  , rotate_angle(rotate_angle)
{
  x_shift = 0.;
  y_shift = 0.;
  z_shift = 0.;
}

const PointValueVector& ToyDepositor::depositions(int frame_number) const{
  if (frame_number >=0 && frame_number < fds->size()){
    if (fds->jump(frame_number) < 0) {
      std::cerr << "Failed to go to frame " << frame_number << std::endl;
    }
  }

  WCP::SimDataSource* sim = dynamic_cast<WCP::SimDataSource*>(fds);
  WCP::SimTruthSelection sts = sim->truth();
  
  std::cout<<"("<<x_center<<","<<y_center<<","<<z_center<<"), angle = "<<rotate_angle*180/TMath::Pi()<<std::endl;
  mchits.clear();
  time_offset.clear();
  for (size_t itruth = 0; itruth < sts.size(); ++itruth) {
    const WCP::SimTruth* st = sts[itruth];
    PointValue p;
    p.first.x = (st->x() * units::cm - x_center*units::cm)*cos(rotate_angle) - (st->z() * units::cm - z_center*units::cm)*sin(rotate_angle) + x_center*units::cm + x_shift*units::cm;
    p.first.y = st->y() * units::cm + y_shift*units::cm;
    p.first.z = (st->x() * units::cm - x_center*units::cm)*sin(rotate_angle) + (st->z() * units::cm - z_center*units::cm)*cos(rotate_angle) + z_center*units::cm + z_shift*units::cm;

    p.second = st->charge()/3.;

    // p.second = 10000;
    //p.first.x /= units::cm;
    //p.first.y /= units::cm;
    //p.first.z /= units::cm;
    //std::cout<< itruth <<": ("<<st->x()<<","<<st->y()<<","<<st->z()<<") --> ("<<p.first.x<<","<<p.first.y<<","<<p.first.z<<")"<<std::endl;
    
    //hack for test, to be removed
    //   if (p.first.x > 0*units::cm && p.first.x < 5.*units::cm ) {
    // && p.first.y > -200 * units::cm && p.first.y < 0 *units::cm &&
    // p.first.z > 920 * units::cm && p.first.z < 1000*units::cm ){
    
    mchits.push_back(p);
    
    int offset = st->tdc() - toffset - (st->x() * units::cm)/(0.5*unit_dis*units::millimeter);
    if (flag==0)
      offset = toffset;
    //std::cout << offset << std::endl;
    
    time_offset.push_back(offset);
    //  }
  }

  return mchits;
}

void ToyDepositor::translation(float xshift, float yshift, float zshift)
{
  x_shift = xshift;
  y_shift = yshift;
  z_shift = zshift;
}
