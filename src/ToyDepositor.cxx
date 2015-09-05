#include "WireCell2dToy/ToyDepositor.h"
#include "WireCellData/Units.h"

using namespace::WireCell;

ToyDepositor::ToyDepositor(WireCell::FrameDataSource* fds1, int flag, float unit_dis, int toffset)
  : fds(fds1)
  , flag(flag)
  , unit_dis(unit_dis)
  , toffset(toffset)
{
}


const PointValueVector& ToyDepositor::depositions(int frame_number) const{
  if (frame_number >=0 && frame_number < fds->size()){
    if (fds->jump(frame_number) < 0) {
      std::cerr << "Failed to go to frame " << frame_number << std::endl;
    }
  }

  WireCell::SimDataSource* sim = dynamic_cast<WireCell::SimDataSource*>(fds);
  WireCell::SimTruthSelection sts = sim->truth();
  

  mchits.clear();
  time_offset.clear();
  for (size_t itruth = 0; itruth < sts.size(); ++itruth) {
     const WireCell::SimTruth* st = sts[itruth];
     PointValue p;
     p.first.x = st->x() * units::cm;
     p.first.y = st->y() * units::cm;
     p.first.z = st->z() * units::cm;
     p.second = st->charge()/3.;
     
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


