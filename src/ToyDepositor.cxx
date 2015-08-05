#include "WireCell2dToy/ToyDepositor.h"
#include "WireCellData/Units.h"

using namespace::WireCell;

ToyDepositor::ToyDepositor(WireCell::FrameDataSource* fds1, int flag)
  : fds(fds1)
  , flag(flag)
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
     mchits.push_back(p);

     int offset = st->tdc() - 3200 - (st->x() * units::cm)/(0.5*1.605723*units::millimeter);
     if (flag==0)
       offset = 3200;
     //std::cout << offset << std::endl;

     time_offset.push_back(offset);
  }

  return mchits;
}


