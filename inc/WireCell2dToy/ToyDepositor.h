#ifndef WIRECELL_TOYDEPOSITOR_H
#define WIRECELL_TOYDEPOSITOR_H

#include "WireCellNav/Depositor.h"
#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/SimDataSource.h"

namespace WireCell {
  class ToyDepositor : public WireCell::Depositor {
  public:
    ToyDepositor(WireCell::FrameDataSource* fds, int flag = 0, float unit_dis = 1.6);
    virtual const PointValueVector& depositions(int frame_number) const;
    virtual std::vector<int>& timeoffset() const{return time_offset;};
    void clear(){mchits.clear();time_offset.clear();};

  private:
    WireCell::FrameDataSource* fds;
    mutable PointValueVector mchits;
    mutable std::vector<int> time_offset;
    int flag;
    float unit_dis;
  };
}

#endif
