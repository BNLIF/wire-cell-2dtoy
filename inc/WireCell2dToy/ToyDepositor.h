#ifndef WIRECELL_TOYDEPOSITOR_H
#define WIRECELL_TOYDEPOSITOR_H

#include "WireCellNav/Depositor.h"
#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/SimDataSource.h"

namespace WireCell {
  class ToyDepositor : public WireCell::Depositor {
  public:
    ToyDepositor(WireCell::FrameDataSource* fds);
    virtual const PointValueVector& depositions(int frame_number) const;

  private:
    WireCell::FrameDataSource* fds;
    mutable PointValueVector mchits;
  };
}

#endif
