#ifndef WIRECELL_TOYDEPOSITOR_H
#define WIRECELL_TOYDEPOSITOR_H

#include "WireCellNav/Depositor.h"
#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/SimDataSource.h"

namespace WireCell {
  class ToyDepositor : public WireCell::Depositor {
  public:
    ToyDepositor(WireCell::FrameDataSource* fds, int flag = 0, float unit_dis = 1.6, int toffset = 3200, float x_center = 0, float y_center = 0, float z_center = 0, float rotate_angle = 0);
    virtual const PointValueVector& depositions(int frame_number) const;
    virtual std::vector<int>& timeoffset() const{return time_offset;};
    void clear(){mchits.clear();time_offset.clear();};
    void translation(float x_shift = 0, float y_shift = 0, float z_shift = 0);
    
  private:
    WireCell::FrameDataSource* fds;
    mutable PointValueVector mchits;
    mutable std::vector<int> time_offset;

    float x_center;
    float y_center;
    float z_center;
    float rotate_angle;
    float x_shift;
    float y_shift;
    float z_shift;
    
    int flag;
    int toffset;
    float unit_dis;
  };
}

#endif
