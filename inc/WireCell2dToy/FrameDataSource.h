#ifndef WIRECELLSST_FRAMEDATASOURCE_H
#define  WIRECELLSST_FRAMEDATASOURCE_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GeomDataSource.h"

namespace WireCell2dToy {
  
    /** A FrameDataSource which generates toy charge on the fly.  This
     * FDS also extends the base to provide access to "mctruth"
     * values. */ 
  class FrameDataSource : public WireCell::FrameDataSource {
  public:
      /// Produce a WireCell2dToy::FrameDataSource that will generate nevents.
    FrameDataSource(int nevents, const WireCell::GeomDataSource& gds);
    virtual ~FrameDataSource(); 
    
    virtual int size() const;
    virtual int jump(int frame_number); 

    /// Access generated cell charges
    const WireCell::PointValueVector& cell_charges() { return mctruth; }
      // fixme: we should talk about this interface addition and come
      // up with a general way to handle MC truth info.

  private:

    int Nevent;
    const WireCell::GeomDataSource& gds;
    WireCell::PointValueVector mctruth;
  };
  
}

#endif
