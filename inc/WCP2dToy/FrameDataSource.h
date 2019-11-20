#ifndef WIRECELLSST_FRAMEDATASOURCE_H
#define  WIRECELLSST_FRAMEDATASOURCE_H

#include "WCPNav/FrameDataSource.h"
#include "WCPNav/GeomDataSource.h"

namespace WCP2dToy {
  
    /** A FrameDataSource which generates toy charge on the fly.  This
     * FDS also extends the base to provide access to "mctruth"
     * values. */ 
  class FrameDataSource : public WCP::FrameDataSource {
  public:
      /// Produce a WCP2dToy::FrameDataSource that will generate nevents.
    FrameDataSource(int nevents, const WCP::GeomDataSource& gds);
    virtual ~FrameDataSource(); 
    
    virtual int size() const;
    virtual int jump(int frame_number); 

    /// Access generated cell charges
    const WCP::PointValueVector& cell_charges() { return mctruth; }
      // fixme: we should talk about this interface addition and come
      // up with a general way to handle MC truth info.

  private:

    int Nevent;
    const WCP::GeomDataSource& gds;
    WCP::PointValueVector mctruth;
  };
  
}

#endif
