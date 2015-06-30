#include "WireCell2dToy/ToySignalSimu.h"

using namespace WireCell;

WireCell2dToy::ToySignalSimuFDS::ToySignalSimuFDS(const WireCell::FrameDataSource& fds, const WireCell::GeomDataSource& gds,
						  int bins_per_frame, int nframes_total)
  : fds(fds)
  , gds(gds)
  , bins_per_frame(bins_per_frame)
  , max_frames(nframes_total)
{
}
