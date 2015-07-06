#include "WireCell2dToy/FrameDataSource.h"
#include "WireCellData/Point.h"
#include "TRandom.h"

WireCell2dToy::FrameDataSource::FrameDataSource(int nevents,
						const WireCell::GeomDataSource& gds)
    : WireCell::FrameDataSource()
    , Nevent(nevents)
    , gds(gds)
{
}
WireCell2dToy::FrameDataSource::~FrameDataSource()
{
}

int WireCell2dToy::FrameDataSource::size() const
{
  return Nevent;
}

int WireCell2dToy::FrameDataSource::jump(int frame_number)
{
  if (frame_number >= Nevent) frame_number = Nevent;

  // fixme: we need to centralize random number seeding!
  gRandom->SetSeed(frame_number);
  frame.clear();
  mctruth.clear();
  
  //main code to construct a frame
  // first create many traces and initialize the channel number
  std::vector <WireCell::Trace> traces(8400);
  for (int i=0;i!=8400;i++){	// fixme: use the gds?
    traces[i].chid = i;
    traces[i].tbin = 0;
    traces[i].charge.push_back(0);
  }
  // start the simulation of random some points and put in charge
  //int npoint = gRandom->Uniform(5,20);
  int npoint = 1;
  float charge;

  // fixme: use the gds?
  //generate point inside a 1m by 1m region
  Double_t y_low = -1 * units::m, y_high = 1 * units::m;
  Double_t z_low = 5 * units::m, z_high = 6 * units::m;


  const WireCell::GeomWire * wire_closest;
  WireCell::WirePlaneType_t plane;

  for (int i=0;i!=npoint;i++){
    WireCell::Point p;
    p.x = -1;
    // p.y = -0.672 * units::m;//gRandom->Uniform(y_low, y_high);
    // p.z = 5.718 * units::m;//gRandom->Uniform(z_low, z_high);

    p.y = gRandom->Uniform(y_low, y_high);
    p.z = gRandom->Uniform(z_low, z_high);
    
    charge = gRandom->Uniform(5,30);
    mctruth.push_back(WireCell::PointValue(p,charge));
    

    for (int j=0;j!=3;j++){
      plane = static_cast<WireCell::WirePlaneType_t>(j);
      wire_closest = gds.closest(p,plane);
      traces[wire_closest->channel()].charge[0] += charge;
    }
  }
  

  // only save the ones with non-zero charge into a frame
  for (int i=0;i!=8400;i++){	// fixme: where does 8400 come from?  use the gds?
    if (traces[i].charge[0] >0){
      // add 5% uncertainties 
      traces[i].charge[0] = gRandom->Gaus(traces[i].charge[0],0.05*traces[i].charge[0]);
      frame.traces.push_back(traces[i]);
    }
  }

  frame.index = frame_number;
  return frame.index;
}
