#include "WireCell2dToy/ToyTracking.h"

using namespace WireCell;
WireCell2dToy::ToyTracking::ToyTracking(WireCell2dToy::ToyCrawler& toycrawler){
  for (int i=0;i!=toycrawler.Get_allMCT().size();i++){
    WCTrack *track = new WCTrack(*toycrawler.Get_allMCT().at(i));
    tracks.push_back(track);
    // std::cout << "Xin " << track->get_end_scells().size()<<std::endl;
  }

  


}

WireCell2dToy::ToyTracking::~ToyTracking(){
  for (int i=0;i!=tracks.size();i++){
    delete tracks.at(i);
  }
  for (int i=0;i!=vertices.size();i++){
    delete vertices.at(i);
  }
}
