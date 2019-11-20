#ifndef WCP2dToy_TOYLIBPE_H
#define WCP2dToy_TOYLIBPE_H

#include "WCPData/Point.h"

#include "TFile.h"
#include "TTree.h"

namespace WCP2dToy{

  typedef std::vector<WCP::Point> pointVec;

  class ToyLibPe{
  public:
    ToyLibPe(const char* root_file);
    ~ToyLibPe();
    
    std::vector<double> getXYZQ(std::vector<float> *x, std::vector<float> *y, std::vector<float> *z, std::vector<float> *e);
    std::vector<double> getShiftedXYZQ(std::vector<float> *x, std::vector<float> *y, std::vector<float> *z, std::vector<float> *e, float xOffset, float yOffset, float zOffset);

    std::pair<float,float> getXminXmax(std::vector<float> *x);
    std::pair<float,float> getShiftedXminXmax(std::vector<float> *x, float xOffset);

    bool xInFidVol(std::pair<float,float> &x);
    bool xAtBoundary(std::pair<float,float> &x);
    std::vector<double> fromClusterFindPeDist(std::vector<float> *x, std::vector<float> *y, std::vector<float> *z, std::vector<float> *e, float xOffset, float yOffset, float zOffset);
    
  protected:
    int convert_to_voxel_id(WCP::Point& p);
    std::map<int,int> getMapLibPmt();
    std::map<int,int> getMapPmtLib();
    std::vector<std::vector<std::pair<int,float> > > getPhotonLibrary();
    
    TFile *file;
    TTree *T;
    std::vector<int> *trackId;
    std::vector<float> *energy;
    std::vector<float> *numElectrons;
    std::vector<float> *x;
    std::vector<float> *y;
    std::vector<float> *z;
  };
}
#endif
