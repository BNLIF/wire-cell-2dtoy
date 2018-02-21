#include "WireCell2dToy/ToyLibPe.h"

#include "WireCellRess/LassoModel.h"
#include "WireCellRess/ElasticNetModel.h"
#include <Eigen/Dense>
#include "TChain.h"

#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"
#include "WireCellData/Units.h"

using namespace std;
using namespace Eigen;
using namespace WireCell;

WireCell2dToy::ToyLibPe::ToyLibPe(const char* root_file){
  file = new TFile(root_file);
  T = (TTree*)file->Get("tout");

  trackId = new std::vector<int>;
  energy = new std::vector<float>;
  numElectrons = new std::vector<float>;
  x = new std::vector<float>;
  y = new std::vector<float>;
  z = new std::vector<float>;
}

WireCell2dToy::ToyLibPe::~ToyLibPe(){
  delete trackId;
  delete energy;
  delete numElectrons;
  delete x;
  delete y;
  delete z;
  delete T;
  delete file;
}

int WireCell2dToy::ToyLibPe::convert_to_voxel_id(WireCell::Point &p){
  int voxel_x_id = round((p.x/units::cm+63.435-5.1096/2.)/5.1096);
  int voxel_y_id = round((p.y/units::cm+191.61-5.1096/2.)/5.1096);
  int voxel_z_id = round((p.z/units::cm+92.375-3.05437/2.)/3.05437);
  if(voxel_x_id<0) voxel_x_id = 0;
  if(voxel_x_id>=75) voxel_x_id = 74;
  if(voxel_y_id<0) voxel_y_id = 0;
  if(voxel_y_id>=75) voxel_y_id = 74;
  if(voxel_z_id<0) voxel_z_id = 0;
  if(voxel_z_id>=400) voxel_z_id = 399;
  int voxel_id = voxel_z_id*75*75 + voxel_y_id*75 + voxel_x_id;
  return voxel_id;
}

std::map<int,int> WireCell2dToy::ToyLibPe::getMapLibPmt(){
  map<int,int> map_lib_pmt;
  map_lib_pmt[1]=2;
  map_lib_pmt[0]=4;
  map_lib_pmt[3]=0;
  map_lib_pmt[2]=5;
  map_lib_pmt[5]=1;
  map_lib_pmt[4]=6;
  map_lib_pmt[6]=3;

  map_lib_pmt[9]=7;
  map_lib_pmt[7]=9;
  map_lib_pmt[8]=11;
  map_lib_pmt[11]=8;
  map_lib_pmt[10]=12;
  map_lib_pmt[12]=10;

  map_lib_pmt[14]=13;
  map_lib_pmt[13]=15;
  map_lib_pmt[15]=17;
  map_lib_pmt[17]=14;
  map_lib_pmt[16]=18;
  map_lib_pmt[18]=16;

  map_lib_pmt[21]=19;
  map_lib_pmt[22]=20;
  map_lib_pmt[19]=21;
  map_lib_pmt[20]=23;
  map_lib_pmt[23]=24;
  map_lib_pmt[24]=22;

  map_lib_pmt[26]=25;
  map_lib_pmt[25]=28;
  map_lib_pmt[27]=30;
  map_lib_pmt[28]=31;
  map_lib_pmt[31]=29;
  map_lib_pmt[30]=27;
  map_lib_pmt[29]=26;
  
  return map_lib_pmt;
}

std::map<int,int> WireCell2dToy::ToyLibPe::getMapPmtLib(){
  map<int,int> map_pmt_lib;
  map_pmt_lib[2]=1;
  map_pmt_lib[4]=0;
  map_pmt_lib[0]=3;
  map_pmt_lib[5]=2;
  map_pmt_lib[1]=5;
  map_pmt_lib[6]=4;
  map_pmt_lib[3]=6;

  map_pmt_lib[7]=9;
  map_pmt_lib[9]=7;
  map_pmt_lib[11]=8;
  map_pmt_lib[8]=11;
  map_pmt_lib[12]=10;
  map_pmt_lib[10]=12;

  map_pmt_lib[13]=14;
  map_pmt_lib[15]=13;
  map_pmt_lib[17]=15;
  map_pmt_lib[14]=17;
  map_pmt_lib[18]=16;
  map_pmt_lib[16]=18;

  map_pmt_lib[19]=21;
  map_pmt_lib[20]=22;
  map_pmt_lib[21]=19;
  map_pmt_lib[23]=20;
  map_pmt_lib[24]=23;
  map_pmt_lib[22]=24;

  map_pmt_lib[25]=26;
  map_pmt_lib[28]=25;
  map_pmt_lib[30]=27;
  map_pmt_lib[31]=28;
  map_pmt_lib[29]=31;
  map_pmt_lib[27]=30;
  map_pmt_lib[26]=29;
  
  return map_pmt_lib;
}

std::vector<std::vector<std::pair<int,float> > > WireCell2dToy::ToyLibPe::getPhotonLibrary(){
  TChain *Tlib = new TChain("/pmtresponse/PhotonLibraryData","/pmtresponse/PhotonLibraryData");
  Tlib->AddFile("./uboone_photon_library.root");
  Int_t Voxel;
  Int_t OpChannel;
  Float_t Visibility;
  Tlib->SetBranchAddress("Voxel",&Voxel);
  Tlib->SetBranchAddress("OpChannel",&OpChannel);
  Tlib->SetBranchAddress("Visibility",&Visibility);
  
  std::vector<std::vector<std::pair<int,float> > > photon_library;
  photon_library.resize(400*75*75);
  
  for (int i=0;i!=Tlib->GetEntries();i++){
    Tlib->GetEntry(i);
    photon_library.at(Voxel).push_back(make_pair(OpChannel,Visibility));
  }
  return photon_library;
}

std::vector<double> WireCell2dToy::ToyLibPe::getXYZQ(std::vector<float> *x, std::vector<float> *y, std::vector<float> *z, std::vector<float> *e){
  double sumX = 0, sumY = 0, sumZ = 0, sumE = 0;

  int nDeps = (int)x->size();
  for(int i=0; i<nDeps; i++){
    sumX += x->at(i)*e->at(i);
    sumY += y->at(i)*e->at(i);
    sumZ += z->at(i)*e->at(i);
    sumE += e->at(i);
  }
  vector<double> xyzq = { sumX/sumE, sumY/sumE, sumZ/sumE, sumE};
  return xyzq;
}

std::vector<double> WireCell2dToy::ToyLibPe::getShiftedXYZQ(std::vector<float> *x, std::vector<float> *y, std::vector<float> *z, std::vector<float> *e, float xOffset, float yOffset, float zOffset){
  double sumX = 0, sumY = 0, sumZ = 0, sumE = 0;

  int nDeps = (int)x->size();
  for(int i=0; i<nDeps; i++){
    sumX += (x->at(i)-xOffset)*e->at(i);
    sumY += (y->at(i)-yOffset)*e->at(i);
    sumZ += (z->at(i)-zOffset)*e->at(i);
    sumE += e->at(i);
  }
  vector<double> xyzq = { sumX/sumE, sumY/sumE, sumZ/sumE, sumE};
  return xyzq;
}

std::pair<float,float> WireCell2dToy::ToyLibPe::getXminXmax(std::vector<float> *x){
  std::pair<float,float> result;
  float min = *std::min_element(x->begin(), x->end());
  float max = *std::max_element(x->begin(), x->end());
  result = std::make_pair(min, max);
  return result;
}

std::pair<float,float> WireCell2dToy::ToyLibPe::getShiftedXminXmax(std::vector<float> *x, float xOffset){
  std::pair<float,float> result;

  float min = *std::min_element(x->begin(), x->end())-xOffset;
  float max = *std::max_element(x->begin(), x->end())-xOffset;
  result = std::make_pair(min, max);
  return result;
}

bool WireCell2dToy::ToyLibPe::xInFidVol(std::pair<float,float> &x){
  bool flag = false;
  double high_x_cut = 256 * units::cm;
  double high_x_cut_ext1 = + 1*units::cm;
  double low_x_cut = 0*units::cm;
  double low_x_cut_ext1 = - 2*units::cm;
  
  if(x.first > low_x_cut + low_x_cut_ext1 &&
     x.second > low_x_cut &&
     x.second < high_x_cut + high_x_cut_ext1 &&
     x.first < high_x_cut){
    flag = true;
  }
  
  return flag;
}

bool WireCell2dToy::ToyLibPe::xAtBoundary(std::pair<float,float> &x){
  bool flag = false;
  double high_x_cut = 256 * units::cm;
  double high_x_cut_ext1 = + 1*units::cm;
  double low_x_cut = 0*units::cm;
  double low_x_cut_ext1 = - 2*units::cm;
  double low_x_cut_ext2 = +0.5*units::cm;
  
  if( (x.first <=low_x_cut + low_x_cut_ext2 &&
       x.first > low_x_cut + low_x_cut_ext1 ) ||
      (x.second >= high_x_cut - high_x_cut_ext1 &&
       x.second < high_x_cut + high_x_cut_ext1) ){
    flag = true;
  }
  return flag;
}

std::vector<double> WireCell2dToy::ToyLibPe::fromClusterFindPeDist(std::vector<float> *x, std::vector<float> *y, std::vector<float> *z, std::vector<float> *e, float xOffset, float yOffset, float zOffset){
  std::vector<double> pred_pmt_light;
  pred_pmt_light.resize(32,0);

  std::vector<std::vector<std::pair<int,float> > > photon_library = getPhotonLibrary();
  std::map<int,int> map_lib_pmt = getMapLibPmt();
  
  int nDep = (int)x->size();
  for(int i=0; i<nDep; i++){
    Point p;
    p.x = x->at(i) - xOffset;
    p.y = y->at(i) - yOffset;
    p.z = z->at(i) - zOffset;
    int voxel_id = convert_to_voxel_id(p);

    std::vector<std::pair<int,float>>& pmt_list = photon_library.at(voxel_id);

    for(auto itr = pmt_list.begin(); itr!=pmt_list.end(); itr++){
      pred_pmt_light.at(map_lib_pmt[itr->first]) += e->at(i) * itr->second;

    }
  }

  double scaling_light_mag = 0.01*1.5;
  for(int i=0; i<32; i++){
    pred_pmt_light.at(i) *= scaling_light_mag;
  }  

    return pred_pmt_light;
}

