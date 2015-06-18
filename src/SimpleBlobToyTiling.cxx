#include "WireCell2dToy/SimpleBlobToyTiling.h"
using namespace WireCell;
WireCell2dToy::SimpleBlobToyTiling::SimpleBlobToyTiling(WireCell2dToy::ToyTiling& toytiling1, WireCell2dToy::MergeToyTiling& mergetiling1, WireCell2dToy::ToyMatrix& toymatrix1){
  toytiling = &toytiling1;
  mergetiling = &mergetiling1;
  toymatrix = &toymatrix1;
  
}

WireCell2dToy::SimpleBlobToyTiling::~SimpleBlobToyTiling(){
}
