#include "WCPSst/GeomDataSource.h"
#include "WCP2dToy/FrameDataSource.h"
#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"

#include "WCPNav/SliceDataSource.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TH1F.h"
#include <iostream>
using namespace WCP;
using namespace std;




int main(int argc, char* argv[])
{
 

  std::map<int,float> abc;
  abc[1] = 10;
  abc[2] = 20;
  
  cout << abc[3] << endl;
  return 0;
}
