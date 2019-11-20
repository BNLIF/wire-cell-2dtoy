#include "WCPNav/ExampleGDS.h"
#include "WCPNav/GenerativeFDS.h"
#include "WCPNav/PepperDepositor.h"
#include "WCPNav/SliceDataSource.h"

#include "WCP2dToy/ToyTiling.h"

#include "WCPData/Units.h"
#include "WCPData/Slice.h"

#include <iostream>
#include <string>

using namespace std;
using namespace WCP;


template<typename OK>
void assert(const OK& ok, string msg="error")
{
    if (ok) { return; }
    cerr << msg << endl;
    exit(-1);
}


int main () {
    using namespace std;
    const float width = 0.5*units::m;
    const PepperDepositor::MinMax drift_dim(0, 2.0*width), trans_dim(-width,width);
    const PepperDepositor::MinMax charge(1,100);
    PepperDepositor dep(drift_dim, trans_dim, trans_dim, charge, 5);
    
    GeomDataSource* gds = make_example_gds(10*units::mm);

    GenerativeFDS fds(dep, *gds);

    fds.jump(0);    
    SliceDataSource sds(fds);
    
    sds.jump(0);

    const Slice& slice = sds.get();
    
    WCP2dToy::ToyTiling *toytiling = new WCP2dToy::ToyTiling(slice,*gds);
    
    GeomCellSelection allcell = toytiling->get_allcell();
    
    cout << allcell.size() << endl;
    
    return 0;
}
