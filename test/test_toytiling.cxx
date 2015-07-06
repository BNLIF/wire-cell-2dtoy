#include "WireCellNav/ExampleWires.h"
#include "WireCellNav/GeomDataSource.h"
#include "WireCellNav/GenerativeFDS.h"
#include "WireCellNav/PepperDepositor.h"
#include "WireCellNav/SliceDataSource.h"

#include "WireCell2dToy/ToyTiling.h"

#include "WireCellData/Units.h"
#include "WireCellData/Slice.h"

#include "WireCellIface/IWireGeometry.h"

#include <iostream>
#include <string>

using namespace std;
using namespace WireCell;


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
    
    IWireGeometry* wires = make_example_wires(10*units::mm);
    GeomDataSource gds;
    gds.use_wires(*wires);

    GenerativeFDS fds(dep, gds);

    fds.jump(0);    
    SliceDataSource sds(fds);
    
    sds.jump(0);

    const Slice& slice = sds.get();
    
    WireCell2dToy::ToyTiling *toytiling = new WireCell2dToy::ToyTiling(slice,gds);
    
    GeomCellSelection allcell = toytiling->get_allcell();
    
    cout << allcell.size() << endl;
    
    return 0;
}
