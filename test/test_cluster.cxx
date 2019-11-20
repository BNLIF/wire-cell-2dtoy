#include "WCPData/Slim3DCluster.h"
#include "WCPSst/GeomDataSource.h"

using namespace WCP;
using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 2) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt" << std::endl;
    return 1;
  }
  WCPSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;

  const GeomWireSelection& uwires = gds.wires_in_plane(WirePlaneType_t(0));
  const GeomWireSelection& vwires = gds.wires_in_plane(WirePlaneType_t(0));
  const GeomWireSelection& wwires = gds.wires_in_plane(WirePlaneType_t(0));
  
  
  SlimMergeGeomCell *mcell1 = new SlimMergeGeomCell(1);
  mcell1->SetTimeSlice(1);
  for (int i=0;i!=100;i++){
    mcell1->AddWire(uwires.at(i),WirePlaneType_t(0));
    mcell1->AddWire(vwires.at(i),WirePlaneType_t(1));
    mcell1->AddWire(wwires.at(i),WirePlaneType_t(2));
  }
  SlimMergeGeomCell *mcell2 = new SlimMergeGeomCell(2);
  mcell2->SetTimeSlice(2);
  for (int i=25; i!=35;i++){
    mcell2->AddWire(uwires.at(i),WirePlaneType_t(0));
    mcell2->AddWire(vwires.at(i),WirePlaneType_t(1));
    mcell2->AddWire(wwires.at(i),WirePlaneType_t(2));
  }
  SlimMergeGeomCell *mcell3 = new SlimMergeGeomCell(2);
  mcell3->SetTimeSlice(2);
  for (int i=0; i!=10;i++){
    mcell3->AddWire(uwires.at(i),WirePlaneType_t(0));
    mcell3->AddWire(vwires.at(i),WirePlaneType_t(1));
    mcell3->AddWire(wwires.at(i),WirePlaneType_t(2));
  }

  SlimMergeGeomCell *mcell4 = new SlimMergeGeomCell(2);
  mcell4->SetTimeSlice(3);
  mcell4->add_bad_planes(WirePlaneType_t(2));
  for (int i=8; i!=23;i++){
    mcell4->AddWire(uwires.at(i),WirePlaneType_t(0));
    mcell4->AddWire(vwires.at(i),WirePlaneType_t(1));
    mcell4->AddWire(wwires.at(i),WirePlaneType_t(2));
  }
  
  SlimMergeGeomCell *mcell4p = new SlimMergeGeomCell(2);
  mcell4p->SetTimeSlice(3);
  mcell4p->add_bad_planes(WirePlaneType_t(2));
  for (int i=6; i!=25;i++){
    mcell4p->AddWire(uwires.at(i),WirePlaneType_t(0));
    mcell4p->AddWire(vwires.at(i),WirePlaneType_t(1));
    mcell4p->AddWire(wwires.at(i),WirePlaneType_t(2));
  }
  

  SlimMergeGeomCell *mcell5 = new SlimMergeGeomCell(2);
  mcell5->SetTimeSlice(4);
  for (int i=8; i!=23;i++){
    mcell5->AddWire(uwires.at(i),WirePlaneType_t(0));
    mcell5->AddWire(vwires.at(i),WirePlaneType_t(1));
    mcell5->AddWire(wwires.at(i),WirePlaneType_t(2));
  }
 
  
  Slim3DCluster *cluster = new Slim3DCluster(*mcell1);
  cluster->AddCell(*mcell2);
  cluster->AddCell(*mcell3);
  cluster->AddCell(*mcell4);
  cluster->AddCell(*mcell5);
  cluster->Calc_Projection();

  
  Slim3DCluster *cluster1 = new Slim3DCluster(*mcell1);
  cluster1->AddCell(*mcell2);
  // cluster1->AddCell(*mcell3);
  cluster1->AddCell(*mcell4p);
  cluster1->AddCell(*mcell5);
  cluster1->Calc_Projection();

  
  Projected2DCluster *u_proj = cluster->get_projection(WirePlaneType_t(0));
  Projected2DCluster *v_proj = cluster->get_projection(WirePlaneType_t(1));
  Projected2DCluster *w_proj = cluster->get_projection(WirePlaneType_t(2));

  Projected2DCluster *u1_proj = cluster1->get_projection(WirePlaneType_t(0));
  Projected2DCluster *v1_proj = cluster1->get_projection(WirePlaneType_t(1));
  Projected2DCluster *w1_proj = cluster1->get_projection(WirePlaneType_t(2));
  
  std::cout << "Status: " << u_proj->judge_coverage(u1_proj) << std::endl;
  std::cout << "Status: " << u1_proj->judge_coverage(u_proj) << std::endl;

  std::vector<int> comp_results = u_proj->calc_coverage(u1_proj);
  std::cout << "Xin: " << " " << comp_results.at(0) << " " << comp_results.at(1) << " " <<
	  comp_results.at(2) << " " << comp_results.at(3) << std::endl;
	

  
}
