#include "WireCell2dToy/TruthToyTiling.h"

using namespace WireCell;

WireCell2dToy::TruthToyTiling::TruthToyTiling(WireCell2dToy::ToyTiling& tiling, const WireCell::PointValueVector &pvv, int tbin, const GeomDataSource& gds, int offset1, float unit_dis){
  offset = offset1;

  
  float sum = 0;
  for (int itruth = 0; itruth < pvv.size(); ++itruth){
    if (int(pvv[itruth].first.x/2.0/unit_dis/units::mm + offset1) == tbin){
      const Point& p = pvv[itruth].first; // get the point
      float charge = pvv[itruth].second;
      
      

      if (gds.contained_yz(p)){

	Point p1 = p;// hack for now, gap in the tiling ... 
	//gds.avoid_gap(p1);
	
	sum += charge;
	const GeomWire* wire_u = gds.closest(p1, static_cast<WirePlaneType_t>(0));
	const GeomWire* wire_v = gds.closest(p1, static_cast<WirePlaneType_t>(1));
	const GeomWire* wire_w = gds.closest(p1, static_cast<WirePlaneType_t>(2));
	
	if (wire_u!=0&&wire_v!=0&&wire_w!=0){
	  
	  GeomWireSelection wires;
	  wires.push_back(wire_u);
	  wires.push_back(wire_v);
	  wires.push_back(wire_w);
	
	  const GeomCell* cell = tiling.cell(wires);
	  
	  if (cell!=0){
	    
	    //Point pp = cell->center();
	    //	std::cout << pp.x << " " << pp.y << " " << pp.z << std::endl;
	    
	    
	    if (cellchargemap.find(cell) == cellchargemap.end()){
	      //not found
	      cellchargemap[cell] = charge;
	    }else{
	      cellchargemap[cell] += charge;
	    }
	  }else{
	    // Point pp = cell->center();
	    // gds.avoid_gap(p1);
	    // std::cout << p.x << " " << p.y << " " << p.z << " " << std::endl;
	    //std::cout << itruth << " " << p1.x << " " << p1.y << " " << p1.z << " " << charge << std::endl;
	    // std::cout << charge << std::endl;
	  }
	}
      }
    }
  }
  
  //  std::cout << "Xin " << sum << std::endl;

}


WireCell2dToy::TruthToyTiling::TruthToyTiling(WireCell2dToy::ToyTiling& tiling, const WireCell::PointValueVector &pvv, const std::vector<int> &time_offset, int tbin, const GeomDataSource& gds, float unit_dis){
  float sum = 0;
  for (int itruth = 0; itruth < pvv.size(); ++itruth){
    if (int(pvv[itruth].first.x/2.0/unit_dis/units::mm + time_offset[itruth]/4. ) == tbin){
      const Point& p = pvv[itruth].first; // get the point
      float charge = pvv[itruth].second;
      
      

      if (gds.contained_yz(p)){

	Point p1 = p;// hack for now, gap in the tiling ... 
	//gds.avoid_gap(p1);
	
	sum += charge;
	const GeomWire* wire_u = gds.closest(p1, static_cast<WirePlaneType_t>(0));
	const GeomWire* wire_v = gds.closest(p1, static_cast<WirePlaneType_t>(1));
	const GeomWire* wire_w = gds.closest(p1, static_cast<WirePlaneType_t>(2));
	
	if (wire_u!=0&&wire_v!=0&&wire_w!=0){
	  
	  GeomWireSelection wires;
	  wires.push_back(wire_u);
	  wires.push_back(wire_v);
	  wires.push_back(wire_w);
	
	  const GeomCell* cell = tiling.cell(wires);
	  
	  if (cell!=0){
	    if (cellchargemap.find(cell) == cellchargemap.end()){
	      //not found
	      cellchargemap[cell] = charge;
	    }else{
	      cellchargemap[cell] += charge;
	    }
	  }
	}
      }
    }
  }
}



WireCell2dToy::TruthToyTiling::TruthToyTiling(WireCell2dToy::ToyTiling& tiling, const WireCell::PointValueVector &pvv, int tbin, const DetectorGDS& gds, int offset1, float unit_dis){
  offset = offset1;

  
  float sum = 0;
  for (int itruth = 0; itruth < pvv.size(); ++itruth){
    const Point& pt = pvv[itruth].first; // get the point
    float charge = pvv[itruth].second;

    short which_cryo = gds.in_which_cryo(pt);
    short which_apa = gds.in_which_apa(pt);
    
    const WrappedGDS *apa_gds = gds.get_apaGDS(gds.in_which_cryo(pt), gds.in_which_apa(pt));

    if (apa_gds == NULL) continue;
    
    std::pair<double, double> xmm = apa_gds->minmax(0); 
    
    if (pt.x>xmm.first && pt.x<xmm.second) continue;// space in between wires in an APA are treated as dead region
    
    int face = 0; // "A face -x"
    float drift_dist;
    
    if (TMath::Abs(pt.x-xmm.first) > TMath::Abs(pt.x-xmm.second)) {
      drift_dist = TMath::Abs(pt.x-xmm.second);
      face = 1; //"B face +x"
    }else{
      drift_dist = TMath::Abs(pt.x-xmm.first);
    }

    int tbin1 = int(drift_dist/2.0/unit_dis/units::mm + offset1);
    tbin1 = TMath::Abs(tbin1);

    //std::cout << tbin1 << " " << tbin << std::endl;

    if (tbin1 == tbin){

      if (apa_gds->contained_yz(pt)){

     	Point p1 = pt;// hack for now, gap in the tiling ... 
	sum += charge;
     	const GeomWire* wire_u = apa_gds->closest(p1, static_cast<WirePlaneType_t>(0),face);
     	const GeomWire* wire_v = apa_gds->closest(p1, static_cast<WirePlaneType_t>(1),face);
     	const GeomWire* wire_w = apa_gds->closest(p1, static_cast<WirePlaneType_t>(2),face);
	
     	if (wire_u!=0&&wire_v!=0&&wire_w!=0){
	  
     	  GeomWireSelection wires;
     	  wires.push_back(wire_u);
     	  wires.push_back(wire_v);
     	  wires.push_back(wire_w);
	
     	  const GeomCell* cell = tiling.cell(wires);
	  
     	  if (cell!=0){
	    
    	    //Point pp = cell->center();
    	    //	std::cout << pp.x << " " << pp.y << " " << pp.z << std::endl;
	    
	    
    	    if (cellchargemap.find(cell) == cellchargemap.end()){
    	      //not found
    	      cellchargemap[cell] = charge;
    	    }else{
    	      cellchargemap[cell] += charge;
    	    }
    	  }else{
    	    // Point pp = cell->center();
    	    // gds.avoid_gap(p1);
    	    // std::cout << p.x << " " << p.y << " " << p.z << " " << std::endl;
    	    //std::cout << itruth << " " << p1.x << " " << p1.y << " " << p1.z << " " << charge << std::endl;
    	    // std::cout << charge << std::endl;
    	  }
	}
      }
    }
  }
  
  //  std::cout << "Xin " << sum << std::endl;

}


WireCell2dToy::TruthToyTiling::TruthToyTiling(WireCell2dToy::ToyTiling& tiling, const WireCell::PointValueVector &pvv, const std::vector<int> &time_offset, int tbin, const DetectorGDS& gds, float unit_dis){
  float sum = 0;
  for (int itruth = 0; itruth < pvv.size(); ++itruth){
    const Point& pt = pvv[itruth].first; // get the point
    float charge = pvv[itruth].second;

    short which_cryo = gds.in_which_cryo(pt);
    short which_apa = gds.in_which_apa(pt);
    
    const WrappedGDS *apa_gds = gds.get_apaGDS(gds.in_which_cryo(pt), gds.in_which_apa(pt));

    if (apa_gds == NULL) continue;
    
    std::pair<double, double> xmm = apa_gds->minmax(0); 
    
    if (pt.x>xmm.first && pt.x<xmm.second) continue;// space in between wires in an APA are treated as dead region
    
    int face = 0; // "A face -x"
    float drift_dist;
    
    if (TMath::Abs(pt.x-xmm.first) > TMath::Abs(pt.x-xmm.second)) {
      drift_dist = TMath::Abs(pt.x-xmm.second);
      face = 1; //"B face +x"
    }else{
      drift_dist = TMath::Abs(pt.x-xmm.first);
    }

    int tbin1 = int(drift_dist/2.0/unit_dis/units::mm + time_offset[itruth]/4.);
    tbin1 = TMath::Abs(tbin1);

    if (tbin1 == tbin){
      
      
      if (apa_gds->contained_yz(pt)){
	
     	Point p1 = pt;// hack for now, gap in the tiling ... 
	sum += charge;
     	const GeomWire* wire_u = apa_gds->closest(p1, static_cast<WirePlaneType_t>(0),face);
     	const GeomWire* wire_v = apa_gds->closest(p1, static_cast<WirePlaneType_t>(1),face);
     	const GeomWire* wire_w = apa_gds->closest(p1, static_cast<WirePlaneType_t>(2),face);
	
     	if (wire_u!=0&&wire_v!=0&&wire_w!=0){
	  
     	  GeomWireSelection wires;
     	  wires.push_back(wire_u);
     	  wires.push_back(wire_v);
     	  wires.push_back(wire_w);
	  
     	  const GeomCell* cell = tiling.cell(wires);
	  
     	  if (cell!=0){
	    
    	    //Point pp = cell->center();
    	    //	std::cout << pp.x << " " << pp.y << " " << pp.z << std::endl;
	    
	    
    	    if (cellchargemap.find(cell) == cellchargemap.end()){
    	      //not found
    	      cellchargemap[cell] = charge;
    	    }else{
    	      cellchargemap[cell] += charge;
    	    }
    	  }else{
    	    // Point pp = cell->center();
    	    // gds.avoid_gap(p1);
    	    // std::cout << p.x << " " << p.y << " " << p.z << " " << std::endl;
    	    //std::cout << itruth << " " << p1.x << " " << p1.y << " " << p1.z << " " << charge << std::endl;
    	    // std::cout << charge << std::endl;
    	  }
	}
      }
            
    }
  }
}



ClassImp(WireCell2dToy::TruthToyTiling);

