#include "WireCell2dToy/LowmemTiling.h"

using namespace WireCell;

WireCell2dToy::LowmemTiling::LowmemTiling(int time_slice, int nrebin, const WireCell::Slice& slice,WireCell::GeomDataSource& gds,
					  WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map,
					  std::vector<float>& uplane_rms, std::vector<float>& vplane_rms, std::vector<float>& wplane_rms)
  : gds(gds)
  , nrebin(nrebin)
  , time_slice(time_slice)
{
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  // form bad wires group
  form_bad_merge_wires(uplane_map, vplane_map, wplane_map);
  
  // form good wires group
  form_fired_merge_wires(slice);

  // create two bad wire cells // these are special ones ... 
  
  // create three good wire cells & two good wire + one bad wire cells
  
  
  
}

void WireCell2dToy::LowmemTiling::form_fired_merge_wires(const WireCell::Slice& slice){
  WireCell::Channel::Group group = slice.group();
  
  
  //double sum = 0;
  for (int i=0;i!=group.size();i++){
    const WireCell::GeomWire *wire = gds.by_channel_segment(group.at(i).first,0);
    float charge = group.at(i).second;
    if (charge < 11) charge = 11;
    if (wirechargemap.find(wire) == wirechargemap.end()){
      //not found
      wirechargemap[wire] = charge;
    }else{
      //wirechargemap[wire] += charge;
    }
  }
  
  // do U
  MergeGeomWire *mwire = 0;
  int ident = 0;
  int last_wire = -1;
  for (int i = 0; i!= nwire_u; i++){
    const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)0,i);
    if (wirechargemap.find(wire)!=wirechargemap.end()){
      if (mwire == 0){
	mwire = new MergeGeomWire(ident,*wire);
	mwire->SetTimeSlice(time_slice);
	fired_wire_u.push_back(mwire);
	ident ++;
      }else{
	if (i==last_wire+1){
	  mwire->AddWire(*wire);
	}else{
	  // create a new wire
	  mwire = new MergeGeomWire(ident,*wire);
	  mwire->SetTimeSlice(time_slice);
	  fired_wire_u.push_back(mwire);
	  ident ++;
	}
      }
      last_wire = i;
    }
  }
  
  // V
  mwire = 0;
  last_wire = -1;
  for (int i = 0; i!= nwire_v; i++){
    const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)1,i);
    if (wirechargemap.find(wire)!=wirechargemap.end()){
      if (mwire == 0){
	mwire = new MergeGeomWire(ident,*wire);
	mwire->SetTimeSlice(time_slice);
	fired_wire_v.push_back(mwire);
	ident ++;
      }else{
	if (i==last_wire+1){
	  mwire->AddWire(*wire);
	}else{
	  // create a new wire
	  mwire = new MergeGeomWire(ident,*wire);
	  mwire->SetTimeSlice(time_slice);
	  fired_wire_v.push_back(mwire);
	  ident ++;
	}
      }
      last_wire = i;
    }
  }
  
  // W
  mwire = 0;
  last_wire = -1;
  for (int i = 0; i!= nwire_w; i++){
    const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)2,i);
    if (wirechargemap.find(wire)!=wirechargemap.end()){
      if (mwire == 0){
	mwire = new MergeGeomWire(ident,*wire);
	mwire->SetTimeSlice(time_slice);
	fired_wire_w.push_back(mwire);
	ident ++;
      }else{
	if (i==last_wire+1){
	  mwire->AddWire(*wire);
	}else{
	  // create a new wire
	  mwire = new MergeGeomWire(ident,*wire);
	  mwire->SetTimeSlice(time_slice);
	  fired_wire_w.push_back(mwire);
	  ident ++;
	}
      }
      last_wire = i;
    }
  }

  
  // int sum = 0;
  // for (int i=0;i!=fired_wire_u.size();i++){
  //   MergeGeomWire *mwire = (MergeGeomWire*)fired_wire_u.at(i);
  //   sum += mwire->get_allwire().size();
  //   for (Int_t j=0;j!=mwire->get_allwire().size();j++){
  //     std::cout << i << " " << mwire->get_allwire().at(j)->index() << std::endl;
  //   }
  // }
  // for (int i=0;i!=fired_wire_v.size();i++){
  //   MergeGeomWire *mwire = (MergeGeomWire*)fired_wire_v.at(i);
  //   sum += mwire->get_allwire().size();
  // }
  // for (int i=0;i!=fired_wire_w.size();i++){
  //   MergeGeomWire *mwire = (MergeGeomWire*)fired_wire_w.at(i);
  //   sum += mwire->get_allwire().size();
  // }

  
  //std::cout << fired_wire_u.size() << " " << fired_wire_v.size() << " " << fired_wire_w.size() << " " << group.size() << " " << sum << std::endl;
  
}

void WireCell2dToy::LowmemTiling::form_bad_merge_wires(WireCell::ChirpMap& uplane_map, WireCell::ChirpMap& vplane_map, WireCell::ChirpMap& wplane_map){
  
  // do U
  MergeGeomWire *mwire = 0;
  int ident = 0;
  int last_wire = -1;

  for (int i = 0; i!= nwire_u; i++){
    if (uplane_map.find(i)!=uplane_map.end()){
      if (time_slice >= uplane_map[i].first/nrebin && 
	  time_slice <= uplane_map[i].second/nrebin){
	const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)0,i);
	// first one 
	if (mwire == 0){
	  mwire = new MergeGeomWire(ident,*wire);
	  mwire->SetTimeSlice(time_slice);
	  bad_wire_u.push_back(mwire);
	  ident ++;
	}else{
	  if (i==last_wire+1){
	    mwire->AddWire(*wire);
	  }else{
	    // create a new wire
	    mwire = new MergeGeomWire(ident,*wire);
	    mwire->SetTimeSlice(time_slice);
	    bad_wire_u.push_back(mwire);
	    ident ++;
	  }
	}
	last_wire = i;
      }
    }
  }

  // now do V plane
  mwire = 0;
  last_wire = -1;

  for (int i = 0; i!= nwire_v; i++){
    if (vplane_map.find(i)!=vplane_map.end()){
      if (time_slice >= vplane_map[i].first/nrebin && 
	  time_slice <= vplane_map[i].second/nrebin){
	const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)1,i);
	// first one 
	if (mwire == 0){
	  mwire = new MergeGeomWire(ident,*wire);
	  mwire->SetTimeSlice(time_slice);
	  bad_wire_v.push_back(mwire);
	  ident ++;
	}else{
	  if (i==last_wire+1){
	    mwire->AddWire(*wire);
	  }else{
	    // create a new wire
	    mwire = new MergeGeomWire(ident,*wire);
	    mwire->SetTimeSlice(time_slice);
	    bad_wire_v.push_back(mwire);
	    ident ++;
	  }
	}
	last_wire = i;
      }
    }
  }
  
  // W plane
  mwire = 0;
  last_wire = -1;
  
  for (int i = 0; i!= nwire_w; i++){
    if (wplane_map.find(i)!=wplane_map.end()){
      if (time_slice >= wplane_map[i].first/nrebin && 
	  time_slice <= wplane_map[i].second/nrebin){
	const GeomWire *wire = gds.by_planeindex((WirePlaneType_t)2,i);
	// first one 
	if (mwire == 0){
	  mwire = new MergeGeomWire(ident,*wire);
	  mwire->SetTimeSlice(time_slice);
	  bad_wire_w.push_back(mwire);
	  ident ++;
	}else{
	  if (i==last_wire+1){
	    mwire->AddWire(*wire);
	  }else{
	    // create a new wire
	    mwire = new MergeGeomWire(ident,*wire);
	    mwire->SetTimeSlice(time_slice);
	    bad_wire_w.push_back(mwire);
	    ident ++;
	  }
	}
	last_wire = i;
      }
    }
  }
  


  // int sum = 0;
  // for (int i=0;i!=bad_wire_w.size();i++){
  //   MergeGeomWire *mwire = (MergeGeomWire*)bad_wire_w.at(i);
  //   sum += mwire->get_allwire().size();
  //   for (int j=0;j!=mwire->get_allwire().size();j++){
  //     GeomWire *wire = (GeomWire*)mwire->get_allwire().at(j);
  //     std::cout << i << " " << wire->index() << std::endl;
  //   }
  // }
  // std::cout << uplane_map.size() << " " << bad_wire_u.size() << " " << sum << std::endl;
  //std::cout << vplane_map.size() << " " << bad_wire_v.size() << " " << sum << std::endl;
  //std::cout << wplane_map.size() << " " << bad_wire_w.size() << " " << sum << std::endl;

}

WireCell2dToy::LowmemTiling::~LowmemTiling(){
  for (int i=0;i!=bad_wire_u.size();i++){
    delete bad_wire_u.at(i);
  }
  bad_wire_u.clear();
  
  for (int i=0;i!=bad_wire_v.size();i++){
    delete bad_wire_v.at(i);
  }
  bad_wire_v.clear();

  for (int i=0;i!=bad_wire_w.size();i++){
    delete bad_wire_w.at(i);
  }
  bad_wire_w.clear();
  
  for (int i=0;i!=fired_wire_u.size();i++){
    delete fired_wire_u.at(i);
  }
  fired_wire_u.clear();
  
  for (int i=0;i!=fired_wire_v.size();i++){
    delete fired_wire_v.at(i);
  }
  fired_wire_v.clear();

  for (int i=0;i!=fired_wire_w.size();i++){
    delete fired_wire_w.at(i);
  }
  fired_wire_w.clear();
}

WireCell::GeomWireSelection WireCell2dToy::LowmemTiling::wires(const WireCell::GeomCell& cell) const{
  GeomWireSelection wires;
  return wires;
}
WireCell::GeomCellSelection WireCell2dToy::LowmemTiling::cells(const WireCell::GeomWire& wire) const{
  GeomCellSelection cells;
  return cells;
}

const WireCell::GeomCell* WireCell2dToy::LowmemTiling::cell(const WireCell::GeomWireSelection& wires) const{
  return 0;
}
