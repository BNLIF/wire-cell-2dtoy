#include "WireCell2dToy/ToyMatrix.h"

using namespace WireCell;

WireCell2dToy::ToyMatrix::ToyMatrix(WireCell2dToy::ToyTiling& toytiling, WireCell2dToy::MergeToyTiling& mergetiling){
  solve_flag = 0;
  chi2 = -1;
  
  // build up index
  Buildup_index(mergetiling);
  
  // std::cout << mcindex << " " << mwindex << " " << swindex << std::endl;

  //Construct Matrix
  MA = new TMatrixD(mwindex,mcindex);
  MB = new TMatrixD(mwindex,swindex);
  MAT = new TMatrixD(mcindex,mwindex);
  MBT = new TMatrixD(swindex,mwindex);
  
  Vy = new TMatrixD(swindex,swindex);
  VBy = new TMatrixD(mwindex,mwindex);
  VBy_inv = new TMatrixD(mwindex,mwindex);
  Vx = new TMatrixD(mcindex,mcindex);
  Vx_inv = new TMatrixD(mcindex,mcindex);

  MC = new TMatrixD(mcindex,mcindex);
  MC_inv = new TMatrixD(mcindex,mcindex);

  //Construct Vector
  Wy = new TVectorD(swindex);
  Cx = new TVectorD(mcindex);
  dCx = new TVectorD(mcindex);
  
  GeomWireSelection allwire = toytiling.get_allwire();
  WireChargeMap wcmap = toytiling.wcmap();
  for (int j=0;j!=allwire.size();j++){
    int index = swimap[allwire[j]];
    float charge = wcmap[allwire[j]];
    (*Wy)[index] =charge;
    // fill covariance matrix as well
    // Vy ...
    // need to differentiate which plane it is
    WirePlaneType_t plane = allwire[j]->plane();
    Double_t charge_noise;
    if (plane == 0){
      charge_noise = 14000*0.05;
    }else if (plane==1){
      charge_noise = 14000*0.03;
    }else if (plane==2){
      charge_noise = 14000*0.02;
    }
    (*Vy)(index,index) = charge_noise*charge_noise + 0.05*0.05 * charge*charge;
  }

  
  GeomWireSelection allmwire = mergetiling.get_allwire();
  for (int j=0;j!=allmwire.size();j++){
    int index = mwimap[allmwire[j]];
    //construct MA
    for (int k=0; k!=mergetiling.cells(*allmwire[j]).size();k++){
      int index1 = mcimap[mergetiling.cells(*allmwire[j]).at(k)];
      (*MA)(index,index1) = 1;
    }
    
    //construct MB
    for (int k=0;k!=((MergeGeomWire*)allmwire[j])->get_allwire().size();k++){
      int index1 = swimap[((MergeGeomWire*)allmwire[j])->get_allwire().at(k)];
      (*MB)(index,index1) = 1;
    }
  }

  // // construct the rest of matrix
  MBT->Transpose(*MB);
  MAT->Transpose(*MA);

  *VBy = (*MB) * (*Vy) * (*MBT);

  *VBy_inv = *VBy;
  VBy_inv->Invert();

  *MC = (*MAT) * (*VBy_inv) * (*MA);

  
}

WireCell2dToy::ToyMatrix::~ToyMatrix(){
  
  
  delete MA, MB, MAT, MBT;
  delete Vy, VBy, Vx, VBy_inv, Vx_inv;
  delete MC, MC_inv;
  delete Wy, Cx, dCx;
}



void WireCell2dToy::ToyMatrix::Buildup_index(WireCell2dToy::MergeToyTiling& mergetiling){

  mcindex = 0;
  mwindex = 0;
  swindex = 0;

  GeomCellSelection allmcell = mergetiling.get_allcell();
 
  for (int j=0;j!=allmcell.size();j++){
    //construct merged cell index
    const MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
    if (mcimap.find(mcell) == mcimap.end()){
      mcimap[mcell] = mcindex;
      mcindex ++;
      
      const GeomWireSelection wires = mergetiling.wires(*allmcell[j]);
      // cout << wires.size() << endl;
      for (int k=0;k!=wires.size();k++){
	//construct merged wire index
	const MergeGeomWire *mwire = (MergeGeomWire*)wires[k];
	if (mwimap.find(mwire) == mwimap.end()){
	  mwimap[mwire] = mwindex;
	  mwindex ++;
	  
	  //construct single wire index
	  GeomWireSelection swires = mwire->get_allwire();
	  for (int kk = 0; kk!=swires.size(); kk++){
	    const GeomWire* wire1 = swires[kk];
	    if (swimap.find(wire1) == swimap.end()){
	      swimap[wire1] = swindex;
	      swindex ++;
	    }
	  }
	  
	}
      }
      
    }
  }
  
}
