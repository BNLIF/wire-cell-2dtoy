#include "WireCell2dToy/ToyMatrix.h"

using namespace WireCell;

WireCell2dToy::ToyMatrix::ToyMatrix(WireCell2dToy::ToyTiling& toytiling, WireCell2dToy::MergeToyTiling& mergetiling){
  solve_flag = -1;
  chi2 = -1;
  
  // build up index
  Buildup_index(mergetiling);
  
  //std::cout << mcindex << " " << mwindex << " " << swindex << std::endl;
  if (mcindex >0){
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
      if (swimap.find(allwire[j])!=swimap.end()){
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
	(*Vy)(index,index) = (charge_noise*charge_noise + 0.05*0.05 * charge*charge)/1e6;
      }
    }
    
    
    GeomWireSelection allmwire = mergetiling.get_allwire();
    for (int j=0;j!=allmwire.size();j++){
      if (mwimap.find(allmwire[j])!=mwimap.end()){
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
    }
    
    
    // // construct the rest of matrix
    MBT->Transpose(*MB);
    MAT->Transpose(*MA);
    
    *VBy = (*MB) * (*Vy) * (*MBT);
    
    *VBy_inv = *VBy;
    //    if (fabs(VBy->Determinant())<0.01) {
    // std::cout << "Problem " << VBy->Determinant() << std::endl;
    //}
    VBy_inv->Invert();
    
    *MC = (*MAT) * (*VBy_inv) * (*MA);
    solve_flag = 0;
    Solve();
  }
  // MA->Print();
  // MAT->Print();
  // MB->Print();
  // MBT->Print();
  // Vy->Print();

  
}

WireCell2dToy::ToyMatrix::~ToyMatrix(){
  
  
  delete MA, MB, MAT, MBT;
  delete Vy, VBy, Vx, VBy_inv, Vx_inv;
  delete MC, MC_inv;
  delete Wy, Cx, dCx;
}

int WireCell2dToy::ToyMatrix::Solve(){
  
  double det = MC->Determinant();

  //  std::cout << det << std::endl;
  
  if (fabs(det)>1e-5){
    *MC_inv = *MC;
    MC_inv->Invert();
    *Cx = (*MC_inv) * (*MAT) * (*VBy_inv) * (*MB) * (*Wy);
    
    *Vx_inv = (*MAT) * (*VBy_inv) * (*MA);
    *Vx = *Vx_inv;
    Vx->Invert();
    
    for (int i=0;i!=mcindex;i++){
      (*dCx)[i] = sqrt( (*Vx)(i,i)) * 1000.;
    }
    
    TVectorD sol = (*MB) * (*Wy) - (*MA) * (*Cx);
    TVectorD sol1 =  (*VBy_inv) * sol;
    chi2 = sol * (sol1)/1e6;
    
    //std::cout << chi2 << std::endl;
    //      for (int i=0;i!=mcindex;i++){
    //	std::cout << (*Cx)[i] << " " << (*dCx)[i]*1000. << std::endl;
    //}
    
    solve_flag = 1;
  }
  
  return solve_flag;
}

double WireCell2dToy::ToyMatrix::Get_Cell_Charge( const WireCell::GeomCell *cell, int flag )  {
  // flag == 1 charge
  // flag == 2 uncertainty
  int index = mcimap[cell];
  if (flag==1){
    return (*Cx)[index];
  }else if (flag==2){
    return (*dCx)[index];
  }
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
