#include "WireCell2dToy/MatrixSolver.h"

using namespace WireCell;

WireCell2dToy::MatrixSolver::MatrixSolver(GeomCellSelection& allmcell, GeomWireSelection& mwires, WireCell::GeomCellMap& cell_wire_map, WireCell::GeomWireMap& wire_cell_map, WireCell::WireChargeMap& wire_charge_map, WireCell::WireChargeMap& wire_charge_error_map)
{
  //
  solve_flag = -1;
  
  // build up the index
  mcindex = 0;
  mwindex = 0;
  swindex = 0;

  for (int j=0;j!=allmcell.size();j++){
    const MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
    if (mcimap.find(mcell) == mcimap.end()){
      mcimap[mcell] = mcindex;
      mcindex ++;
      const GeomWireSelection wires = cell_wire_map[mcell];
      // std::cout << wires.size() << std::endl;
      for (int k=0;k!=wires.size();k++){
       	//construct merged wire index
       	const MergeGeomWire *mwire = (MergeGeomWire*)wires[k];
       	// require all the wire must be good to be used in the matrix solving 
       	if (mwimap.find(mwire) == mwimap.end() 
	    && wire_charge_map[mwire] >10){
       	  //if (mwimap.find(mwire) == mwimap.end() ){
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

  //  std::cout << swindex << " " << mwindex << " " << mcindex << std::endl;

  if (mcindex >0 ){
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
    
    MWy = new TVectorD(mwindex);
    MWy_pred = new TVectorD(mwindex);
    
    Cx = new TVectorD(mcindex);
    dCx = new TVectorD(mcindex);
    
    // fill in the matrix content ...
    for ( auto it = wire_charge_map.begin(); it != wire_charge_map.end(); it++){
      if (swimap.find(it->first) != swimap.end()){
	int index = swimap[it->first];
	float charge = wire_charge_map[it->first];
	float charge_err = wire_charge_error_map[it->first];
	(*Wy)[index] =charge;
	(*Vy)(index,index) = charge_err*charge_err/1e6;
	//std::cout << index << " " << charge << " " << charge_err << std::endl;
      }
    }
    

    for (int j=0;j!=mwires.size();j++){
      if (mwimap.find(mwires[j])!=mwimap.end()){
       	int index = mwimap[mwires[j]];
       	//construct MA
       	for (int k=0; k!=wire_cell_map[mwires[j]].size();k++){
	  int index1 = mcimap[wire_cell_map[mwires[j]].at(k)];
       	  (*MA)(index,index1) = 1;
       	}
	
       	//construct MB
       	for (int k=0;k!=((MergeGeomWire*)mwires[j])->get_allwire().size();k++){
       	  int index1 = swimap[((MergeGeomWire*)mwires[j])->get_allwire().at(k)];
      	  (*MB)(index,index1) = 1;
       	}
      }
    }

    //MA->Print();
    // // construct the rest of matrix
    MBT->Transpose(*MB);
    MAT->Transpose(*MA);
    
    *VBy = (*MB) * (*Vy) * (*MBT);
    
    //    std::cout << "fill in the data" << std::endl;
    
    *VBy_inv = *VBy;
    //    if (fabs(VBy->Determinant())<0.01) {
    // std::cout << "Problem " << VBy->Determinant() << std::endl;
    //}
    VBy_inv->Invert();
    
    //std::cout << "invert the matrix " << std::endl; 
    *MC = (*MAT) * (*VBy_inv) * (*MA);

    solve_flag = 0;
    double det;
    if (mcindex <= mwindex){
      det = MC->Determinant();
    }else{
      det = 0;
    }

    if (det > 1e-5){
      Direct_Solve();
    }else{
      // first time L1 solve ...
      L1_Solve();
    }
    
    //std::cout << det << std::endl;
    
  }
  
}
void WireCell2dToy::MatrixSolver::L1_Solve(){
  
  
  solve_flag = 2;
}


void WireCell2dToy::MatrixSolver::Direct_Solve(){
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
  direct_chi2 = sol * (sol1)/1e6;
  
  direct_ndf = mwindex - mcindex;
  
  // std::cout << chi2 << " " << ndf << std::endl;
  // for (int i=0;i!=mcindex;i++){
  //   std::cout << (*Cx)[i] << " " << (*dCx)[i] << std::endl;
  // }

  // MA->Print();
  *MWy = (*MB) * (*Wy);
  *MWy_pred = (*MA)*(*Cx);
  // for (int i=0;i!=mwindex;i++){
  //   std::cout << (*MWy)[i] << std::endl;
  // }
  solve_flag = 1;
}


WireCell2dToy::MatrixSolver::~MatrixSolver(){
  delete MA;
  delete MB;
  delete MAT;
  delete MBT;
  delete Vy;
  delete VBy;
  delete Vx;
  delete VBy_inv;
  delete Vx_inv;
  delete MC;
  delete MC_inv;
  delete Wy;
  delete Cx;
  delete dCx;
  delete MWy_pred;
  delete MWy;
}
					   
