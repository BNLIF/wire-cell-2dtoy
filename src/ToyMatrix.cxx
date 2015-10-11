#include "WireCell2dToy/ToyMatrix.h"

#include "TDecompSVD.h"

using namespace WireCell;

void WireCell2dToy::ToyMatrix::JudgeSimpleBlob(WireCell2dToy::ToyTiling& toytiling, WireCell2dToy::MergeToyTiling& mergetiling){
   GeomCellSelection allmcell = mergetiling.get_allcell();
   num_blob = 0;
   for (int j=0;j!=allmcell.size();j++){
     MergeGeomCell *mcell =(MergeGeomCell*)allmcell.at(j);
     
     double charge =Get_Cell_Charge(mcell);
     if (charge>recon_threshold){
       mcell->FindEdges();
       // std::cout << mcell->get_allcell().size() << " " << mcell->get_edgecells().size() << std::endl;
       int nwire = 0, max_wire = 0;
       int ncell = 0;
       if(mcell->IsBlob()) {
	 mcell->FindCorners(toytiling.cmap(), toytiling.wmap());
	 num_blob ++; 
	 // GeomWireSelection n_mwires = mergetiling[i]->wires(*mcell);
	 // for (int k=0;k!=n_mwires.size();k++){
	 //   int ncells = 0;
	 //   for (int kk= 0; kk!=mergetiling[i]->cells(*n_mwires.at(k)).size();kk++){
	 //     if (toymatrix[i]->Get_Cell_Charge(mergetiling[i]->cells(*n_mwires.at(k)).at(kk))>recon_threshold){
	 // 	ncells ++;
	 //     }
	 //   }
	 //   std::cout << ncells << std::endl;
	 
	 // find simple blob ... 
	 // find the merged wire regarding the merged blob
	 GeomWireSelection wires = mergetiling.wires(*mcell);
	 // find the number of wires 
	 for (int k =0; k!=wires.size();k++){
	   MergeGeomWire* mwire = (MergeGeomWire*)wires.at(k);
	   //std::cout << mwire->get_allwire().size() << std::endl;
	   nwire += mwire->get_allwire().size();
	   if (mwire->get_allwire().size() > max_wire) 
	     max_wire = mwire->get_allwire().size();
	   GeomCellSelection cells = mergetiling.cells(*mwire);
	   for (int kk=0; kk!= cells.size(); kk++){
	     if (cells.at(kk) != mcell && Get_Cell_Charge(cells.at(kk))>recon_threshold){
	       ncell += ((MergeGeomCell*)cells.at(kk))->get_allcell().size();
	     }
	   }
	 }
	 nwire -= max_wire;
	 // for the merged blob passed the cut, find the number of cells
	 //std::cout << "Xin: " << nwire << " " << max_wire << " " << ncell << std::endl;
	 if (ncell < nwire){
	   mcell->SetSimpleBlob(true);
	   simple_blob_reduction = true;
	 }else{
	   mcell->SetSimpleBlob(false);
	 }
       }
     }
     
     
    // 	GeomCellSelection corners = mcell->get_cornercells();
    // 	total_corner_cells.insert(total_corner_cells.end(),corners.begin(),corners.end());
    // 	for (int k=0;k!=mcell->get_allcell().size();k++){
    // 	  total_recon_cells.push_back(mcell->get_allcell().at(k));
    // 	  //get charge
    // 	  double sc_charge = 0;
    // 	  GeomCellMap scmap = toytiling[i]->cmap();
    // 	  WireChargeMap wcmap = toytiling[i]->wcmap();
    // 	  GeomWireSelection wires = scmap[mcell->get_allcell().at(k)];
    // 	  double aa[3];
    // 	  for (int q=0;q!=wires.size();q++){
    // 	    sc_charge += wcmap[wires.at(q)];
    // 	    aa[q] = wcmap[wires.at(q)];
    // 	  }
    // 	  total_scmap[mcell->get_allcell().at(k)] = sc_charge/3;
    // 	  total_scrms[mcell->get_allcell().at(k)] = rms(aa[0],aa[1],aa[2])*3./sc_charge;
    // 	  if (sc_charge/3 > charge_max) charge_max = sc_charge/3;
    // 	  if (sc_charge/3 < charge_min) charge_min = sc_charge/3;
    // 	}
    //   }
   }

    // toymatrix[i]->Set_blob(num_blob);
    
    // cout << "# of blobs " << toymatrix[i]->Get_blob() << endl;

}

WireCell2dToy::ToyMatrix::ToyMatrix(){
  solve_flag = -1;
  chi2 = -1;
  svd_flag = 0;
  num_blob = 0;
  simple_blob_reduction = false;
  recon_threshold = 2000;
  
  mwindex = 1;
  mcindex = 1;

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

}


WireCell2dToy::ToyMatrix::ToyMatrix(WireCell2dToy::ToyTiling& toytiling, WireCell2dToy::MergeToyTiling& mergetiling, int svd_flag1, int recon_t){
  solve_flag = -1;
  chi2 = -1;
  svd_flag = svd_flag1;
  num_blob = 0;
  simple_blob_reduction = false;
  recon_threshold = recon_t;
  
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
    
    MWy = new TVectorD(mwindex);
    MWy_pred = new TVectorD(mwindex);

    Cx = new TVectorD(mcindex);
    dCx = new TVectorD(mcindex);
    
    GeomWireSelection allwire = toytiling.get_allwire();
    WireChargeMap wcmap = toytiling.wcmap();
    WireChargeMap wcemap = toytiling.wcemap();
    
    for (int j=0;j!=allwire.size();j++){
      if (swimap.find(allwire[j])!=swimap.end()){
	int index = swimap[allwire[j]];
	float charge = wcmap[allwire[j]];
	float charge_err = wcemap[allwire[j]];
	(*Wy)[index] =charge;
	(*Vy)(index,index) = charge_err*charge_err/1e6;
	// fill covariance matrix as well
	// Vy ...
	// need to differentiate which plane it is
	// WirePlaneType_t plane = allwire[j]->plane();
	// Double_t charge_noise;
	// if (plane == 0){
	//   charge_noise = 14000*0.05;
	// }else if (plane==1){
	//   charge_noise = 14000*0.03;
	// }else if (plane==2){
	//   charge_noise = 14000*0.02;
	// }
	// (*Vy)(index,index) = (charge_noise*charge_noise + 0.05*0.05 * charge*charge)/1e6;
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
    //double det = MC->Determinant();
    double det;
    if (mcindex <= mwindex){
      det = MC->Determinant();
    }else{
      det = 0;
    }
    solve_flag = 0;
    if (svd_flag==0){
      if (det > 1e-5)
	Solve();
    }else{
      if (det > 1e-5){
	Solve();
      }else{
	Solve_SVD();
      }
      //Solve_SVD();
    }
  }
  //MA->Print();
  // MAT->Print();
  // MB->Print();
  // MBT->Print();
  // Vy->Print();

  
}



WireCell2dToy::ToyMatrix::ToyMatrix(WireCell2dToy::ToyTiling& toytiling, WireCell2dToy::MergeToyTiling& mergetiling, int recon_t){
  solve_flag = -1;
  chi2 = -1;
  svd_flag = 0;
  num_blob = 0;
  simple_blob_reduction = false;
  // build up index
  Buildup_index(mergetiling);
  recon_threshold = recon_t;
  
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
    
    MWy = new TVectorD(mwindex);
    MWy_pred = new TVectorD(mwindex);

    Cx = new TVectorD(mcindex);
    dCx = new TVectorD(mcindex);
    
    GeomWireSelection allwire = toytiling.get_allwire();
    WireChargeMap wcmap = toytiling.wcmap();
    WireChargeMap wcemap = toytiling.wcemap();
    for (int j=0;j!=allwire.size();j++){
      if (swimap.find(allwire[j])!=swimap.end()){
	int index = swimap[allwire[j]];
	float charge = wcmap[allwire[j]];
	float charge_err = wcemap[allwire[j]];
	(*Wy)[index] =charge;
	(*Vy)(index,index) = charge_err*charge_err/1e6;
	// fill covariance matrix as well
	// Vy ...
	// need to differentiate which plane it is
	// WirePlaneType_t plane = allwire[j]->plane();
	// Double_t charge_noise;
	// if (plane == 0){
	//   charge_noise = 14000*0.05;
	// }else if (plane==1){
	//   charge_noise = 14000*0.03;
	// }else if (plane==2){
	//   charge_noise = 14000*0.02;
	// }
	// (*Vy)(index,index) = (charge_noise*charge_noise + 0.05*0.05 * charge*charge)/1e6;
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
      Solve();
    }
  }
  //MA->Print();
  // MAT->Print();
  // MB->Print();
  // MBT->Print();
  // Vy->Print();

  
}



void WireCell2dToy::ToyMatrix::Update_pred(){
  *MWy = (*MB) * (*Wy);
  *MWy_pred = (*MA) * (*Cx);
}

void WireCell2dToy::ToyMatrix::Print(){

  MA->Print();
  std::cout << std::endl;
  Cx->Print();

  for (int i=0;i!=mwindex;i++){
    std::cout << "wire: " << i << " " << (*MWy)[i] << " " << (*MWy_pred)[i] << " " << sqrt((*VBy)[i][i]) << std::endl;
  }
}

double WireCell2dToy::ToyMatrix::Get_residual(const WireCell::GeomCell *cell){
  double res=0;
  int index = mcimap[cell];

  for (int i=0;i!=mwindex;i++){
    res += (*MA)(i,index) * fabs((*MWy)[i] - (*MWy_pred)[i]);
  }
  
  return res;
}

WireCell2dToy::ToyMatrix::~ToyMatrix(){
  
  
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

int WireCell2dToy::ToyMatrix::Solve(){
  
  // double det = MC->Determinant();

  //  std::cout << det << std::endl;
  
  //if (fabs(det)>1e-5){
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
    
    ndf = mwindex - mcindex;
    
    //std::cout << chi2 << std::endl;
    //      for (int i=0;i!=mcindex;i++){
    //	std::cout << (*Cx)[i] << " " << (*dCx)[i]*1000. << std::endl;
    //}
    
    solve_flag = 1;
    //}
  
  return solve_flag;
}

int WireCell2dToy::ToyMatrix::Solve_SVD(){
  
  // if (svd_removed.size()!=0){
  //   //update *MC
  //   int n_removed = svd_removed.size();
  //   TMatrixD MA2(mwindex,mcindex - n_removed);
  //   TMatrixD MA2T(mcindex - n_removed,mwindex);
    
  //   int index2 = 0;
  //   for (int k=0;k!=mcindex;k++){ //loop all possibility
  //     auto it3 = find(svd_removed.begin(),svd_removed.end(),k);
  //     if (it3 != svd_removed.end()){ // not in removed
  // 	for (int j=0;j!=mwindex;j++){
  // 	  MA2(j,index2) = (*MA)(j,k);
  // 	}
  // 	index2 ++;
  //     }
  //   }
  //   MA2T.Transpose(MA2);
  //   TMatrixD MC2(mcindex-svd_removed.size(),mcindex-svd_removed.size());
  //   TMatrixD MC2_inv(mcindex-svd_removed.size(),mcindex-svd_removed.size());

  //   MC2 = MA2T * (*VBy_inv) * MA2;
  //   TDecompSVD svd(MC2);
  //   MC2_inv = svd.Invert();
     
        
  //   TVectorD Cxt(mcindex-n_removed);
  //   TVectorD dCxt(mcindex-n_removed);
  //   TMatrixD Vxt(mcindex-n_removed,mcindex-n_removed);
  //   TMatrixD Vxt_inv(mcindex-n_removed,mcindex-n_removed);
    
  //   Cxt = (MC2_inv) * (MA2T) * (*VBy_inv) * (*MB) * (*Wy);
  //   Vxt_inv = (MA2T) * (*VBy_inv) * (MA2);
  //   Vxt = Vxt_inv;
  //   Vxt.Invert();
    
  //   for (int i=0;i!=mcindex-n_removed;i++){
  //     dCxt[i] = sqrt( Vxt(i,i)) * 1000.;
  //   }

  //   TVectorD sol = (*MB) * (*Wy) - (MA2) * (Cxt);
  //   TVectorD sol1 =  (*VBy_inv) * sol;
  //   chi2 = sol * (sol1)/1e6;
  //   ndf = mwindex - mcindex + n_removed;
  //   for (int i=0;i!=mcindex-n_removed;i++){
  //     if (Cxt[i] <0){
  // 	chi2 += 10*pow(Cxt[i]/dCxt[i],2);
  //     }
  //   }

  //   int index = 0;
  //   for (int i=0;i!=mcindex;i++){
  //     auto it = find(svd_removed.begin(),svd_removed.end(),i);
  //     if (it==svd_removed.end()){
  // 	(*Cx)[i] = Cxt[index];
  // 	(*dCx)[i] = dCxt[index];
  // 	index ++;
  //     }else{
  // 	(*Cx)[i] = 0.;
  // 	(*dCx)[i] = 0.;
  //     }
  //   }


  // }else{
  TDecompSVD svd(*MC);
  *MC_inv = svd.Invert();
  *Cx = (*MC_inv) * (*MAT) * (*VBy_inv) * (*MB) * (*Wy);
  *Vx_inv = (*MAT) * (*VBy_inv) * (*MA);
  *Vx = *Vx_inv;
  Vx->Invert();
  for (int i=0;i!=mcindex;i++){
    (*dCx)[i] = sqrt( (*Vx)(i,i)) * 1000.;
    if ((*Cx)[i]<0) svd_removed.push_back(i);
  }
  TVectorD sol = (*MB) * (*Wy) - (*MA) * (*Cx);
  TVectorD sol1 =  (*VBy_inv) * sol;
  //chi2 = sol * (sol1)/1e6;
  //ndf = mwindex - mcindex;
  solve_flag = 0;
  //}
    
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
	// require all the wire must be good to be used in the matrix solving 
	if (mwimap.find(mwire) == mwimap.end() && mergetiling.wcmap()[mwire] >10){
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
  
  //std::cout << mcindex << " " << mwindex << " " << swindex << std::endl;
}


