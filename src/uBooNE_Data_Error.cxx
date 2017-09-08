#include "WireCell2dToy/uBooNE_Data_Error.h"
#include "TGraph.h"

using namespace WireCell;

WireCell2dToy::uBooNEDataError::uBooNEDataError(const WireCell::GeomDataSource& gds, TH2I *hu_decon, TH2I *hv_decon, TH2I *hw_decon, int eve_num, int nrebin)
  : gds(gds)
  , nrebin(nrebin)
{
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();


#include "charge_err.txt"
  TGraph *gu = new TGraph(201,xu,yu);
  TGraph *gv = new TGraph(201,xv,yv);
  TGraph *gw = new TGraph(201,xw,yw);

  
  frame.clear();		// win or lose, we start anew

  frame.index =eve_num;
  bins_per_frame = hu_decon->GetNbinsY();

  int nbin = bins_per_frame;

  double fudge_factor = 1.1;//
  double fudge_factor_ind = 2.1;
  // U plane
  for (size_t ind=0; ind < hu_decon->GetNbinsX(); ++ind) {

    // do ROI ...
    std::vector<bool> signalsBool;
    signalsBool.resize(nbin, false);
    for (int i = 0; i < nbin; ++ i){
      if (hu_decon->GetBinContent(ind+1,i+1)!=0)
	signalsBool.at(i) = true;
    }
    std::vector<std::vector<int>> rois;
    bool inside = false;
    for (int i=0; i<nbin; ++i) {
      if (inside) {
	if (signalsBool[i]) { // still inside
	  rois.back().push_back(i);
	}else{
	  inside = false;
	}
      }
      else {                  // outside the Rio
	if (signalsBool[i]) { // just entered ROI
	  std::vector<int> roi;
	  roi.push_back(i);
	  rois.push_back(roi);
	  inside = true;
	}
      }
    }
    //    std::cout << rois.size() << std::endl;


    
    
    WireCell::Trace trace;
    trace.chid = ind;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    for (int i=0;i!=rois.size();i++){
      int time = rois.at(i).size() * nrebin;

      if (time <0) time = 12;
      if (time > 800) {
	std::cout << "too long time U " << ind << " " << time << std::endl;
	time = 800;
      }
      for (int j=0; j!=rois.at(i).size();j++){
       	trace.charge.at(rois.at(i).at(j)) = gu->Eval(time) * fudge_factor *fudge_factor_ind * nrebin / 4.;
	// if (trace.charge.at(rois.at(i).at(j))<0)
	//   std::cout << ind << " " << time << " " << gu->Eval(time) << std::endl; 
      }
      //      std::cout << ind << " " << time << std::endl;
    }

    // for (int ibin=0; ibin != bins_per_frame; ibin++) {
    //   if (hu_decon->GetBinContent(ind+1,ibin+1)!=0){
    // 	std::cout << ind << " " << trace.charge.at(ibin) << " " << hu_decon->GetBinContent(ind+1,ibin+1) << std::endl;
    //   }else{
    // 	if (trace.charge.at(ibin) !=0)
    // 	  std::cout << "wrong! " << std::endl;
    //   }
    //   //trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    // }

    
    frame.traces.push_back(trace);
    
  }
  
  // V plane
  for (size_t ind=0; ind < hv_decon->GetNbinsX(); ++ind) {
    // do ROI ...
    std::vector<bool> signalsBool;
    signalsBool.resize(nbin, false);
    for (int i = 0; i < nbin; ++ i){
      if (hv_decon->GetBinContent(ind+1,i+1)!=0)
	signalsBool.at(i) = true;
    }
    std::vector<std::vector<int>> rois;
    bool inside = false;
    for (int i=0; i<nbin; ++i) {
      if (inside) {
	if (signalsBool[i]) { // still inside
	  rois.back().push_back(i);
	}else{
	  inside = false;
	}
      }
      else {                  // outside the Rio
	if (signalsBool[i]) { // just entered ROI
	  std::vector<int> roi;
	  roi.push_back(i);
	  rois.push_back(roi);
	  inside = true;
	}
      }
    }



    WireCell::Trace trace;
    trace.chid = ind + nwire_u;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);

    for (int i=0;i!=rois.size();i++){
      int time = rois.at(i).size() * nrebin;

      if (time <0) time = 12;
      if (time > 800) {
	std::cout << "too long time V " << ind << " " << time << std::endl;
	time = 800;
      }
      
      for (int j=0; j!=rois.at(i).size();j++){
	trace.charge.at(rois.at(i).at(j)) = gv->Eval(time) * fudge_factor * fudge_factor_ind * nrebin / 4.; 
      }
      //      std::cout << ind << " " << time << std::endl;
    }
    
    // for (int ibin=0; ibin != bins_per_frame; ibin++) {
    //   if (hv_decon->GetBinContent(ind+1,ibin+1)!=0){
    // 	std::cout << ind << " " << trace.charge.at(ibin) << " " << hv_decon->GetBinContent(ind+1,ibin+1) << std::endl;
    //   }else{
    // 	if (trace.charge.at(ibin) !=0)
    // 	  std::cout << "wrong! " << std::endl;
    //   }
    //   //trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    // }
    frame.traces.push_back(trace);
  }

  // W plane
  for (size_t ind=0; ind < hw_decon->GetNbinsX(); ++ind) {
    // do ROI ...
    std::vector<bool> signalsBool;
    signalsBool.resize(nbin, false);
    for (int i = 0; i < nbin; ++ i){
      if (hw_decon->GetBinContent(ind+1,i+1)!=0)
	signalsBool.at(i) = true;
    }
    std::vector<std::vector<int>> rois;
    bool inside = false;
    for (int i=0; i<nbin; ++i) {
      if (inside) {
	if (signalsBool[i]) { // still inside
	  rois.back().push_back(i);
	}else{
	  inside = false;
	}
      }
      else {                  // outside the Rio
	if (signalsBool[i]) { // just entered ROI
	  std::vector<int> roi;
	  roi.push_back(i);
	  rois.push_back(roi);
	  inside = true;
	}
      }
    }

    WireCell::Trace trace;
    trace.chid = ind + nwire_u + nwire_v;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);

    for (int i=0;i!=rois.size();i++){
      int time = rois.at(i).size() * nrebin;

      if (time <0) time = 12;
      if (time > 800) {
	//std::cout << "too long time W " << ind << " " << time << std::endl;
	time = 800;
      }
      
      for (int j=0; j!=rois.at(i).size();j++){
	trace.charge.at(rois.at(i).at(j)) = gw->Eval(time) * fudge_factor * nrebin / 4.; 
      }
      //      std::cout << ind << " " << time << std::endl;
    }

    // for (int ibin=0; ibin != bins_per_frame; ibin++) {
    //   if (hw_decon->GetBinContent(ind+1,ibin+1)!=0){
    // 	std::cout << ind << " " << trace.charge.at(ibin) << " " << hw_decon->GetBinContent(ind+1,ibin+1) << std::endl;
    //   }else{
    // 	if (trace.charge.at(ibin) !=0)
    // 	  std::cout << "wrong! " << std::endl;
    //   }
    //   //trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    // }

    
    // for (int ibin=0; ibin != bins_per_frame; ibin++) {
    //   trace.charge.at(ibin) = hw_decon->GetBinContent(ind+1,ibin+1);
    // }
    frame.traces.push_back(trace);
  }

  //std::cout << frame.traces.size() << " " << bins_per_frame << std::endl;
}

void WireCell2dToy::uBooNEDataError::refresh(TH2F *hu_decon, TH2F *hv_decon, TH2F *hw_decon, int eve_num){
   GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();


#include "charge_err.txt"
  TGraph *gu = new TGraph(201,xu,yu);
  TGraph *gv = new TGraph(201,xv,yv);
  TGraph *gw = new TGraph(201,xw,yw);

  
  frame.clear();		// win or lose, we start anew

  frame.index =eve_num;
  bins_per_frame = hu_decon->GetNbinsY();

  int nbin = bins_per_frame;

  double fudge_factor = 1.1;//
  double fudge_factor_ind = 2.1;
  // U plane
  for (size_t ind=0; ind < hu_decon->GetNbinsX(); ++ind) {

    // do ROI ...
    std::vector<bool> signalsBool;
    signalsBool.resize(nbin, false);
    for (int i = 0; i < nbin; ++ i){
      if (hu_decon->GetBinContent(ind+1,i+1)!=0)
	signalsBool.at(i) = true;
    }
    std::vector<std::vector<int>> rois;
    bool inside = false;
    for (int i=0; i<nbin; ++i) {
      if (inside) {
	if (signalsBool[i]) { // still inside
	  rois.back().push_back(i);
	}else{
	  inside = false;
	}
      }
      else {                  // outside the Rio
	if (signalsBool[i]) { // just entered ROI
	  std::vector<int> roi;
	  roi.push_back(i);
	  rois.push_back(roi);
	  inside = true;
	}
      }
    }
    //    std::cout << rois.size() << std::endl;


    
    
    WireCell::Trace trace;
    trace.chid = ind;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    for (int i=0;i!=rois.size();i++){
      int time = rois.at(i).size() * nrebin;

      if (time <0) time = 12;
      if (time > 800) {
	std::cout << "too long time U " << ind << " " << time << std::endl;
	time = 800;
      }
      for (int j=0; j!=rois.at(i).size();j++){
       	trace.charge.at(rois.at(i).at(j)) = gu->Eval(time) * fudge_factor *fudge_factor_ind * nrebin / 4.;
	// if (trace.charge.at(rois.at(i).at(j))<0)
	//   std::cout << ind << " " << time << " " << gu->Eval(time) << std::endl; 
      }
      //      std::cout << ind << " " << time << std::endl;
    }

    // for (int ibin=0; ibin != bins_per_frame; ibin++) {
    //   if (hu_decon->GetBinContent(ind+1,ibin+1)!=0){
    // 	std::cout << ind << " " << trace.charge.at(ibin) << " " << hu_decon->GetBinContent(ind+1,ibin+1) << std::endl;
    //   }else{
    // 	if (trace.charge.at(ibin) !=0)
    // 	  std::cout << "wrong! " << std::endl;
    //   }
    //   //trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    // }

    
    frame.traces.push_back(trace);
    
  }
  
  // V plane
  for (size_t ind=0; ind < hv_decon->GetNbinsX(); ++ind) {
    // do ROI ...
    std::vector<bool> signalsBool;
    signalsBool.resize(nbin, false);
    for (int i = 0; i < nbin; ++ i){
      if (hv_decon->GetBinContent(ind+1,i+1)!=0)
	signalsBool.at(i) = true;
    }
    std::vector<std::vector<int>> rois;
    bool inside = false;
    for (int i=0; i<nbin; ++i) {
      if (inside) {
	if (signalsBool[i]) { // still inside
	  rois.back().push_back(i);
	}else{
	  inside = false;
	}
      }
      else {                  // outside the Rio
	if (signalsBool[i]) { // just entered ROI
	  std::vector<int> roi;
	  roi.push_back(i);
	  rois.push_back(roi);
	  inside = true;
	}
      }
    }



    WireCell::Trace trace;
    trace.chid = ind + nwire_u;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);

    for (int i=0;i!=rois.size();i++){
      int time = rois.at(i).size() * nrebin;

      if (time <0) time = 12;
      if (time > 800) {
	std::cout << "too long time V " << ind << " " << time << std::endl;
	time = 800;
      }
      
      for (int j=0; j!=rois.at(i).size();j++){
	trace.charge.at(rois.at(i).at(j)) = gv->Eval(time) * fudge_factor * fudge_factor_ind * nrebin / 4.; 
      }
      //      std::cout << ind << " " << time << std::endl;
    }
    
    // for (int ibin=0; ibin != bins_per_frame; ibin++) {
    //   if (hv_decon->GetBinContent(ind+1,ibin+1)!=0){
    // 	std::cout << ind << " " << trace.charge.at(ibin) << " " << hv_decon->GetBinContent(ind+1,ibin+1) << std::endl;
    //   }else{
    // 	if (trace.charge.at(ibin) !=0)
    // 	  std::cout << "wrong! " << std::endl;
    //   }
    //   //trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    // }
    frame.traces.push_back(trace);
  }

  // W plane
  for (size_t ind=0; ind < hw_decon->GetNbinsX(); ++ind) {
    // do ROI ...
    std::vector<bool> signalsBool;
    signalsBool.resize(nbin, false);
    for (int i = 0; i < nbin; ++ i){
      if (hw_decon->GetBinContent(ind+1,i+1)!=0)
	signalsBool.at(i) = true;
    }
    std::vector<std::vector<int>> rois;
    bool inside = false;
    for (int i=0; i<nbin; ++i) {
      if (inside) {
	if (signalsBool[i]) { // still inside
	  rois.back().push_back(i);
	}else{
	  inside = false;
	}
      }
      else {                  // outside the Rio
	if (signalsBool[i]) { // just entered ROI
	  std::vector<int> roi;
	  roi.push_back(i);
	  rois.push_back(roi);
	  inside = true;
	}
      }
    }

    WireCell::Trace trace;
    trace.chid = ind + nwire_u + nwire_v;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);

    for (int i=0;i!=rois.size();i++){
      int time = rois.at(i).size() * nrebin;

      if (time <0) time = 12;
      if (time > 800) {
	//std::cout << "too long time W " << ind << " " << time << std::endl;
	time = 800;
      }
      
      for (int j=0; j!=rois.at(i).size();j++){
	trace.charge.at(rois.at(i).at(j)) = gw->Eval(time) * fudge_factor * nrebin / 4.; 
      }
      //      std::cout << ind << " " << time << std::endl;
    }

    // for (int ibin=0; ibin != bins_per_frame; ibin++) {
    //   if (hw_decon->GetBinContent(ind+1,ibin+1)!=0){
    // 	std::cout << ind << " " << trace.charge.at(ibin) << " " << hw_decon->GetBinContent(ind+1,ibin+1) << std::endl;
    //   }else{
    // 	if (trace.charge.at(ibin) !=0)
    // 	  std::cout << "wrong! " << std::endl;
    //   }
    //   //trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    // }

    
    // for (int ibin=0; ibin != bins_per_frame; ibin++) {
    //   trace.charge.at(ibin) = hw_decon->GetBinContent(ind+1,ibin+1);
    // }
    frame.traces.push_back(trace);
  }

}

WireCell2dToy::uBooNEDataError::uBooNEDataError(const WireCell::GeomDataSource& gds, TH2F *hu_decon, TH2F *hv_decon, TH2F *hw_decon, int eve_num, int nrebin)
  : gds(gds)
  , nrebin(nrebin)
{
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();


#include "charge_err.txt"
  TGraph *gu = new TGraph(201,xu,yu);
  TGraph *gv = new TGraph(201,xv,yv);
  TGraph *gw = new TGraph(201,xw,yw);

  
  frame.clear();		// win or lose, we start anew

  frame.index =eve_num;
  bins_per_frame = hu_decon->GetNbinsY();

  int nbin = bins_per_frame;

  double fudge_factor = 1.1;//
  double fudge_factor_ind = 2.1;
  // U plane
  for (size_t ind=0; ind < hu_decon->GetNbinsX(); ++ind) {

    // do ROI ...
    std::vector<bool> signalsBool;
    signalsBool.resize(nbin, false);
    for (int i = 0; i < nbin; ++ i){
      if (hu_decon->GetBinContent(ind+1,i+1)!=0)
	signalsBool.at(i) = true;
    }
    std::vector<std::vector<int>> rois;
    bool inside = false;
    for (int i=0; i<nbin; ++i) {
      if (inside) {
	if (signalsBool[i]) { // still inside
	  rois.back().push_back(i);
	}else{
	  inside = false;
	}
      }
      else {                  // outside the Rio
	if (signalsBool[i]) { // just entered ROI
	  std::vector<int> roi;
	  roi.push_back(i);
	  rois.push_back(roi);
	  inside = true;
	}
      }
    }
    //    std::cout << rois.size() << std::endl;


    
    
    WireCell::Trace trace;
    trace.chid = ind;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    for (int i=0;i!=rois.size();i++){
      int time = rois.at(i).size() * nrebin;

      if (time <0) time = 12;
      if (time > 800) {
	std::cout << "too long time U " << ind << " " << time << std::endl;
	time = 800;
      }
      for (int j=0; j!=rois.at(i).size();j++){
       	trace.charge.at(rois.at(i).at(j)) = gu->Eval(time) * fudge_factor *fudge_factor_ind * nrebin / 4.;
	// if (trace.charge.at(rois.at(i).at(j))<0)
	//   std::cout << ind << " " << time << " " << gu->Eval(time) << std::endl; 
      }
      //      std::cout << ind << " " << time << std::endl;
    }

    // for (int ibin=0; ibin != bins_per_frame; ibin++) {
    //   if (hu_decon->GetBinContent(ind+1,ibin+1)!=0){
    // 	std::cout << ind << " " << trace.charge.at(ibin) << " " << hu_decon->GetBinContent(ind+1,ibin+1) << std::endl;
    //   }else{
    // 	if (trace.charge.at(ibin) !=0)
    // 	  std::cout << "wrong! " << std::endl;
    //   }
    //   //trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    // }

    
    frame.traces.push_back(trace);
    
  }
  
  // V plane
  for (size_t ind=0; ind < hv_decon->GetNbinsX(); ++ind) {
    // do ROI ...
    std::vector<bool> signalsBool;
    signalsBool.resize(nbin, false);
    for (int i = 0; i < nbin; ++ i){
      if (hv_decon->GetBinContent(ind+1,i+1)!=0)
	signalsBool.at(i) = true;
    }
    std::vector<std::vector<int>> rois;
    bool inside = false;
    for (int i=0; i<nbin; ++i) {
      if (inside) {
	if (signalsBool[i]) { // still inside
	  rois.back().push_back(i);
	}else{
	  inside = false;
	}
      }
      else {                  // outside the Rio
	if (signalsBool[i]) { // just entered ROI
	  std::vector<int> roi;
	  roi.push_back(i);
	  rois.push_back(roi);
	  inside = true;
	}
      }
    }



    WireCell::Trace trace;
    trace.chid = ind + nwire_u;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);

    for (int i=0;i!=rois.size();i++){
      int time = rois.at(i).size() * nrebin;

      if (time <0) time = 12;
      if (time > 800) {
	std::cout << "too long time V " << ind << " " << time << std::endl;
	time = 800;
      }
      
      for (int j=0; j!=rois.at(i).size();j++){
	trace.charge.at(rois.at(i).at(j)) = gv->Eval(time) * fudge_factor * fudge_factor_ind * nrebin / 4.; 
      }
      //      std::cout << ind << " " << time << std::endl;
    }
    
    // for (int ibin=0; ibin != bins_per_frame; ibin++) {
    //   if (hv_decon->GetBinContent(ind+1,ibin+1)!=0){
    // 	std::cout << ind << " " << trace.charge.at(ibin) << " " << hv_decon->GetBinContent(ind+1,ibin+1) << std::endl;
    //   }else{
    // 	if (trace.charge.at(ibin) !=0)
    // 	  std::cout << "wrong! " << std::endl;
    //   }
    //   //trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    // }
    frame.traces.push_back(trace);
  }

  // W plane
  for (size_t ind=0; ind < hw_decon->GetNbinsX(); ++ind) {
    // do ROI ...
    std::vector<bool> signalsBool;
    signalsBool.resize(nbin, false);
    for (int i = 0; i < nbin; ++ i){
      if (hw_decon->GetBinContent(ind+1,i+1)!=0)
	signalsBool.at(i) = true;
    }
    std::vector<std::vector<int>> rois;
    bool inside = false;
    for (int i=0; i<nbin; ++i) {
      if (inside) {
	if (signalsBool[i]) { // still inside
	  rois.back().push_back(i);
	}else{
	  inside = false;
	}
      }
      else {                  // outside the Rio
	if (signalsBool[i]) { // just entered ROI
	  std::vector<int> roi;
	  roi.push_back(i);
	  rois.push_back(roi);
	  inside = true;
	}
      }
    }

    WireCell::Trace trace;
    trace.chid = ind + nwire_u + nwire_v;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);

    for (int i=0;i!=rois.size();i++){
      int time = rois.at(i).size() * nrebin;

      if (time <0) time = 12;
      if (time > 800) {
	//std::cout << "too long time W " << ind << " " << time << std::endl;
	time = 800;
      }
      
      for (int j=0; j!=rois.at(i).size();j++){
	trace.charge.at(rois.at(i).at(j)) = gw->Eval(time) * fudge_factor * nrebin / 4.; 
      }
      //      std::cout << ind << " " << time << std::endl;
    }

    // for (int ibin=0; ibin != bins_per_frame; ibin++) {
    //   if (hw_decon->GetBinContent(ind+1,ibin+1)!=0){
    // 	std::cout << ind << " " << trace.charge.at(ibin) << " " << hw_decon->GetBinContent(ind+1,ibin+1) << std::endl;
    //   }else{
    // 	if (trace.charge.at(ibin) !=0)
    // 	  std::cout << "wrong! " << std::endl;
    //   }
    //   //trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    // }

    
    // for (int ibin=0; ibin != bins_per_frame; ibin++) {
    //   trace.charge.at(ibin) = hw_decon->GetBinContent(ind+1,ibin+1);
    // }
    frame.traces.push_back(trace);
  }

  //std::cout << frame.traces.size() << " " << bins_per_frame << std::endl;
}

WireCell2dToy::uBooNEDataError::~uBooNEDataError(){
}

int WireCell2dToy::uBooNEDataError::jump(int frame_number){
  return frame.index;
}

int WireCell2dToy::uBooNEDataError::size() const{
  return 1;
}

