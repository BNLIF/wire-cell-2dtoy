#include "WireCellSst/GeomDataSource.h"
#include "WireCellSst/DatauBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/BadTiling.h"

#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCell2dToy/SimpleBlobToyTiling.h"

#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/ToyMatrixExclusive.h"
#include "WireCell2dToy/ToyMatrixKalman.h"
#include "WireCell2dToy/ToyMatrixIterate.h"
#include "WireCell2dToy/ToyMatrixIterate_SingleWire.h"
#include "WireCell2dToy/ToyMatrixIterate_Only.h"


#include "WireCell2dToy/ToyMatrixMarkov.h"
#include "WireCell2dToy/ToyMetric.h"
#include "WireCell2dToy/BlobMetric.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"

#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"

#include "WireCellData/GeomCluster.h"
//#include "WireCellNav/SliceDataSource.h"


#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/SimDataSource.h"
#include "WireCellNav/SliceDataSource.h"
#include "WireCellSst/Util.h"
#include "WireCellData/SimTruth.h"
#include "WireCell2dToy/ToyDepositor.h"
#include "WireCellNav/GenerativeFDS.h"
#include "WireCell2dToy/ToySignalSimu.h"
#include "WireCell2dToy/ToySignalSimuTrue.h"
#include "WireCell2dToy/DataSignalGaus.h"
#include "WireCell2dToy/DataSignalWien.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include <iostream>

using namespace WireCell;
using namespace std;




int main(int argc, char* argv[])
{
  if (argc < 4) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num -t[0,1] -s[0,1,2]" << endl;
    return 1;
  }

  int two_plane = 0;
  int save_file = 0;
  int nt_off1 = 0;
  int nt_off2 = 0;
  int deghost = 1;
  for(Int_t i = 1; i != argc; i++){
     switch(argv[i][1]){
     case 't':
       two_plane = atoi(&argv[i][2]); 
       break;
     case 's':
       save_file = atoi(&argv[i][2]); 
       break;
     case 'd':
       deghost = atoi(&argv[i][2]); 
       break;
     case 'a':
       nt_off1 = atoi(&argv[i][2]);
       break;
     case 'b':
       nt_off2 = atoi(&argv[i][2]);
       break;
     }
  }
  
  if (two_plane)
    cout << "Enable Two Plane Reconstruction " << endl; 
  
  WireCellSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;

  cout << "Pitch: " << gds.pitch(WirePlaneType_t(0)) 
       << " " << gds.pitch(WirePlaneType_t(1)) 
       << " " << gds.pitch(WirePlaneType_t(2))
       << endl;
  cout << "Angle: " << gds.angle(WirePlaneType_t(0)) 
       << " " << gds.angle(WirePlaneType_t(1)) 
       << " " << gds.angle(WirePlaneType_t(2))
       << endl;
  
  
  
  
  // WireCell::FrameDataSource* fds = 0;
  // fds = WireCellSst::make_fds(root_file);
  // if (!fds) {
  //   cerr << "ERROR: failed to get FDS from " << root_file << endl;
  //   return 1;
  // }
  
  //float unit_dis = 1.01483;  // 58KV @ 226.5 V/cm
  float unit_dis = 1.14753;  // 70 KV @ 226.5 V/cm

  //**** time offset for 58kV ****// 
  //float toffset_1=-(1.834-1.647) -0.1;//+ 0.1 + (nt_off1 * 0.1 - 0.5 );  // time offset between u/v 
  //float toffset_2=-(1.834+1.555-1.539-1.647) +0.1;//+ 0.3 + (nt_off2*0.2 - 1); // time offset between u/w
  //float toffset_3=-0.5; //overall time shift

  //final offset after time scan (70kV)
  float toffset_1=-(1.834-1.647) -0.1 - 0.5;
  float toffset_2=-(1.834+1.555-1.539-1.647) +0.1 - 0.5;
  float toffset_3=0.0;
  
  int save_image_outline_flag = 0; // prescale flag 
  

  int total_time_bin=9592;
  int recon_threshold = 2000;
  int frame_length = 3200;
  int max_events = 100;
  int eve_num  = atoi(argv[3]);
  int nrebin = 4;

  TPCParams& mp = Singleton<TPCParams>::Instance();
  
  double pitch_u = gds.pitch(WirePlaneType_t(0));
  double pitch_v = gds.pitch(WirePlaneType_t(1));
  double pitch_w = gds.pitch(WirePlaneType_t(2));
  double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;

  mp.set_pitch_u(pitch_u);
  mp.set_pitch_v(pitch_v);
  mp.set_pitch_w(pitch_w);
  mp.set_ts_width(time_slice_width);
  
  std::cout << "Singleton: " << mp.get_pitch_u() << " " << mp.get_pitch_v() << " " << mp.get_pitch_w() << " " << mp.get_ts_width() << std::endl;
  


  float threshold_u = 5.87819e+02 * 4.0;
  float threshold_v = 8.36644e+02 * 4.0;
  float threshold_w = 5.67974e+02 * 4.0;

  float threshold_ug = 755.96;
  float threshold_vg = 822.81;
  float threshold_wg = 510.84;
  
  int time_offset = -52.;
  


  const char* root_file = argv[2];
 
  
  int run_no, subrun_no, event_no;
  // sst->SetBranchAddress("eventNo",&event_no);
  // sst->SetBranchAddress("runNo",&run_no);
  // sst->SetBranchAddress("subRunNo",&subrun_no);
  // sst->GetEntry(eve_num);
  
  
  WireCellSst::DatauBooNEFrameDataSource data_fds(root_file,gds,total_time_bin);
  if (save_file != 2){
    data_fds.jump(eve_num);
    if (save_file == 1)
      data_fds.Save();
  }
  
  run_no = data_fds.get_run_no();
  subrun_no = data_fds.get_subrun_no();
  event_no = data_fds.get_event_no();
  
  cout << "Run No: " << run_no << " " << subrun_no << " " << eve_num << endl;


  // WireMap& uplane_map = data_fds.get_u_map();
  // WireMap& vplane_map = data_fds.get_v_map();
  // WireMap& wplane_map = data_fds.get_w_map();

  ChirpMap& uplane_map = data_fds.get_u_cmap();
  ChirpMap& vplane_map = data_fds.get_v_cmap();
  ChirpMap& wplane_map = data_fds.get_w_cmap();

  // std::cout << uplane_map.size() << " " << vplane_map.size() << " " << wplane_map.size() << std::endl;
    
  cout << "Deconvolution with Gaussian filter" << endl;
  WireCell2dToy::DataSignalGausFDS gaus_fds(data_fds,gds,total_time_bin/nrebin,max_events,toffset_1,toffset_2,toffset_3); // gaussian smearing for charge estimation
  if (save_file != 2){
    gaus_fds.jump(eve_num);
    if (save_file == 1)
      gaus_fds.Save();
  }else{
  
  }

  cout << "Deconvolution with Wiener filter" << endl; 
  WireCell2dToy::DataSignalWienFDS wien_fds(data_fds,gds,uplane_map, vplane_map, wplane_map, total_time_bin/nrebin,max_events,toffset_1,toffset_2,toffset_3); // weiner smearing for hit identification
  if (save_file !=2 ){
    wien_fds.jump(eve_num);
    if (save_file == 1)
      wien_fds.Save();
  }else{
    
  }


  data_fds.Clear();

  std::vector<float>& uplane_rms = wien_fds.get_uplane_rms();
  std::vector<float>& vplane_rms = wien_fds.get_vplane_rms();
  std::vector<float>& wplane_rms = wien_fds.get_wplane_rms();

  // // hack for now ...  remove the very busy wires ... 
  // for (int i=0;i!=uplane_rms.size();i++){
  //   //cout << "U " << i << " " << uplane_rms.at(i) << endl;
  //   if (uplane_rms.at(i) > 1500) {
  //     uplane_rms.at(i) *=2;
  //     uplane_map.erase(i);
  //   }
  // }
  // for (int i=0;i!=vplane_rms.size();i++){
  //   //cout << "V " << i << " " << vplane_rms.at(i) << endl;
  //   if (vplane_rms.at(i) > 2000 && vplane_rms.at(i)<3000){
  //     vplane_rms.at(i) *=2;
  //     vplane_map.erase(i);
  //   }else if (vplane_rms.at(i)>=3000){
  //     vplane_rms.at(i) *=10;
  //     vplane_map.erase(i);
  //   }
  // }
  // for (int i=0;i!=wplane_rms.size();i++){
  //   //cout << "W " << i << " " << wplane_rms.at(i) << endl;
  //   if (wplane_rms.at(i) > 1000) {
  //     wplane_rms.at(i) *=2;
  //     wplane_map.erase(i);
  //   }
  // }
  

  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  int nwire_u = wires_u.size();
  int nwire_v = wires_v.size();
  int nwire_w = wires_w.size();
  
  
 
  

 

  // WireCellSst::ToyuBooNESliceDataSource sds(wien_fds,gaus_fds,threshold_u, 
  // 					    threshold_v, threshold_w, 
  // 					    threshold_ug, 
  // 					    threshold_vg, threshold_wg, 
  // 					    nwire_u, 
  // 					    nwire_v, nwire_w); 

  WireCellSst::ToyuBooNESliceDataSource sds(wien_fds,gaus_fds,threshold_u, 
  					    threshold_v, threshold_w, 
  					    threshold_ug, 
  					    threshold_vg, threshold_wg, 
  					    nwire_u, 
  					    nwire_v, nwire_w,
  					    &uplane_rms, &vplane_rms, &wplane_rms); 
    
  
  
  int ncount = 0;
  int ncount1 = 0;
  int ncount2 = 0;

  int ncount_t = 0;
  

  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
  WireCell2dToy::BadTiling **badtiling = new WireCell2dToy::BadTiling*[2400];
  WireCell2dToy::MergeToyTiling **mergetiling = new WireCell2dToy::MergeToyTiling*[2400];
  WireCell2dToy::ToyMatrix **toymatrix = new WireCell2dToy::ToyMatrix*[2400];
    
  //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  
  int ncount_mcell = 0;

  //  delete fds;

  int start_num = 0 ;
  int end_num = sds.size()-1;

  // int start_num = 61 ;
  // int end_num = 61;

  for (int i=start_num;i!=end_num+1;i++){
 
    sds.jump(i);
    WireCell::Slice slice = sds.get();

    //cout << i << " " << slice.group().size() << std::endl;
    toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds,0.15,0.2,0.1,threshold_ug,threshold_vg, threshold_wg, &uplane_rms, &vplane_rms, &wplane_rms);

    if (two_plane)
      toytiling[i]->twoplane_tiling(i,nrebin,gds,uplane_rms,vplane_rms,wplane_rms, uplane_map, vplane_map, wplane_map);


    GeomCellSelection allcell = toytiling[i]->get_allcell();
    GeomWireSelection allwire = toytiling[i]->get_allwire();

    cout << i << " " << allcell.size() << " " << allwire.size() << endl;

    mergetiling[i] = new WireCell2dToy::MergeToyTiling(*toytiling[i],i,3);
    
   
    
    if (i==0){
      badtiling[i] = new WireCell2dToy::BadTiling(i,nrebin,uplane_map,vplane_map,wplane_map,gds,1);
    }

    //badtiling[i] = new WireCell2dToy::BadTiling(i,nrebin,uplane_map,vplane_map,wplane_map,gds);

    // // if (toymatrix[i]->Get_Solve_Flag()!=0){
    // //   toymatrix[i]->Update_pred();
    // //   toymatrix[i]->Print();
    // // }

    // //draw ... 
    // TApplication theApp("theApp",&argc,argv);
    // theApp.SetReturnFromRun(true);
    
    // TCanvas c1("ToyMC","ToyMC",800,600);
    // c1.Draw();
    
    // WireCell2dToy::ToyEventDisplay display(c1, gds);
    // display.charge_min = 0;
    // display.charge_max = 5e4;


    // gStyle->SetOptStat(0);
    
    // const Int_t NRGBs = 5;
    // const Int_t NCont = 255;
    // Int_t MyPalette[NCont];
    // Double_t stops[NRGBs] = {0.0, 0.34, 0.61, 0.84, 1.0};
    // Double_t red[NRGBs] = {0.0, 0.0, 0.87 ,1.0, 0.51};
    // Double_t green[NRGBs] = {0.0, 0.81, 1.0, 0.2 ,0.0};
    // Double_t blue[NRGBs] = {0.51, 1.0, 0.12, 0.0, 0.0};
    // Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    // gStyle->SetNumberContours(NCont);
    // for (int kk=0;kk!=NCont;kk++) MyPalette[kk] = FI+kk;
    // gStyle->SetPalette(NCont,MyPalette);

    

    // display.init(0,10.3698,-2.33/2.,2.33/2.);
    // display.draw_mc(1,WireCell::PointValueVector(),"colz");
    
    // // display.draw_bad_region(uplane_map,i,nrebin,0,"same");
    // // display.draw_bad_region(vplane_map,i,nrebin,1,"same");
    // // display.draw_bad_region(wplane_map,i,nrebin,2,"same");
    // display.draw_bad_cell(badtiling[i]->get_cell_all());
  
    // display.draw_cells(toytiling[i]->get_allcell(),"*same");
    // display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell
    
    // display.draw_wires_charge(toytiling[i]->wcmap(),"Fsame",FI);
    // display.draw_cells_charge(toytiling[i]->get_allcell(),"Fsame");
  
    
    
    // theApp.Run();
  }
  


  //initial clustering ... 
  // form a map to illustrate connectivities 
  GeomCellCellsMap cell_prev_map;
  GeomCellCellsMap cell_next_map;
  
  for (int i=start_num;i!=end_num;i++){
    GeomCellSelection curr_mcell = mergetiling[i]->get_allcell();
    GeomCellSelection next_mcell = mergetiling[i+1]->get_allcell();
    for (int j=0;j!=curr_mcell.size();j++){
      const MergeGeomCell *curr_cell = (MergeGeomCell*)curr_mcell.at(j);
      for (int k=0;k!=next_mcell.size();k++){
	const MergeGeomCell *next_cell = (MergeGeomCell*)next_mcell.at(k);
	if (curr_cell->Overlap(*next_cell)){
	  if (cell_next_map.find(curr_cell)==cell_next_map.end()){
	    GeomCellSelection cells;
	    cell_next_map[curr_cell] = cells;
	  }
	  cell_next_map[curr_cell].push_back(next_cell);
	  
	  if (cell_prev_map.find(next_cell)==cell_prev_map.end()){
	    GeomCellSelection cells;
	    cell_prev_map[next_cell].push_back(curr_cell);
	    cell_prev_map[next_cell] = cells;
	  }
	  cell_prev_map[next_cell].push_back(curr_cell);
	}
      }
    }
  }
    
  //do clustering ... 
  for (int i=start_num;i!=end_num+1;i++){
    GeomCellSelection pallmcell = mergetiling[i]->get_allcell();
    GeomCellSelection allmcell;
    for (int j=0;j!=pallmcell.size();j++){
      const GeomCell* mcell = pallmcell[j];
      int flag_save_cell = 1;
            
      if (flag_save_cell == 1)
      	allmcell.push_back(mcell);
    }
    
    if (cluster_set.empty()){
      // if cluster is empty, just insert all the mcell, each as a cluster
      for (int j=0;j!=allmcell.size();j++){
	GeomCluster *cluster = new GeomCluster(*((MergeGeomCell*)allmcell[j]));
	cluster_set.insert(cluster);
      }
    }else{
      for (int j=0;j!=allmcell.size();j++){
	int flag = 0;
	int flag_save = 0;
	GeomCluster *cluster_save = 0;
	cluster_delset.clear();
	// loop through merged cell
	for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
	  //loop through clusters
	  
	  flag += (*it)->AddCell(*((MergeGeomCell*)allmcell[j]));
	  if (flag==1 && flag != flag_save){
	    cluster_save = *it;
	  }else if (flag>1 && flag != flag_save){
	    cluster_save->MergeCluster(*(*it));
	    cluster_delset.insert(*it);
	  }
	  flag_save = flag;
  	  
	}
	
	for (auto it = cluster_delset.begin();it!=cluster_delset.end();it++){
	  cluster_set.erase(*it);
	  delete (*it);
	}
	
	if (flag==0){
	  GeomCluster *cluster = new GeomCluster(*((MergeGeomCell*)allmcell[j]));
	  cluster_set.insert(cluster);
	}
	
      }
    }
  }
  
  

  // create a vector of array of merged cells
  // identify the good cells, # of cells in cluster > ???
  // also need to have previous and next merged cells
  std::vector<GeomCellSelection> good_mcells;
  for (int i=start_num;i!=end_num+1;i++){
    GeomCellSelection mcells;
    good_mcells.push_back(mcells);
  }
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    int number_mcells = (*it)->get_allcell().size();
    int number_time = (*it)->get_ordercell().size();

    
    // if very big cluster, then no need to be both side ...

    if (number_time >=5 && number_mcells >=6){
      int flag = 0;
      // if cluster contains a three-wire cell?
      for (int i=0;i!=(*it)->get_allcell().size();i++){
	const MergeGeomCell* mcell = (MergeGeomCell*)((*it)->get_allcell().at(i));
	int time_slice = mcell->GetTimeSlice();
	GeomCellSelection& three_wires_cells = mergetiling[time_slice]->get_three_wires_cells();
	auto it = find(three_wires_cells.begin(),three_wires_cells.end(),mcell);
	if (it != three_wires_cells.end()){
	  flag = 1;
	  break;
	}
      }
      
      //two wire cluster, need longer ... 
      if (flag == 0 && (number_time < 10 && number_mcells < 13)) continue;
      
      for (int i=0;i!=(*it)->get_allcell().size();i++){
	const MergeGeomCell* mcell = (MergeGeomCell*)((*it)->get_allcell().at(i));
	// see front or back
	int time_slice = mcell->GetTimeSlice();
	if (cell_prev_map[mcell].size()>0 && cell_next_map.size()>0){
	  //std::cout << "Xin1: " << time_slice << std::endl;
	  good_mcells.at(time_slice).push_back(mcell);
	}
      }
    }
  }



  //separate the deghosting ... 
  for (int i=start_num;i!=end_num+1;i++){
     if ( deghost) 
       mergetiling[i]->deghost(good_mcells.at(i));
   
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    GeomWireSelection allmwire = mergetiling[i]->get_allwire();
    
    cout << i << " " << allmcell.size() << " " << allmwire.size() << endl;
    

    toymatrix[i] = new WireCell2dToy::ToyMatrix(*toytiling[i],*mergetiling[i]);
    
    GeomCellSelection& two_wires_cells = mergetiling[i]->get_two_wires_cells();
    GeomCellSelection& three_wires_cells = mergetiling[i]->get_three_wires_cells();

    // if (two_plane){
    //   if (toymatrix[i]->Get_Solve_Flag()==0)
    //    	WireCell2dToy::ToyMatrixIterate_SingleWire toymatrix_single_wire_it(*toymatrix[i],mergetiling[i]);
      
    //   // double chi2_3p = 0;
    //   // double ndf_3p = 0;
    //   // std::vector<double> Cxt, dCxt;
    //   // Cxt.resize(toymatrix[i]->Get_mcindex());
    //   // dCxt.resize(toymatrix[i]->Get_mcindex());
      
    //   // //clean up the results first ... 
    //   // toymatrix[i]->Set_Solve_Flag(0);
    //   // toymatrix[i]->Set_chi2(-1);
      
    //   // // deal with three planes
    //   // std::vector<int> already_removed;
    //   // for (int j=0;j!=two_wires_cells.size();j++){
    //   // 	int index = toymatrix[i]->Get_mcindex(two_wires_cells.at(j));
    //   // 	already_removed.push_back(index);
    //   // }
    //   // if (three_wires_cells.size() > 0){
    //   // 	WireCell2dToy::ToyMatrixIterate_Only toymatrix_it_only(*toymatrix[i],mergetiling[i]);
    //   // 	// WireCell2dToy::ToyMatrixIterate toymatrix_it(*toymatrix[i],already_removed, 2000 , 1e5);
    //   // 	// if (toymatrix[i]->Get_Solve_Flag()==0 ){
    //   // 	//   // if not solved
    //   // 	//   // deal everything together ... 
    //   // 	//   WireCell2dToy::ToyMatrixIterate_Only toymatrix_it_only(*toymatrix[i],mergetiling[i]);
    //   // 	// }else{ 
    //   // 	//   // if solved?
    //   // 	//   chi2_3p = toymatrix[i]->Get_Chi2();
    //   // 	//   ndf_3p = toymatrix[i]->Get_ndf();
    //   // 	//   for (int j=0;j!=toymatrix[i]->Get_mcindex();j++){
    //   // 	//     Cxt.at(j) = toymatrix[i]->Get_value(j);
    //   // 	//     dCxt.at(j) = toymatrix[i]->Get_error(j);
    //   // 	//   }
    //   // 	//   // deal with two planes
	
    //   // 	//   //combine them together ... 
    //   // 	// }
    //   // }else{
    //   //   // deal with two planes without any constraints, like to deal everything together
    //   // 	WireCell2dToy::ToyMatrixIterate_Only toymatrix_it_only(*toymatrix[i],mergetiling[i]);
    //   // }
    // }else{
    //   if (toymatrix[i]->Get_Solve_Flag()==0)
    // 	WireCell2dToy::ToyMatrixIterate toymatrix_it(*toymatrix[i],2000,1e5);
    // }
    // 
    //   // deal with three-planes cells, 
    //   // try iterative method first
    //   
    // }
    // // deal with two-planes cells
    // if (two_wires_cells.size() > 0){
    // }
    //comebine them 
    
    cout << "chi2: " << toymatrix[i]->Get_Chi2() << endl;
    cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
  }


  
  // double penalty = 6;
  // std::cout << "Starting to use connectivitiy" << std::endl;
  //  std::list<int> solve_list;
   
  //  if (start_num != end_num){
  //    int first_solve=-1;
  //    for (int i=start_num; i!=end_num+1;i++){
  //      if (toymatrix[i]->Get_Solve_Flag()!=0){
  //  	 first_solve = i;
  //  	 break;
  //      }
  //    }
  //    if (first_solve == -1) first_solve = start_num;

    
  //    for (int i=first_solve+1;i<=end_num-1;i++){
  //      if (toymatrix[i]->Get_Solve_Flag()==0 && toymatrix[i]->Get_mcindex()>0){ 
  // 	 if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
  // 	   if (two_plane){
  // 	     WireCell2dToy::ToyMatrixIterate_SingleWire toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty);
  // 	   }else{
  // 	     WireCell2dToy::ToyMatrixIterate toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty);
  // 	   }
  // 	   GeomCellSelection allmcell = mergetiling[i]->get_allcell();
	 	   
  // 	   cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
  // 	   cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
  // 	 }else{
  // 	   solve_list.push_back(i); 
  // 	 }
  //      }
  //    }
     
  //    for (int i=first_solve-1;i>=start_num+1;i--){
  //      if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
  // 	 if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
  // 	   if (two_plane){
  // 	     WireCell2dToy::ToyMatrixIterate_SingleWire toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty);
  // 	   }else{
  // 	     WireCell2dToy::ToyMatrixIterate toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty);
  // 	   }
	   
  // 	   GeomCellSelection allmcell = mergetiling[i]->get_allcell();
	   
	 
  // 	   cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
  // 	   cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
  // 	 }else{
  // 	   solve_list.push_back(i);
  // 	 }
  //      }
  //    }
  //  }
   
  //  // start second round ...
  //  // std::cout << "Connectivitiy rest " << solve_list.size() << std::endl;
  //  int prev_count = 0;
  //  while (solve_list.size() >0){
  //    int curr_count = solve_list.size();
  //    std::cout << "Connectivitiy rest " << solve_list.size() << std::endl;

  //    if (curr_count == prev_count){
  //      int i = solve_list.front(); // pick the first element ... 
  //      if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
  // 	 if (two_plane){
  // 	   WireCell2dToy::ToyMatrixIterate_SingleWire toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty);
  // 	 }else{
  // 	   WireCell2dToy::ToyMatrixIterate toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty);
  // 	 }
  // 	 GeomCellSelection allmcell = mergetiling[i]->get_allcell();
	 
  // 	 cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
  // 	 cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
  // 	 solve_list.erase(solve_list.begin());
  //      }
  //    }else{
  //      for (auto it = solve_list.begin(); it!= solve_list.end(); it++){
  // 	 int i = *it;
  // 	 if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
  // 	   if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
  // 	     if (two_plane){
  // 	       WireCell2dToy::ToyMatrixIterate_SingleWire toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty);
  // 	     }else{
  // 	       WireCell2dToy::ToyMatrixIterate toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty);
  // 	     }
  // 	     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
  // 	     cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
  // 	     cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
  // 	     it = solve_list.erase(it);
  // 	   }
  // 	 }
  //      }
  //    }
     
  //    prev_count = curr_count;
  //  }


   
  //  // by the end do the final two
  //  if (toymatrix[end_num]->Get_Solve_Flag()==0&& toymatrix[end_num]->Get_mcindex()>0){
  //    if (two_plane){
  //      WireCell2dToy::ToyMatrixIterate_SingleWire toymatrix_it(toymatrix[end_num-1],toymatrix[end_num],0,mergetiling[end_num-1],mergetiling[end_num],0,2000,1e5,penalty);
  //    }else{
  //      WireCell2dToy::ToyMatrixIterate toymatrix_it(toymatrix[end_num-1],toymatrix[end_num],0,mergetiling[end_num-1],mergetiling[end_num],0,2000,1e5,penalty);
  //    }
     
  //    GeomCellSelection allmcell = mergetiling[end_num]->get_allcell();
          
  //    cout << "chi2: " << end_num << " " << toymatrix[end_num]->Get_Chi2() << endl;
  //    cout << "NDF: " << toymatrix[end_num]->Get_ndf() << endl;
  //  }
   
  //  if (toymatrix[start_num]->Get_Solve_Flag()==0){
  //    if (two_plane){
  //      WireCell2dToy::ToyMatrixIterate_SingleWire toymatrix_it(0,toymatrix[start_num],toymatrix[start_num+1],0,mergetiling[start_num],mergetiling[start_num+1],2000,1e5,penalty);
  //    }else{
  //      WireCell2dToy::ToyMatrixIterate toymatrix_it(0,toymatrix[start_num],toymatrix[start_num+1],0,mergetiling[start_num],mergetiling[start_num+1],2000,1e5,penalty);
  //    }
     
     
  //    GeomCellSelection allmcell = mergetiling[start_num]->get_allcell();
          
  //    cout << "chi2: " << start_num << " " << toymatrix[start_num]->Get_Chi2() << endl;
  //    cout << "NDF: " << toymatrix[start_num]->Get_ndf() << endl;
  //  }

   
   // std::cout << "Starting MCMC" << std::endl;
   // solve_list.clear();
   
   // //without  time information
   // // for (int i=start_num;i!=end_num+1;i++){
   // //   if (toymatrix[i]->Get_Solve_Flag()==0){
   // //     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   // //     WireCell2dToy::ToyMatrixMarkov toymatrix_markov(*toymatrix[i],*mergetiling[i],&allmcell,Good_MCells.at(i-start_num));
   // //     //WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i],&allmcell);
   // //     CellChargeMap ccmap = truthtiling[i]->ccmap();
   // //     if (toymatrix[i]->Get_Solve_Flag()!=0)
   // // 	toymetric.Add(allmcell,*toymatrix[i],ccmap);
   // //     toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
   // //     cout << " chi2: " << i << " " << toymatrix[i]->Get_Chi2() 
   // // 	   << " NDF: " << toymatrix[i]->Get_ndf() << endl;
   // //   }
   // // }
   
   
   // //with time information
   // if (start_num != end_num){
   //   int first_solve=-1;
   //   for (int i=start_num; i!=end_num+1;i++){
   //     if (toymatrix[i]->Get_Solve_Flag()!=0){
   // 	 first_solve = i;
   // 	 break;
   //     }
   //   }
   //   if (first_solve == -1) first_solve = start_num;
     
     
   //   for (int i=first_solve+1;i<=end_num-1;i++){
   //     if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
   // 	 if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
   // 	   GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   // 	   WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell,1500,2000,penalty);
   	  	   
   // 	   cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
   // 	   cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
   // 	 }else{
   // 	   solve_list.push_back(i);
   // 	 }
   // 	 //toymetric.Print();
   //     }
   //   }
     
   //   // go to early ones 
   //   for (int i=first_solve-1;i>=start_num+1;i--){
   //     if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
   // 	 if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
   // 	   GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   // 	   WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell,1500,2000,penalty);
	   

	   
   // 	   cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
   // 	   cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
   // 	 }else{
   // 	   solve_list.push_back(i);
   // 	 }
   //     }
   //   }
      
   //   // do the while ... 
   //   int prev_count = 0;
   // while (solve_list.size() >0){
   //   int curr_count = solve_list.size();
   //   std::cout << "Connectivitiy rest " << solve_list.size() << std::endl;

   //   if (curr_count == prev_count){
   //     int i = solve_list.front(); // pick the first element ... 
   //     if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
   // 	 GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   // 	 WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell,1500,2000,penalty);
	 

   // 	 cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
   // 	 cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
   // 	 solve_list.erase(solve_list.begin());
   //     }
   //   }else{
   //     for (auto it = solve_list.begin(); it!= solve_list.end(); it++){
   // 	 int i = *it;
   // 	 if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
   // 	   if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
   // 	     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   // 	     WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell,1500,2000,penalty);
	     

   // 	     cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
   // 	     cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
   // 	     it = solve_list.erase(it);
   // 	   }
   // 	 }
   //     }
   //   }
     
   //   prev_count = curr_count;
   // }


   //   //deal with the start/end ones ... 
   //   if (toymatrix[end_num]->Get_Solve_Flag()==0&& toymatrix[end_num]->Get_mcindex()>0){
   //     GeomCellSelection allmcell = mergetiling[end_num]->get_allcell();
   //     WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[end_num-1],toymatrix[end_num],0,mergetiling[end_num-1],mergetiling[end_num],0,&allmcell,1500,2000,penalty);
       
       

       
   //     cout << "chi2: " << end_num << " " << toymatrix[end_num]->Get_Chi2() << endl;
   //     cout << "NDF: " << toymatrix[end_num]->Get_ndf() << endl;
   //   }
     
   //   if (toymatrix[start_num]->Get_Solve_Flag()==0&& toymatrix[start_num]->Get_mcindex()>0){
   //     GeomCellSelection allmcell = mergetiling[start_num]->get_allcell();
   //     WireCell2dToy::ToyMatrixMarkov toymatrix_markov(0,toymatrix[start_num],toymatrix[start_num+1],0,mergetiling[start_num],mergetiling[start_num+1],&allmcell,1500,2000,penalty);
       
       
       
   //     cout << "chi2: " << start_num << " " << toymatrix[start_num]->Get_Chi2() << endl;
   //     cout << "NDF: " << toymatrix[start_num]->Get_ndf() << endl;
   //   }
   // }




  //clear the previous cluster ... 
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    delete (*it);
  }
  cluster_set.clear();
  cluster_delset.clear();


  // form a map to illustrate connectivities 
  cell_prev_map.clear();
  cell_next_map.clear();
   
  for (int i=start_num;i!=end_num;i++){
    GeomCellSelection curr_mcell = mergetiling[i]->get_allcell();
    GeomCellSelection next_mcell = mergetiling[i+1]->get_allcell();
    for (int j=0;j!=curr_mcell.size();j++){
      const MergeGeomCell *curr_cell = (MergeGeomCell*)curr_mcell.at(j);
      for (int k=0;k!=next_mcell.size();k++){
	const MergeGeomCell *next_cell = (MergeGeomCell*)next_mcell.at(k);
	if (curr_cell->Overlap(*next_cell)){
	  
	  if (cell_next_map.find(curr_cell)==cell_next_map.end()){
	    GeomCellSelection cells;
	    cell_next_map[curr_cell] = cells;
	  }
	  cell_next_map[curr_cell].push_back(next_cell);
	  
	  if (cell_prev_map.find(next_cell)==cell_prev_map.end()){
	    GeomCellSelection cells;
	    cell_prev_map[next_cell].push_back(curr_cell);
	    cell_prev_map[next_cell] = cells;
	  }
	  cell_prev_map[next_cell].push_back(curr_cell);
	}
      }
    }
  }
  // save good cluster_cells;
  GeomCellSelection good_cluster_cells;
  for (int i=start_num;i!=end_num+1;i++){
    GeomCellSelection curr_mcell = mergetiling[i]->get_allcell();
    for (int j=0;j!=curr_mcell.size();j++){
      const MergeGeomCell *curr_cell = (MergeGeomCell*)curr_mcell.at(j);
      int flag_good_cluster_cell = 0;
      
      if (cell_next_map.find(curr_cell) != cell_next_map.end()){
	for (int k=0;k!=cell_next_map[curr_cell].size();k++){
	  if (toymatrix[i+1]->Get_Cell_Charge(cell_next_map[curr_cell].at(k))> recon_threshold){
	    flag_good_cluster_cell = 1;
	    break;
	  }
	}
      }
      if (flag_good_cluster_cell == 0){
	if (cell_prev_map.find(curr_cell) != cell_prev_map.end()){
	  for (int k=0;k!=cell_prev_map[curr_cell].size();k++){
	    if (toymatrix[i-1]->Get_Cell_Charge(cell_prev_map[curr_cell].at(k))> recon_threshold){
	      flag_good_cluster_cell = 1;
	      break;
	    }
	  }
	}
      }
      
      if (flag_good_cluster_cell == 1)
	good_cluster_cells.push_back(curr_cell);
    }
  }
  
    
  //do clustering ... 
  for (int i=start_num;i!=end_num+1;i++){
    GeomCellSelection pallmcell = mergetiling[i]->get_allcell();
    GeomCellSelection allmcell;
    for (int j=0;j!=pallmcell.size();j++){
      const GeomCell* mcell = pallmcell[j];
      // if (toymatrix[i]->Get_Cell_Charge(mcell)> recon_threshold || toymatrix[i]->Get_Solve_Flag()==0){
      // 	allmcell.push_back(mcell);
      // }

       int flag_save_cell = 0;

      if (toymatrix[i]->Get_Solve_Flag()==0){
      	flag_save_cell = 1;
      }else{
      	if (toymatrix[i]->Get_Cell_Charge(mcell)> recon_threshold){
      	  flag_save_cell = 1;
      	}else{
      	  if (i == start_num || i == end_num + 1) continue;
	  
      	  // there are good cells from the prev and next
      	  flag_save_cell = 0;
      	  //	  std::cout << "Xin: " << cell_next_map[mcell].size() <<  " " << cell_prev_map[mcell].size() << std::endl;

      	  for (int k=0;k!=cell_next_map[mcell].size();k++){
	    if (find(good_cluster_cells.begin(),good_cluster_cells.end(),cell_next_map[mcell].at(k))!=good_cluster_cells.end()){
	      //if (toymatrix[i+1]->Get_Cell_Charge(cell_next_map[mcell].at(k))> recon_threshold){
      	      flag_save_cell = 1;
      	      break;
      	    }
      	  }
      	  if (flag_save_cell==1){
      	    flag_save_cell = 0;
      	    for (int k=0;k!=cell_prev_map[mcell].size();k++){
	      if (find(good_cluster_cells.begin(),good_cluster_cells.end(),cell_prev_map[mcell].at(k))!=good_cluster_cells.end()){
		//if (toymatrix[i-1]->Get_Cell_Charge(cell_prev_map[mcell].at(k))> recon_threshold){
      		flag_save_cell = 1;
      		break;
      	      }
      	    }
      	  }
      	}
      }
      
      
      if (flag_save_cell == 1)
      	allmcell.push_back(mcell);

    }
    
    
    if (cluster_set.empty()){
      // if cluster is empty, just insert all the mcell, each as a cluster
      for (int j=0;j!=allmcell.size();j++){
	GeomCluster *cluster = new GeomCluster(*((MergeGeomCell*)allmcell[j]));
	cluster_set.insert(cluster);
      }
    }else{
      for (int j=0;j!=allmcell.size();j++){
	int flag = 0;
	int flag_save = 0;
	GeomCluster *cluster_save = 0;
	
	cluster_delset.clear();
	
	// loop through merged cell
	for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
	  //loop through clusters
	  
	  flag += (*it)->AddCell(*((MergeGeomCell*)allmcell[j]));
	  if (flag==1 && flag != flag_save){
	    cluster_save = *it;
	  }else if (flag>1 && flag != flag_save){
	    cluster_save->MergeCluster(*(*it));
	    cluster_delset.insert(*it);
	  }
	  flag_save = flag;
  	  
	}
	
	for (auto it = cluster_delset.begin();it!=cluster_delset.end();it++){
	  cluster_set.erase(*it);
	  delete (*it);
	}
	
	if (flag==0){
	  GeomCluster *cluster = new GeomCluster(*((MergeGeomCell*)allmcell[j]));
	  cluster_set.insert(cluster);
	}
	
      }
    }
    
    int ncount_mcell_cluster = 0;
    for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
      ncount_mcell_cluster += (*it)->get_allcell().size();
    }
    ncount_mcell += allmcell.size();
    cout << i << " " << allmcell.size()  << " " << cluster_set.size()  << endl;
  }
  
  int ncount_mcell_cluster = 0;
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    ncount_mcell_cluster += (*it)->get_allcell().size();
  }
  cout << "Summary: " << ncount << " " << ncount_mcell << " " << ncount_mcell_cluster << endl;
  

  TFile *file = new TFile(Form("%d_%d.root",run_no,eve_num),"RECREATE");
  TTree *t_rec = new TTree("T_rec","T_rec");
  TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
  TTree *t_rec_charge_blob = new TTree("T_rec_charge_blob","T_rec_charge_blob");

  Double_t x_save, y_save, z_save;
  Double_t charge_save;
  Double_t ncharge_save;
  Double_t chi2_save;
  Double_t ndf_save;
  Double_t type_save;
 
  t_rec->SetDirectory(file);
  t_rec->Branch("x",&x_save,"x/D");
  t_rec->Branch("y",&y_save,"y/D");
  t_rec->Branch("z",&z_save,"z/D");
  
  t_rec_charge->SetDirectory(file);
  t_rec_charge->Branch("x",&x_save,"x/D");
  t_rec_charge->Branch("y",&y_save,"y/D");
  t_rec_charge->Branch("z",&z_save,"z/D");
  t_rec_charge->Branch("q",&charge_save,"q/D");
  t_rec_charge->Branch("nq",&ncharge_save,"nq/D");
  t_rec_charge->Branch("chi2",&chi2_save,"chi2/D");
  t_rec_charge->Branch("ndf",&ndf_save,"ndf/D");

  TTree *t_bad = new TTree("T_bad","T_bad");
  t_bad->SetDirectory(file);
  Int_t bad_npoints;
  Double_t bad_y[100],bad_z[100];
  t_bad->Branch("bad_npoints",&bad_npoints,"bad_npoints/I");
  t_bad->Branch("bad_y",bad_y,"bad_y[bad_npoints]/D");
  t_bad->Branch("bad_z",bad_z,"bad_z[bad_npoints]/D");
  
  for (int i=0; i!=badtiling[0]->get_cell_all().size();i++){
    const GeomCell *cell = badtiling[0]->get_cell_all().at(i);
    PointVector ps = cell->boundary();
    bad_npoints = ps.size();
    for (int j=0;j!=bad_npoints;j++){
      bad_y[j] = ps.at(j).y/units::cm;
      bad_z[j] = ps.at(j).z/units::cm;
    }
    t_bad->Fill();
  }

  //blob stuff
  t_rec_charge_blob->SetDirectory(file);
  t_rec_charge_blob->Branch("x",&x_save,"x/D");
  t_rec_charge_blob->Branch("y",&y_save,"y/D");
  t_rec_charge_blob->Branch("z",&z_save,"z/D");
  t_rec_charge_blob->Branch("q",&charge_save,"q/D");
  t_rec_charge_blob->Branch("nq",&ncharge_save,"nq/D");
  
  if (save_image_outline_flag==1){
    t_rec->Branch("type",&type_save,"type/D");
    t_rec_charge->Branch("type",&type_save,"type/D");
    t_rec_charge_blob->Branch("type",&type_save,"type/D");
  }


  TGraph2D *g = new TGraph2D();
  TGraph2D *g_rec = new TGraph2D();
  TGraph2D *g_rec_blob = new TGraph2D();

 
  //save results 
  for (int i=start_num;i!=end_num+1;i++){
    //recon 1
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];

      if (save_image_outline_flag==1){
	Point p = mcell->get_allcell().at(0)->center();
	x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.-time_offset*unit_dis/10.;
	y_save = p.y/units::cm;
	z_save = p.z/units::cm;
	g->SetPoint(ncount,x_save,y_save,z_save);
	ncount ++;
	type_save = 1; //center
	t_rec->Fill();
	
	for (int k=0;k!=mcell->get_edgecells().size();k++){
	  Point p = mcell->get_edgecells().at(k)->center();
	  y_save = p.y/units::cm;
    	  z_save = p.z/units::cm;
	  g->SetPoint(ncount,x_save,y_save,z_save);
	  ncount ++;
	  type_save = 2; //boundary ...
	  t_rec->Fill();
	}
	
	for (int k=0;k!=mcell->get_allcell().size();k++){
    	  Point p = mcell->get_allcell().at(k)->center();
	  x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.-time_offset*unit_dis/10.;
	  y_save = p.y/units::cm;
    	  z_save = p.z/units::cm;
	  type_save = 0;
	  t_rec->Fill();
	 
	}
      }else{
	for (int k=0;k!=mcell->get_allcell().size();k++){
    	  Point p = mcell->get_allcell().at(k)->center();
	  x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.-time_offset*unit_dis/10.;
	  y_save = p.y/units::cm;
    	  z_save = p.z/units::cm;
	  g->SetPoint(ncount,x_save,y_save,z_save);
	  ncount ++;
	  t_rec->Fill();
	 
	}
      }
    }
    // GeomCellSelection allcell = toytiling[i]->get_allcell();
    // for (int j=0;j!=allcell.size();j++){
    //   Point p = allcell[j]->center();
    //   x_save = i*unit_dis/10.*nrebin/2. - unit_dis/10.*frame_length/2. -time_offset*unit_dis/10.;
    //   y_save = p.y/units::cm;
    //   z_save = p.z/units::cm;
      

    //   g->SetPoint(ncount,x_save,y_save,z_save);
    //   t_rec->Fill();

    //   ncount ++;
    // }

    //recon 2 with charge
    //GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
      if (charge> recon_threshold || toymatrix[i]->Get_Solve_Flag()==0){
		if (toymatrix[i]->Get_Solve_Flag()==0)
	  charge = toytiling[i]->get_ave_charge();

	if (save_image_outline_flag==1){
	  Point p = mcell->get_allcell().at(0)->center();
	  x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.-time_offset*unit_dis/10.;
	  y_save = p.y/units::cm;
	  z_save = p.z/units::cm;
	  charge_save = charge;
	  ncharge_save = mcell->get_allcell().size();
	  chi2_save = toymatrix[i]->Get_Chi2();
	  ndf_save = toymatrix[i]->Get_ndf();
	  g_rec->SetPoint(ncount1,x_save,y_save,z_save);
	  ncount1 ++;
	  type_save = 1; //center
	  t_rec_charge->Fill();

	  charge_save = 0;
	  ncharge_save = 0;
	  chi2_save= 0;
	  ndf_save = 0;

	  for (int k=0;k!=mcell->get_edgecells().size();k++){
	    Point p = mcell->get_edgecells().at(k)->center();
	    y_save = p.y/units::cm;
	    z_save = p.z/units::cm;
	    g_rec->SetPoint(ncount1,x_save,y_save,z_save);
	    ncount1 ++;
	    type_save = 2; //boundary ...
	    t_rec_charge->Fill();
	  }

	   for (int k=0;k!=mcell->get_allcell().size();k++){
	    Point p = mcell->get_allcell().at(k)->center();
	    x_save = i*unit_dis/10.*nrebin/2. - unit_dis/10.*frame_length/2. -time_offset*unit_dis/10.;
	    y_save = p.y/units::cm;
	    z_save = p.z/units::cm;
	    charge_save = charge/mcell->get_allcell().size();
	    ncharge_save = mcell->get_allcell().size();
	    chi2_save = toymatrix[i]->Get_Chi2();
	    ndf_save = toymatrix[i]->Get_ndf();
	    type_save = 0;
	    t_rec_charge->Fill();
	    
	  }
	}else{	
	  for (int k=0;k!=mcell->get_allcell().size();k++){
	    Point p = mcell->get_allcell().at(k)->center();
	    x_save = i*unit_dis/10.*nrebin/2. - unit_dis/10.*frame_length/2. -time_offset*unit_dis/10.;
	    y_save = p.y/units::cm;
	    z_save = p.z/units::cm;
	    charge_save = charge/mcell->get_allcell().size();
	    ncharge_save = mcell->get_allcell().size();
	    chi2_save = toymatrix[i]->Get_Chi2();
	    ndf_save = toymatrix[i]->Get_ndf();
	    
	    g_rec->SetPoint(ncount1,x_save,y_save,z_save);
	    t_rec_charge->Fill();
	    
	    ncount1 ++;
	  }
	}
      }
    }

    
  //   //recon 3 with charge and deblob
  //   for (int j=0;j!=allmcell.size();j++){
  //     MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
  //     double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
  //     //if (charge> recon_threshold && !(mcell->IsSimpleBlob() && mcell->IsBlob())){
  //     if (charge> recon_threshold ){
  // 	for (int k=0;k!=mcell->get_allcell().size();k++){
  // 	  Point p = mcell->get_allcell().at(k)->center();
  // 	  x_save = i*unit_dis/10.*nrebin/2. - unit_dis/10.*frame_length/2. -time_offset*unit_dis/10.;
  // 	  y_save = p.y/units::cm;
  // 	  z_save = p.z/units::cm;
  // 	  charge_save = charge/mcell->get_allcell().size();
  // 	  ncharge_save = mcell->get_allcell().size();
	  
  // 	  g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
  // 	  t_rec_charge_blob->Fill();
	  
  // 	  ncount2 ++;
  // 	}
  //     }
  // }
  //    // if (toymatrix[i]->GetSimpleBlobReduction()){
  //    //   for (int j=0;j!=blobtiling[i]->Get_Cells().size();j++){
  //    // 	 const GeomCell *cell = blobtiling[i]->Get_Cells().at(j);
  //    // 	 Point p = cell->center();
  //    // 	 x_save = i*unit_dis/10.*nrebin/2. - unit_dis/10.*frame_length/2.-time_shift*unit_dis/10.;
  //    // 	 y_save = p.y/units::cm;
  //    // 	 z_save = p.z/units::cm;
  //    // 	 charge_save = blobtiling[i]->Get_Cell_Charge(cell,1);
  //    // 	 ncharge_save = 1;
	 
  //    // 	 g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
  //    // 	 t_rec_charge_blob->Fill();
	  
  //    // 	 ncount2 ++;
  //    //   }
  //    // }
    
  //   //save all results
  //   // file->Write(Form("toytiling_%d",i),toytiling[i]);
  //   // file->Write(Form("mergetiling_%d",i),mergetiling[i]);
  //   // file->Write(Form("truthtiling_%d",i),truthtiling[i]);
  //   // file->Write(Form("toymatrix_%d",i),toymatrix[i]);

  }
 
  

 
  
  Double_t x,y,z;
  //save cluster
  // int ncluster = 0;
  // for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  //   ncount = 0;
  //   TGraph2D *g1 = new TGraph2D();
  //   for (int i=0; i!=(*it)->get_allcell().size();i++){
  //     const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
  //     for (int j=0; j!=mcell->get_allcell().size();j++){
  // 	Point p = mcell->get_allcell().at(j)->center();
  // 	x = mcell->GetTimeSlice()*unit_dis/10.*nrebin/2. - unit_dis/10.0*frame_length/2.-time_offset*unit_dis/10.;
  // 	y = p.y/units::cm;
  // 	z = p.z/units::cm;
  // 	g1->SetPoint(ncount,x,y,z);
  // 	ncount ++;
  //     }
  //   }
  //   // cout << ncount << endl;
  //   g1->Write(Form("cluster_%d",ncluster));
  //   ncluster ++;
  // }

  

  // save all the toy tiling stuff
  // WireCell2dToy::ToyTiling* tt1 = 0;
  int time_slice;
  
  // TTree* ttree = new TTree("T","T");
  // ttree->Branch("time_slice",&time_slice,"time_slice/I");
  // ttree->Branch("toytiling",&tt1);
  // ttree->SetDirectory(file);
  // for (int i=start_num;i!=end_num+1;i++){
  //   tt1 = toytiling[i];
  //   time_slice = i;
  //   ttree->Fill();
  // }
  // ttree->Write();


  // TTree *ttree2 = new TTree("TB","TB");
  // // save center, and boundar points
  // Double_t x_center, y_center, z_center;
  // Double_t x_b[20], y_b[20], z_b[20];
  // Int_t npoints_b;
  // // save all the wires 
  // Int_t ub_index[2400], vb_index[2400], wb_index[3500];
  // Int_t unc, vnc, wnc;
  
  // ttree2->SetDirectory(file);

  // ttree2->Branch("time_slice",&time_slice,"time_slice/I"); // done
  
  // ttree2->Branch("x_center",&x_center,"x_center/D");
  // ttree2->Branch("y_center",&y_center,"y_center/D");
  // ttree2->Branch("z_center",&z_center,"z_center/D");
  
  // ttree2->Branch("npoints_b",&npoints_b,"npoints_b/I");
  // ttree2->Branch("x_b",x_b,"x_b[npoints_b]/D");
  // ttree2->Branch("y_b",y_b,"y_b[npoints_b]/D");
  // ttree2->Branch("z_b",z_b,"z_b[npoints_b]/D");

  // ttree2->Branch("unc",&unc,"unc/I");
  // ttree2->Branch("vnc",&vnc,"vnc/I");
  // ttree2->Branch("wnc",&wnc,"wnc/I");

  // ttree2->Branch("ub_index",ub_index,"ub_index[unc]");
  // ttree2->Branch("vb_index",vb_index,"vb_index[vnc]");
  // ttree2->Branch("wb_index",wb_index,"wb_index[wnc]");
  
  
  

  // for (int i=start_num;i!=end_num+1;i++){
  //   time_slice = i;
  //   for (int j=0;j!= badtiling[i]->get_cell_all().size();j++){
  //     const GeomCell *cell = badtiling[i]->get_cell_all().at(j);
  //     x_center = i*unit_dis/10.*nrebin/2. - unit_dis/10.*frame_length/2. -time_offset*unit_dis/10.;
  //     y_center = cell->center().y/units::cm;
  //     z_center = cell->center().z/units::cm;

  //     npoints_b = cell->boundary().size();
  //     for (int k=0;k!=npoints_b;k++){
  // 	x_b[k] = i*unit_dis/10.*nrebin/2. - unit_dis/10.*frame_length/2. -time_offset*unit_dis/10.;
  // 	y_b[k] = cell->boundary().at(k).y/units::cm;
  // 	z_b[k] = cell->boundary().at(k).z/units::cm;
  //     }

  //     unc = 0;
  //     vnc = 0;
  //     wnc = 0;
      
  //     for (int k=0;k!=badtiling[i]->cmap()[cell].size();k++){
  // 	MergeGeomWire *mwire = (MergeGeomWire*)badtiling[i]->cmap()[cell].at(k);
  // 	for (int k1 = 0; k1 != mwire->get_allwire().size(); k1++){
  // 	  GeomWire *wire = (GeomWire*) mwire->get_allwire().at(k1);
  // 	  if (wire->plane()==0){
  // 	    ub_index[unc] = wire->index();
  // 	    unc ++;
  // 	  }else if (wire->plane()==1){
  // 	    vb_index[vnc] = wire->index();
  // 	    vnc++;
  // 	  }else if (wire->plane()==2){
  // 	    wb_index[wnc] = wire->index();
  // 	    wnc++;
  // 	  }
  // 	}
  //     }
      
  //     ttree2->Fill();
      
  //   }
  // }
  

  

  TTree *ttree1 = new TTree("TC","TC");
  // To save cluster, we need to save
  // 1. time slice
  // 2. single cell
  // 3. charge
  // 4. cluster number
  const GeomCell* cell_save = 0;
  int cluster_num = -1;
  int mcell_id = -1;
  
  ttree1->Branch("time_slice",&time_slice,"time_slice/I"); // done
  ttree1->Branch("cell",&cell_save);
  ttree1->Branch("ncluster",&cluster_num,"cluster_num/I"); //done
  ttree1->Branch("mcell_id",&mcell_id,"mcell_id/I");
  ttree1->Branch("charge",&charge_save,"charge/D"); 
  
  double xx;
  //ttree1->Branch("xx",&xx,"xx/D");    //don
  ttree1->Branch("xx",&x,"xx/D");    //done
  ttree1->Branch("yy",&y,"yy/D");
  ttree1->Branch("zz",&z,"zz/D");

  // save information to reconstruct the toytiling
  int u_index, v_index, w_index;
  double u_charge, v_charge, w_charge;
  double u_charge_err, v_charge_err, w_charge_err;
  
  int apa_no=0, cryostat_no=0;
  ttree1->Branch("apa_no",&apa_no,"apa_no/I");
  int face = 0;
  ttree1->Branch("face",&face,"face/I");
  ttree1->Branch("cryostat_no",&cryostat_no,"cryostat_no/I");

  ttree1->Branch("u_index",&u_index,"u_index/I");
  ttree1->Branch("v_index",&v_index,"v_index/I");
  ttree1->Branch("w_index",&w_index,"w_index/I");
  
  ttree1->Branch("u_charge",&u_charge,"u_charge/D");
  ttree1->Branch("v_charge",&v_charge,"v_charge/D");
  ttree1->Branch("w_charge",&w_charge,"w_charge/D");

  ttree1->Branch("u_charge_err",&u_charge_err,"u_charge_err/D");
  ttree1->Branch("v_charge_err",&v_charge_err,"v_charge_err/D");
  ttree1->Branch("w_charge_err",&w_charge_err,"w_charge_err/D");

  //end save 

  ttree1->SetDirectory(file);
  
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    cluster_num ++;
    //loop merged cell
    for (int i=0; i!=(*it)->get_allcell().size();i++){
      const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
      mcell_id ++;
      time_slice = mcell->GetTimeSlice();
      x = time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.-time_offset*unit_dis/10.;
      xx =x;
      //loop single cell
      for (int j=0; j!=mcell->get_allcell().size();j++){
	cell_save = mcell->get_allcell().at(j);

	//fill the information needed for toytiling
	GeomWireSelection wires = toytiling[time_slice]->wires(*cell_save);
	
	//	if (i==0 && j==0) cout << "abc: " << time_slice << " " << toytiling[time_slice]->get_allcell().size() << " " << wires.size() << endl;

	for (int k=0;k!=wires.size();k++){
	  const GeomWire *wire = wires.at(k);
	  WirePlaneType_t plane = wire->plane();
	  if (plane==0){
	    u_index = wire->index();
	    u_charge = toytiling[time_slice]->wcmap()[wire];
	    u_charge_err = toytiling[time_slice]->wcemap()[wire];
	  }else if (plane==1){
	    v_index = wires.at(k)->index();
	    v_charge = toytiling[time_slice]->wcmap()[wire];
	    v_charge_err = toytiling[time_slice]->wcemap()[wire];
	  }else if (plane==2){
	    w_index = wire->index();
	    w_charge = toytiling[time_slice]->wcmap()[wire];
	    w_charge_err = toytiling[time_slice]->wcemap()[wire];
	  }
	}
	//end fill

	Point p = mcell->get_allcell().at(j)->center();
	charge_save = toymatrix[time_slice]->Get_Cell_Charge(mcell,1)/mcell->cross_section() * cell_save->cross_section();
	y = p.y/units::cm;
  	z = p.z/units::cm;
	ttree1->Fill();

	if (save_image_outline_flag==0){
	  //save the g_rec_blob tree ... 
	  x_save = x;
	  y_save = y;
	  z_save = z;
	  ncharge_save = mcell->get_allcell().size();
	  g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
	  ncount2++;
	  t_rec_charge_blob->Fill();
	}else if (save_image_outline_flag==1){
	  //save the g_rec_blob tree ... 
	  x_save = x;
	  y_save = y;
	  z_save = z;
	  ncharge_save = mcell->get_allcell().size();
	  type_save = 0;
	  t_rec_charge_blob->Fill();
	}
      }

       if (save_image_outline_flag==1){
	Point p = mcell->get_allcell().at(0)->center();
	x_save = time_slice *nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
	y_save = mcell->center().y/units::cm;
	z_save = mcell->center().z/units::cm;
	charge_save = toymatrix[time_slice]->Get_Cell_Charge(mcell,1);
	ncharge_save = mcell->get_allcell().size();
	g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
	ncount2++;
	type_save = 1;
	t_rec_charge_blob->Fill();
	
	charge_save = 0;
	ncharge_save = 0;
	
	for (int k=0;k!=mcell->get_edgecells().size();k++){
	  Point p = mcell->get_edgecells().at(k)->center();
	  y_save = p.y/units::cm;
	  z_save = p.z/units::cm;
	  g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
	  ncount2 ++;
	  type_save = 2; //boundary ...
	  t_rec_charge_blob->Fill();
	}
      }

    }
  }
  ttree1->Write();
   g->Write("shower3D");
  g_rec->Write("shower3D_charge");
  g_rec_blob->Write("shower3D_charge_blob");
  
  TTree *Trun = new TTree("Trun","Trun");
  Trun->SetDirectory(file);

  int detector = 0; // MicroBooNE
  Trun->Branch("detector",&detector,"detector/I");

  Trun->Branch("eventNo",&event_no,"eventNo/I");
  Trun->Branch("runNo",&run_no,"runNo/I");
  Trun->Branch("subRunNo",&subrun_no,"runRunNo/I");
  
  Trun->Branch("unit_dis",&unit_dis,"unit_dis/F");
  Trun->Branch("toffset_uv",&toffset_1,"toffset_uv/F");
  Trun->Branch("toffset_uw",&toffset_2,"toffset_uw/F");
  Trun->Branch("toffset_u",&toffset_3,"toffset_u/F");
  Trun->Branch("total_time_bin",&total_time_bin,"total_time_bin/I");
  Trun->Branch("recon_threshold",&recon_threshold,"recon_threshold/I");
  Trun->Branch("frame_length",&frame_length,"frame_length/I");
  Trun->Branch("max_events",&max_events,"max_events/I");
  Trun->Branch("eve_num",&eve_num,"eve_num/I");
  Trun->Branch("nrebin",&nrebin,"nrebin/I");
  Trun->Branch("threshold_u",&threshold_u,"threshold_u/F");
  Trun->Branch("threshold_v",&threshold_v,"threshold_v/F");
  Trun->Branch("threshold_w",&threshold_w,"threshold_w/F");
  Trun->Branch("time_offset",&time_offset,"time_offset/I");

   pitch_u = pitch_u/units::cm;
  pitch_v = pitch_v/units::cm;
  pitch_w = pitch_w/units::cm;
  Trun->Branch("pitch_u",&pitch_u,"pitch_u/D");
  Trun->Branch("pitch_v",&pitch_v,"pitch_v/D");
  Trun->Branch("pitch_w",&pitch_w,"pitch_w/D");



  Trun->Fill();


  file->Write();
  file->Close();

  

  return 0;
  
} // main()
