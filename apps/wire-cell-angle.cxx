#include "WireCellSst/GeomDataSource.h"
#include "WireCellSst/ToyuBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ClusterDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCell2dToy/SimpleBlobToyTiling.h"

#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/ToyMatrixExclusive.h"
#include "WireCell2dToy/ToyMatrixKalman.h"
#include "WireCell2dToy/ToyMatrixIterate.h"
#include "WireCell2dToy/ToyMatrixMarkov.h"
#include "WireCell2dToy/ToyMetric.h"
#include "WireCell2dToy/BlobMetric.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"

#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"

#include "WireCellData/GeomCluster.h"

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/SimDataSource.h"
#include "WireCellNav/SliceDataSource.h"
#include "WireCellSst/Util.h"
#include "WireCellData/SimTruth.h"
#include "WireCell2dToy/ToyDepositor.h"
#include "WireCellNav/GenerativeFDS.h"
#include "WireCell2dToy/ToySignalSimu.h"
#include "WireCell2dToy/ToySignalSimuTrue.h"
#include "WireCell2dToy/ToySignalGaus.h"
#include "WireCell2dToy/ToySignalWien.h"


#include "WireCellData/SpaceCell.h"
#include "WireCellData/MergeSpaceCell.h"
#include "WireCell2dToy/ToyCrawler.h"


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
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num -x[x_center] -y[y_center] -z[z_center] -a[rotate_angle]" << endl;
    return 1;
  }
  
  float x_center = 0;
  float y_center = 0;
  float z_center = 0;
  float rotate_angle = 0;
  
  for(Int_t i = 1; i != argc; i++){
     switch(argv[i][1]){
     case 'x':
       x_center = atof(&argv[i][2]); //m
       break;
     case 'y':
       y_center = atof(&argv[i][2]); //m
       break;
     case 'z':
       z_center = atof(&argv[i][2]); // m
       break;
     case 'a':
       rotate_angle = atof(&argv[i][2]); //degrees 
       break;
     }
  }



  x_center *= units::m;
  y_center *= units::m;
  z_center *= units::m;
  rotate_angle = rotate_angle/180.*3.1415926;

  
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
  
  
  const char* root_file = argv[2];
  const char* tpath = "/Event/Sim";
  
  int recon_threshold = 2000;
  int max_events = 100;
  int eve_num = atoi(argv[3]);

  float unit_dis = 1.6;

  float toffset_1=1.647;
  float toffset_2=1.539+1.647;
  float toffset_3=0;

  int total_time_bin=9600;
  int frame_length = 3200;
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

  int time_offset = 0;

  TFile *tfile = TFile::Open(root_file);

  TTree* sst = dynamic_cast<TTree*>(tfile->Get(tpath));

  int run_no, subrun_no, event_no;
  sst->SetBranchAddress("eventNo",&event_no);
  sst->SetBranchAddress("runNo",&run_no);
  sst->SetBranchAddress("subRunNo",&subrun_no);
  sst->GetEntry(eve_num);

  cout << "Run No: " << run_no << " " << subrun_no << " " << eve_num << endl;

  WireCell::FrameDataSource* fds = 0;
  fds = WireCellSst::make_fds(*tfile);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  
  TH1::AddDirectory(kFALSE);
  
  
 
  
  WireCell::ToyDepositor *toydep = new WireCell::ToyDepositor(fds,0,unit_dis,frame_length,x_center,y_center,z_center,rotate_angle);
  const PointValueVector& pvv = toydep->depositions(eve_num);
  
  WireCell::GenerativeFDS *gfds = new WireCell::GenerativeFDS(*toydep,gds,total_time_bin,max_events,0.5*unit_dis*units::millimeter); // 87 K at 0.5 kV/cm
  gfds->jump(eve_num);

  
  int abc;
  
  tfile->Close("R");
  delete tfile;

  
  cout << "Put in Truth " << endl; 
  WireCell2dToy::ToySignalSimuTrueFDS *st_fds = new WireCell2dToy::ToySignalSimuTrueFDS(*gfds,gds,total_time_bin/nrebin,max_events,0); //truth
  st_fds->jump(eve_num);
  
  cout << "Simulate Raw WaveForm " << endl; 
  WireCell2dToy::ToySignalSimuFDS *simu_fds = new WireCell2dToy::ToySignalSimuFDS(*gfds,gds,total_time_bin,max_events,toffset_1,toffset_2,1); // time offset among different planes for the time electrons travel among different planes
  simu_fds->jump(eve_num);
  
  cout << "Deconvolution with Gaussian filter" << endl;
  WireCell2dToy::ToySignalGausFDS *gaus_fds = new WireCell2dToy::ToySignalGausFDS(*simu_fds,gds,total_time_bin/nrebin,max_events,toffset_1,toffset_2); // gaussian smearing for charge estimation
  gaus_fds->jump(eve_num);
  
  cout << "Deconvolution with Wiener filter" << endl;
   WireCell2dToy::ToySignalWienFDS *wien_fds = new WireCell2dToy::ToySignalWienFDS(*simu_fds,gds,total_time_bin/nrebin,max_events,toffset_1,toffset_2); // weiner smearing for hit identification
  wien_fds->jump(eve_num);
  
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  int nwire_u = wires_u.size();
  int nwire_v = wires_v.size();
  int nwire_w = wires_w.size();
  
 
  
  // cin >> abc;

  WireCellSst::ToyuBooNESliceDataSource *sds = new WireCellSst::ToyuBooNESliceDataSource(*wien_fds,*gaus_fds,threshold_u, 
  					    threshold_v, threshold_w, 
  					    threshold_ug, 
  					    threshold_vg, threshold_wg, 
  					    nwire_u, 
  					    nwire_v, nwire_w); 

  WireCellSst::ToyuBooNESliceDataSource *sds_th = new WireCellSst::ToyuBooNESliceDataSource(*st_fds,*st_fds,500, 
  					    500, 500, 
  					    threshold_ug, 
  					    threshold_vg, threshold_wg, 
  					    nwire_u, 
  					    nwire_v, nwire_w); 

 
  int ncount = 0;
  int ncount1 = 0;  
  int ncount2 = 0;

  int ncount_t = 0;
  
  
  

  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
  WireCell2dToy::MergeToyTiling **mergetiling = new WireCell2dToy::MergeToyTiling*[2400];
  WireCell2dToy::TruthToyTiling **truthtiling = new WireCell2dToy::TruthToyTiling*[2400];
  WireCell2dToy::SimpleBlobToyTiling **blobtiling = new WireCell2dToy::SimpleBlobToyTiling*[2400];
  WireCell2dToy::ToyMatrix **toymatrix = new WireCell2dToy::ToyMatrix*[2400];
  //save truth ...
  WireCell2dToy::ToyTiling **toytiling_th = new WireCell2dToy::ToyTiling*[2400];
  WireCell2dToy::TruthToyTiling **truthtiling_th = new WireCell2dToy::TruthToyTiling*[2400];
  WireCell2dToy::ToyMetric toymetric;
  WireCell2dToy::BlobMetric blobmetric;



  
  int start_num = 0 ;
  int end_num = sds->size()-1;
  // int start_num = 800 ;
  // int end_num = 810;
  // GeomCellSelection total_cells;
  // GeomCellSelection total_edge_cells;

  cout << "Start the Tiling " << endl; 
  WireCell::Slice slice;
  for (int i=start_num;i!=end_num+1;i++){
    
    sds->jump(i);
    sds_th->jump(i);
    slice = sds->get();
    WireCell::Slice slice_th = sds_th->get();
    // cout << i << " " << slice.group().size() << " " << slice_th.group().size() << endl;
    
    toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
    GeomCellSelection allcell = toytiling[i]->get_allcell();
    
    GeomWireSelection allwire = toytiling[i]->get_allwire();
    cout << "Single Cell: " << i << " "  << allcell.size() << " " << allwire.size() << endl;
    //    mergetiling[i] = new WireCell2dToy::MergeToyTiling(*toytiling[i],i,3,1);
    mergetiling[i] = new WireCell2dToy::MergeToyTiling(*toytiling[i],i,3);
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    
    GeomWireSelection allmwire = mergetiling[i]->get_allwire();
    cout <<"Blob: " << i << " " << allmcell.size() << " " << allmwire.size() << endl;
    truthtiling[i] = new WireCell2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds,frame_length/nrebin,unit_dis);
    
    
    toytiling_th[i] = new WireCell2dToy::ToyTiling(slice_th,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
    truthtiling_th[i] = new WireCell2dToy::TruthToyTiling(*toytiling_th[i],pvv,i,gds,frame_length/nrebin,unit_dis);
    CellChargeMap ccmap = truthtiling[i]->ccmap();
        
    // for (int j=0;j!=allcell.size();j++){
    //   total_cells.push_back(allcell.at(j));
    // }
    // for (int j=0;j!=allmcell.size();j++){
    //   const MergeGeomCell *mcell = (const MergeGeomCell*)allmcell.at(j);
    //   GeomCellSelection edgecells = mcell->get_edgecells();
    //   for (int k=0;k!=edgecells.size();k++){
    // 	total_edge_cells.push_back(edgecells.at(k));
    //   }
    // }

  }

  delete sds;
  delete sds_th;
  
  delete simu_fds;
  delete gaus_fds;
  delete wien_fds;

  delete st_fds;
  delete gfds;
  delete toydep;
  delete fds;


  
    //do display
  
  // TApplication theApp("theApp",&argc,argv);
  // theApp.SetReturnFromRun(true);
  
  // TCanvas c1("ToyMC","ToyMC",800,600);
  // WireCell2dToy::ToyEventDisplay display(c1, gds);
  
  // display.charge_min = 0.;
  // display.charge_max = 1.;
  
  
  
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
  // display.draw_slice(slice,"");
  // display.draw_cells(total_cells,"*same");
  // display.draw_cells(total_edge_cells,"*same",2);
  // theApp.Run();



  


   for (int i=start_num;i!=end_num+1;i++){
     toymatrix[i] = new WireCell2dToy::ToyMatrix(*toytiling[i],*mergetiling[i]);
    // cout << "start the iterate " << endl; 
    if (toymatrix[i]->Get_Solve_Flag()==0){
      WireCell2dToy::ToyMatrixIterate toymatrix_it(*toymatrix[i],2000,1e5);
    }
    
    cout << i << " chi2: " << toymatrix[i]->Get_Chi2() <<
      " NDF: " << toymatrix[i]->Get_ndf() << endl;
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    CellChargeMap ccmap = truthtiling[i]->ccmap();
    if (toymatrix[i]->Get_Solve_Flag()!=0)
      toymetric.Add(allmcell,*toymatrix[i],ccmap);
    
    toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
   }

   toymetric.Print();




  //add in cluster
  std::cout << "Start Clustering " << std::endl;
  GeomClusterList cluster_list, cluster_dellist;
  int ncount_mcell = 0;
  //Now do cluster
  for (int i=start_num;i!=end_num+1;i++){
    //for (int i=800;i!=810;i++){
    GeomCellSelection pallmcell = mergetiling[i]->get_allcell();
    GeomCellSelection allmcell;
    for (int j=0;j!=pallmcell.size();j++){
      const GeomCell* mcell = pallmcell[j];
      
      //hack for now
      if (toymatrix[i]->Get_Solve_Flag()!=0){
      	if (toymatrix[i]->Get_Cell_Charge(mcell)> recon_threshold){
      	  allmcell.push_back(mcell);
      	}
      }else{
      	allmcell.push_back(mcell);
      }
      
      //allmcell.push_back(mcell);
    }
    
    
    if (cluster_list.empty()){
      // if cluster is empty, just insert all the mcell, each as a cluster
      for (int j=0;j!=allmcell.size();j++){
	MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
	if (mcell->get_allcell().size()>0){
	  GeomCluster *cluster = new GeomCluster(*mcell);
	  cluster_list.push_back(cluster);
	}
      }
    }else{
      for (int j=0;j!=allmcell.size();j++){
	MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
	if (mcell->get_allcell().size()>0){
	  int flag = 0;
	  int flag_save = 0;
	  GeomCluster *cluster_save = 0;
	  cluster_dellist.clear();
	  
	  // loop through merged cell
	  int tmp_num = 0;
	  for (auto it = cluster_list.begin();it!=cluster_list.end();it++){
	    //loop through clusters
	    
	    flag += (*it)->AddCell(*mcell);
	    // std::cout << i << " " << j << " " << tmp_num << " " << flag << std::endl;
	    tmp_num ++;
	    
	    if (flag==1 && flag != flag_save){
	      cluster_save = *it;
	    }else if (flag>1 && flag != flag_save){
	      cluster_save->MergeCluster(*(*it));
	      cluster_dellist.push_back(*it);
	    }
	    flag_save = flag;
	    
	  }
	  
	  for (auto it = cluster_dellist.begin();it!=cluster_dellist.end();it++){
	    auto it1 = find(cluster_list.begin(),cluster_list.end(),*it);
	    cluster_list.erase(it1);
	    delete (*it);
	  }
	  
	  if (flag==0){
	    GeomCluster *cluster = new GeomCluster(*mcell);
	    cluster_list.push_back(cluster);
	  }
	}
	
      }
    }
    
    // int ncount_mcell_cluster = 0;
    // for (auto it = cluster_list.begin();it!=cluster_list.end();it++){
    //   ncount_mcell_cluster += (*it)->get_allcell().size();
    // }
    ncount_mcell += allmcell.size();
    //cout << i << " " << allmcell.size()  << " " << cluster_set.size()  << endl;
  }
  
  int ncount_mcell_cluster = 0;
  for (auto it = cluster_list.begin();it!=cluster_list.end();it++){
    ncount_mcell_cluster += (*it)->get_allcell().size();
  }
  cout << "Summary: " << ncount << " " << ncount_mcell << " " << ncount_mcell_cluster << endl;
  
  int ncluster;

   

  // start crawler
  cout << "Start Crawling " << endl;
  std::vector<WireCell2dToy::ToyCrawler*> crawlers;
  MergeSpaceCellSelection all_msc_cells;
  SpaceCellSelection all_sc_cells;
  
  ncluster = 0;
  for (auto it = cluster_list.begin();it!=cluster_list.end();it++){
    
    MergeSpaceCellSelection mscells;
    for (int i=0; i!=(*it)->get_allcell().size();i++){
      const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
      MergeSpaceCell *mscell = new MergeSpaceCell();
      all_msc_cells.push_back(mscell);
      mscell->set_mcell(mcell);
      for (int j=0;j!=mcell->get_allcell().size();j++){
  	const GeomCell *cell = mcell->get_allcell().at(j);
  	SpaceCell *space_cell = new SpaceCell(ncluster,*cell,(mcell->GetTimeSlice()*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.)*units::cm,1,unit_dis/10.*nrebin/2.*units::cm);
  	all_sc_cells.push_back(space_cell);
	mscell->AddSpaceCell(space_cell);
      }
      mscells.push_back(mscell);
    }
    WireCell2dToy::ToyCrawler* toycrawler = new WireCell2dToy::ToyCrawler(mscells);
    crawlers.push_back(toycrawler);
    
    // std::cout << ncluster << " " << toycrawler->Get_mcells_map().size() << " " << toycrawler->Get_allCT().size() << " " << toycrawler->Get_allMCT().size()  << std::endl;

    // if (toycrawler->Get_mcells_map().size()>200){
    //   TApplication theApp("theApp",&argc,argv);
    //   theApp.SetReturnFromRun(true);
      
    //   TCanvas c1("ToyMC","ToyMC",800,600);
    //   c1.Draw();
      
    //   WireCell2dToy::ClusterDisplay display(c1);
    //   //display.DrawCluster(cells);
    //   display.DrawCluster(mscells);
      
      
    //   display.DrawCrawler(*toycrawler,"psame",1);
      
    //   theApp.Run();
    // }


    ncluster ++;
  }
  
  
  //check # of clusters 
  int sum = 0 ;
  for (int i=0;i!=crawlers.size();i++){
    for (auto it = crawlers.at(i)->Get_mcells_map().begin(); it!= crawlers.at(i)->Get_mcells_map().end();it++){
      MergeSpaceCell *mcell1 = it->first;
      sum += mcell1->Get_all_spacecell().size();
    }
  }
  
  std::cout << "Check: " << crawlers.size() << " "  << sum << std::endl;

  

  //start to prepare the important merge cell vectors
  std::vector<GeomCellSelection> Good_MCells;
  for (int i=start_num;i!=end_num+1;i++){
    GeomCellSelection cells;
    Good_MCells.push_back(cells);
  }
 
  for (int i=0;i!=crawlers.size();i++){
    WireCell2dToy::ToyCrawler *toycrawler = crawlers.at(i);
    for (int j=0;j!=toycrawler->Get_allMCT().size();j++){
      MergeClusterTrack *mct = toycrawler->Get_allMCT().at(j);
      int ntime = mct->Get_TimeLength();
      int ntime_flag = 0;

      if (ntime < 6 && ntime >= 2){
	ntime_flag = 1;
	for (int k=0;k!=ntime;k++){
	  MergeSpaceCellSelection cells = mct->Get_MSCS(k);
	  ntime_flag = 0;
	  for (int i1 = 0; i1!=cells.size(); i1++){
	    if (cells.at(i1)->Get_all_spacecell().size() > 300/ntime){
	      ntime_flag = 1;
	      break;
	    }
	  }
	  if (ntime_flag == 0)
	    break;
	}
      }else{
	ntime_flag = 1;
      }
      

      if (ntime_flag){
	// do something	
	for (int k=0;k!=ntime;k++){
	  MergeSpaceCellSelection cells = mct->Get_MSCS(k);
	  int time = mct->Get_Time(k);
	  if (cells.size()==1){
	    Good_MCells.at(time).push_back(cells.at(0)->get_mcell());  
	  }else if (cells.size()>1){
	    
	    MergeSpaceCell *cell = cells.at(0);
	    MergeSpaceCell *next_cell = cells.at(1);
	    for (int i1 = 1; i1!=cells.size();i1++){
	      if (cell->Get_all_spacecell().size() < cells.at(i1)->Get_all_spacecell().size()){
	   	cell = cells.at(i1);
		next_cell = cell;
	      }
	    }
	    if (cell->Get_all_spacecell().size() > 3* next_cell->Get_all_spacecell().size())
	      Good_MCells.at(time).push_back(cell->get_mcell());
	  }
	}
      }
    }
  }
  

  //start to solve matrix ...
  for (int i=start_num;i!=end_num+1;i++){
    //   //toymatrix[i] = new WireCell2dToy::ToyMatrix(*toytiling[i],*mergetiling[i]);
    //  // cout << "start the iterate " << endl; 
    if (toymatrix[i]->Get_Solve_Flag()==0){
      WireCell2dToy::ToyMatrixIterate toymatrix_it(*toymatrix[i], *mergetiling[i], Good_MCells.at(i-start_num));
      cout << i << " chi2: " << toymatrix[i]->Get_Chi2() <<
      " NDF: " << toymatrix[i]->Get_ndf() << endl;
    }
    //  GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    //  CellChargeMap ccmap = truthtiling[i]->ccmap();
    //  if (toymatrix[i]->Get_Solve_Flag()!=0)
    //    toymetric.Add(allmcell,*toymatrix[i],ccmap);
    
    //  toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
  }

  toymetric.Print();
  // std::cout << "Starting MCMC" << std::endl;
  
//     //without  time information
//   // for (int i=start_num;i!=end_num+1;i++){
//   //   if (toymatrix[i]->Get_Solve_Flag()==0){
//   //     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
//   //     WireCell2dToy::ToyMatrixMarkov toymatrix_markov(*toymatrix[i],*mergetiling[i],&allmcell,Good_MCells.at(i-start_num));
//   //     //WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i],&allmcell);
//   //     CellChargeMap ccmap = truthtiling[i]->ccmap();
//   //     if (toymatrix[i]->Get_Solve_Flag()!=0)
//   // 	toymetric.Add(allmcell,*toymatrix[i],ccmap);
//   //     toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
//   //     cout << " chi2: " << i << " " << toymatrix[i]->Get_Chi2() 
//   // 	   << " NDF: " << toymatrix[i]->Get_ndf() << endl;
//   //   }
//   // }


// //with time information
//   if (start_num != end_num){
//     int first_solve=-1;
//     for (int i=start_num; i!=end_num+1;i++){
//       if (toymatrix[i]->Get_Solve_Flag()!=0){
//   	first_solve = i;
//   	break;
//       }
//     }
//     if (first_solve <0){
//       for (int i=start_num;i!=end_num+1;i++){
//   	if (toymatrix[i]->Get_Solve_Flag()==0){
//   	  GeomCellSelection allmcell = mergetiling[i]->get_allcell();
//   	  WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i],&allmcell);
//   	  CellChargeMap ccmap = truthtiling[i]->ccmap();
//   	  if (toymatrix[i]->Get_Solve_Flag()!=0)
//   	    toymetric.Add(allmcell,*toymatrix[i],ccmap);
//   	  toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
//   	  cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
//   	  cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
//   	}
//       }
//     }else{
//       for (int i=first_solve+1;i<=end_num-1;i++){
//   	if (toymatrix[i]->Get_Solve_Flag()==0){
//   	  GeomCellSelection allmcell = mergetiling[i]->get_allcell();
//   	  WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell);
//   	  CellChargeMap ccmap = truthtiling[i]->ccmap();
//   	  if (toymatrix[i]->Get_Solve_Flag()!=0)
//   	    toymetric.Add(allmcell,*toymatrix[i],ccmap);
//   	  toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	  
//   	  cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
//   	  cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
	  
//   	  //toymetric.Print();
//   	}
//       }
      
//       if (toymatrix[end_num]->Get_Solve_Flag()==0){
//   	GeomCellSelection allmcell = mergetiling[end_num]->get_allcell();
//   	WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[end_num-1],toymatrix[end_num],0,mergetiling[end_num-1],mergetiling[end_num],0,&allmcell);
	
	
//   	CellChargeMap ccmap = truthtiling[end_num]->ccmap();
//   	if (toymatrix[end_num]->Get_Solve_Flag()!=0)
//   	  toymetric.Add(allmcell,*toymatrix[end_num],ccmap);
//   	toymetric.AddSolve(toymatrix[end_num]->Get_Solve_Flag());
	
//   	cout << "chi2: " << end_num << " " << toymatrix[end_num]->Get_Chi2() << endl;
//   	cout << "NDF: " << toymatrix[end_num]->Get_ndf() << endl;
//       }
      
//       // go to early ones 
//       for (int i=first_solve-1;i>=start_num+1;i--){
//   	if (toymatrix[i]->Get_Solve_Flag()==0){
//   	  GeomCellSelection allmcell = mergetiling[i]->get_allcell();
//   	  WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell);
	  
//   	  CellChargeMap ccmap = truthtiling[i]->ccmap();
//   	  if (toymatrix[i]->Get_Solve_Flag()!=0)
//   	    toymetric.Add(allmcell,*toymatrix[i],ccmap);
//   	  toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	  
//   	  cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
//   	  cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
//   	}
//       }
      
//       if (toymatrix[start_num]->Get_Solve_Flag()==0){
//   	GeomCellSelection allmcell = mergetiling[start_num]->get_allcell();
//   	WireCell2dToy::ToyMatrixMarkov toymatrix_markov(0,toymatrix[start_num],toymatrix[start_num+1],0,mergetiling[start_num],mergetiling[start_num+1],&allmcell);
	
	
//   	CellChargeMap ccmap = truthtiling[start_num]->ccmap();
//   	if (toymatrix[start_num]->Get_Solve_Flag()!=0)
//   	  toymetric.Add(allmcell,*toymatrix[start_num],ccmap);
//   	toymetric.AddSolve(toymatrix[start_num]->Get_Solve_Flag());
	
//   	cout << "chi2: " << start_num << " " << toymatrix[start_num]->Get_Chi2() << endl;
//   	cout << "NDF: " << toymatrix[start_num]->Get_ndf() << endl;
//       }
//     }
//   }


  // reset existing cluster and redo it ...
  std::cout << "Start Clustering after solving equations " << std::endl;
  cluster_dellist.clear();
  for (auto it = cluster_list.begin();it!=cluster_list.end();it++){
    delete (*it);
  }
  cluster_list.clear();
  
  ncount_mcell = 0;

  //Now do cluster
  for (int i=start_num;i!=end_num+1;i++){
    GeomCellSelection pallmcell = mergetiling[i]->get_allcell();
    GeomCellSelection allmcell;
    for (int j=0;j!=pallmcell.size();j++){
      const GeomCell* mcell = pallmcell[j];
      if (toymatrix[i]->Get_Cell_Charge(mcell)> recon_threshold || toymatrix[i]->Get_Solve_Flag()==0){
	allmcell.push_back(mcell);
      }
    }
    
    
    if (cluster_list.empty()){
      // if cluster is empty, just insert all the mcell, each as a cluster
      for (int j=0;j!=allmcell.size();j++){
	MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
	if (mcell->get_allcell().size()>0){
	  GeomCluster *cluster = new GeomCluster(*mcell);
	  cluster_list.push_back(cluster);
	}
      }
    }else{
      for (int j=0;j!=allmcell.size();j++){
	MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
	if (mcell->get_allcell().size()>0){
	  int flag = 0;
	  int flag_save = 0;
	  GeomCluster *cluster_save = 0;
	  cluster_dellist.clear();
	  
	  // loop through merged cell
	  int tmp_num = 0;
	  for (auto it = cluster_list.begin();it!=cluster_list.end();it++){
	    //loop through clusters
	    
	    flag += (*it)->AddCell(*mcell);
	    // std::cout << i << " " << j << " " << tmp_num << " " << flag << std::endl;
	    tmp_num ++;
	    
	    if (flag==1 && flag != flag_save){
	      cluster_save = *it;
	    }else if (flag>1 && flag != flag_save){
	      cluster_save->MergeCluster(*(*it));
	      cluster_dellist.push_back(*it);
	    }
	    flag_save = flag;
	    
	  }
	  
	  for (auto it = cluster_dellist.begin();it!=cluster_dellist.end();it++){
	    auto it1 = find(cluster_list.begin(),cluster_list.end(),*it);
	    cluster_list.erase(it1);
	    delete (*it);
	  }
	  
	  if (flag==0){
	    GeomCluster *cluster = new GeomCluster(*mcell);
	    cluster_list.push_back(cluster);
	  }
	}
	
      }
    }
    
    // int ncount_mcell_cluster = 0;
    // for (auto it = cluster_list.begin();it!=cluster_list.end();it++){
    //   ncount_mcell_cluster += (*it)->get_allcell().size();
    // }
    ncount_mcell += allmcell.size();
    //cout << i << " " << allmcell.size()  << " " << cluster_set.size()  << endl;
  }
  
  ncount_mcell_cluster = 0;
  for (auto it = cluster_list.begin();it!=cluster_list.end();it++){
    ncount_mcell_cluster += (*it)->get_allcell().size();
  }
  cout << "Summary: " << ncount << " " << ncount_mcell << " " << ncount_mcell_cluster << endl;
  
  // reset crawler ... 
  for (int i =0; i!=crawlers.size();i++){
    delete crawlers.at(i);
  }
  crawlers.clear();
  for (int i=0;i!=all_msc_cells.size();i++){
    delete all_msc_cells.at(i);
  }
  all_msc_cells.clear();
  for (int i=0;i!=all_sc_cells.size();i++){
    delete all_sc_cells.at(i);
  }
  all_sc_cells.clear();

  // start crawler ... 

  
  // start crawler
  cout << "Start Crawling after solving equations " << endl;
  
  ncluster = 0;
  for (auto it = cluster_list.begin();it!=cluster_list.end();it++){
    
    MergeSpaceCellSelection mscells;
    for (int i=0; i!=(*it)->get_allcell().size();i++){
      const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
      MergeSpaceCell *mscell = new MergeSpaceCell();
      all_msc_cells.push_back(mscell);
      mscell->set_mcell(mcell);
      for (int j=0;j!=mcell->get_allcell().size();j++){
  	const GeomCell *cell = mcell->get_allcell().at(j);
  	SpaceCell *space_cell = new SpaceCell(ncluster,*cell,(mcell->GetTimeSlice()*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.)*units::cm,1,unit_dis/10.*nrebin/2.*units::cm);
	all_sc_cells.push_back(space_cell);
  	mscell->AddSpaceCell(space_cell);
      }
      mscells.push_back(mscell);
    }
    WireCell2dToy::ToyCrawler* toycrawler = new WireCell2dToy::ToyCrawler(mscells,2);
    crawlers.push_back(toycrawler);
    
    // std::cout << ncluster << " " << toycrawler->Get_mcells_map().size() << " " << toycrawler->Get_allCT().size() << " " << toycrawler->Get_allMCT().size()  << std::endl;

    // if (toycrawler->Get_mcells_map().size()>200){
    //   TApplication theApp("theApp",&argc,argv);
    //   theApp.SetReturnFromRun(true);
      
    //   TCanvas c1("ToyMC","ToyMC",800,600);
    //   c1.Draw();
      
    //   WireCell2dToy::ClusterDisplay display(c1);
    //   //display.DrawCluster(cells);
    //   display.DrawCluster(mscells);
      
      
    //   display.DrawCrawler(*toycrawler,"psame",1);
      
    //   theApp.Run();
    // }


    ncluster ++;
  }
  
  
  //check # of clusters 
  sum = 0 ;
  for (int i=0;i!=crawlers.size();i++){
    for (auto it = crawlers.at(i)->Get_mcells_map().begin(); it!= crawlers.at(i)->Get_mcells_map().end();it++){
      MergeSpaceCell *mcell1 = it->first;
      sum += mcell1->Get_all_spacecell().size();
    }
  }
  
  std::cout << "Check: " << crawlers.size() << " "  << sum << std::endl;


  //save files
  TFile *file = new TFile(Form("shower3D_cluster_%d_%d.root",eve_num, int(round(rotate_angle/3.1415926*180.))),"RECREATE");
  TTree *t_true = new TTree("T_true","T_true");
  TTree *t_rec = new TTree("T_rec","T_rec");
  TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
  TTree *t_rec_charge_blob = new TTree("T_rec_charge_blob","T_rec_charge_blob");

  Double_t x_save, y_save, z_save;
  Double_t charge_save;
  Double_t ncharge_save;
  Double_t chi2_save;
  Double_t ndf_save;
  
  t_true->SetDirectory(file);
  t_true->Branch("x",&x_save,"x/D");
  t_true->Branch("y",&y_save,"y/D");
  t_true->Branch("z",&z_save,"z/D");
  t_true->Branch("q",&charge_save,"q/D");
  
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

  //blob stuff
  t_rec_charge_blob->SetDirectory(file);
  t_rec_charge_blob->Branch("x",&x_save,"x/D");
  t_rec_charge_blob->Branch("y",&y_save,"y/D");
  t_rec_charge_blob->Branch("z",&z_save,"z/D");
  t_rec_charge_blob->Branch("q",&charge_save,"q/D");
  t_rec_charge_blob->Branch("nq",&ncharge_save,"nq/D");

  
  TGraph2D *g = new TGraph2D();
  TGraph2D *gt = new TGraph2D();
  TGraph2D *g_rec = new TGraph2D();
  TGraph2D *g_rec_blob = new TGraph2D();

  //save results 
  for (int i=start_num;i!=end_num+1;i++){
    //truth
    CellChargeMap ccmap = truthtiling_th[i]->ccmap();
    for (auto it = ccmap.begin();it!=ccmap.end(); it++){
      Point p = it->first->center();
      x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
      y_save = p.y/units::cm;
      z_save = p.z/units::cm;
      charge_save = it->second;
      
      gt->SetPoint(ncount_t,x_save,y_save,z_save);
      t_true->Fill();
      
      ncount_t ++;
    }
    
    //recon 1
    GeomCellSelection allcell = toytiling[i]->get_allcell();
    for (int j=0;j!=allcell.size();j++){
      Point p = allcell[j]->center();
      x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
      y_save = p.y/units::cm;
      z_save = p.z/units::cm;
      

      g->SetPoint(ncount,x_save,y_save,z_save);
      t_rec->Fill();

      ncount ++;
    }
    
    //recon 2 with charge
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
      if (charge> recon_threshold || toymatrix[i]->Get_Solve_Flag()==0){

    	if (toymatrix[i]->Get_Solve_Flag()==0)
	  charge = toytiling[i]->get_ave_charge();

	//truth
    	for (int k=0;k!=mcell->get_allcell().size();k++){
    	  Point p = mcell->get_allcell().at(k)->center();
    	  x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
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
    
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
      if (charge> recon_threshold || toymatrix[i]->Get_Solve_Flag()==0){
	if (toymatrix[i]->Get_Solve_Flag()==0)
	  charge = toytiling[i]->get_ave_charge();

    	for (int k=0;k!=mcell->get_allcell().size();k++){
    	  Point p = mcell->get_allcell().at(k)->center();
    	  x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
    	  y_save = p.y/units::cm;
    	  z_save = p.z/units::cm;
    	  charge_save = charge/mcell->get_allcell().size();
    	  ncharge_save = mcell->get_allcell().size();
	  
    	  g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
    	  t_rec_charge_blob->Fill();
	  
    	  ncount2 ++;
    	}
      }
    }


  }
 
  

  g->Write("shower3D");
  gt->Write("shower3D_truth");
  g_rec->Write("shower3D_charge");
  g_rec_blob->Write("shower3D_charge_blob");
  
  // const int N = 100000;
  // Double_t x[N],y[N],z[N];
  Double_t x,y,z;
  //save cluster
  ncluster = 0;
  for (auto it = cluster_list.begin();it!=cluster_list.end();it++){
    ncount = 0;
    TGraph2D *g1 = new TGraph2D();
    for (int i=0; i!=(*it)->get_allcell().size();i++){
      const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
      for (int j=0; j!=mcell->get_allcell().size();j++){
  	Point p = mcell->get_allcell().at(j)->center();
  	x = mcell->GetTimeSlice()*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
  	y = p.y/units::cm;
  	z = p.z/units::cm;
  	g1->SetPoint(ncount,x,y,z);
  	ncount ++;
      }
    }
    // cout << ncount << endl;
    g1->Write(Form("cluster_%d",ncluster));
    ncluster ++;
  }


  // save all the toy tiling stuff
  WireCell2dToy::ToyTiling* tt1 = 0;
  int time_slice;
  
  TTree* ttree = new TTree("T","T");
  ttree->Branch("time_slice",&time_slice,"time_slice/I");
  ttree->Branch("toytiling",&tt1);
  ttree->SetDirectory(file);
  for (int i=start_num;i!=end_num+1;i++){
    tt1 = toytiling[i];
    time_slice = i;
    ttree->Fill();
  }
  ttree->Write();

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
  double xx,yy,zz;
  ttree1->Branch("xx",&xx,"xx/D");    //don
  ttree1->Branch("yy",&yy,"yy/D");    //don
  ttree1->Branch("zz",&zz,"zz/D");    //don

  // ttree1->Branch("x",&x_save,"x/D");    //done
  //ttree1->Branch("y",&y_save,"y/D");
  //ttree1->Branch("z",&z_save,"z/D");
  

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
  
  for (auto it = cluster_list.begin();it!=cluster_list.end();it++){
    cluster_num ++;
    //loop merged cell
    for (int i=0; i!=(*it)->get_allcell().size();i++){
      const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
      mcell_id ++;
      time_slice = mcell->GetTimeSlice();
      x_save = time_slice *nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
      xx = x_save;
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
	//hack for now
	charge_save = toymatrix[time_slice]->Get_Cell_Charge(mcell,1)/mcell->cross_section() * cell_save->cross_section();
	//charge_save = 1;
	
	//std::cout << time_slice << " " << x_save << std::endl;
	yy = p.y/units::cm;
  	zz = p.z/units::cm;
	ttree1->Fill();
	
      }
    }
  }
  ttree1->Write();
  
 
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

  toymetric.Print();
  blobmetric.Print();

  // cin >> abc;
  
  return 0;
  
} // main()
