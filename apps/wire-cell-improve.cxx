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
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num" << endl;
    return 1;
  }
  
  
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
  
  TFile *tfile = TFile::Open(root_file);

  WireCell::FrameDataSource* fds = 0;
  fds = WireCellSst::make_fds(*tfile);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  
  TH1::AddDirectory(kFALSE);
  

  int recon_threshold = 2000;
  int max_events = 100;
  int eve_num = atoi(argv[3]);
  
  WireCell::ToyDepositor *toydep = new WireCell::ToyDepositor(fds);
  const PointValueVector& pvv = toydep->depositions(eve_num);
  
  WireCell::GenerativeFDS *gfds = new WireCell::GenerativeFDS(*toydep,gds,9600,max_events,0.5*1.60*units::millimeter); // 87 K at 0.5 kV/cm
  gfds->jump(eve_num);

  
  int abc;
  
  tfile->Close("R");
  delete tfile;

  
  cout << "Put in Truth " << endl; 
  WireCell2dToy::ToySignalSimuTrueFDS *st_fds = new WireCell2dToy::ToySignalSimuTrueFDS(*gfds,gds,9600/4,max_events,0); //truth
  st_fds->jump(eve_num);
  
  cout << "Simulate Raw WaveForm " << endl; 
  WireCell2dToy::ToySignalSimuFDS *simu_fds = new WireCell2dToy::ToySignalSimuFDS(*gfds,gds,9600,max_events,1.647,1.539+1.647,1); // time offset among different planes for the time electrons travel among different planes
  simu_fds->jump(eve_num);
  
  cout << "Deconvolution with Gaussian filter" << endl;
  WireCell2dToy::ToySignalGausFDS *gaus_fds = new WireCell2dToy::ToySignalGausFDS(*simu_fds,gds,9600/4,max_events,1.647,1.539+1.647); // gaussian smearing for charge estimation
  gaus_fds->jump(eve_num);
  
  cout << "Deconvolution with Wiener filter" << endl;
   WireCell2dToy::ToySignalWienFDS *wien_fds = new WireCell2dToy::ToySignalWienFDS(*simu_fds,gds,9600/4,max_events,1.647,1.539+1.647); // weiner smearing for hit identification
  wien_fds->jump(eve_num);
  
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  int nwire_u = wires_u.size();
  int nwire_v = wires_v.size();
  int nwire_w = wires_w.size();
  
  float threshold_u = 5.87819e+02 * 4.0;
  float threshold_v = 8.36644e+02 * 4.0;
  float threshold_w = 5.67974e+02 * 4.0;

  float threshold_ug = 755.96;
  float threshold_vg = 822.81;
  float threshold_wg = 510.84;
  
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
    cout << i << " " << slice.group().size() << " " << slice_th.group().size() << endl;
    
    toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
    GeomCellSelection allcell = toytiling[i]->get_allcell();
    
    GeomWireSelection allwire = toytiling[i]->get_allwire();
    cout << "Single Cell: " << i << " "  << allcell.size() << " " << allwire.size() << endl;
    mergetiling[i] = new WireCell2dToy::MergeToyTiling(*toytiling[i],i,3,1);
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    
    GeomWireSelection allmwire = mergetiling[i]->get_allwire();
    cout <<"Blob: " << i << " " << allmcell.size() << " " << allmwire.size() << endl;
    truthtiling[i] = new WireCell2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds,800);
    
    
    toytiling_th[i] = new WireCell2dToy::ToyTiling(slice_th,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
    truthtiling_th[i] = new WireCell2dToy::TruthToyTiling(*toytiling_th[i],pvv,i,gds,800);
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
      //     if (toymatrix[i]->Get_Cell_Charge(mcell)> recon_threshold){
      allmcell.push_back(mcell);
      //}
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
  

   

  // start crawler
  cout << "Start Crawling " << endl;
  std::vector<WireCell2dToy::ToyCrawler*> crawlers;
  
  int ncluster = 0;
  for (auto it = cluster_list.begin();it!=cluster_list.end();it++){
    
    MergeSpaceCellSelection mscells;
    for (int i=0; i!=(*it)->get_allcell().size();i++){
      const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
      MergeSpaceCell *mscell = new MergeSpaceCell();
      mscell->set_mcell(mcell);
      for (int j=0;j!=mcell->get_allcell().size();j++){
  	const GeomCell *cell = mcell->get_allcell().at(j);
  	SpaceCell *space_cell = new SpaceCell(ncluster,*cell,(mcell->GetTimeSlice()*0.32-256)*units::cm,1,0.32*units::cm);
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

  


  //save files
  TFile *file = new TFile(Form("shower3D_cluster_%d.root",eve_num),"RECREATE");
  
  Double_t x_save, y_save, z_save;
  Double_t charge_save;
  Double_t ncharge_save;
  Double_t chi2_save;
  Double_t ndf_save;


  
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
      x_save = i*0.32 - 256;
      y_save = p.y/units::cm;
      z_save = p.z/units::cm;
      charge_save = it->second;
      
      gt->SetPoint(ncount_t,x_save,y_save,z_save);
      
      
      ncount_t ++;
    }
    
    //recon 1
    GeomCellSelection allcell = toytiling[i]->get_allcell();
    for (int j=0;j!=allcell.size();j++){
      Point p = allcell[j]->center();
      x_save = i*0.32- 256;
      y_save = p.y/units::cm;
      z_save = p.z/units::cm;
      

      g->SetPoint(ncount,x_save,y_save,z_save);
      

      ncount ++;
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
  	x = mcell->GetTimeSlice()*0.32- 256;
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
  ttree1->Branch("x",&x,"x/D");    //done
  ttree1->Branch("y",&y,"y/D");
  ttree1->Branch("z",&z,"z/D");
  
  ttree1->SetDirectory(file);
  
  for (auto it = cluster_list.begin();it!=cluster_list.end();it++){
    cluster_num ++;
    //loop merged cell
    for (int i=0; i!=(*it)->get_allcell().size();i++){
      const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
      mcell_id ++;
      time_slice = mcell->GetTimeSlice();
      x = time_slice *0.32- 256;
      //loop single cell
      for (int j=0; j!=mcell->get_allcell().size();j++){
	cell_save = mcell->get_allcell().at(j);
	Point p = mcell->get_allcell().at(j)->center();
	charge_save = 1;
	y = p.y/units::cm;
  	z = p.z/units::cm;
	ttree1->Fill();
	
      }
    }
  }
  ttree1->Write();
  
 


  file->Write();
  file->Close();

  // cin >> abc;
  
  return 0;
  
} // main()
