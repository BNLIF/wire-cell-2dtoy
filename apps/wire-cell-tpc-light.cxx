#include "WireCellSst/GeomDataSource.h"
#include "WireCellData/PR3DCluster.h"
#include "WireCellData/SlimMergeGeomCell.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"
#include "WireCell2dToy/ExecMon.h"
#include "WireCell2dToy/CalcPoints.h"
#include "WireCell2dToy/ToyClustering.h"
#include "WireCell2dToy/uBooNE_light_reco.h"

#include "WireCell2dToy/ToyMatching.h"


#include "TH1F.h"
#include "TH2F.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace WireCell;
using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 4) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/imaging.root light.root " << endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);

  int flag_pos_corr = 1; // correct X position after matching ... 
  for(Int_t i = 1; i != argc; i++){
     switch(argv[i][1]){
     case 'c':
       flag_pos_corr = atoi(&argv[i][2]); 
       break;
     }
  }
  
  
  ExecMon em("starting");
  cerr << em("load geometry") << endl;

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



  // test geometry ...
  const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),0);
  const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),0);
  const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),0);
  double first_u_dis = gds.wire_dist(*uwire) ;
  double first_v_dis = gds.wire_dist(*vwire) ;
  double first_w_dis = gds.wire_dist(*wwire) ; 
  
  
  TString filename = argv[2];
  TFile *file = new TFile(filename);
  TTree *Trun = (TTree*)file->Get("Trun");

  int run_no, subrun_no, event_no;
  int time_offset;
  int nrebin;
  int frame_length;
  int eve_num;
  float unit_dis;
  
  Trun->SetBranchAddress("eventNo",&event_no);
  Trun->SetBranchAddress("runNo",&run_no);
  Trun->SetBranchAddress("subRunNo",&subrun_no);
  Trun->SetBranchAddress("unit_dis",&unit_dis);
  Trun->SetBranchAddress("frame_length",&frame_length);
  Trun->SetBranchAddress("eve_num",&eve_num);
  Trun->SetBranchAddress("nrebin",&nrebin);
  Trun->SetBranchAddress("time_offset",&time_offset);
  Trun->GetEntry(0);

  //std::cout << nrebin << " " << time_offset << std::endl;
  
  // define singleton ... 
  TPCParams& mp = Singleton<TPCParams>::Instance();
  
  double pitch_u = gds.pitch(WirePlaneType_t(0));
  double pitch_v = gds.pitch(WirePlaneType_t(1));
  double pitch_w = gds.pitch(WirePlaneType_t(2));
  double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;

  double angle_u = gds.angle(WirePlaneType_t(0));
  double angle_v = gds.angle(WirePlaneType_t(1));
  double angle_w = gds.angle(WirePlaneType_t(2));

  //std::cout << angle_u << " " << angle_v << " " << angle_w << std::endl;
  
  mp.set_pitch_u(pitch_u);
  mp.set_pitch_v(pitch_v);
  mp.set_pitch_w(pitch_w);
  mp.set_angle_u(angle_u);
  mp.set_angle_v(angle_v);
  mp.set_angle_w(angle_w);
  mp.set_ts_width(time_slice_width);
  mp.set_first_u_dis(first_u_dis);
  mp.set_first_v_dis(first_v_dis);
  mp.set_first_w_dis(first_w_dis);
  
  std::set<int> dead_u_index;
  std::set<int> dead_v_index;
  std::set<int> dead_w_index;
  // load mcell
  
  TTree *TC = (TTree*)file->Get("TC");
  Int_t cluster_id;
  Int_t time_slice;
  Double_t q, uq, vq, wq, udq, vdq, wdq;
  Int_t nwire_u=0, flag_u; //number of wires, dead?
  Int_t nwire_v=0, flag_v;
  Int_t nwire_w=0, flag_w;
  Int_t wire_index_u[2400];
  Int_t wire_index_v[2400];
  Int_t wire_index_w[3256];
  Double_t wire_charge_u[2400];
  Double_t wire_charge_v[2400];
  Double_t wire_charge_w[2400];
  Double_t wire_charge_err_u[2400];
  Double_t wire_charge_err_v[2400];
  Double_t wire_charge_err_w[2400];

  TC->SetBranchAddress("cluster_id",&cluster_id);
  TC->SetBranchAddress("time_slice",&time_slice);
  TC->SetBranchAddress("q",&q);
  TC->SetBranchAddress("uq",&uq);
  TC->SetBranchAddress("vq",&vq);
  TC->SetBranchAddress("wq",&wq);
  TC->SetBranchAddress("udq",&udq);
  TC->SetBranchAddress("vdq",&vdq);
  TC->SetBranchAddress("wdq",&wdq);
  TC->SetBranchAddress("nwire_u",&nwire_u);
  TC->SetBranchAddress("nwire_v",&nwire_v);
  TC->SetBranchAddress("nwire_w",&nwire_w);
  TC->SetBranchAddress("flag_u",&flag_u);
  TC->SetBranchAddress("flag_v",&flag_v);
  TC->SetBranchAddress("flag_w",&flag_w);
  TC->SetBranchAddress("wire_index_u",wire_index_u);
  TC->SetBranchAddress("wire_index_v",wire_index_v);
  TC->SetBranchAddress("wire_index_w",wire_index_w);
  TC->SetBranchAddress("wire_charge_u",wire_charge_u);
  TC->SetBranchAddress("wire_charge_v",wire_charge_v);
  TC->SetBranchAddress("wire_charge_w",wire_charge_w);
  TC->SetBranchAddress("wire_charge_err_u",wire_charge_err_u);
  TC->SetBranchAddress("wire_charge_err_v",wire_charge_err_v);
  TC->SetBranchAddress("wire_charge_err_w",wire_charge_err_w);


  //load mcell
  
  TTree *TDC = (TTree*)file->Get("TDC");
  int ntime_slice = 0,time_slices[2400];
  TDC->SetBranchAddress("cluster_id",&cluster_id);
  TDC->SetBranchAddress("ntime_slice",&ntime_slice);
  TDC->SetBranchAddress("time_slice",time_slices);
  
  TDC->SetBranchAddress("nwire_u",&nwire_u);
  TDC->SetBranchAddress("nwire_v",&nwire_v);
  TDC->SetBranchAddress("nwire_w",&nwire_w);
  TDC->SetBranchAddress("flag_u",&flag_u);
  TDC->SetBranchAddress("flag_v",&flag_v);
  TDC->SetBranchAddress("flag_w",&flag_w);
  TDC->SetBranchAddress("wire_index_u",wire_index_u);
  TDC->SetBranchAddress("wire_index_v",wire_index_v);
  TDC->SetBranchAddress("wire_index_w",wire_index_w);

  // load cells ... 
  GeomCellSelection mcells;
  PR3DClusterSelection live_clusters;
  PR3DClusterSelection dead_clusters;
  PR3DCluster *cluster;
  int prev_cluster_id=-1;
  int ident = 0;
  for (int i=0;i!=TC->GetEntries();i++){
    TC->GetEntry(i);
    SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
    mcell->SetTimeSlice(time_slice);

    mcell->set_uq(uq);
    mcell->set_vq(vq);
    mcell->set_wq(wq);

    mcell->set_udq(udq);
    mcell->set_vdq(vdq);
    mcell->set_wdq(wdq);

    mcell->set_q(q);
    if (flag_u==0){
      mcell->add_bad_planes(WirePlaneType_t(0));
      for (int i=0;i!=nwire_u;i++){
	dead_u_index.insert(wire_index_u[i]);
      }
    }
    if (flag_v==0){
      mcell->add_bad_planes(WirePlaneType_t(1));
      for (int i=0;i!=nwire_v;i++){
	dead_v_index.insert(wire_index_v[i]);
      }
    }
    if (flag_w==0){
      mcell->add_bad_planes(WirePlaneType_t(2));
      for (int i=0;i!=nwire_w;i++){
	dead_w_index.insert(wire_index_w[i]);
      }
    }
    for (int i=0;i!=nwire_u;i++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(0),wire_index_u[i]);
      mcell->AddWire(wire,WirePlaneType_t(0),wire_charge_u[i],wire_charge_err_u[i]);
    }
    for (int i=0;i!=nwire_v;i++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(1),wire_index_v[i]);
      mcell->AddWire(wire,WirePlaneType_t(1),wire_charge_v[i],wire_charge_err_v[i]);
    }
    for (int i=0;i!=nwire_w;i++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(2),wire_index_w[i]);
      mcell->AddWire(wire,WirePlaneType_t(2),wire_charge_w[i],wire_charge_err_w[i]);
    }
    mcells.push_back(mcell);

    if (cluster_id != prev_cluster_id){
      cluster = new PR3DCluster(cluster_id);
      live_clusters.push_back(cluster);
    }
    cluster->AddCell(mcell,time_slice);

    prev_cluster_id = cluster_id;
    ident++;
  }
  //  std::cout << live_clusters.size() << std::endl;

  prev_cluster_id = -1;
  // TDC
   for (int i=0;i!=TDC->GetEntries();i++){
    TDC->GetEntry(i);

    SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
    mcell->SetTimeSlice(time_slices[0]);

    if (flag_u==0){
      mcell->add_bad_planes(WirePlaneType_t(0));
      for (int i=0;i!=nwire_u;i++){
	dead_u_index.insert(wire_index_u[i]);
      }
    }
    if (flag_v==0){
      mcell->add_bad_planes(WirePlaneType_t(1));
      for (int i=0;i!=nwire_v;i++){
	dead_v_index.insert(wire_index_v[i]);
      }
    }
    if (flag_w==0){
      mcell->add_bad_planes(WirePlaneType_t(2));
      for (int i=0;i!=nwire_w;i++){
	dead_w_index.insert(wire_index_w[i]);
      }
    }
    for (int i=0;i!=nwire_u;i++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(0),wire_index_u[i]);
      mcell->AddWire(wire,WirePlaneType_t(0));
    }
    for (int i=0;i!=nwire_v;i++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(1),wire_index_v[i]);
      mcell->AddWire(wire,WirePlaneType_t(1));
    }
    for (int i=0;i!=nwire_w;i++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(2),wire_index_w[i]);
      mcell->AddWire(wire,WirePlaneType_t(2));
    }
    mcells.push_back(mcell);

    if (cluster_id!=prev_cluster_id){
      cluster = new PR3DCluster(cluster_id);
      dead_clusters.push_back(cluster);
    }
    for (int i=0;i!=ntime_slice;i++){
      cluster->AddCell(mcell,time_slices[i]);
    }
    
    prev_cluster_id=cluster_id;
    ident++;
  }

   // std::cout << live_clusters.size() << std::endl;
   // for (size_t i=0;i!=live_clusters.size();i++){
   //   std::cout << live_clusters.at(i)->get_cluster_id() << " " 
   // 	       << live_clusters.at(i)->get_num_mcells() << " "
   // 	       << live_clusters.at(i)->get_num_time_slices() << std::endl;
   // }
   // std::cout << dead_clusters.size() << std::endl;
   for (size_t i=0;i!=dead_clusters.size();i++){
     dead_clusters.at(i)->Remove_duplicated_mcells();
     // std::cout << dead_clusters.at(i)->get_cluster_id() << " " 
     // 	       << dead_clusters.at(i)->get_num_mcells() << " "
     // 	       << dead_clusters.at(i)->get_num_time_slices() << std::endl;
   }
  
   
   cerr << em("load clusters from file") << endl;

  

   // Start to add X, Y, Z points
   // form boundaries for the bad cells ... 
   for (size_t j = 0; j!= dead_clusters.size(); j++){
     WireCell2dToy::calc_boundary_points_dead(gds,dead_clusters.at(j));
   }
   // form sampling points for the normal cells ...
   for (size_t i=0; i!=live_clusters.size();i++){
     WireCell2dToy::calc_sampling_points(gds,live_clusters.at(i),nrebin, frame_length, unit_dis);
     live_clusters.at(i)->Calc_PCA();
   }
   cerr << em("Add X, Y, Z points") << std::endl;

   
   // WireCell2dToy::Clustering_live_dead(live_clusters, dead_clusters);
   // cerr << em("Clustering live and dead clusters") << std::endl;

   WireCell2dToy::Clustering_jump_gap_cosmics(live_clusters, dead_clusters,dead_u_index, dead_v_index, dead_w_index);
   cerr << em("Clustering to jump gap in cosmics") << std::endl;

   // need to further cluster things ...
   std::map<PR3DCluster*,std::vector<std::pair<PR3DCluster*,double>>> group_clusters =  WireCell2dToy::Clustering_isolated(live_clusters);
   cerr << em("Clustering isolated") << std::endl;
   
   
   
   // for (size_t i=0;i!=live_clusters.size();i++){
   //   //    std::cout << live_clusters.at(i)->get_mcells().size() << " " << live_clusters.at(i)->get_num_time_slices() << std::endl;
   //   live_clusters.at(i)->Create_graph();
   //   std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = live_clusters.at(i)->get_highest_lowest_wcps();
   //   live_clusters.at(i)->dijkstra_shortest_paths(wcps.first);
   //   live_clusters.at(i)->cal_shortest_path(wcps.second);
   //   live_clusters.at(i)->fine_tracking();
   // }
   // cerr << em("Create Graph in all clusters") << std::endl;
   
   // processing light information
   const char* root_file = argv[3];
   WireCell2dToy::uBooNE_light_reco uboone_flash(root_file);
   uboone_flash.load_event_raw(0);
   TFile *file2 = new TFile(root_file);
   TTree *T = (TTree*)file2->Get("/Event/Sim");
   
   TClonesArray* cosmic_hg_wf = new TClonesArray;
   TClonesArray* cosmic_lg_wf = new TClonesArray;
   TClonesArray* beam_hg_wf = new TClonesArray;
   TClonesArray* beam_lg_wf = new TClonesArray;
   vector<short> *cosmic_hg_opch = new vector<short>;
   vector<short> *cosmic_lg_opch = new vector<short>;
   vector<short> *beam_hg_opch = new vector<short>;
   vector<short> *beam_lg_opch = new vector<short>;
   vector<double> *cosmic_hg_timestamp = new vector<double>;
   vector<double> *cosmic_lg_timestamp = new vector<double>;
   vector<double> *beam_hg_timestamp = new vector<double>;
   vector<double> *beam_lg_timestamp = new vector<double>;
   vector<short> *opch_to_opdet = new vector<short>; 
   std::vector<float> *op_gain = new std::vector<float>;
   std::vector<float> *op_gainerror = new std::vector<float>; 
   double triggerTime;
   
   // for "saturation" waveforms
   TClonesArray* op_wf = new TClonesArray("TH1S");
   std::vector<short> *op_femch = new std::vector<short>;
   //std::vector<int> *op_femch = new std::vector<int>;
   std::vector<double> *op_timestamp = new std::vector<double>;
   /*
     std::vector<float> *op_gain = new std::vector<float>;
     std::vector<float> *op_gainerror = new std::vector<float>; 
     double triggerTime;
   */
   
   T->SetBranchAddress("cosmic_hg_wf",&cosmic_hg_wf);
   T->SetBranchAddress("cosmic_lg_wf",&cosmic_lg_wf);
   T->SetBranchAddress("beam_hg_wf",&beam_hg_wf);
   T->SetBranchAddress("beam_lg_wf",&beam_lg_wf);
   T->SetBranchAddress("cosmic_hg_opch",&cosmic_hg_opch);
   T->SetBranchAddress("cosmic_lg_opch",&cosmic_lg_opch);
   T->SetBranchAddress("beam_hg_opch",&beam_hg_opch);
   T->SetBranchAddress("beam_lg_opch",&beam_lg_opch);
   T->SetBranchAddress("cosmic_hg_timestamp",&cosmic_hg_timestamp);
   T->SetBranchAddress("cosmic_lg_timestamp",&cosmic_lg_timestamp);
   T->SetBranchAddress("beam_hg_timestamp",&beam_hg_timestamp);
   T->SetBranchAddress("beam_lg_timestamp",&beam_lg_timestamp);
   T->SetBranchAddress("opch_to_opdet",&opch_to_opdet); 
   T->SetBranchAddress("op_gain",&op_gain);
   T->SetBranchAddress("op_gainerror",&op_gainerror);
   /*
     T->SetBranchAddress("op_femch",&op_femch);
     T->SetBranchAddress("op_gain",&op_gain);
     T->SetBranchAddress("op_gainerror",&op_gainerror);
     T->SetBranchAddress("op_timestamp",&op_timestamp);
     T->SetBranchAddress("op_wf",&op_wf);
   */
   T->SetBranchAddress("triggerTime",&triggerTime);
   T->SetBranchStatus("eventNo",1);
   T->SetBranchAddress("eventNo" , &event_no);
   T->SetBranchStatus("runNo",1);
   T->SetBranchAddress("runNo"   , &run_no);
   T->SetBranchStatus("subRunNo",1);
   T->SetBranchAddress("subRunNo", &subrun_no);
   
   T->GetEntry(eve_num);
   

   // prepare light matching ....
   WireCell::OpflashSelection& flashes = uboone_flash.get_flashes();
   
   std::vector<std::tuple<WireCell::PR3DCluster*, WireCell::Opflash*, double, std::vector<double>>> matched_results = WireCell2dToy::tpc_light_match(time_offset,nrebin,group_clusters,flashes);
   cerr << em("TPC Light Matching") << std::endl;
   //


   
   TFile *file1 = new TFile(Form("match_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
   TTree *T_match = new TTree("T_match","T_match");
   T_match->SetDirectory(file1);
   Int_t ncluster=0;
   Int_t flash_id;
   Double_t strength;
   Double_t pe_pred[32];
   Double_t pe_meas[32];
   Double_t pe_meas_err[32];
   T_match->Branch("tpc_cluster_id",&ncluster,"tpc_cluster_id/I");
   T_match->Branch("flash_id",&flash_id,"flash_id/I");
   T_match->Branch("strength",&strength,"strength/D");
   T_match->Branch("pe_pred",pe_pred,"pe_pred[32]/D");
   T_match->Branch("pe_meas",pe_meas,"pe_meas[32]/D");
   T_match->Branch("pe_meas_err",pe_meas_err,"pe_meas_err[32]/D");

   for (auto it = matched_results.begin(); it!=matched_results.end(); it++){
     Opflash *flash = std::get<1>(*it);
     if (flash!=0){
       auto it1 = find(flashes.begin(),flashes.end(),flash);
       flash_id = it1 - flashes.begin();
       strength = std::get<2>(*it);
       std::vector<double> temp = std::get<3>(*it);
       for (int i=0;i!=32;i++){
     	 pe_pred[i] = temp.at(i);
     	 pe_meas[i] = flash->get_PE(i);
     	 pe_meas_err[i] = flash->get_PE_err(i);
       }
     }else{
       flash_id = -1;
       strength = 0;
       for (int i=0;i!=32;i++){
     	 pe_pred[i] = 0;
     	 pe_meas[i] = 0;
     	 pe_meas_err[i] = 0.;
       }
     }
     T_match->Fill();
     ncluster++;
   }
   
   
   
   
   TTree *t_bad = new TTree("T_bad","T_bad");
   t_bad->SetDirectory(file1);
   Int_t bad_npoints;
   
   Double_t bad_y[100],bad_z[100];
   t_bad->Branch("cluster_id",&ncluster,"cluster_id/I");
   t_bad->Branch("bad_npoints",&bad_npoints,"bad_npoints/I");
   t_bad->Branch("bad_y",bad_y,"bad_y[bad_npoints]/D");
   t_bad->Branch("bad_z",bad_z,"bad_z[bad_npoints]/D");
   
   for (size_t j = 0; j!= dead_clusters.size(); j++){
     SMGCSelection& mcells = dead_clusters.at(j)->get_mcells();
     ncluster = dead_clusters.at(j)->get_cluster_id();
     for (size_t i=0;i!=mcells.size();i++){
       PointVector ps = mcells.at(i)->boundary();
       bad_npoints = ps.size();
       for (int k=0;k!=bad_npoints;k++){
	 bad_y[k] = ps.at(k).y/units::cm;
	 bad_z[k] = ps.at(k).z/units::cm;
       }
       t_bad->Fill();
     }
   }
   
   TTree *T_cluster ;
   Double_t x,y,z;
   
   T_cluster = new TTree("T_cluster","T_cluster");
   T_cluster->Branch("cluster_id",&ncluster,"cluster_id/I");
   T_cluster->Branch("x",&x,"x/D");
   T_cluster->Branch("y",&y,"y/D");
   T_cluster->Branch("z",&z,"z/D");
   T_cluster->SetDirectory(file1);

   TTree *T_rec = new TTree("T_rec","T_rec");
   T_rec->Branch("x",&x,"x/D");
   T_rec->Branch("y",&y,"y/D");
   T_rec->Branch("z",&z,"z/D");
   T_rec->SetDirectory(file1);

   Double_t charge_save=1, ncharge_save=1, chi2_save=1, ndf_save=1;
   TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
   t_rec_charge->SetDirectory(file1);
   t_rec_charge->Branch("x",&x,"x/D");
   t_rec_charge->Branch("y",&y,"y/D");
   t_rec_charge->Branch("z",&z,"z/D");
   t_rec_charge->Branch("q",&charge_save,"q/D");
   t_rec_charge->Branch("nq",&ncharge_save,"nq/D");
   t_rec_charge->Branch("chi2",&chi2_save,"chi2/D");
   t_rec_charge->Branch("ndf",&ndf_save,"ndf/D");
   
   // test ... 
   // ncluster = 0;
   // x=  0;
   // y=0;
   // z=0;
   // T_cluster->Fill();

   
   // note did not save the unmatched cluster ... 
   ncluster = 0;
   for (auto it = matched_results.begin(); it!= matched_results.end(); it++){
     PR3DCluster *main_cluster = std::get<0>(*it);
     Opflash *flash = std::get<1>(*it);
     double offset_x ;
     if (flash!=0){
       offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
     }else{
       offset_x = 0;
     }
     
     if (flag_pos_corr==0)
       offset_x = 0;
     
   //   for (auto it = group_clusters.begin(); it!= group_clusters.end(); it++){
     PR3DClusterSelection temp_clusters;
     temp_clusters.push_back(main_cluster);
     for (auto it1 = group_clusters[main_cluster].begin(); it1!=group_clusters[main_cluster].end(); it1++){
       temp_clusters.push_back((*it1).first);
       //std::cout << (*it1).second/units::cm << std::endl;
     }
     for (size_t j = 0; j!= temp_clusters.size(); j++){
       SMGCSelection& mcells = temp_clusters.at(j)->get_mcells();
       //ncluster = temp_clusters.at(0)->get_cluster_id();
       for (size_t i=0;i!=mcells.size();i++){
	 PointVector ps = mcells.at(i)->get_sampling_points();
	 int time_slice = mcells.at(i)->GetTimeSlice();
	 if (ps.size()==0) std::cout << "zero sampling points!" << std::endl;
	 for (int k=0;k!=ps.size();k++){
	   x = (ps.at(k).x- offset_x)/units::cm ;
	   y = ps.at(k).y/units::cm;
	   z = ps.at(k).z/units::cm;
	   T_cluster->Fill();
	 }
       }
     }
     ncluster ++;
   }
   
   for (size_t j = 0; j!= live_clusters.size(); j++){
     // SMGCSelection& mcells = live_clusters.at(j)->get_mcells();
     // ncluster = live_clusters.at(j)->get_cluster_id();
     // for (size_t i=0;i!=mcells.size();i++){
     //   PointVector ps = mcells.at(i)->get_sampling_points();
     //   int time_slice = mcells.at(i)->GetTimeSlice();
     //   if (ps.size()==0) std::cout << "zero sampling points!" << std::endl;
     //   for (int k=0;k!=ps.size();k++){
     // 	 x = ps.at(k).x/units::cm;//time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
     // 	 y = ps.at(k).y/units::cm;
     // 	 z = ps.at(k).z/units::cm;
     // 	 T_cluster->Fill();
     //   }
     // }

     // save wcps
     std::list<WCPointCloud<double>::WCPoint>& wcps_list = live_clusters.at(j)->get_path_wcps();
     //ncluster = -1 * ncluster-100;
     for (auto it = wcps_list.begin(); it!=wcps_list.end(); it++){
       x = (*it).x/units::cm;
       y = (*it).y/units::cm;
       z = (*it).z/units::cm;
       T_rec->Fill();
     }

     PointVector& pts = live_clusters.at(j)->get_fine_tracking_path();
     for (size_t i=0; i!=pts.size(); i++){
       x = pts.at(i).x/units::cm;
       y = pts.at(i).y/units::cm;
       z = pts.at(i).z/units::cm;
       t_rec_charge->Fill();
     }

     

     // // save mcells
     // std::list<SlimMergeGeomCell*>& mcells_list = live_clusters.at(j)->get_path_mcells();
     // ncluster = -1 * ncluster-100;
     // for (auto it = mcells_list.begin(); it!=mcells_list.end(); it++){
     //   Point p = (*it)->center();
     //   x = p.x/units::cm;
     //   y = p.y/units::cm;
     //   z = p.z/units::cm;
     //   T_cluster->Fill();
     // }
     
     
     // if (live_clusters.at(j)->get_num_mcells()>30){
     //   // add PCA axis point
     //   Vector center = live_clusters.at(j)->get_center();
     //   Vector dir = live_clusters.at(j)->get_PCA_axis(0);
     //   for (int i=-200;i!=200;i++){
     // 	 x = (center.x + dir.x *(i*units::cm) )/units::cm;
     // 	 y = (center.y + dir.y *(i*units::cm) )/units::cm;
     // 	 z = (center.z + dir.z *(i*units::cm) )/units::cm;
     // 	 T_cluster->Fill();
     //   }
     // }
    }

   // ncluster = 0;
   // for (auto it = dead_live_cluster_mapping.begin(); it!= dead_live_cluster_mapping.end(); it++){
   //   std::vector<PR3DCluster*> clusters = (*it).second;
   //   if (clusters.size()>1){
   //     //std::cout << clusters.size() << std::endl;
   //     for (auto it1 = clusters.begin(); it1!=clusters.end(); it1++){
   // 	 PR3DCluster* cluster = (*it1);
   // 	 ncluster = cluster->get_cluster_id();
   // 	 SMGCSelection& mcells = cluster->get_mcells();
   // 	 for (size_t i=0;i!=mcells.size();i++){
   // 	   PointVector ps = mcells.at(i)->get_sampling_points();
   // 	   for (int k=0;k!=ps.size();k++){
   // 	     x = ps.at(k).x/units::cm;//time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
   // 	     y = ps.at(k).y/units::cm;
   // 	     z = ps.at(k).z/units::cm;
   // 	     T_cluster->Fill();
   // 	   }
   // 	 }
   //     }
   //   }
   //   // ncluster++;
   // }
   // for (auto it = dead_live_mcells_mapping.begin(); it!= dead_live_mcells_mapping.end(); it++){
   //   std::vector<std::vector<SlimMergeGeomCell*>> mcellss = (*it).second;
   //   // std::cout << mcellss.size() << std::endl;
   //   if (mcellss.size()>1){
   //     for (auto it1 = mcellss.begin(); it1!=mcellss.end(); it1++){
   // 	 std::vector<SlimMergeGeomCell*> mcells = (*it1);
   // 	 for (size_t i=0;i!=mcells.size();i++){
   // 	   PointVector ps = mcells.at(i)->get_sampling_points();
   // 	   for (int k=0;k!=ps.size();k++){
   // 	     x = ps.at(k).x/units::cm;//time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
   // 	     y = ps.at(k).y/units::cm;
   // 	     z = ps.at(k).z/units::cm;
   // 	     T_cluster->Fill();
   // 	   }
   // 	 }
   //     }
   //   }
   //   ncluster++;
   // }
   
   Trun->CloneTree()->Write();


   
   TTree *t1 = new TTree("T_data","T_data");
   t1->SetDirectory(file1);

  t1->Branch("op_gain",&op_gain);
  t1->Branch("op_gainerror",&op_gainerror);
  t1->Branch("op_femch",&op_femch);
  t1->Branch("op_timestamp",&op_timestamp);
  t1->Branch("op_wf",&op_wf,256000,0);
  t1->Branch("triggerTime",&triggerTime);
  
  t1->Branch("runNo",&run_no);
  t1->Branch("subRunNo",&subrun_no);
  t1->Branch("eventNo",&event_no);

  op_wf = uboone_flash.get_rawWfm();
  op_femch = uboone_flash.get_rawChan();
  op_timestamp = uboone_flash.get_rawTimestamp();

  t1->Fill();

  TH2F *h1 = new TH2F("hraw","hraw",1500,0,1500,32,0,32);
  TH2F *h2 = new TH2F("hdecon","hdecon",250,0,250,32,0,32);
  TH2F *h3 = new TH2F("hl1","hl1",250,0,250,32,0,32);
  h1->SetDirectory(file1);
  h2->SetDirectory(file1);
  h3->SetDirectory(file1);
  for (int i=0;i!=32;i++){
    TH1F *h10 = uboone_flash.get_raw_hist(i);
    TH1F *h20 = uboone_flash.get_decon_hist(i);
    TH1F *h30 = uboone_flash.get_l1_hist(i);
    for (int j=0;j!=1500;j++){
      h1->SetBinContent(j+1,i+1,h10->GetBinContent(j+1));
    }
    for (int j=0;j!=250;j++){
      h2->SetBinContent(j+1,i+1,h20->GetBinContent(j+1));
      h3->SetBinContent(j+1,i+1,h30->GetBinContent(j+1));
    }
  }

  TH1F *h_totPE = (TH1F*)uboone_flash.get_totPE()->Clone("totPE");
  TH1F *h_mult = (TH1F*)uboone_flash.get_mult()->Clone("mult");
  TH1F *h_l1_mult = (TH1F*)uboone_flash.get_l1_mult()->Clone("l1_mult");
  TH1F *h_l1_totPE = (TH1F*)uboone_flash.get_l1_totPE()->Clone("l1_totPE");

  h_totPE->SetDirectory(file1);
  h_mult->SetDirectory(file1);
  h_l1_mult->SetDirectory(file1);
  h_l1_totPE->SetDirectory(file1);

  TTree *T_flash = new TTree("T_flash","T_flash");
  T_flash->SetDirectory(file1);
  int type;
  flash_id;
  double low_time, high_time;
  double time;
  double total_PE;
  double PE[32],PE_err[32];
  std::vector<int> fired_channels;
  std::vector<double> l1_fired_time;
  std::vector<double> l1_fired_pe;

  T_flash->Branch("type",&type);
  T_flash->Branch("flash_id",&flash_id);
  T_flash->Branch("low_time",&low_time);
  T_flash->Branch("high_time",&high_time);
  T_flash->Branch("time",&time);
  T_flash->Branch("total_PE",&total_PE);
  T_flash->Branch("PE",PE,"PE[32]/D");
  T_flash->Branch("PE_err",PE_err,"PE_err[32]/D");
  T_flash->Branch("fired_channels",&fired_channels);
  T_flash->Branch("l1_fired_time",&l1_fired_time);
  T_flash->Branch("l1_fired_pe",&l1_fired_pe);

  
  for (auto it = flashes.begin(); it!=flashes.end(); it++){
    fired_channels.clear();
    flash_id = it - flashes.begin();
    Opflash *flash = (*it);
    type = flash->get_type();
    low_time = flash->get_low_time();
    high_time = flash->get_high_time();
    time = flash->get_time();
    total_PE = flash->get_total_PE();
    for (int i=0;i!=32;i++){
      PE[i] = flash->get_PE(i);
      PE_err[i] = flash->get_PE_err(i);
      if (flash->get_fired(i))
	fired_channels.push_back(i);
    }
    l1_fired_time = flash->get_l1_fired_time();
    l1_fired_pe = flash->get_l1_fired_pe();
    T_flash->Fill();
    
  }






   
   file1->Write();
   file1->Close();
}
