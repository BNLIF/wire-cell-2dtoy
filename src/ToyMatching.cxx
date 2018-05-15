#include "WireCell2dToy/ToyMatching.h"

#include "WireCellRess/LassoModel.h"
#include "WireCellRess/ElasticNetModel.h"
#include <Eigen/Dense>

#include "TChain.h"

#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"
#include "WireCellData/FlashTPCBundle.h"


using namespace Eigen;
using namespace WireCell;

int WireCell2dToy::convert_xyz_voxel_id(WireCell::Point &p){
  // int voxel_x_id = std::round((p.x/units::cm+64.825-5.14667/2.)/5.14667);
  // int voxel_y_id = std::round((p.y/units::cm+193-5.14667/2.)/5.14667);
  // int voxel_z_id = std::round((p.z/units::cm+128.243-3.23122/2.)/3.23122);
  
  int voxel_x_id = std::round((p.x/units::cm+63.435-5.1096/2.)/5.1095);
  int voxel_y_id = std::round((p.y/units::cm+191.61-5.1096/2.)/5.1096);
  int voxel_z_id = std::round((p.z/units::cm+92.375-3.05437/2.)/3.05437);
  //fit
  // int voxel_y_id = std::round((p.y/units::cm+188.416-5.11698/2.)/5.11698);
  // int voxel_z_id = std::round((p.z/units::cm+88.3554-3.05768/2.)/3.05768);
  
  if (voxel_x_id<0) voxel_x_id=0;
  if (voxel_x_id>=75) voxel_x_id=74;
  if (voxel_y_id<0) voxel_y_id=0;
  if (voxel_y_id>=75) voxel_y_id=74;
  if (voxel_z_id<0) voxel_z_id=0;
  if (voxel_z_id>=400) voxel_z_id=399;

  int voxel_id = voxel_z_id*75*75 + voxel_y_id*75 + voxel_x_id;
  return voxel_id;
}


//std::vector<std::tuple<PR3DCluster*, Opflash*, double, std::vector<double>>> WireCell2dToy::tpc_light_match(int time_offset, int nrebin, std::map<WireCell::PR3DCluster*,std::vector<std::pair<WireCell::PR3DCluster*,double>>>& group_clusters, WireCell::OpflashSelection& flashes){
FlashTPCBundleSelection WireCell2dToy::tpc_light_match(int time_offset, int nrebin, std::map<WireCell::PR3DCluster*,std::vector<std::pair<WireCell::PR3DCluster*,double>>>& group_clusters, WireCell::OpflashSelection& flashes){
  TChain *T = new TChain("/pmtresponse/PhotonLibraryData","/pmtresponse/PhotonLibraryData");
  T->AddFile("./uboone_photon_library.root");
  //  std::cout << T->GetEntries();
  Int_t Voxel;
  Int_t OpChannel;
  Float_t Visibility;
  T->SetBranchAddress("Voxel",&Voxel);
  T->SetBranchAddress("OpChannel",&OpChannel);
  T->SetBranchAddress("Visibility",&Visibility);

  // PMT data is based on OpChannel
  // Library is NOT ...
  // these are PE numbers ... 
  // Double_t cos_pe_low[32]={11,11,11,11,10,
  // 			   7,8,8,10,7,
  // 			   11,11,11,11,10,
  // 			   9,11,11,7,9,
  // 			   11,10,11,11,11,
  // 			   11,10,11,11,9,
  // 			   10,11};
  // Double_t cos_pe_mid[32]={34,32,28,35,22,
  // 			   23,22,24,33,30,
  // 			   35,35,33,36,33,
  // 			   33,36,33,19,27,
  // 			   32,23,42,32,33,
  // 			   34,24,33,35,25,
  // 			   32,34};

  Double_t cos_pe_low[32]={11,11,11,11,10,
			   7,8,8,10,7,
			   11,11,11,11,10,
			   9,11,11,7,9,
			   11,10,11,11,11,
			   11,11,10,11,11,
			   9,10};
  Double_t cos_pe_mid[32]={34,32,28,35,22,
			   23,22,24,33,30,
			   35,35,33,36,33,
			   33,36,33,19,27,
			   32,23,42,32,33,
			   34,34,24,33,35,
			   25,32};
  
  std::map<int,int> map_lib_pmt,map_pmt_lib;
  map_lib_pmt[1]=2; map_pmt_lib[2]=1; 
  map_lib_pmt[0]=4; map_pmt_lib[4]=0; 
  map_lib_pmt[3]=0; map_pmt_lib[0]=3; 
  map_lib_pmt[2]=5; map_pmt_lib[5]=2; 
  map_lib_pmt[5]=1; map_pmt_lib[1]=5; 
  map_lib_pmt[4]=6; map_pmt_lib[6]=4; 
  map_lib_pmt[6]=3; map_pmt_lib[3]=6; 
  
  map_lib_pmt[9]=7; map_pmt_lib[7]=9; 
  map_lib_pmt[7]=9; map_pmt_lib[9]=7; 
  map_lib_pmt[8]=11; map_pmt_lib[11]=8; 
  map_lib_pmt[11]=8; map_pmt_lib[8]=11;  
  map_lib_pmt[10]=12; map_pmt_lib[12]=10; 
  map_lib_pmt[12]=10; map_pmt_lib[10]=12; 

  map_lib_pmt[14]=13; map_pmt_lib[13]=14;  
  map_lib_pmt[13]=15; map_pmt_lib[15]=13; 
  map_lib_pmt[15]=17; map_pmt_lib[17]=15; 
  map_lib_pmt[17]=14; map_pmt_lib[14]=17; 
  map_lib_pmt[16]=18; map_pmt_lib[18]=16; 
  map_lib_pmt[18]=16; map_pmt_lib[16]=18; 

  map_lib_pmt[21]=19; map_pmt_lib[19]=21; 
  map_lib_pmt[22]=20; map_pmt_lib[20]=22; 
  map_lib_pmt[19]=21; map_pmt_lib[21]=19; 
  map_lib_pmt[20]=23; map_pmt_lib[23]=20; 
  map_lib_pmt[23]=24; map_pmt_lib[24]=23; 
  map_lib_pmt[24]=22; map_pmt_lib[22]=24; 

  map_lib_pmt[26]=25; map_pmt_lib[25]=26; 
  map_lib_pmt[27]=30; map_pmt_lib[30]=27; 
  map_lib_pmt[28]=31; map_pmt_lib[31]=28; 
  map_lib_pmt[31]=29; map_pmt_lib[29]=31;
  // original map
  map_lib_pmt[25]=28; map_pmt_lib[28]=25; 
  map_lib_pmt[30]=27; map_pmt_lib[27]=30; 
  map_lib_pmt[29]=26; map_pmt_lib[26]=29;

  // fixed map ... (if not swap in the flash reconstruction ...)
  // map_lib_pmt[25]=27; map_pmt_lib[27]=25; 
  // map_lib_pmt[30]=26; map_pmt_lib[26]=30; 
  // map_lib_pmt[29]=28; map_pmt_lib[28]=29;
  
  std::vector<std::list<std::pair<int,float>>> photon_library;
  photon_library.resize(400*75*75);
  for (int i=0;i!=T->GetEntries();i++){
    T->GetEntry(i);
    photon_library.at(Voxel).push_back(std::make_pair(OpChannel,Visibility));
  }

  TPCParams& mp = Singleton<TPCParams>::Instance();
  double time_slice_width = mp.get_ts_width();
  
  //Figure out how to use library ... 
  //flashes measurement = 32 * num_flashes 
  const int num_flashes = flashes.size();
  const int num_tpc_objs = group_clusters.size();

  //std::cout << num_flashes << " " << num_tpc_objs << std::endl;
  double high_x_cut = 256 * units::cm;
  double high_x_cut_ext1 = + 1.2*units::cm;
  double high_x_cut_ext2 = - 2.0*units::cm;
  double low_x_cut = 0*units::cm;
  double low_x_cut_ext1 = - 2*units::cm;
  double low_x_cut_ext2 = + 2.0*units::cm;
  double scaling_light_mag = 0.01 * 1.5;

  int solv_type = 1; // new matching code ... 
  
  if (solv_type==1){
    FlashTPCBundleSet all_bundles;
    Flash_bundles_map flash_bundles_map;
    Cluster_bundles_map cluster_bundles_map;
    std::map<std::pair<Opflash*,PR3DCluster*>,FlashTPCBundle*> fc_bundles_map;
      
    int flash_index_id = 0;
    for (auto it1 = flashes.begin(); it1!=flashes.end(); it1++){
      Opflash *flash = (*it1);
      double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
      int cluster_index_id = 0;
      for (auto it2 = group_clusters.begin(); it2!=group_clusters.end(); it2++){
  	PR3DCluster* main_cluster = it2->first;
  	std::vector<std::pair<WireCell::PR3DCluster*,double>>& additional_clusters = it2->second;
  	FlashTPCBundle *bundle = new FlashTPCBundle(flash, main_cluster, flash_index_id, cluster_index_id);
  	bool flag_good_bundle = false;
	
  	std::vector<double>& pred_pmt_light = bundle->get_pred_pmt_light();
  	PR3DClusterSelection& other_clusters = bundle->get_other_clusters();
  	PR3DClusterSelection& more_clusters = bundle->get_more_clusters();

  	double first_pos_x = (*((main_cluster->get_time_cells_set_map().begin())->second.begin()))->get_sampling_points().front().x;
  	double last_pos_x = (*((main_cluster->get_time_cells_set_map().rbegin())->second.begin()))->get_sampling_points().front().x;

	bool flag_spec_end = false;
	
	// improve the position code ... 
	if (first_pos_x - offset_x <= low_x_cut + low_x_cut_ext1 &&
	    first_pos_x - offset_x > low_x_cut - 120*units::cm ){
	  
	  std::map<int,SMGCSet>& time_cells_set_map = main_cluster->get_time_cells_set_map();
	  int num_mcells_outside = 0;
	  int num_time_slices_outside = 0;

	  int num_mcells_def_outside = 0;
	  
	  double prev_pos_x= first_pos_x;
	  double current_pos_x = first_pos_x;
	  for (auto it3 = time_cells_set_map.begin(); it3 != time_cells_set_map.end(); it3++){
	    current_pos_x = (*(it3->second.begin()))->get_sampling_points().front().x;


	      // if (flash->get_flash_id()==35&&abs(main_cluster->get_cluster_id()-5)<=0)
	      // std::cout << num_time_slices_outside<< " " << num_mcells_outside << "  " << (first_pos_x - offset_x)/units::cm << " " <<
	      // 	(prev_pos_x-offset_x)/units::cm << " " << (current_pos_x-offset_x)/units::cm << std::endl;
	    
	    if (current_pos_x -offset_x > low_x_cut + low_x_cut_ext1 && current_pos_x - prev_pos_x > 0.75*units::cm)
	      break;
	    if (num_time_slices_outside > 60) break;
	    
	    if (current_pos_x -offset_x < low_x_cut + low_x_cut_ext1) num_mcells_def_outside += it3->second.size();

	    num_time_slices_outside += 1;
	    num_mcells_outside += it3->second.size();
	    prev_pos_x = current_pos_x;
	  }
	  if (num_time_slices_outside <=36 && num_mcells_outside < 0.05*main_cluster->get_num_mcells()){
	    first_pos_x = current_pos_x;
	    if (num_time_slices_outside > 10 && fabs(current_pos_x - prev_pos_x)<10*units::cm)
	      flag_spec_end = true;
	  }else if (num_time_slices_outside <=60 && num_mcells_outside < 0.06*main_cluster->get_num_mcells() && fabs(current_pos_x - prev_pos_x)>10*units::cm){
	    first_pos_x = current_pos_x;
	  }else if (num_time_slices_outside <=25 && num_mcells_outside < 0.12 * main_cluster->get_num_mcells() && fabs(current_pos_x - prev_pos_x)>20*units::cm){
	    first_pos_x = current_pos_x;
	  }

	  if (num_mcells_def_outside < 0.0015 * main_cluster->get_num_mcells()&&num_mcells_def_outside>0)
	    first_pos_x = offset_x;
	  
	  //  if (flash->get_flash_id()==61&&main_cluster->get_cluster_id()==1)
	  // std::cout << num_mcells_outside << " " << main_cluster->get_num_mcells() << "  A " << (first_pos_x - offset_x)/units::cm << " " <<
	  //  (prev_pos_x-offset_x)/units::cm << " " << (current_pos_x-offset_x)/units::cm << " " << num_mcells_def_outside << std::endl;
	  
	}
	if (last_pos_x - offset_x >= high_x_cut + high_x_cut_ext1 &&
	    last_pos_x - offset_x < high_x_cut + 120*units::cm){
	  std::map<int,SMGCSet>& time_cells_set_map = main_cluster->get_time_cells_set_map();
	  int num_mcells_outside = 0;
	  int num_time_slices_outside = 0;
	  int num_mcells_def_outside = 0;
	  double prev_pos_x= last_pos_x;
	  double current_pos_x = last_pos_x;

	  for (auto it3 = time_cells_set_map.rbegin(); it3 != time_cells_set_map.rend(); it3++){
	    current_pos_x = (*(it3->second.begin()))->get_sampling_points().front().x;
	    if (current_pos_x -offset_x<high_x_cut + high_x_cut_ext1 && fabs(current_pos_x - prev_pos_x) > 0.75*units::cm)
	      break;
	    if (num_time_slices_outside > 60) break;


	    if (current_pos_x -offset_x>high_x_cut + high_x_cut_ext1) num_mcells_def_outside +=it3->second.size();
	    
	    num_time_slices_outside += 1;
	    num_mcells_outside += it3->second.size();
	    prev_pos_x = current_pos_x;
	  }
	  if (num_time_slices_outside <=36 && num_mcells_outside < 0.05*main_cluster->get_num_mcells()){
	    last_pos_x = current_pos_x;
	    if (num_time_slices_outside > 10 && fabs(current_pos_x - prev_pos_x)<10*units::cm)
	      flag_spec_end = true;
	  }else if (num_time_slices_outside <=60 && num_mcells_outside < 0.06*main_cluster->get_num_mcells() && fabs(current_pos_x - prev_pos_x)>10*units::cm){
	    last_pos_x = current_pos_x;
	  }else if (num_time_slices_outside <=25 && num_mcells_outside < 0.12 * main_cluster->get_num_mcells() && fabs(current_pos_x - prev_pos_x)>20*units::cm){
	    last_pos_x = current_pos_x;
	  }
	  
	  if (num_mcells_def_outside < 0.0015 * main_cluster->get_num_mcells()&&num_mcells_def_outside>0)
	    last_pos_x = offset_x+high_x_cut;

	  // if (flash->get_flash_id()==19&&main_cluster->get_cluster_id()==19)
	  //   std::cout << flash->get_flash_id() << " "<< main_cluster->get_cluster_id() << " " << (first_pos_x-offset_x)/units::cm << " " << (last_pos_x-offset_x)/units::cm << " " << num_time_slices_outside << " " << num_mcells_outside << " " << main_cluster->get_num_mcells() << " " << fabs(current_pos_x - prev_pos_x)/units::cm << std::endl;
	  
	}
	
	// if (flash->get_flash_id()==19 && main_cluster->get_cluster_id()==19 )
	//  std::cout << flash->get_flash_id() << " "<< main_cluster->get_cluster_id() << " " << (first_pos_x-offset_x)/units::cm << " " << (last_pos_x-offset_x)/units::cm << std::endl;

	//	if (flash->get_flash_id()==14)
	//      std::cout << flash->get_flash_id() << " " << main_cluster->get_cluster_id() << " " << offset_x/units::cm << " " << (first_pos_x-offset_x)/units::cm << " " << (last_pos_x-offset_x)/units::cm << " " << std::endl;
	
  	if (first_pos_x-offset_x > low_x_cut + low_x_cut_ext1 -1.0*units::cm &&
  	    last_pos_x-offset_x > low_x_cut &&
  	    last_pos_x-offset_x < high_x_cut + high_x_cut_ext1 &&
  	    first_pos_x-offset_x < high_x_cut){
	  
	  bundle->set_spec_end_flag(flag_spec_end);
	  
  	  if (first_pos_x-offset_x <=low_x_cut + low_x_cut_ext2 && first_pos_x-offset_x > low_x_cut + low_x_cut_ext1 ){
  	    bundle->set_flag_close_to_PMT(true);
  	    bundle->set_flag_at_x_boundary(true);
  	  }
  	  if (last_pos_x-offset_x >= high_x_cut + high_x_cut_ext2 && last_pos_x-offset_x < high_x_cut + high_x_cut_ext1){
  	    bundle->set_flag_at_x_boundary(true);
  	  }

  	  PR3DClusterSelection temp_clusters;
  	  temp_clusters.push_back(main_cluster);
  	  for (auto it3 = additional_clusters.begin(); it3!=additional_clusters.end(); it3++){
  	    temp_clusters.push_back(it3->first);
  	    other_clusters.push_back(it3->first);
  	  }

  	  for (auto it3 = temp_clusters.begin(); it3!=temp_clusters.end(); it3++){
  	    SMGCSelection& mcells = (*it3)->get_mcells();
  	    bool flag_save = true;

	    if ((*it3) == main_cluster)  flag_save = false;
	    
  	    for (auto it4 = mcells.begin(); it4!=mcells.end(); it4++){
  	      SlimMergeGeomCell *mcell = (*it4);
  	      if (mcell->get_q()>0){
  		PointVector& pts = mcell->get_sampling_points();
  		if (pts.at(0).x-offset_x < low_x_cut+low_x_cut_ext1 ||
  		    pts.at(0).x-offset_x > high_x_cut+high_x_cut_ext1){
  		  flag_save = false;
  		  continue;
  		}
		
  		float charge = mcell->get_q()/pts.size();
  		Point p;
  		for (size_t i=0;i!=pts.size();i++){
  		  p.x = pts.at(i).x - offset_x;
  		  p.y = pts.at(i).y;
  		  p.z = pts.at(i).z;
		  
		  int voxel_id = WireCell2dToy::convert_xyz_voxel_id(p);
		  std::list<std::pair<int,float>>& pmt_list = photon_library.at(voxel_id);
		  
		  for (auto it5 = pmt_list.begin(); it5!=pmt_list.end(); it5++){
		    pred_pmt_light.at(map_lib_pmt[it5->first]) += charge * it5->second;
		  }
		}
  	      }
  	    }
  	    if (flag_save)
  	      more_clusters.push_back(*it3);
	    
  	  } // loop over all clusters within this TPC object...

  	  double sum1 = 0, sum2 = 0, max_pe = 0;
  	  for (size_t i=0;i!=32;i++){
  	    pred_pmt_light.at(i) *= scaling_light_mag;

  	    sum1 += flash->get_PE(i);
  	    sum2 += pred_pmt_light.at(i);
  	    if (pred_pmt_light.at(i) > max_pe)
  	      max_pe = pred_pmt_light.at(i);
  	  }
	  
	  // if (sum2 < sum1 * 3){ // three times allowrance ... 
	  flag_good_bundle = true;
	  // }
  	}

  	if (flag_good_bundle){
  	  all_bundles.insert(bundle);
  	  if (flash_bundles_map.find(flash)==flash_bundles_map.end()){
  	    FlashTPCBundleSelection bundles;
  	    bundles.push_back(bundle);
  	    flash_bundles_map[flash] = bundles;
  	  }else{
  	    flash_bundles_map[flash].push_back(bundle);
  	  }
  	  if (cluster_bundles_map.find(main_cluster)==cluster_bundles_map.end()){
  	    FlashTPCBundleSelection bundles;
  	    bundles.push_back(bundle);
  	    cluster_bundles_map[main_cluster] = bundles;
  	  }else{
  	    cluster_bundles_map[main_cluster].push_back(bundle);
  	  }
  	  fc_bundles_map[std::make_pair(flash,main_cluster)] = bundle;

	  // if (flash->get_flash_id()==14)
	  //   std::cout << flash->get_flash_id() << " " << main_cluster->get_cluster_id() << " " << bundle << std::endl;
	  
  	}else{
  	  delete bundle;
  	}
	cluster_index_id++;
      }
      flash_index_id++;
    }

    
    // examine the bundles ... 
    //std::cout << "Starting: " << cluster_bundles_map.size() << " A " << flash_bundles_map.size() << " " << all_bundles.size() << std::endl;
    {
      FlashTPCBundleSelection to_be_removed;
      for (auto it = all_bundles.begin(); it!=all_bundles.end();it++){
	FlashTPCBundle *bundle = *it;
	if (!bundle->examine_bundle(cos_pe_low, cos_pe_mid)){
	  to_be_removed.push_back(bundle);
	}
      }
      
      for (auto it1 = to_be_removed.begin(); it1!=to_be_removed.end(); it1++){
	FlashTPCBundle *bundle = *it1;
	all_bundles.erase(bundle);
	
	Opflash *flash = bundle->get_flash();
	PR3DCluster *cluster = bundle->get_main_cluster();
	
	fc_bundles_map.erase(std::make_pair(flash,cluster));
	
	{
	  FlashTPCBundleSelection& temp_bundles = flash_bundles_map[flash];
	  temp_bundles.erase(find(temp_bundles.begin(), temp_bundles.end(), bundle));
	  if (temp_bundles.size()==0)
	    flash_bundles_map.erase(flash);
	}
	{
	  FlashTPCBundleSelection& temp_bundles = cluster_bundles_map[cluster];
	  temp_bundles.erase(find(temp_bundles.begin(), temp_bundles.end(), bundle));
	  if (temp_bundles.size()==0)
	    cluster_bundles_map.erase(cluster);
	}
	
	delete bundle;
      }
      to_be_removed.clear();


      for (auto it = cluster_bundles_map.begin(); it!= cluster_bundles_map.end(); it++){
       	PR3DCluster *main_cluster = it->first;
       	FlashTPCBundleSelection& bundles = it->second;
       	bool flag_consist = false;
	FlashTPCBundleSelection temp_removed;
	
       	for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	  FlashTPCBundle *bundle = *it1;
	  if (bundle->get_potential_bad_match_flag())
	    temp_removed.push_back(bundle);
	  if (bundle->get_consistent_flag() ||
	       (bundle->get_ks_dis() < 0.12 || bundle->get_chi2() < 3 * bundle->get_ndf() )&& bundle->get_ndf()>=3 ||
	       bundle->get_ks_dis()<0.33 && bundle->get_chi2() < 50 * bundle->get_ndf()&&bundle->get_ndf()>=5 && bundle->get_flag_close_to_PMT())
	    flag_consist = true;

	  // if (main_cluster->get_cluster_id()==30)
	  //   std::cout << bundle->get_flash()->get_flash_id() << " " << bundle->get_main_cluster()->get_cluster_id() << " " << bundle->get_flag_at_x_boundary() << " " << bundle->get_ks_dis() << " " << bundle->get_chi2() << " " << bundle->get_ndf() << " " << bundle->get_consistent_flag() << " " << bundle->get_flag_close_to_PMT() << " " << bundle->get_potential_bad_match_flag() << " " << flag_consist << std::endl;
	  
	}
	if (flag_consist)
	  to_be_removed.insert(to_be_removed.end(),temp_removed.begin(),temp_removed.end());
      }

      
      for (auto it1 = to_be_removed.begin(); it1!=to_be_removed.end(); it1++){
	FlashTPCBundle *bundle = *it1;
	all_bundles.erase(bundle);
	
	Opflash *flash = bundle->get_flash();
	PR3DCluster *cluster = bundle->get_main_cluster();
	
	fc_bundles_map.erase(std::make_pair(flash,cluster));
	
	{
	  FlashTPCBundleSelection& temp_bundles = flash_bundles_map[flash];
	  temp_bundles.erase(find(temp_bundles.begin(), temp_bundles.end(), bundle));
	  if (temp_bundles.size()==0)
	    flash_bundles_map.erase(flash);
	}
	{
	  FlashTPCBundleSelection& temp_bundles = cluster_bundles_map[cluster];
	  temp_bundles.erase(find(temp_bundles.begin(), temp_bundles.end(), bundle));
	  if (temp_bundles.size()==0)
	    cluster_bundles_map.erase(cluster);
	}
	
	delete bundle;
      }
      
      
    }

    

    
    //    std::cout << "After Cleaning 1 : " << cluster_bundles_map.size() << " A " << flash_bundles_map.size() << " " << all_bundles.size() << std::endl;
    
    for (auto it = cluster_bundles_map.begin(); it!= cluster_bundles_map.end(); it++){
      PR3DCluster *main_cluster = it->first;
      FlashTPCBundleSelection& bundles = it->second;
      bool flag_tight_bundle = false;

      bool flag_highly_consistent_bundle = false;

      
      for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	FlashTPCBundle *bundle = *it1;

	//	std::cout << bundle->get_flash()->get_flash_id() << " " << bundle->get_main_cluster()->get_cluster_id() << " " << bundle->get_flag_at_x_boundary() << " " << bundle->get_ks_dis() << " " << bundle->get_chi2() << " " << bundle->get_ndf() << " " << bundle->get_consistent_flag() << " " << bundle->get_flag_close_to_PMT() << " " << bundle->get_potential_bad_match_flag() << std::endl;
	
	if (bundle->get_consistent_flag()){
	  flag_tight_bundle = true;
	  
	  if (bundle->get_ks_dis()<0.05 && bundle->get_ndf() >= 10 && bundle->get_chi2() < bundle->get_ndf()  * 12  ||
	      bundle->get_ks_dis()<0.07 && bundle->get_ndf() >= 10 && bundle->get_chi2() < bundle->get_ndf()  * 9  ||
	      bundle->get_ks_dis()<0.1 && bundle->get_ndf() >=5 && bundle->get_chi2() < bundle->get_ndf()  * 4  ||
	      bundle->get_ks_dis()<0.15 && bundle->get_ndf()>=5 && bundle->get_chi2() < bundle->get_ndf()  * 3 	      )
	    flag_highly_consistent_bundle = true;
	  //  break;
	}
      }
      
      if (!flag_tight_bundle){
	for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	  FlashTPCBundle *bundle = *it1;
	  if (bundle->get_ks_dis()<0.07 && bundle->get_ndf()>=10 && bundle->get_chi2() < bundle->get_ndf() * 60){
	    bundle->set_consistent_flag(true);
	    flag_tight_bundle = true;
	  }else if (bundle->get_ks_dis()<0.33 && bundle->get_ndf()>=3 && bundle->get_chi2() < bundle->get_ndf() * 10){
	    bundle->set_consistent_flag(true);
	    flag_tight_bundle = true;
	  }else if (bundle->get_ks_dis()<0.33 && bundle->get_ndf()>=3 && (bundle->get_chi2() < bundle->get_ndf() * 50 && bundle->get_flag_close_to_PMT() || bundle->get_chi2() < bundle->get_ndf() * 16 && bundle->get_flag_at_x_boundary())){
	    bundle->set_consistent_flag(true);
	    flag_tight_bundle = true;
	  }else if (bundle->get_ks_dis()<0.22 && bundle->get_ndf()>=3 && bundle->get_chi2() < bundle->get_ndf() * 16){
	    bundle->set_consistent_flag(true);
	    flag_tight_bundle = true;
	  }else if (bundle->get_ks_dis()<0.16 && bundle->get_ndf()>=6 && bundle->get_chi2() < bundle->get_ndf() * 20){
	    bundle->set_consistent_flag(true);
	    flag_tight_bundle = true;
	  }
	}

	
	FlashTPCBundle *min_bundle = *bundles.begin();
	
	for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	  FlashTPCBundle *bundle = *it1;
	  if (bundle->get_ks_dis()<min_bundle->get_ks_dis()){
	    min_bundle = bundle;
	  }
	}
	
	FlashTPCBundle *min_bundle1 = 0;
	for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	  FlashTPCBundle *bundle = *it1;
	  if (bundle==min_bundle) continue;
	  if (min_bundle1==0) {
	    min_bundle1 = bundle;
	  }else if (bundle->get_ks_dis()<min_bundle1->get_ks_dis()){
	    min_bundle1 = bundle;
	  }
	}
	
	
	// if (main_cluster->get_cluster_id()==12)
	//   std::cout << "Xin: " << min_bundle->get_flash()->get_flash_id() << " " << min_bundle1->get_flash()->get_flash_id() << std::endl;
	

	bool flag_set = false;
	if (min_bundle->get_ks_dis()<0.15 && min_bundle->get_ndf()>=6 && min_bundle->get_chi2() < min_bundle->get_ndf() * 40){
	  min_bundle->set_consistent_flag(true);
	  flag_tight_bundle = true;
	  flag_set = true;
	}else if ( min_bundle->get_ks_dis()<0.11 && min_bundle->get_ndf() >= 3 && min_bundle->get_chi2() < min_bundle->get_ndf() * 36){
	  min_bundle->set_consistent_flag(true);
	  flag_tight_bundle = true;
	  flag_set = true;
	}else if (min_bundle->get_ks_dis()<0.075 && min_bundle->get_ndf()>=10 && min_bundle->get_chi2() < min_bundle->get_ndf() * 60){
	  min_bundle->set_consistent_flag(true);
	  flag_tight_bundle = true;
	  flag_set = true;
	}else if (min_bundle->get_ks_dis()<0.11 && min_bundle->get_ndf()>=10 && min_bundle->get_chi2() < min_bundle->get_ndf() * 120){
	  min_bundle->set_consistent_flag(true);
	  flag_tight_bundle = true;
	  flag_set = true;
	}else if (min_bundle->get_ks_dis()<0.17 && min_bundle->get_ndf()>=10 && min_bundle->get_chi2() < min_bundle->get_ndf() * 40){
	  min_bundle->set_consistent_flag(true);
	  flag_tight_bundle = true;
	  flag_set = true;
	}else if (min_bundle->get_ks_dis()<0.22 && min_bundle->get_ndf()>=20 && min_bundle->get_chi2() < min_bundle->get_ndf() * 40){
	  min_bundle->set_consistent_flag(true);
	  flag_tight_bundle = true;
	  flag_set = true;
	}else if (min_bundle->get_flag_at_x_boundary() && min_bundle->get_ks_dis()<0.3 && min_bundle->get_ndf()>=10 && min_bundle->get_chi2() < min_bundle->get_ndf() * 30){
	  min_bundle->set_consistent_flag(true);
	  flag_tight_bundle = true;
	  flag_set = true;
	}
	
	if (!flag_set && min_bundle1!=0){
	  if (min_bundle1->get_ks_dis()<0.15 && min_bundle1->get_ndf()>=6 && min_bundle1->get_chi2() < min_bundle1->get_ndf() * 40){
	    min_bundle1->set_consistent_flag(true);
	    flag_tight_bundle = true;
	  }else if ( min_bundle1->get_ks_dis()<0.11 && min_bundle1->get_ndf() >= 3 && min_bundle1->get_chi2() < min_bundle1->get_ndf() * 36){
	    min_bundle1->set_consistent_flag(true);
	    flag_tight_bundle = true;
	  }else if (min_bundle1->get_ks_dis()<0.075 && min_bundle1->get_ndf()>=10 && min_bundle1->get_chi2() < min_bundle1->get_ndf() * 60){
	    min_bundle1->set_consistent_flag(true);
	    flag_tight_bundle = true;
	  }else if (min_bundle1->get_ks_dis()<0.11 && min_bundle1->get_ndf()>=10 && min_bundle1->get_chi2() < min_bundle1->get_ndf() * 120){
	    min_bundle1->set_consistent_flag(true);
	    flag_tight_bundle = true;
	  }else if (min_bundle1->get_ks_dis()<0.17 && min_bundle1->get_ndf()>=10 && min_bundle1->get_chi2() < min_bundle1->get_ndf() * 40){
	    min_bundle1->set_consistent_flag(true);
	    flag_tight_bundle = true;
	  }else if (min_bundle1->get_ks_dis()<0.22 && min_bundle1->get_ndf()>=20 && min_bundle1->get_chi2() < min_bundle1->get_ndf() * 40){
	    min_bundle1->set_consistent_flag(true);
	    flag_tight_bundle = true;
	  }else if (min_bundle1->get_flag_at_x_boundary() && min_bundle1->get_ks_dis()<0.3 && min_bundle1->get_ndf()>=10 && min_bundle1->get_chi2() < min_bundle1->get_ndf() * 30){
	    min_bundle1->set_consistent_flag(true);
	    flag_tight_bundle = true;
	  }
	  
	  //else if (min_bundle1->get_ks_dis()<0.17 && min_bundle1->get_ndf()>=10 && min_bundle1->get_chi2() < min_bundle1->get_ndf() * 40){
	  // min_bundle1->set_consistent_flag(true);
	  // flag_tight_bundle = true;
	  //}
	}
      }

      
       if (!flag_tight_bundle){
	// last pieces ... 
	for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	  FlashTPCBundle *bundle = *it1;
	  if (bundle->get_flag_close_to_PMT() && bundle->get_ks_dis()<0.7 && bundle->get_ndf()>5 && bundle->get_chi2() < bundle->get_ndf() * 40){
	    bundle->set_consistent_flag(true);
	    flag_tight_bundle = true;
	  }
	}
      }
     


      if (flag_tight_bundle){
	if (flag_highly_consistent_bundle){
	  for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	    FlashTPCBundle *bundle = *it1;
	    if (bundle->get_consistent_flag()){
	      if (bundle->get_ks_dis()<0.15 && bundle->get_ndf() >=5 && bundle->get_chi2() < bundle->get_ndf()  * 9 ){
	      }else if (bundle->get_ks_dis()<0.075 && bundle->get_ndf() >=8 && bundle->get_chi2() < bundle->get_ndf()  * 12){
	      }else if (bundle->get_ks_dis()<0.07 && bundle->get_ndf() >= 16 && bundle->get_chi2() < bundle->get_ndf()  * 20){
	      }else{
		bundle->set_consistent_flag(false);
	      }
	      
	      
	    }
	  }

	 
	  
	  
	  
	}

	FlashTPCBundleSelection ndf1_bundles;
	FlashTPCBundleSelection other_bundles;
	for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	  FlashTPCBundle *bundle = *it1;
	  if (bundle->get_consistent_flag()){
	    if (bundle->get_ndf()==1){
	      ndf1_bundles.push_back(bundle);
	      	}else{
	      other_bundles.push_back(bundle);
	    }
	  }
	}

	// now remove ... 
	bool flag_temp = false;
	// if (main_cluster->get_cluster_id()==15)
	//   std::cout << "Xin: " << other_bundles.size() << " " << ndf1_bundles.size() << std::endl;
	
	for (auto it1 = other_bundles.begin(); it1!=other_bundles.end(); it1++){
	  FlashTPCBundle *bundle = *it1;
	  if (flash_bundles_map.find(bundle->get_flash())!=flash_bundles_map.end()){
	    FlashTPCBundleSelection bundles1 = flash_bundles_map[bundle->get_flash()];
	    int temp_num = 0;
	    for (auto it2 = bundles1.begin(); it2!=bundles1.end(); it2++){
	      if ((*it2)->get_consistent_flag())
		temp_num++;
	    }
	    // if (bundle->get_main_cluster()->get_cluster_id()==15)
	    //   std::cout << bundle->get_flash()->get_flash_id() << " " << temp_num << std::endl;
	    
	    if (temp_num==1){
	      flag_temp = true;
	      break;
	    }
	    
	  }
	}
	
	if (flag_temp){
	  for (auto it1 = ndf1_bundles.begin(); it1!=ndf1_bundles.end(); it1++){
	    FlashTPCBundle *bundle = *it1;
	    bundle->set_consistent_flag(false);
	  }
	}
	
	
      }
    }

    //
    for (auto it = cluster_bundles_map.begin(); it!= cluster_bundles_map.end(); it++){
      PR3DCluster *main_cluster = it->first;
      FlashTPCBundleSelection& bundles = it->second;

      FlashTPCBundleSelection temp_bundles;
      bool flag_remove = false;

      for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	FlashTPCBundle *bundle = *it1;
	if (bundle->get_consistent_flag() && bundle->get_spec_end_flag()){
	  temp_bundles.push_back(bundle);
	}else if (bundle->get_consistent_flag()){
	  flag_remove = true;
	}
      }
      
      if (flag_remove){
	for (auto it1 = temp_bundles.begin(); it1!=temp_bundles.end(); it1++){
	  FlashTPCBundle *bundle = *it1;
	  bundle->set_consistent_flag(false);
	}
      }

      // examine the again
      flag_remove = false;
      temp_bundles.clear();
      for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	FlashTPCBundle *bundle = *it1;
	if (bundle->get_consistent_flag()){
	  if (bundle->get_ks_dis() < 0.06 && bundle->get_chi2() < 3.*bundle->get_ndf() && bundle->get_flag_at_x_boundary()){
	    flag_remove = true;
	  }else{
	    temp_bundles.push_back(bundle);
	  }
	}
      }
      if (flag_remove){
	for (auto it1 = temp_bundles.begin(); it1!=temp_bundles.end(); it1++){
	  FlashTPCBundle *bundle = *it1;
	  bundle->set_consistent_flag(false);
	}
      }
      
    }
    
    

    // further examine the map ...
    {
      Flash_bundles_map flash_good_bundles_map;
      Cluster_bundles_map cluster_good_bundles_map;
      FlashTPCBundleSelection good_bundles;
            
      Flash_bundles_map flash_other_bundles_map;
      Cluster_bundles_map cluster_other_bundles_map;
      
      for (auto it = all_bundles.begin(); it!=all_bundles.end();it++){
	FlashTPCBundle *bundle = *it;
	Opflash *flash = bundle->get_flash();
	PR3DCluster *cluster = bundle->get_main_cluster();

	if (bundle->get_consistent_flag()){
	  if (flash_good_bundles_map.find(flash)==flash_good_bundles_map.end()){
	    FlashTPCBundleSelection bundles;
	    bundles.push_back(bundle);
	    flash_good_bundles_map[flash] = bundles;
	  }else{
	    flash_good_bundles_map[flash].push_back(bundle);
	  }
	  if (cluster_good_bundles_map.find(cluster)==cluster_good_bundles_map.end()){
	    FlashTPCBundleSelection bundles;
	    bundles.push_back(bundle);
	    cluster_good_bundles_map[cluster] = bundles;
	  }else{
	    cluster_good_bundles_map[cluster].push_back(bundle);
	  }
	  good_bundles.push_back(bundle);
	}else if (bundle->get_flag_at_x_boundary()){
	  if (flash_other_bundles_map.find(flash)==flash_other_bundles_map.end()){
	    FlashTPCBundleSelection bundles;
	    bundles.push_back(bundle);
	    flash_other_bundles_map[flash] = bundles;
	  }else{
	    flash_other_bundles_map[flash].push_back(bundle);
	  }
	  if (cluster_other_bundles_map.find(cluster)==cluster_other_bundles_map.end()){
	    FlashTPCBundleSelection bundles;
	    bundles.push_back(bundle);
	    cluster_other_bundles_map[cluster] = bundles;
	    }else{
	    cluster_other_bundles_map[cluster].push_back(bundle);
	  }
	}
	// if ( bundle->get_consistent_flag() || bundle->get_flag_at_x_boundary())
	// 	std::cout << bundle->get_flash()->get_flash_id() << " " << bundle->get_main_cluster()->get_cluster_id() << " " << bundle->get_flag_at_x_boundary() << " " << bundle->get_ks_dis() << " " << bundle->get_chi2() << " " << bundle->get_ndf() << std::endl;
      }

      for (auto it = cluster_good_bundles_map.begin(); it!= cluster_good_bundles_map.end(); it++){
	PR3DCluster *cluster =  it->first;
	FlashTPCBundleSelection& bundles = it->second;

	
	if (bundles.size()>1){ // more than one flash

	  // find the min bundle ... 
	  FlashTPCBundle *min_bundle = bundles.at(0);
	  for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	    FlashTPCBundle *bundle = *it1;
	    if (bundle->get_ks_dis()+0.003 * flash_good_bundles_map[bundle->get_flash()].size() < min_bundle->get_ks_dis() + 0.003 * flash_good_bundles_map[min_bundle->get_flash()].size())
	      min_bundle = bundle;
	  }

	  if (min_bundle->get_ks_dis() < 0.18){
	    //  std::cout << min_bundle->get_flash()->get_flash_id() << " Xin: " << min_bundle->get_main_cluster()->get_cluster_id() << std::endl;

	    FlashTPCBundleSelection temp_bundles;
	    
	    //examine the rest of bundles;
	    for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	      FlashTPCBundle *bundle = *it1;
	      if (bundle==min_bundle) continue;
	      bool flag_remove = false;
	      
	      if (min_bundle->get_ks_dis() +0.015 < bundle->get_ks_dis() &&
		  min_bundle->get_chi2()/min_bundle->get_ndf() *1.1 < bundle->get_chi2()/bundle->get_ndf()){
		// prepare remove bundle from the list ...
		flag_remove = true;
	      }else if(min_bundle->get_ks_dis() < bundle->get_ks_dis() &&
		       min_bundle->get_chi2()/min_bundle->get_ndf() * 3 < bundle->get_chi2()/bundle->get_ndf()){
		flag_remove = true;
	      }else if (min_bundle->get_ks_dis() +0.025 < bundle->get_ks_dis() &&
			min_bundle->get_chi2()/min_bundle->get_ndf() *0.85 < bundle->get_chi2()/bundle->get_ndf()){
		flag_remove = true;
		if (bundle->get_ks_dis() < 0.075&&bundle->get_chi2()/bundle->get_ndf()<1.6)
		  flag_remove = false;
	      }else if (min_bundle->get_ks_dis() +0.03 < bundle->get_ks_dis() &&
			min_bundle->get_chi2()/min_bundle->get_ndf() *0.75 < bundle->get_chi2()/bundle->get_ndf() &&
			min_bundle->get_ks_dis() <0.06){
		flag_remove = true;
	      }else if (min_bundle->get_ks_dis()+0.075 < bundle->get_ks_dis() &&
			min_bundle->get_ks_dis()<0.075 &&
			min_bundle->get_chi2()/min_bundle->get_ndf() /15. < bundle->get_chi2()/bundle->get_ndf()){
		flag_remove = true;
	      }else if (flash_good_bundles_map[min_bundle->get_flash()].size()==1){
		if (min_bundle->get_ks_dis()<0.06){
		  if (min_bundle->get_chi2()/min_bundle->get_ndf() < bundle->get_chi2()/bundle->get_ndf() * 2.2){
		    flag_remove = true;
		  }else if (min_bundle->get_ks_dis() + 0.06 < bundle->get_ks_dis() &&
			    min_bundle->get_chi2()/min_bundle->get_ndf() < bundle->get_chi2()/bundle->get_ndf() * 3.0){
		    flag_remove = true;
		  }else if (min_bundle->get_chi2()/min_bundle->get_ndf()<9 && min_bundle->get_ks_dis() + 0.01 < bundle->get_ks_dis()){
		    flag_remove = true;
		  }else if (min_bundle->get_ks_dis()<0.05 && min_bundle->get_chi2()/min_bundle->get_ndf()<38 && min_bundle->get_ks_dis() + 0.05 < bundle->get_ks_dis()){
		    flag_remove = true;
		  }
		}else{
		  if (min_bundle->get_ks_dis() + 0.015 < bundle->get_ks_dis() &&
		      min_bundle->get_chi2()/min_bundle->get_ndf()*0.75 < bundle->get_chi2()/bundle->get_ndf() &&flash_good_bundles_map[bundle->get_flash()].size()>1
		      ){
		    flag_remove = true;
		  }else if (min_bundle->get_ks_dis() + 0.08 < bundle->get_ks_dis() &&
			    flash_good_bundles_map[bundle->get_flash()].size()>1 &&
			    min_bundle->get_chi2()/min_bundle->get_ndf() < 25
			    ){
		    flag_remove = true;
		  }
		}
	      }
	      
	      // if (min_bundle->get_main_cluster()->get_cluster_id()==22)
	      //  	std::cout << min_bundle->get_flash()->get_flash_id() << " A " << min_bundle->get_main_cluster()->get_cluster_id() << " " << flag_remove << " " << flash_good_bundles_map[min_bundle->get_flash()].size()<< std::endl;
	      
	      if (flag_remove){
		if (flash_good_bundles_map[bundle->get_flash()].size()>1){
		  
		  FlashTPCBundle *min_bundle1 = flash_good_bundles_map[bundle->get_flash()].at(0);
		  FlashTPCBundle *min_bundle2 = flash_good_bundles_map[bundle->get_flash()].at(0);
		  for (auto it2 = flash_good_bundles_map[bundle->get_flash()].begin(); it2!=flash_good_bundles_map[bundle->get_flash()].end(); it2++){
		    FlashTPCBundle *bundle1 = *it2;
		    if (bundle1->get_ks_dis() < min_bundle1->get_ks_dis()){
		      min_bundle1 = bundle1;
		      min_bundle2 = min_bundle1;
		    }
		  }
		  if (bundle!=min_bundle1){
		    temp_bundles.push_back(bundle);
		  }else{
		    if (min_bundle1->get_ks_dis() + 0.01 > min_bundle2->get_ks_dis() &&
			min_bundle1->get_chi2() > min_bundle2->get_chi2()*0.85){
		      temp_bundles.push_back(bundle);
		    }else if (min_bundle1->get_chi2()/min_bundle1->get_ndf() > min_bundle2->get_chi2()/min_bundle2->get_ndf() * 3 &&
			      min_bundle1->get_ndf()>=5 && 
			      min_bundle2->get_ks_dis() < 0.12){
		      temp_bundles.push_back(bundle);
		    }
		  }
		}
	      }
	    }
	  

	    for (auto it1 = temp_bundles.begin(); it1!=temp_bundles.end(); it1++){
	      FlashTPCBundle *bundle = *it1;
	      PR3DCluster *cluster =  bundle->get_main_cluster();
	      Opflash *flash = bundle->get_flash();
	      
	      bundle->set_consistent_flag(false);
	      flash_good_bundles_map[flash].erase(find(flash_good_bundles_map[flash].begin(),flash_good_bundles_map[flash].end(),bundle));
	      cluster_good_bundles_map[cluster].erase(find(cluster_good_bundles_map[cluster].begin(),cluster_good_bundles_map[cluster].end(),bundle));
	      // remove from 
	    }
	  }
	}
      }

      // new round according to chi2 .. .
      for (auto it = cluster_good_bundles_map.begin(); it!= cluster_good_bundles_map.end(); it++){
	PR3DCluster *cluster =  it->first;
	FlashTPCBundleSelection& bundles = it->second;

	
	if (bundles.size()>1){ // more than one flash
	  // find the min bundle ... 
	  FlashTPCBundle *min_bundle = bundles.at(0);
	  for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	    FlashTPCBundle *bundle = *it1;
	    if (bundle->get_chi2()/bundle->get_ndf() < min_bundle->get_chi2()/min_bundle->get_ndf())
	      min_bundle = bundle;
	  }

	  if (min_bundle->get_ks_dis() < 0.15){
	    FlashTPCBundleSelection temp_bundles;
	    
	    //examine the rest of bundles;
	    for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	      FlashTPCBundle *bundle = *it1;
	      if (bundle==min_bundle) continue;
	      bool flag_remove = false;
	      
	      if (flash_good_bundles_map[min_bundle->get_flash()].size()==1){
		if (min_bundle->get_ks_dis()<0.08){
		  if (min_bundle->get_chi2()/min_bundle->get_ndf() *1.1 < bundle->get_chi2()/bundle->get_ndf() ){
		    flag_remove = true;
		  }
		}else if (min_bundle->get_ks_dis()<0.1 && min_bundle->get_chi2()/min_bundle->get_ndf() *1.33< bundle->get_chi2()/bundle->get_ndf()){
		  flag_remove = true;
		}else if (min_bundle->get_ks_dis()<0.15 && min_bundle->get_chi2()/min_bundle->get_ndf() * 3 < bundle->get_chi2()/bundle->get_ndf()){
		  flag_remove = true;
		}
	      }else if (flash_good_bundles_map[min_bundle->get_flash()].size()==2){
		if (min_bundle->get_ks_dis()<0.1 && min_bundle->get_chi2()/min_bundle->get_ndf() * 6 < bundle->get_chi2()/bundle->get_ndf()){
		  flag_remove = true;
		}
	      }
	      
	      if (flag_remove){
		if (flash_good_bundles_map[bundle->get_flash()].size()>flash_good_bundles_map[min_bundle->get_flash()].size()){
		  if (bundle!=min_bundle){
		    temp_bundles.push_back(bundle);
		  }
		}
	      }
	    }
	  

	    for (auto it1 = temp_bundles.begin(); it1!=temp_bundles.end(); it1++){
	      FlashTPCBundle *bundle = *it1;
	      PR3DCluster *cluster =  bundle->get_main_cluster();
	      Opflash *flash = bundle->get_flash();
	      
	      bundle->set_consistent_flag(false);
	      flash_good_bundles_map[flash].erase(find(flash_good_bundles_map[flash].begin(),flash_good_bundles_map[flash].end(),bundle));
	      cluster_good_bundles_map[cluster].erase(find(cluster_good_bundles_map[cluster].begin(),cluster_good_bundles_map[cluster].end(),bundle));
	      // remove from 
	    }
	  }
	}
      }
      //finish chi2 ...
      
      

      //

      
      FlashTPCBundleSelection to_be_removed;
      for (auto it = cluster_good_bundles_map.begin(); it!= cluster_good_bundles_map.end(); it++){
	PR3DCluster *cluster =  it->first;
	FlashTPCBundleSelection& bundles = it->second;
	bool flag_remove_other = false;
	for (auto it1 = bundles.begin(); it1!= bundles.end(); it1++){
	  FlashTPCBundle *bundle = *it1;
	  Opflash *flash = bundle->get_flash();
	  
	  if (flash_good_bundles_map[flash].size()==1){
	    flag_remove_other = true;
	    break;
	  }
	}
	if (flag_remove_other){
	  if (cluster_other_bundles_map.find(cluster)!=cluster_other_bundles_map.end()){
	    std::copy(cluster_other_bundles_map[cluster].begin(),
		      cluster_other_bundles_map[cluster].end(),
		      std::back_inserter(to_be_removed));
	  }
	}else{
	  if (cluster_other_bundles_map.find(cluster)!=cluster_other_bundles_map.end()){
	    for (auto it1 = cluster_other_bundles_map[cluster].begin(); it1!= cluster_other_bundles_map[cluster].end(); it1++){
	      FlashTPCBundle *bundle1 = *it1;
	      if (bundle1->get_chi2() > bundle1->get_ndf()*25 && (!bundle1->get_flag_close_to_PMT())){
		to_be_removed.push_back(bundle1);
	      }
	    }
	  }
	}
	
      }

      for (auto it1 = to_be_removed.begin(); it1!=to_be_removed.end(); it1++){
	FlashTPCBundle *bundle = *it1;
	all_bundles.erase(bundle);
	
	Opflash *flash = bundle->get_flash();
	PR3DCluster *cluster = bundle->get_main_cluster();
	
	fc_bundles_map.erase(std::make_pair(flash,cluster));
	
	//	bundles.erase(find(bundles.begin(),bundles.end(),bundle));
	{
	  FlashTPCBundleSelection& temp_bundles = cluster_bundles_map[cluster];
	  temp_bundles.erase(find(temp_bundles.begin(), temp_bundles.end(), bundle));
	  if (temp_bundles.size()==0)
	    cluster_bundles_map.erase(cluster);
	}

	{
	  FlashTPCBundleSelection& temp_bundles = flash_bundles_map[flash];
	  temp_bundles.erase(find(temp_bundles.begin(), temp_bundles.end(), bundle));
	  if (temp_bundles.size()==0)
	    flash_bundles_map.erase(flash);
	}
	
	delete bundle;
      }
      
    }

    for (auto it = cluster_bundles_map.begin(); it!= cluster_bundles_map.end(); it++){
      PR3DCluster *main_cluster = it->first;
      FlashTPCBundleSelection& bundles = it->second;
      bool flag_tight_bundle = false;

      for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	FlashTPCBundle *bundle = *it1;
	if (bundle->get_consistent_flag()){
	  flag_tight_bundle = true;
	  break;
	}
      }

      // if (main_cluster->get_cluster_id()==5)
      // 	std::cout << "Xin: " << flag_tight_bundle << " " << bundles.size() << std::endl;
      
      
      // clean up the map ...
      if (flag_tight_bundle){
	// all_bundles;
	// flash_bundles_map;
	// cluster_bundles_map;
	// std::map<std::pair<Opflash*,PR3DCluster*>,FlashTPCBundle*> fc_bundles_map;
	FlashTPCBundleSelection to_be_removed;
	for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	  FlashTPCBundle *bundle = *it1;
	  if ((!bundle->get_consistent_flag()) && (!bundle->get_flag_at_x_boundary())){
	    to_be_removed.push_back(bundle);
	  }
	  //	  std::cout << bundle->get_flash()->get_flash_id() << " " << bundle->get_main_cluster()->get_cluster_id() << " " << bundle->get_flag_at_x_boundary() << " " << bundle->get_ks_dis() << " " << bundle->get_chi2() << " " << bundle->get_ndf() << std::endl;
	}
	for (auto it1 = to_be_removed.begin(); it1!=to_be_removed.end(); it1++){
	  FlashTPCBundle *bundle = *it1;
	  all_bundles.erase(bundle);
	  
	  Opflash *flash = bundle->get_flash();
	  PR3DCluster *cluster = bundle->get_main_cluster();
	  
	  fc_bundles_map.erase(std::make_pair(flash,cluster));

	  bundles.erase(find(bundles.begin(),bundles.end(),bundle));

	  FlashTPCBundleSelection& temp_bundles = flash_bundles_map[flash];
	  temp_bundles.erase(find(temp_bundles.begin(), temp_bundles.end(), bundle));
	  if (temp_bundles.size()==0)
	    flash_bundles_map.erase(flash);
	  delete bundle;
	}
      }
    }

    
    // finish further examine the bundle ... 
    
    // std::cout << "After Cleaning 2 : " << cluster_bundles_map.size() << " A " << flash_bundles_map.size() << " " << all_bundles.size() << std::endl;


    // std::cout << std::endl << std::endl;
    // for (auto it = all_bundles.begin(); it!=all_bundles.end();it++){
    //   FlashTPCBundle *bundle = *it;
      
    //   if ( bundle->get_consistent_flag() || bundle->get_flag_at_x_boundary())
    // 	std::cout << bundle->get_flash()->get_flash_id() << " " << bundle->get_main_cluster()->get_cluster_id() << " " << bundle->get_flag_at_x_boundary() << " " << bundle->get_ks_dis() << " " << bundle->get_chi2() << " " << bundle->get_ndf() << " " << bundle->get_consistent_flag()  << std::endl;
    // }
    // std::cout << std::endl << std::endl;

    
    // examining flash ... 
    for (auto it = flash_bundles_map.begin(); it!=flash_bundles_map.end(); it++){
      Opflash *flash = it->first;
      FlashTPCBundleSelection& bundles = it->second;
      FlashTPCBundleSelection consistent_bundles;
      FlashTPCBundleSelection remaining_bundles;
      FlashTPCBundleSelection to_be_removed;
      
      for (auto it1 = bundles.begin(); it1!=bundles.end(); it1++){
	FlashTPCBundle *bundle = *it1;
	if (bundle->get_consistent_flag()){
	  consistent_bundles.push_back(bundle);
	}else{
	  remaining_bundles.push_back(bundle);
	}
      }

      if (consistent_bundles.size()>0){
	for (auto it1 = remaining_bundles.begin(); it1!=remaining_bundles.end(); it1++){
	  FlashTPCBundle *bundle1 = *it1;
	  bool flag_remove = true;
	  for (auto it2 = consistent_bundles.begin(); it2!=consistent_bundles.end(); it2++){
	    FlashTPCBundle *bundle2 = *it2;
	    if (bundle2->examine_bundle(bundle1,cos_pe_low, cos_pe_mid)){
	      flag_remove = false;
	      break;
	    }
	  }
	  if (flag_remove)
	    to_be_removed.push_back(bundle1);
	}
      }

      
      for (auto it1 = to_be_removed.begin(); it1!=to_be_removed.end(); it1++){
	FlashTPCBundle *bundle = *it1;
	all_bundles.erase(bundle);
	
	Opflash *flash = bundle->get_flash();
	PR3DCluster *cluster = bundle->get_main_cluster();
	
	fc_bundles_map.erase(std::make_pair(flash,cluster));

	// remaing flash
	bundles.erase(find(bundles.begin(),bundles.end(),bundle));
	
	FlashTPCBundleSelection& temp_bundles = cluster_bundles_map[cluster];
	temp_bundles.erase(find(temp_bundles.begin(), temp_bundles.end(), bundle));
	if (temp_bundles.size()==0)
	  cluster_bundles_map.erase(cluster);
	delete bundle;
      }
    }


    

    
    // std::cout << "After Cleaning 3 : " << cluster_bundles_map.size() << " A " << flash_bundles_map.size() << " " << all_bundles.size() << std::endl;
    
    // std::cout << std::endl << std::endl;
    // for (auto it = all_bundles.begin(); it!=all_bundles.end();it++){
    //   FlashTPCBundle *bundle = *it;
      
    //   if ( bundle->get_consistent_flag() || bundle->get_flag_at_x_boundary())
    // 	std::cout << bundle->get_flash()->get_flash_id() << " " << bundle->get_main_cluster()->get_cluster_id() << " " << bundle->get_flag_at_x_boundary() << " " << bundle->get_ks_dis() << " " << bundle->get_chi2() << " " << bundle->get_ndf() << " " << bundle->get_consistent_flag() << std::endl;
    // }
    // std::cout << std::endl << std::endl;



    // matching code // first round ... 
    {
      // regularization strength ... 
      double lambda = 0.1; // note the coefficient is all around 1
      //form matrix ...
      double fudge_factor1 = 0.06; // add 6% relative uncertainty for pe
      double fudge_factor2 = 1.0; // increase the original uncertainties by 50% ... 
      int num_unknowns = all_bundles.size();
      std::map<PR3DCluster*, int> map_tpc_index;
      std::map<Opflash*, int> map_flash_index;
    
      int tpc_index = 0;
      
      for(auto it=cluster_bundles_map.begin(); it!=cluster_bundles_map.end(); it++){
	PR3DCluster* main_cluster = it->first;
	map_tpc_index[main_cluster] = tpc_index;
	tpc_index++;
      }
      
      // improve the chisquare definition ...
      double delta_track = 0.01; // track can only be used once
      double delta_flash = 0.025;
      //    double delta_flash1 = 0.1;
      
      double num_unused_flash = flash_bundles_map.size() - cluster_bundles_map.size();
      if (num_unused_flash<0) num_unused_flash = 0;
      
      VectorXd M = VectorXd::Zero(32*flash_bundles_map.size()); // measurement from each PMT from each flash
      MatrixXd R = MatrixXd::Zero(32*flash_bundles_map.size(), num_unknowns+flash_bundles_map.size()); // unknowns to measurement matrix
      VectorXd MF = VectorXd::Zero(map_tpc_index.size() + flash_bundles_map.size());
      MatrixXd RF = MatrixXd::Zero(map_tpc_index.size() + flash_bundles_map.size(), num_unknowns + flash_bundles_map.size()); // penalty matrix term
      std::vector<std::pair<Opflash*,PR3DCluster*>> total_pairs;
      std::vector<double> total_weights;
      
      size_t i=0;
      for (auto it = flash_bundles_map.begin(); it != flash_bundles_map.end(); it++){
	Opflash *flash = it->first;
	FlashTPCBundleSelection& bundles = it->second;
	for (size_t j=0;j!=32;j++){
	  double pe = flash->get_PE(j);
	  double pe_err = sqrt(pow(flash->get_PE_err(j)*fudge_factor2,2) + pow(pe*fudge_factor1,2));
	  
	  M(32*i+j) = pe/pe_err;
	  R(32*i+j,num_unknowns+i) = pe/pe_err; // flash alone term
	}
	
	for (size_t j=0;j!=bundles.size();j++){
	  FlashTPCBundle *bundle = bundles.at(j);
	  PR3DCluster* main_cluster = bundle->get_main_cluster();
	  std::vector<double>& pred_pmt_light = bundle->get_pred_pmt_light();
	  for (size_t k=0;k!=32;k++){
	    double pe = flash->get_PE(k);
	    double pe_err = sqrt(pow(flash->get_PE_err(k)*fudge_factor2,2) + pow(pe*fudge_factor1,2));
	    R(32*i+k,total_pairs.end()-total_pairs.begin()) = 1./pe_err * pred_pmt_light.at(k);
	  }
	  
	  total_pairs.push_back(std::make_pair(flash,main_cluster));
	  
	  if (bundle->get_flag_at_x_boundary()){
	    total_weights.push_back(0.2);
	  }else{
	    total_weights.push_back(1.0);
	  }
	}
	
	i++;
      }
      
      for (auto it = flash_bundles_map.begin(); it != flash_bundles_map.end(); it++){
	Opflash *flash = it->first;
	total_weights.push_back(1);
      }
      
      
      
      // normalization of the flashes 
      i=0;
      for (auto it = flash_bundles_map.begin(); it != flash_bundles_map.end(); it++){
	MF(cluster_bundles_map.size()+i) = 0;//num_unused_flash/delta_flash;
	RF(cluster_bundles_map.size()+i,num_unknowns+i) = 1./delta_flash;
	
	map_flash_index[it->first] = num_unknowns + i;
	//	MF(cluster_bundles_map.size()+flash_bundles_map.size()) = num_unused_flash/delta_flash1;
	//	RF(cluster_bundles_map.size()+flash_bundles_map.size(),num_unknowns+i) = 1./delta_flash1;
	i++;
      }
      
      // normalization of the tracks
      i=0;
      for (auto it=cluster_bundles_map.begin(); it!=cluster_bundles_map.end(); it++){
	MF(i) = 1./delta_track;
	i++;
      }
      {
	for (size_t i=0; i!=total_pairs.size(); i++){
	  Opflash *flash = total_pairs.at(i).first;
	  PR3DCluster *main_cluster = total_pairs.at(i).second;
	  RF(map_tpc_index[main_cluster], i) = 1./delta_track;	
	}
      }
      
      
      MatrixXd RT = R.transpose();
      MatrixXd RFT = RF.transpose();
      
      VectorXd W = RT * M + RFT * MF;
      MatrixXd G = RT * R + RFT * RF;

      // for (size_t i = 0; i!= 32*flash_bundles_map.size(); i++){
      // 	std::cout << i << " " << M(i) << std::endl;
      // }
      
      
      WireCell::LassoModel m2(lambda, 100000, 0.01);
      m2.SetData(G, W);

      std::vector<double> init_values;
      init_values.resize(num_unknowns + flash_bundles_map.size(),0);
     
      for (size_t i=0; i!=total_pairs.size(); i++){
	FlashTPCBundle *bundle = fc_bundles_map[std::make_pair(total_pairs.at(i).first, total_pairs.at(i).second)];
	if (bundle->get_consistent_flag()){
	  init_values.at(i) = 1;
	}else if (bundle->get_flag_at_x_boundary()){
	  init_values.at(i) = 0.5;
	}
      }
      m2.Set_init_values(init_values);
      
      for (size_t i=0; i!=total_weights.size(); i++){
	m2.SetLambdaWeight(i,total_weights.at(i));
      }
      m2.Fit();
      VectorXd beta = m2.Getbeta();

      
      FlashTPCBundleSelection to_be_removed;
      for (size_t i=0;i!=total_pairs.size();i++){
	//	std::cout << i << " " << beta(i)  << std::endl;
	if (beta(i) < 0.05){
	  to_be_removed.push_back(fc_bundles_map[std::make_pair(total_pairs.at(i).first, total_pairs.at(i).second)]);
	}
      }
      
      for (auto it1 = to_be_removed.begin(); it1!=to_be_removed.end(); it1++){
	FlashTPCBundle *bundle = *it1;
	all_bundles.erase(bundle);
	
	Opflash *flash = bundle->get_flash();
	PR3DCluster *cluster = bundle->get_main_cluster();
	
	fc_bundles_map.erase(std::make_pair(flash,cluster));
	
	{
	  FlashTPCBundleSelection& temp_bundles = flash_bundles_map[flash];
	  temp_bundles.erase(find(temp_bundles.begin(), temp_bundles.end(), bundle));
	  if (temp_bundles.size()==0)
	    flash_bundles_map.erase(flash);
	}
	{
	  FlashTPCBundleSelection& temp_bundles = cluster_bundles_map[cluster];
	  temp_bundles.erase(find(temp_bundles.begin(), temp_bundles.end(), bundle));
	  if (temp_bundles.size()==0)
	    cluster_bundles_map.erase(cluster);
	}
	
	delete bundle;
      }
    }
    
    //    std::cout << "After Cleaning 4 : " << cluster_bundles_map.size() << " A " << flash_bundles_map.size() << " " << all_bundles.size() << std::endl;
    
    // regularization strength ... 
    double lambda = 0.1; // note the coefficient is all around 1
    //form matrix ...
    double fudge_factor1 = 0.05; // add 6% relative uncertainty for pe
    double fudge_factor2 = 1.0; // increase the original uncertainties by 50% ... 
    int num_unknowns = all_bundles.size();
    std::map<PR3DCluster*, int> map_tpc_index;
    std::map<Opflash*, int> map_flash_index;
    
    int tpc_index = 0;
    
    for(auto it=cluster_bundles_map.begin(); it!=cluster_bundles_map.end(); it++){
      PR3DCluster* main_cluster = it->first;
      map_tpc_index[main_cluster] = tpc_index;
      tpc_index++;
    }
    
    // improve the chisquare definition ...
    double delta_track = 0.01; // track can only be used once
    
    VectorXd M = VectorXd::Zero(32*flash_bundles_map.size()); // measurement from each PMT from each flash
    MatrixXd R = MatrixXd::Zero(32*flash_bundles_map.size(), num_unknowns); // unknowns to measurement matrix
    VectorXd MF = VectorXd::Zero(map_tpc_index.size() );
    MatrixXd RF = MatrixXd::Zero(map_tpc_index.size() , num_unknowns ); // penalty matrix term
    std::vector<std::pair<Opflash*,PR3DCluster*>> total_pairs;
    std::vector<double> total_weights;
    
    size_t i=0;
    for (auto it = flash_bundles_map.begin(); it != flash_bundles_map.end(); it++){
      Opflash *flash = it->first;
      FlashTPCBundleSelection& bundles = it->second;
      for (size_t j=0;j!=32;j++){
	double pe = flash->get_PE(j);
	double pe_err = sqrt(pow(flash->get_PE_err(j)*fudge_factor2,2) + pow(pe*fudge_factor1,2));
	
	M(32*i+j) = pe/pe_err;
	//	R(32*i+j,num_unknowns+i) = pe/pe_err; // flash alone term
      }
      
      for (size_t j=0;j!=bundles.size();j++){
	FlashTPCBundle *bundle = bundles.at(j);
	PR3DCluster* main_cluster = bundle->get_main_cluster();
	std::vector<double>& pred_pmt_light = bundle->get_pred_pmt_light();
	for (size_t k=0;k!=32;k++){
	  double pe = flash->get_PE(k);
	  double pe_err = sqrt(pow(flash->get_PE_err(k)*fudge_factor2,2) + pow(pe*fudge_factor1,2));
	  R(32*i+k,total_pairs.end()-total_pairs.begin()) = 1./pe_err * pred_pmt_light.at(k);
	}
	
	total_pairs.push_back(std::make_pair(flash,main_cluster));
	
	if (bundle->get_flag_at_x_boundary()){
	  total_weights.push_back(0.2);
	}else{
	  total_weights.push_back(1.0);
	}
      }
      
      i++;
    }
    
    // normalization of the tracks
    i=0;
    for (auto it=cluster_bundles_map.begin(); it!=cluster_bundles_map.end(); it++){
      MF(i) = 1./delta_track;
      i++;
    }
    {
      for (size_t i=0; i!=total_pairs.size(); i++){
	Opflash *flash = total_pairs.at(i).first;
	PR3DCluster *main_cluster = total_pairs.at(i).second;
	RF(map_tpc_index[main_cluster], i) = 1./delta_track;	
      }
    }
    
      
    MatrixXd RT = R.transpose();
    MatrixXd RFT = RF.transpose();
    
    VectorXd W = RT * M + RFT * MF;
    MatrixXd G = RT * R + RFT * RF;
    
    WireCell::LassoModel m2(lambda, 100000, 0.01);
    m2.SetData(G, W);
    for (size_t i=0; i!=total_weights.size(); i++){
      m2.SetLambdaWeight(i,total_weights.at(i));
    }
    m2.Fit();
    VectorXd beta = m2.Getbeta();
    
    // for (auto it = flash_bundles_map.begin(); it != flash_bundles_map.end(); it++){
    //   Opflash *flash = it->first;
    //   std::cout << flash->get_flash_id() << " " << beta(map_flash_index[flash]) << std::endl;
    // }
    
    std::map<int,std::pair<Opflash*,double>> matched_pairs;
    for (size_t i=0;i!=total_pairs.size();i++){
      if(beta(i)!=0){
  	int tpc_index = map_tpc_index[total_pairs.at(i).second];
  	Opflash *flash = total_pairs.at(i).first;
  	if (matched_pairs.find(tpc_index)==matched_pairs.end()){
  	  matched_pairs[tpc_index] = std::make_pair(flash,beta(i));
  	}else{
  	  if (beta(i) > matched_pairs[tpc_index].second){
  	    matched_pairs[tpc_index] = std::make_pair(flash,beta(i));
  	  }
  	}
	//std::cout << i << " Q " <<  tpc_index << " " << flash->get_flash_id() << " " << total_pairs.at(i).second->get_cluster_id() << " " << total_weights.at(i) << " " << beta(i)  << " " << flash->get_time() << std::endl;
      }
    }
   
    // Need some organization ...
    
    
    FlashTPCBundleSelection results_bundles;
    // return bundles ...    
    for (auto it = group_clusters.begin(); it!=group_clusters.end(); it++){
      PR3DCluster* main_cluster = it->first;
      if (map_tpc_index.find(main_cluster)!=map_tpc_index.end()){
	int tpc_index = map_tpc_index[main_cluster];
	if (matched_pairs.find(tpc_index)!=matched_pairs.end()){
	  Opflash* flash = matched_pairs[tpc_index].first;
	  double strength = matched_pairs[tpc_index].second;
	  FlashTPCBundle* bundle = fc_bundles_map[std::make_pair(flash,main_cluster)];
	  bundle->set_strength(strength);
	  results_bundles.push_back(bundle);
	}else{
	  Opflash *flash = 0;
	  double strength  =0;
	  FlashTPCBundle *bundle = new FlashTPCBundle(flash,main_cluster,-1,tpc_index);
	  results_bundles.push_back(bundle);
	}
      }else{
	Opflash *flash = 0;
	double strength  =0;
	FlashTPCBundle *bundle = new FlashTPCBundle(flash,main_cluster,-1,-1);
	results_bundles.push_back(bundle);
      }
    }
    
    return results_bundles;
    
  }
  
}
