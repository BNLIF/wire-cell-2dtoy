#include "WireCell2dToy/ToyMatching.h"

#include "WireCellRess/LassoModel.h"
#include "WireCellRess/ElasticNetModel.h"
#include <Eigen/Dense>

#include "TChain.h"

#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"


using namespace Eigen;
using namespace WireCell;

int WireCell2dToy::convert_xyz_voxel_id(WireCell::Point &p){
  // int voxel_x_id = std::round((p.x/units::cm+64.825-5.14667/2.)/5.14667);
  // int voxel_y_id = std::round((p.y/units::cm+193-5.14667/2.)/5.14667);
  // int voxel_z_id = std::round((p.z/units::cm+128.243-3.23122/2.)/3.23122);
  
  int voxel_x_id = std::round((p.x/units::cm+64.825-5.14667/2.)/5.14667);
  int voxel_y_id = std::round((p.y/units::cm+188.416-5.11698/2.)/5.11698);
  int voxel_z_id = std::round((p.z/units::cm+88.3554-3.05768/2.)/3.05768);
  
  if (voxel_x_id<0) voxel_x_id=0;
  if (voxel_x_id>=75) voxel_x_id=74;
  if (voxel_y_id<0) voxel_y_id=0;
  if (voxel_y_id>=75) voxel_y_id=74;
  if (voxel_z_id<0) voxel_z_id=0;
  if (voxel_z_id>=400) voxel_z_id=399;

  int voxel_id = voxel_z_id*75*75 + voxel_y_id*75 + voxel_x_id;
  return voxel_id;
}


void WireCell2dToy::tpc_light_match(int time_offset, int nrebin, std::map<WireCell::PR3DCluster*,std::vector<std::pair<WireCell::PR3DCluster*,double>>>& group_clusters, WireCell::OpflashSelection& flashes){
  TChain *T = new TChain("/pmtresponse/PhotonLibraryData","/pmtresponse/PhotonLibraryData");
  T->AddFile("./uboone_photon_library.root");
  //  std::cout << T->GetEntries();
  Int_t Voxel;
  Int_t OpChannel;
  Float_t Visibility;
  T->SetBranchAddress("Voxel",&Voxel);
  T->SetBranchAddress("OpChannel",&OpChannel);
  T->SetBranchAddress("Visibility",&Visibility);

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
  map_lib_pmt[25]=28; map_pmt_lib[28]=25; 
  map_lib_pmt[27]=30; map_pmt_lib[30]=27; 
  map_lib_pmt[28]=31; map_pmt_lib[31]=28; 
  map_lib_pmt[31]=29; map_pmt_lib[29]=31; 
  map_lib_pmt[30]=27; map_pmt_lib[27]=30; 
  map_lib_pmt[29]=26; map_pmt_lib[26]=29; 
  
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
  double high_x_cut_ext1 = + 1*units::cm;
  double low_x_cut = 0*units::cm;
  double low_x_cut_ext1 = - 2*units::cm;
  double scaling_light_mag = 0.01;


  std::vector<Opflash*> good_flashes; // save flashes 
  std::map<Opflash*,std::vector<int>> map_flash_tpc_ids; // save tpc ids
  std::map<Opflash*,std::vector<std::vector<double>>> map_flash_tpc_light_preds; // save tpc predictions
  std::map<Opflash*,std::vector<bool>> map_flash_tpc_boundaries; // save if the tpc predictions is at boundaries ... 
  std::map<Opflash*,int> map_flash_index;
  
  
  
  int flash_num = 0;
  for (auto it1 =flashes.begin(); it1!=flashes.end(); it1++){
    Opflash *flash = (*it1);
    double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
    int cluster_id = 0;
    bool flag_good_flash = false;
    
    for (auto it2 = group_clusters.begin(); it2!=group_clusters.end(); it2++){
      PR3DCluster* main_cluster = it2->first;
      std::vector<std::pair<WireCell::PR3DCluster*,double>>& more_clusters = it2->second;
      // judge if the main clusters are within the detector range ...
      double first_pos_x = (*((main_cluster->get_time_cells_set_map().begin())->second.begin()))->get_sampling_points().front().x;
      double last_pos_x = (*((main_cluster->get_time_cells_set_map().rbegin())->second.begin()))->get_sampling_points().front().x;

      bool flag_at_x_boundary = false;
      if (first_pos_x-offset_x > low_x_cut + low_x_cut_ext1 &&
	  last_pos_x-offset_x > low_x_cut &&
	  last_pos_x-offset_x < high_x_cut + high_x_cut_ext1 &&
	  first_pos_x-offset_x < high_x_cut){
	flag_good_flash = true;

	if (map_flash_tpc_ids.find(flash)==map_flash_tpc_ids.end()){
	  std::vector<int> temp_vec_int;
	  map_flash_tpc_ids[flash] = temp_vec_int;
	  std::vector<std::vector<double>> temp_vec_preds;
	  map_flash_tpc_light_preds[flash] = temp_vec_preds;
	  std::vector<bool> temp_vec_flag;
	  map_flash_tpc_boundaries[flash] = temp_vec_flag;
	}
	
	
	// tracks ends at boundary ... 
	if (first_pos_x-offset_x <=low_x_cut && first_pos_x-offset_x > low_x_cut + low_x_cut_ext1 ||
	    last_pos_x-offset_x >= high_x_cut && last_pos_x-offset_x < high_x_cut + high_x_cut_ext1){
	  flag_at_x_boundary = true;
	}
	
	std::vector<double> pred_pmt_light;
	//std::vector<double> pred_pmt_light1;
	pred_pmt_light.resize(32,0);
	//pred_pmt_light1.resize(32,0);
	
	PR3DClusterSelection temp_clusters;
	// fill in stuff for the main cluster
	temp_clusters.push_back(main_cluster);
	// fill in the stuff for additional clusters
	for (auto it3 = more_clusters.begin(); it3!=more_clusters.end(); it3++){
	  temp_clusters.push_back(it3->first);
	}
	for (auto it3 = temp_clusters.begin(); it3!=temp_clusters.end(); it3++){
	  SMGCSelection& mcells = (*it3)->get_mcells();
	  for (auto it4 = mcells.begin(); it4!=mcells.end(); it4++){
	    SlimMergeGeomCell *mcell = (*it4);
	    if (mcell->get_q()>0){
	      PointVector& pts = mcell->get_sampling_points();
	      if (pts.at(0).x-offset_x < low_x_cut+low_x_cut_ext1 ||
		  pts.at(0).x-offset_x > high_x_cut+high_x_cut_ext1)
		continue;
	      float charge = mcell->get_q()/pts.size();
	      Point p;
	      for (size_t i=0;i!=pts.size();i++){
		p.x = pts.at(i).x - offset_x;
		p.y = pts.at(i).y;
		p.z = pts.at(i).z;

		// if (flash_num==27&&flag_at_x_boundary){
		//   std::cout << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << " "  << flash_num << " " << cluster_id << std::endl;
		// }
		
		int voxel_id = WireCell2dToy::convert_xyz_voxel_id(p);
		std::list<std::pair<int,float>>& pmt_list = photon_library.at(voxel_id);
		for (auto it5 = pmt_list.begin(); it5!=pmt_list.end(); it5++){
		  pred_pmt_light.at(it5->first) += charge * it5->second;
		  //pred_pmt_light1.at(it5->first) += charge * it5->second;
		}
	      }
	    }
	  }
	}

	for (size_t i=0;i!=32;i++){
	  pred_pmt_light.at(i) *= scaling_light_mag;
	}

	map_flash_tpc_ids[flash].push_back(cluster_id);
	map_flash_tpc_light_preds[flash].push_back(pred_pmt_light);
	
	// fill in the content ... 
	if (flag_at_x_boundary){
	  map_flash_tpc_boundaries[flash].push_back(true);
	  //	  for (size_t i=0;i!=32;i++){
	  //  std::cout << flash_num << " " << cluster_id << " " << i << " " << flash->get_PE(i) << " " << flash->get_PE_err(i) << " " << pred_pmt_light.at(map_pmt_lib[i]) << " " <<  std::endl;
	  // }
	}else{
	  map_flash_tpc_boundaries[flash].push_back(false);
	}
      }
      cluster_id++;
    }

    if (flag_good_flash){
      good_flashes.push_back(flash);
      map_flash_index[flash] = flash_num;
    }
    flash_num ++;
  }

  
  // regularization strength ... 
  double lambda = 32;
  //form matrix ...
  double fudge_factor = 0.05; // add 5% relative uncertainty for pe
  int num_unknowns = 0;
  
  for (auto it = map_flash_tpc_ids.begin(); it!= map_flash_tpc_ids.end(); it++){
    num_unknowns += it->second.size();
  }
  
  VectorXd W = VectorXd::Zero(32*good_flashes.size());
  MatrixXd G = MatrixXd::Zero(32*good_flashes.size(), num_unknowns);
  std::vector<std::pair<Opflash*,int>> total_pairs;
  std::vector<double> total_weights;
  
  for (size_t i=0; i!= good_flashes.size(); i++){
    Opflash *flash = good_flashes.at(i);
    for (size_t j=0;j!=32;j++){
      double pe = flash->get_PE(j);
      double pe_err = sqrt(pow(flash->get_PE_err(j),2) + pow(pe*fudge_factor,2));
      W(32*i+j) = pe/pe_err;
    }
    for (size_t j=0;j!=map_flash_tpc_ids[flash].size();j++){
   
      for (size_t k=0;k!=32;k++){
	double pe = flash->get_PE(k);
	double pe_err = sqrt(pow(flash->get_PE_err(k),2) + pow(pe*fudge_factor,2));
	G(32*i+k,total_pairs.end()-total_pairs.begin()) = 1./pe_err * map_flash_tpc_light_preds[flash].at(j).at(map_pmt_lib[k]);
      }
      total_pairs.push_back(std::make_pair(flash,map_flash_tpc_ids[flash].at(j)));
      if (map_flash_tpc_boundaries[flash].at(j)){ // add boundary ... 
	total_weights.push_back(0.1);
      }else{
	total_weights.push_back(1);
      }
    }
  }


  // std::cout << 32*good_flashes.size() << " " << num_unknowns << " " << total_weights.size() << " " << total_pairs.size() << std::endl;
  

  WireCell::LassoModel m2(lambda, 100000, 0.05);
  m2.SetData(G, W);
  for (size_t i=0; i!=total_weights.size(); i++){
    m2.SetLambdaWeight(i,total_weights.at(i));
  }
  m2.Fit();
  VectorXd beta = m2.Getbeta();
  for (size_t i=0;i!=total_pairs.size();i++){
   std::cout << i << " " <<  map_flash_index[total_pairs.at(i).first] << " " << total_pairs.at(i).second << " " << total_weights.at(i) << " " << beta(i)  << std::endl;
  }

  // for (size_t i=0;i!=32;i++){
  //   std::cout << i << " " << W(32*4+i) << " " << G(32*4+i,60) << " " << total_weights.at(60) << std::endl;
  // }
  
  
  // number of unknonws are matched pairs #matched_pair < num_flashes * tpc_object
  // Matrix would then be (32*num_flashes, #matched_pair)
  // for each coefficient, we need to save its flash and the tpc_object
  // for each tpc object, we need to save all its coefficients, so we can pick up the latest one
  
  // For each TPC object matched with a particular, need to generate 


  
  
}
