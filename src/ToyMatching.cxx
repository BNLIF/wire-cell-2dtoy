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
  int voxel_x_id = std::round((p.x/units::cm+64.825-5.14667/2.)/5.14667);
  int voxel_y_id = std::round((p.y/units::cm+193-5.14667/2.)/5.14667);
  int voxel_z_id = std::round((p.z/units::cm+128.243-3.23122/2.)/3.23122);
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

  int flash_num = 0;
  for (auto it1 =flashes.begin(); it1!=flashes.end(); it1++){
    Opflash *flash = (*it1);
    double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
    int cluster_id = 0;
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

	// tracks ends at boundary ... 
	if (first_pos_x-offset_x <=low_x_cut && first_pos_x-offset_x > low_x_cut + low_x_cut_ext1 ||
	    last_pos_x-offset_x >= high_x_cut && last_pos_x-offset_x < high_x_cut + high_x_cut_ext1){
	  flag_at_x_boundary = true;
	}

	std::vector<double> pred_pmt_light;
	pred_pmt_light.resize(32,0);

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

		if (flash_num==27&&flag_at_x_boundary){
		  std::cout << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << " "  << flash_num << " " << cluster_id << std::endl;
		}
		
		int voxel_id = WireCell2dToy::convert_xyz_voxel_id(p);
		std::list<std::pair<int,float>>& pmt_list = photon_library.at(voxel_id);
		for (auto it5 = pmt_list.begin(); it5!=pmt_list.end(); it5++){
		  pred_pmt_light.at(it5->first) += charge * it5->second;
		}
	      }
	    }
	  }
	}

	

	if (flag_at_x_boundary){
	  for (size_t i=0;i!=32;i++){
	    std::cout << flash_num << " " << cluster_id << " " << i << " " << flash->get_PE(i) << " " << flash->get_PE_err(i) << " " << pred_pmt_light.at(i) << std::endl;
	  }
	}
	
      }
      cluster_id++;
    }
    flash_num ++;
  }
  // number of unknonws are matched pairs #matched_pair < num_flashes * tpc_object
  // Matrix would then be (32*num_flashes, #matched_pair)
  // for each coefficient, we need to save its flash and the tpc_object
  // for each tpc object, we need to save all its coefficients, so we can pick up the latest one
  
  // For each TPC object matched with a particular, need to generate 


  
  
}
