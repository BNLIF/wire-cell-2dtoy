#include "WCP2dToy/ToyFiducial.h"
#include "WCPData/TPCParams.h"
#include "WCPData/Singleton.h"

using namespace WCP;

int pnpoly(std::vector<double>& vertx, std::vector<double>& verty, double testx, double testy)
{
  int i, j, c = 0;
  for (i = 0, j = int(vertx.size())-1; i < int(vertx.size()); j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
      c = !c;
  }
  return c;
}


//Helper function that returns the number of boundary contacts for a set of extreme points, drift offset, and tolerances
//checks whether a given cluster's extreme points are touching any detector boundaries for a given flash time / drift offest_x
//-2 means short track, -1 means outside boundary, 0 means inside boudnary, 1 means STM, 2 meanst TGM
int WCP2dToy::ToyFiducial::check_boundary(std::vector<std::vector<WCP::WCPointCloud<double>::WCPoint>> extreme_points, double offset_x, std::vector<double>* tol_vec){

	//Check whether the extreme points are contained within the boundary planes (allowing for tolerance), and if they are near a boundary.
	bool front_flag = false;
	bool back_flag = false;
	for(int i=0;i<int(extreme_points.size());i++){
		for(int ii=0;ii<int(extreme_points[i].size());ii++){
			Point p(extreme_points[i][ii].x,extreme_points[i][ii].y,extreme_points[i][ii].z);
			//Check whether the point is inside the extended fiducial volume
			std::vector<double> neg_tol_vec = {-1*tol_vec->at(0),-1*tol_vec->at(1),-1*tol_vec->at(2),-1*tol_vec->at(3)};
			if(!inside_fiducial_volume(p,offset_x,tol_vec)){
				return -1;
			//Now that the point is known to be within the extended fiducial volume, check whether it is near any TGM boundaries, but only if it is near a PCA endpoint
			} else if(!inside_fiducial_volume(p,offset_x,&neg_tol_vec)){
				if(i==0){front_flag = true;}
				if(i==1){back_flag = true;}
			}
		}
	}
	return front_flag + back_flag;
}

//Main cosmic tagger function
void WCP2dToy::ToyFiducial::cosmic_tagger(WCP::OpflashSelection& flashes, FlashTPCBundleSelection *matched_bundles, FlashTPCBundle* main_bundle, WCP::Photon_Library *pl, int time_offset, int nrebin, float unit_dis, WCP::ToyCTPointCloud& ct_point_cloud,
						std::map<WCP::PR3DCluster*, WCP::PR3DCluster*>& old_new_cluster_map, int run_no, int subrun_no, int event_no, bool flag_data, bool debug_tagger){

	//std::cout << "starting cosmic tagger ===================================================" << std::endl;

	//Take in flash and cluster info
	Opflash* main_flash = main_bundle->get_flash();
	PR3DCluster *main_cluster = main_bundle->get_main_cluster();
	std::vector<std::vector<WCPointCloud<double>::WCPoint>> extreme_points = main_cluster->get_extreme_wcps();
	double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;
	double main_offset_x = (main_flash->get_time() - time_offset)*2./nrebin*time_slice_width;

	//set tolerances
	double ks_main_stm_tol = 0.08;					//This is the tolerance below which a flash is considered a very good match and should not be considered
	double ks_main_tgm_tol = 0.00;					//TGM-with-flash is pure already, no need to screen by KS
	double ks_main_tgm_nof_tol = 0.08;
	double stm_flash_tol = 3.0*units::m;
	double tgm_flash_tol = 3.0*units::m;
	double stm_pe_frac_tol = 2.0;
	double tgm_pe_frac_tol = 2.0;
	std::vector<double> stm_tol_vec =     {1.0, 1.0, 1.5, 1.0};	//x, ybot, ytop, z
	std::vector<double> tgm_tol_vec =     {2.0, 2.0, 3.0, 3.0};
	std::vector<double> tgm_nof_tol_vec = {0.8, 1.0, 1.2, 1.0};
	for(int i=0;i<4;i++){
		stm_tol_vec[i] *= units::cm;
		tgm_tol_vec[i] *= units::cm;
		tgm_nof_tol_vec[i] *= units::cm;
	}

	//Get existing match parameters
	//Later they are used to skip cases where the neutrino is already well matched
	//Only KS is used right now
	double chi2_main = main_bundle->get_chi2();
	double ndf_main = main_bundle->get_ndf();
	double ks_main = main_bundle->get_ks_dis();
	bool bad_match_main = main_bundle->get_potential_bad_match_flag();
	std::cout << "chi2/ndf = " << chi2_main/ndf_main << ",  \t  ks = " << ks_main << ",  \t  bm_flag = " << bad_match_main << std::endl;

	//PMT map and geometry info
	std::vector<int> OpDet_to_OpChannel_map = {4,2,5,0,6,1,3,9,11,7,12,8, 10,15,13,17,18,14,16,21,23,19,20,24,22,28,25,30,31,26,27,29};
	std::vector<int> OpChannel_to_OpDet_map = {3,5,1,6,0,2,4,9,11,7,12,8,10,14,17,13,18,15,16,21,22,19,24,20,23,26,29,30,25,31,27,28};
	std::vector<std::vector<double>> opDet_xyz_map = {{2.645, -28.625, 990.356},{2.682,  27.607, 989.712},{2.324, -56.514, 951.865},{2.458,  55.313, 951.861},{2.041, -56.309, 911.939},{2.265,  55.822, 911.066},{1.923, -0.722,  865.598},{1.795, -0.502,  796.208},{1.495, -56.284, 751.905},{1.559,  55.625, 751.884},{1.487, -56.408, 711.274},{1.438,  55.800, 711.073},{1.475, -0.051,  664.203},{1.448, -0.549,  585.284},{1.226,  55.822, 540.929},{1.479, -56.205, 540.616},{1.505, -56.323, 500.221},{1.116,  55.771, 500.134},{1.481, -0.875,  453.096},{1.014, -0.706,  373.839},{1.451, -57.022, 328.341},{0.913,  54.693, 328.212},{0.682,  54.646, 287.976},{1.092, -56.261, 287.639},{0.949, -0.829,  242.014},{0.658, -0.303,  173.743},{0.703,  55.249, 128.355},{0.821, -56.203, 128.179},{0.862, -56.615, 87.8695},{0.558,  55.249, 87.7605},{0.665,  27.431, 51.1015},{0.947, -28.576, 50.4745}};

	//Vectors that store STM / TGM / TGM-noflash cases
	std::vector<int> n_boundary_list;
	std::vector<Opflash*> candidate_flash_list;
	std::vector<double> candidate_offset_x_list;

	//Check whether the track is at least 15 cm long
	Point p0(extreme_points[0][0].x,extreme_points[0][0].y,extreme_points[0][0].z);
	Point p1(extreme_points[1][0].x,extreme_points[1][0].y,extreme_points[1][0].z);
	double distance = sqrt(pow(p0.x-p1.x,2)+pow(p0.y-p1.y,2)+pow(p0.z-p1.z,2));
	if(distance > 10*units::cm){
		//iterate over all the flashes to check each flash for STM/TGM conditions
		for (auto it = flashes.begin(); it!= flashes.end(); it++){
			Opflash *flash = (*it);
			double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;

			//Skip the already matched flash
			if(flash->get_time() == main_flash->get_time()){
				continue;
			}
			//Skip cases where the flash is clearly out of the beam window or on the anode side
			if(extreme_points[0][0].x-offset_x < 128*units::cm-tgm_tol_vec[0] || extreme_points[0][0].x-offset_x > 256*units::cm+tgm_tol_vec[0]
			|| extreme_points[1][0].x-offset_x < 128*units::cm-tgm_tol_vec[0] || extreme_points[1][0].x-offset_x > 256*units::cm+tgm_tol_vec[0]){
				continue;
			}
			//Skip cases where the replacement flash is well matched
			bool flash_well_matched = false;
			for (auto it2 = matched_bundles->begin(); it2!= matched_bundles->end(); it2++){
				FlashTPCBundle *old_bundle = *it2;
				Opflash* old_flash = old_bundle->get_flash();
				if(old_flash){
					//find the previous bundle for the flash being iterated over, if it exists
					if(flash->get_flash_id() == old_bundle->get_flash()->get_flash_id()){
						double chi2 = old_bundle->get_chi2();
						double ndf = old_bundle->get_ndf();
						double ks = old_bundle->get_ks_dis();
						bool bad_match = old_bundle->get_potential_bad_match_flag();
						//std::cout << "flash # " << flash->get_flash_id() << ",  \t  chi2 = " << chi2 << ",  \t  ndf = " << ndf << ",  \t  ks = " << ks << ",  \t  bm_flag = " << bad_match << std::endl;
						//Should only KS be used?  Anything useful at all?
	//					if(chi2/ndf < 1.5 && ks < 0.1 && !bad_match){
	//						flash_well_matched = true;
	//					}
						//don't keep searching after the flash-bundle match is found
						break;
					}
				}
			}
			if(flash_well_matched){
				continue;
			}

			//Flash PE Prediction
			std::vector<std::pair<WCP::PR3DCluster*,double>> additional_clusters;
			PR3DClusterSelection other_clusters;
			PR3DClusterSelection more_clusters;
			bool flag_good_bundle;
			FlashTPCBundle *new_bundle =  new FlashTPCBundle(flash, main_cluster,flash->get_flash_id(),main_cluster->get_cluster_id());
			std::vector<double>& pred_pmt_light = new_bundle->get_pred_pmt_light();
			WCP2dToy::calculate_pred_pe(run_no, time_offset, nrebin, time_slice_width, pl, new_bundle, &pred_pmt_light, &additional_clusters, &other_clusters, &more_clusters, flag_good_bundle, flag_data);
			//std::vector<double>& pred_pmt_light = main_bundle->get_pred_pmt_light();

			//Compute the PE centroids for the flash PE and predicted PE
			double pred_pe_tot = 0;
			double flash_pe_tot = 0;
			double pred_pe_z_centroid = 0;
			double flash_pe_z_centroid = 0;
			for(int i=0;i<int(pred_pmt_light.size());i++){
				double flash_pe = flash->get_PE(i);
				flash_pe_tot += flash_pe;
				pred_pe_tot += pred_pmt_light[i];
				flash_pe_z_centroid += flash_pe*opDet_xyz_map[i][2]*units::cm;
				pred_pe_z_centroid += pred_pmt_light[i]*opDet_xyz_map[i][2]*units::cm;
			}
			flash_pe_z_centroid /= flash_pe_tot;
			pred_pe_z_centroid  /= pred_pe_tot;
			delete new_bundle;

			//Geometry evaluation
			//Call the helper function check_boundary to get the number of boundary intersections for a given flash
			//Only store a STM / TGM if the flash/predicted-PE roughly match
			int boundary_num_tgm = 0;
			int boundary_num_stm = 0;
			int boundary_num = 0;
			//TGM-with-flash case
			if(ks_main > ks_main_tgm_tol && flash_pe_tot < pred_pe_tot*tgm_pe_frac_tol && flash_pe_z_centroid > pred_pe_z_centroid-tgm_flash_tol && flash_pe_z_centroid < pred_pe_z_centroid+tgm_flash_tol){
				boundary_num_tgm = check_boundary(extreme_points, offset_x, &tgm_tol_vec);
			}
			//STM-with-flash case
			if(ks_main > ks_main_stm_tol && flash_pe_tot < pred_pe_tot*stm_pe_frac_tol && flash_pe_z_centroid > pred_pe_z_centroid-stm_flash_tol && flash_pe_z_centroid < pred_pe_z_centroid+stm_flash_tol){
				boundary_num_stm = check_boundary(extreme_points, offset_x, &stm_tol_vec);
			}
			//TGM gets priority over STM
			if(boundary_num_tgm==2)					{boundary_num = 2;}
			else if(boundary_num_stm==1)				{boundary_num = 1;}
			else if(boundary_num_tgm<0 || boundary_num_stm<0)	{boundary_num = -1;}

			//Store a STM / TGM tag if appropriate
			if(boundary_num==1 || (boundary_num==2 && check_tgm(main_bundle,offset_x,ct_point_cloud,old_new_cluster_map,1)) ){
				n_boundary_list.push_back(boundary_num);
				candidate_flash_list.push_back(flash);
			}
		}

		//Check whether a TGM w/o flash hypothesis works
		//Skip if the existing neutrino match is very good
		if(ks_main > ks_main_tgm_nof_tol){
			double max_offset =  1000000*units::cm;		//The max offset allowed in the +x direction.  The starting value is meant to be immediately overriden.
			double offset_x = (main_bundle->get_flash()->get_time() - time_offset)*2./nrebin*time_slice_width;
			double tx = tgm_nof_tol_vec[0];
			double ty_bot = tgm_nof_tol_vec[1];
			double ty_top = tgm_nof_tol_vec[2];
			double tz = tgm_nof_tol_vec[3];
			//Iterate extreme points.  For each, compute the drift distance until the boundary+tolerance is crossed.
			//Set the max_offset variable to the lowest drift distance that causes any extreme point to (nearly) leave the boundary+tolerance.
			for(int i=0;i<int(extreme_points.size());i++){
				for(int ii=0;ii<int(extreme_points[i].size());ii++){
					Point p(extreme_points[i][ii].x,extreme_points[i][ii].y,extreme_points[i][ii].z);

					std::vector<double> boundary_points_xy_x = get_boundary_SCB_xy_x(p);
					std::vector<double> boundary_points_xy_y = get_boundary_SCB_xy_y(p);
					std::vector<double> boundary_points_xz_x = get_boundary_SCB_xz_x(p);
					std::vector<double> boundary_points_xz_z = get_boundary_SCB_xz_z(p);
					//Don't let the offset be larger than the distance to the cathode side
					max_offset = std::min(max_offset,boundary_points_xy_x[2]-p.x+tx);
					max_offset = std::min(max_offset,boundary_points_xz_x[2]-p.x+tx);

					//Conditions in Y and Z for whether the endpoint will drift into the SCB before the cathode side
					if     (p.y > boundary_points_xy_y[1] && p.y < boundary_points_xy_y[2]-ty_bot){
						double x0 = boundary_points_xy_x[1];	double y0 = boundary_points_xy_y[1]-ty_bot;
						double x1 = boundary_points_xy_x[2]+tx;	double y1 = boundary_points_xy_y[2]-ty_bot;
						double x_intersection = x0 + (x1-x0)*(p.y-y0)/(y1-y0);
						max_offset = std::min(max_offset,x_intersection-p.x);
					} else if(p.y > boundary_points_xy_y[3]+ty_top && p.y < boundary_points_xy_y[4]){
						double x0 = boundary_points_xy_x[3]+tx;	double y0 = boundary_points_xy_y[3]+ty_top;
						double x1 = boundary_points_xy_x[4];	double y1 = boundary_points_xy_y[4]+ty_top;
						double x_intersection = x0 + (x1-x0)*(p.y-y0)/(y1-y0);
						max_offset = std::min(max_offset,x_intersection-p.x);
					}
					if     (p.z > boundary_points_xz_z[1] && p.z < boundary_points_xz_z[2]-tz){
						double x0 = boundary_points_xz_x[1];	double z0 = boundary_points_xz_z[1]-tz;
						double x1 = boundary_points_xz_x[2]+tx;	double z1 = boundary_points_xz_z[2]-tz;
						double x_intersection = x0 + (x1-x0)*(p.z-z0)/(z1-z0);
						max_offset = std::min(max_offset,x_intersection-p.x);
					} else if(p.z > boundary_points_xz_z[3]+tz && p.z < boundary_points_xz_z[4]){
						double x0 = boundary_points_xz_x[3]+tx;	double z0 = boundary_points_xz_z[3]+tz;
						double x1 = boundary_points_xz_x[4];	double z1 = boundary_points_xz_z[4]+tz;
						double x_intersection = x0 + (x1-x0)*(p.z-z0)/(z1-z0);
						max_offset = std::min(max_offset,x_intersection-p.x);
					}
				}
			}
			//Keep the offset within the boundary, and then check for TGM conditions.  If so, insert to the list, and use n_boudnary id of 3 to clarify that no flash was used.
			max_offset -= 0.01*units::cm;
			int nboundaries = check_boundary(extreme_points,offset_x-max_offset, &tgm_nof_tol_vec);
			if(nboundaries==2){
				n_boundary_list.push_back(3);
				candidate_offset_x_list.push_back(-offset_x+max_offset);
			}
		}

	}

	//Iterate the list and pick the highest priority tagging (TGM > TGM-noflash > STM > no tag)
	int tag_type = 11;
	for(int n_match = 0; n_match<int(n_boundary_list.size());n_match++){
		if(n_boundary_list[n_match]==2)							{tag_type = 22;}
		else if(n_boundary_list[n_match]==3 && tag_type != 22)				{tag_type = 24;}
		else if(n_boundary_list[n_match]==1 && tag_type != 22 && tag_type != 24)	{tag_type = 23;}

		//print out results for debugging purposes
		if(debug_tagger){
			if(n_match < int(candidate_flash_list.size())){
				std::cout << "match type = " << n_boundary_list[n_match] << ", flash # = " << candidate_flash_list[n_match]->get_flash_id() << std::endl;		
			} else {
				std::cout << "match type = " << n_boundary_list[n_match] << std::endl;
			}
		}
	}

	//Write to file for debugging purposes
	if(debug_tagger){
		srand(time(NULL));
		int random_number = std::rand();
		std::ofstream debug_file;
		std::string path = "/uboone/data/users/lcoopert/cosmic_tagger/data";
//		debug_file.open(path+"/temp_numu_cc_all/numu_cc_tagger_results_"+std::to_string(random_number)+".txt", std::ios_base::app);
//		debug_file.open(path+"/temp_extbnb_all/extbnb_tagger_results_"+std::to_string(random_number)+".txt", std::ios_base::app);

//		debug_file.open(path+"/temp_extbnb_test/extbnb_tagger_results_stm_"    +std::to_string(random_number)+".txt", std::ios_base::app);
//		debug_file.open(path+"/temp_extbnb_test/extbnb_tagger_results_tgm_"    +std::to_string(random_number)+".txt", std::ios_base::app);
//		debug_file.open(path+"/temp_extbnb_test/extbnb_tagger_results_tgm_nof_"+std::to_string(random_number)+".txt", std::ios_base::app);

		debug_file << run_no << " " << subrun_no << " " << event_no << " " << tag_type << std::endl;
		debug_file.close();
		//std::cout << "run subrun event: " << run_no << " " << subrun_no << " " << event_no << ", tag_type = " << tag_type << std::endl;
	}
	//std::cout << "ending cosmic tagger =====================================================" << std::endl;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
   

/*
WCP2dToy::ToyFiducial::ToyFiducial(int dead_region_ch_ext, double offset_t, double offset_u, double offset_v, double offset_w, double slope_t, double slope_u, double slope_v, double slope_w, double angle_u, double angle_v, double angle_w, double boundary_dis_cut, double top, double bottom, double upstream, double downstream, double anode, double cathode, int flag_data)
  : dead_region_ch_ext(dead_region_ch_ext)
  , offset_t(offset_t)
  , offset_u(offset_u)
  , offset_v(offset_v)
  , offset_w(offset_w)
  , slope_t(slope_t)
  , slope_u(slope_u)
  , slope_v(slope_v)
  , slope_w(slope_w)
  , angle_u(angle_u)
  , angle_v(angle_v)
  , angle_w(angle_w)
  , m_top(top)
  , m_bottom(bottom)
  , m_upstream(upstream)
  , m_downstream(downstream)
  , m_anode(anode)
  , m_cathode(cathode)
{

  if (flag_data){
    std::cout << "Data Reco Fiducial Volume! " << std::endl;
    // data 
    m_sc_bottom_1_y=-116*units::cm;
    m_sc_bottom_1_x=80*units::cm;
    
    m_sc_bottom_2_y=-99*units::cm;
    m_sc_bottom_2_x=256*units::cm;
    
    m_sc_top_1_y = 116*units::cm; // used to be 118 cm
    m_sc_top_1_x = 100*units::cm;
    
    m_sc_top_2_y = 102*units::cm; // used to be 103 cm
    m_sc_top_2_x = 256*units::cm;
    
    m_sc_upstream_1_z = 0*units::cm;
    m_sc_upstream_1_x = 120*units::cm;
    
    m_sc_upstream_2_z = 11*units::cm;
    m_sc_upstream_2_x = 256*units::cm;
    
    m_sc_downstream_1_z=1037*units::cm;
    m_sc_downstream_1_x=120*units::cm;
    
    m_sc_downstream_2_z=1026*units::cm;
    m_sc_downstream_2_x=256*units::cm;
  }else{
    // MC
    std::cout << "MC Truth Fiducial Volume! " << std::endl;
    m_sc_bottom_1_y=-116*units::cm;
    m_sc_bottom_1_x=34*units::cm;
    
    m_sc_bottom_2_y=-98*units::cm;
    m_sc_bottom_2_x=256*units::cm;
    
    m_sc_top_1_y = 116*units::cm;
    m_sc_top_1_x = 70*units::cm;
    
    m_sc_top_2_y = 100*units::cm;
    m_sc_top_2_x = 256*units::cm;
    
    m_sc_upstream_1_z = 0*units::cm;
    m_sc_upstream_1_x = 50*units::cm;
    
    m_sc_upstream_2_z = 14*units::cm;
    m_sc_upstream_2_x = 256*units::cm;
    
    m_sc_downstream_1_z=1037*units::cm;
    m_sc_downstream_1_x=40*units::cm;
    
    m_sc_downstream_2_z=1023*units::cm;
    m_sc_downstream_2_x=256*units::cm;
  }
  

  //
  boundary_xy_x.clear(); boundary_xy_y.clear();
  boundary_xy_x.push_back(m_anode + boundary_dis_cut); boundary_xy_y.push_back(m_bottom + boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_bottom_1_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_bottom_1_y+boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_bottom_2_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_bottom_2_y+boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_top_2_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_top_2_y-boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_top_1_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_top_1_y-boundary_dis_cut);
  boundary_xy_x.push_back(m_anode + boundary_dis_cut); boundary_xy_y.push_back(m_top - boundary_dis_cut);
  // boundary_xy_x.push_back(m_anode + boundary_dis_cut); boundary_xy_y.push_back(m_bottom + boundary_dis_cut);
  
  // for (size_t i=0;i!=boundary_xy_x.size();i++){
  //   std::cout << boundary_xy_x.at(i)/units::cm << " XY " << boundary_xy_y.at(i)/units::cm << std::endl;
  // }
  
  boundary_xz_x.clear(); boundary_xz_z.clear();
  boundary_xz_x.push_back(m_anode + boundary_dis_cut); boundary_xz_z.push_back(m_upstream + boundary_dis_cut+1*units::cm);
  boundary_xz_x.push_back(m_sc_upstream_1_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_upstream_1_z + boundary_dis_cut+1*units::cm);
  boundary_xz_x.push_back(m_sc_upstream_2_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_upstream_2_z + boundary_dis_cut+1*units::cm);
  boundary_xz_x.push_back(m_sc_downstream_2_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_downstream_2_z - boundary_dis_cut-1*units::cm);
  boundary_xz_x.push_back(m_sc_downstream_1_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_downstream_1_z - boundary_dis_cut-1*units::cm);
  boundary_xz_x.push_back(m_anode + boundary_dis_cut); boundary_xz_z.push_back(m_downstream - boundary_dis_cut-1*units::cm);
  // boundary_xz_x.push_back(m_anode + boundary_dis_cut); boundary_xz_z.push_back(m_upstream + boundary_dis_cut+2*units::cm);

  // for (size_t i=0;i!=boundary_xz_x.size();i++){
  //   std::cout << boundary_xz_x.at(i)/units::cm << " XZ " << boundary_xz_z.at(i)/units::cm << std::endl;
  // }
}
*/

WCP2dToy::ToyFiducial::ToyFiducial(int dead_region_ch_ext, double offset_t, double offset_u, double offset_v, double offset_w, double slope_t, double slope_u, double slope_v, double slope_w, double angle_u, double angle_v, double angle_w, double boundary_dis_cut, double top, double bottom, double upstream, double downstream, double anode, double cathode, int flag_data)
  : dead_region_ch_ext(dead_region_ch_ext)
  , offset_t(offset_t)
  , offset_u(offset_u)
  , offset_v(offset_v)
  , offset_w(offset_w)
  , slope_t(slope_t)
  , slope_u(slope_u)
  , slope_v(slope_v)
  , slope_w(slope_w)
  , angle_u(angle_u)
  , angle_v(angle_v)
  , angle_w(angle_w)
  , m_top(top)
  , m_bottom(bottom)
  , m_upstream(upstream)
  , m_downstream(downstream)
  , m_anode(anode)
  , m_cathode(cathode)
{

  //old -----------------------------
  if (flag_data){
    std::cout << "Data Fiducial Volume! " << std::endl;
    // data 
    m_sc_bottom_1_y=-116*units::cm;
    m_sc_bottom_1_x=80*units::cm;
    
    m_sc_bottom_2_y=-99*units::cm;
    m_sc_bottom_2_x=256*units::cm;
    
    m_sc_top_1_y = 116*units::cm; // used to be 118 cm
    m_sc_top_1_x = 100*units::cm;
    
    m_sc_top_2_y = 102*units::cm; // used to be 103 cm
    m_sc_top_2_x = 256*units::cm;
    
    m_sc_upstream_1_z = 0*units::cm;
    m_sc_upstream_1_x = 120*units::cm;
    
    m_sc_upstream_2_z = 11*units::cm;
    m_sc_upstream_2_x = 256*units::cm;
    
    m_sc_downstream_1_z=1037*units::cm;
    m_sc_downstream_1_x=120*units::cm;
    
    m_sc_downstream_2_z=1026*units::cm;
    m_sc_downstream_2_x=256*units::cm;
  }else{
    // MC
    std::cout << "MC Fiducial Volume! " << std::endl;
    m_sc_bottom_1_y=-116*units::cm;
    m_sc_bottom_1_x=34*units::cm;
    
    m_sc_bottom_2_y=-98*units::cm;
    m_sc_bottom_2_x=256*units::cm;
    
    m_sc_top_1_y = 116*units::cm;
    m_sc_top_1_x = 70*units::cm;
    
    m_sc_top_2_y = 100*units::cm;
    m_sc_top_2_x = 256*units::cm;
    
    m_sc_upstream_1_z = 0*units::cm;
    m_sc_upstream_1_x = 50*units::cm;
    
    m_sc_upstream_2_z = 14*units::cm;
    m_sc_upstream_2_x = 256*units::cm;
    
    m_sc_downstream_1_z=1037*units::cm;
    m_sc_downstream_1_x=40*units::cm;
    
    m_sc_downstream_2_z=1023*units::cm;
    m_sc_downstream_2_x=256*units::cm;
  }
  
  //
  boundary_xy_x.clear(); boundary_xy_y.clear();
  boundary_xy_x.push_back(m_anode + boundary_dis_cut); boundary_xy_y.push_back(m_bottom + boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_bottom_1_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_bottom_1_y+boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_bottom_2_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_bottom_2_y+boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_top_2_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_top_2_y-boundary_dis_cut);
  boundary_xy_x.push_back(m_sc_top_1_x-boundary_dis_cut); boundary_xy_y.push_back(m_sc_top_1_y-boundary_dis_cut);
  boundary_xy_x.push_back(m_anode + boundary_dis_cut); boundary_xy_y.push_back(m_top - boundary_dis_cut);
  
  // for (size_t i=0;i!=boundary_xy_x.size();i++){
  //   std::cout << boundary_xy_x.at(i)/units::cm << " XY " << boundary_xy_y.at(i)/units::cm << std::endl;
  // }
  
  boundary_xz_x.clear(); boundary_xz_z.clear();
  boundary_xz_x.push_back(m_anode + boundary_dis_cut); boundary_xz_z.push_back(m_upstream + boundary_dis_cut+1*units::cm);
  boundary_xz_x.push_back(m_sc_upstream_1_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_upstream_1_z + boundary_dis_cut+1*units::cm);
  boundary_xz_x.push_back(m_sc_upstream_2_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_upstream_2_z + boundary_dis_cut+1*units::cm);
  boundary_xz_x.push_back(m_sc_downstream_2_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_downstream_2_z - boundary_dis_cut-1*units::cm);
  boundary_xz_x.push_back(m_sc_downstream_1_x - boundary_dis_cut); boundary_xz_z.push_back(m_sc_downstream_1_z - boundary_dis_cut-1*units::cm);
  boundary_xz_x.push_back(m_anode + boundary_dis_cut); boundary_xz_z.push_back(m_downstream - boundary_dis_cut-1*units::cm);
  // boundary_xz_x.push_back(m_anode + boundary_dis_cut); boundary_xz_z.push_back(m_upstream + boundary_dis_cut+2*units::cm);

  boundary_SCB_xy_x.clear(); boundary_SCB_xy_y.clear();
  boundary_SCB_xz_x.clear(); boundary_SCB_xz_z.clear();

  //define boundary_SCB_xy
  boundary_SCB_xy_x.push_back(m_anode);				boundary_SCB_xy_y.push_back(m_bottom);
  boundary_SCB_xy_x.push_back(m_sc_bottom_1_x);			boundary_SCB_xy_y.push_back(m_sc_bottom_1_y);
  boundary_SCB_xy_x.push_back(m_sc_bottom_2_x);			boundary_SCB_xy_y.push_back(m_sc_bottom_2_y);
  boundary_SCB_xy_x.push_back(m_sc_top_2_x);			boundary_SCB_xy_y.push_back(m_sc_top_2_y);
  boundary_SCB_xy_x.push_back(m_sc_top_1_x);			boundary_SCB_xy_y.push_back(m_sc_top_1_y);
  boundary_SCB_xy_x.push_back(m_anode);				boundary_SCB_xy_y.push_back(m_top);

//why the extra 1cm on z???
  //define boundary_SCB_xz
  boundary_SCB_xz_x.push_back(m_anode);				boundary_SCB_xz_z.push_back(m_upstream+1*units::cm);
  boundary_SCB_xz_x.push_back(m_sc_upstream_1_x);		boundary_SCB_xz_z.push_back(m_sc_upstream_1_z+1*units::cm);
  boundary_SCB_xz_x.push_back(m_sc_upstream_2_x);		boundary_SCB_xz_z.push_back(m_sc_upstream_2_z+1*units::cm);
  boundary_SCB_xz_x.push_back(m_sc_downstream_2_x);		boundary_SCB_xz_z.push_back(m_sc_downstream_2_z-1*units::cm);
  boundary_SCB_xz_x.push_back(m_sc_downstream_1_x);		boundary_SCB_xz_z.push_back(m_sc_downstream_1_z-1*units::cm);
  boundary_SCB_xz_x.push_back(m_anode);				boundary_SCB_xz_z.push_back(m_downstream-1*units::cm);

  // for (size_t i=0;i!=boundary_xz_x.size();i++){
  //   std::cout << boundary_xz_x.at(i)/units::cm << " XZ " << boundary_xz_z.at(i)/units::cm << std::endl;
  // }

//new -----------------------------------------

  double YX_TOP_y1_array     = 116*units::cm;
  double YX_TOP_x1_array[10] = {150.00, 132.56, 122.86, 119.46, 114.22, 110.90, 115.85, 113.48, 126.36, 144.21};
  double YX_TOP_y2_array[10] = {110.00, 108.14, 106.77, 105.30, 103.40, 102.18, 101.76, 102.27, 102.75, 105.10};
//  double YX_TOP_x1_array[10] = {100.00, 100.00, 100.00, 100.00, 100.00, 100.00, 100.00, 100.00, 100.00, 100.00};
//  double YX_TOP_y2_array[10] = {102.00, 102.00, 102.00, 102.00, 102.00, 102.00, 102.00, 102.00, 102.00, 102.00};
  double YX_TOP_x2_array = 256*units::cm;
    
  double YX_BOT_y1_array     = -115*units::cm;
  double YX_BOT_x1_array[10] = {115.71, 98.05, 92.42, 91.14, 92.25, 85.38, 78.19, 74.46, 78.86, 108.90};
  double YX_BOT_y2_array[10] = {-101.72, -99.46, -99.51, -100.43, -99.55, -98.56, -98.00, -98.30, -99.32, -104.20};
//  double YX_BOT_x1_array[10] = {80.00,  80.00,  80.00,  80.00,  80.00,  80.00,  80.00,  80.00,  80.00,  80.00};
//  double YX_BOT_y2_array[10] = {-99.00, -99.00, -99.00, -99.00, -99.00, -99.00, -99.00, -99.00, -99.00, -99.00};
  double YX_BOT_x2_array = 256*units::cm;

  double ZX_Up_z1_array = 0*units::cm;
  double ZX_Up_x1_array = 120*units::cm;
  double ZX_Up_z2_array = 11*units::cm;
  double ZX_Up_x2_array = 256*units::cm;
    
  double ZX_Dw_z1_array     = 1037*units::cm;
  double ZX_Dw_x1_array[10] = {120.00, 115.24, 108.50, 110.67, 120.90, 126.43, 140.51, 157.15, 120.00, 120.00};
  double ZX_Dw_z2_array[10] = {1029.00, 1029.12, 1027.21, 1026.01, 1024.91, 1025.27, 1025.32, 1027.61, 1026.00, 1026.00};
//  double ZX_Dw_x1_array[10] = {120.00,  120.00,  120.00,  120.00,  120.00,  120.00,  120.00,  120.00,  120.00,  120.00};
//  double ZX_Dw_z2_array[10] = {1026.00, 1026.00, 1026.00, 1026.00, 1026.00, 1026.00, 1026.00, 1026.00, 1026.00, 1026.00};
  double ZX_Dw_x2_array     = 256*units::cm;

  for(int idx=0; idx<=9; idx++) {
    YX_BOT_x1_array[idx] *= units::cm;
    YX_BOT_y2_array[idx] *= units::cm;

    YX_TOP_x1_array[idx] *= units::cm;
    YX_TOP_y2_array[idx] *= units::cm;

    ZX_Dw_x1_array[idx] *= units::cm;
    ZX_Dw_z2_array[idx] *= units::cm;
  }
    boundary_xy_x_array.clear();						boundary_xy_y_array.clear();
    boundary_xz_x_array.clear();						boundary_xz_z_array.clear();
    boundary_SCB_xy_x_array.clear();						boundary_SCB_xy_y_array.clear();
    boundary_SCB_xz_x_array.clear();						boundary_SCB_xz_z_array.clear();

  for(int idx=0; idx<=9; idx++) {
    boundary_xy_x_array.push_back({0,0,0,0,0,0});				boundary_xy_y_array.push_back({0,0,0,0,0,0});
    boundary_xz_x_array.push_back({0,0,0,0,0,0});				boundary_xz_z_array.push_back({0,0,0,0,0,0});
    boundary_SCB_xy_x_array.push_back({0,0,0,0,0,0});				boundary_SCB_xy_y_array.push_back({0,0,0,0,0,0});
    boundary_SCB_xz_x_array.push_back({0,0,0,0,0,0});				boundary_SCB_xz_z_array.push_back({0,0,0,0,0,0});

    boundary_xy_x_array[idx][0] = m_anode              + boundary_dis_cut;	boundary_xy_y_array[idx][0] = m_bottom             + boundary_dis_cut;
    boundary_xy_x_array[idx][1] = YX_BOT_x1_array[idx] - boundary_dis_cut;	boundary_xy_y_array[idx][1] = YX_BOT_y1_array      + boundary_dis_cut;
    boundary_xy_x_array[idx][2] = YX_BOT_x2_array      - boundary_dis_cut;	boundary_xy_y_array[idx][2] = YX_BOT_y2_array[idx] + boundary_dis_cut;
    boundary_xy_x_array[idx][3] = YX_TOP_x2_array      - boundary_dis_cut;	boundary_xy_y_array[idx][3] = YX_TOP_y2_array[idx] - boundary_dis_cut;
    boundary_xy_x_array[idx][4] = YX_TOP_x1_array[idx] - boundary_dis_cut;	boundary_xy_y_array[idx][4] = YX_TOP_y1_array      - boundary_dis_cut;
    boundary_xy_x_array[idx][5] = m_anode              + boundary_dis_cut;	boundary_xy_y_array[idx][5] = m_top                - boundary_dis_cut;

    boundary_xz_x_array[idx][0] = m_anode             + boundary_dis_cut;	boundary_xz_z_array[idx][0] = m_upstream          + boundary_dis_cut+1*units::cm;
    boundary_xz_x_array[idx][1] = ZX_Up_x1_array      - boundary_dis_cut;	boundary_xz_z_array[idx][1] = ZX_Up_z1_array      + boundary_dis_cut+1*units::cm;
    boundary_xz_x_array[idx][2] = ZX_Up_x2_array      - boundary_dis_cut;	boundary_xz_z_array[idx][2] = ZX_Up_z2_array      + boundary_dis_cut+1*units::cm;
    boundary_xz_x_array[idx][3] = ZX_Dw_x2_array      - boundary_dis_cut;	boundary_xz_z_array[idx][3] = ZX_Dw_z2_array[idx] - boundary_dis_cut-1*units::cm;
    boundary_xz_x_array[idx][4] = ZX_Dw_x1_array[idx] - boundary_dis_cut;	boundary_xz_z_array[idx][4] = ZX_Dw_z1_array      - boundary_dis_cut-1*units::cm;
    boundary_xz_x_array[idx][5] = m_anode             + boundary_dis_cut;	boundary_xz_z_array[idx][5] = m_downstream        - boundary_dis_cut-1*units::cm;

    boundary_SCB_xy_x_array[idx][0] = m_anode             ;			boundary_SCB_xy_y_array[idx][0] = m_bottom            ;
    boundary_SCB_xy_x_array[idx][1] = YX_BOT_x1_array[idx];			boundary_SCB_xy_y_array[idx][1] = YX_BOT_y1_array     ;
    boundary_SCB_xy_x_array[idx][2] = YX_BOT_x2_array     ;			boundary_SCB_xy_y_array[idx][2] = YX_BOT_y2_array[idx];
    boundary_SCB_xy_x_array[idx][3] = YX_TOP_x2_array     ;			boundary_SCB_xy_y_array[idx][3] = YX_TOP_y2_array[idx];
    boundary_SCB_xy_x_array[idx][4] = YX_TOP_x1_array[idx];			boundary_SCB_xy_y_array[idx][4] = YX_TOP_y1_array     ;
    boundary_SCB_xy_x_array[idx][5] = m_anode             ;			boundary_SCB_xy_y_array[idx][5] = m_top               ;

    boundary_SCB_xz_x_array[idx][0] = m_anode            ;			boundary_SCB_xz_z_array[idx][0] = m_upstream          +1*units::cm;
    boundary_SCB_xz_x_array[idx][1] = ZX_Up_x1_array     ;			boundary_SCB_xz_z_array[idx][1] = ZX_Up_z1_array      +1*units::cm;
    boundary_SCB_xz_x_array[idx][2] = ZX_Up_x2_array     ;			boundary_SCB_xz_z_array[idx][2] = ZX_Up_z2_array      +1*units::cm;
    boundary_SCB_xz_x_array[idx][3] = ZX_Dw_x2_array     ;			boundary_SCB_xz_z_array[idx][3] = ZX_Dw_z2_array[idx] -1*units::cm;
    boundary_SCB_xz_x_array[idx][4] = ZX_Dw_x1_array[idx];			boundary_SCB_xz_z_array[idx][4] = ZX_Dw_z1_array      -1*units::cm;
    boundary_SCB_xz_x_array[idx][5] = m_anode            ;			boundary_SCB_xz_z_array[idx][5] = m_downstream        -1*units::cm;
    
  }  
}

int WCP2dToy::ToyFiducial::check_LM(WCP::FlashTPCBundle *bundle, double& cluster_length){

  PR3DCluster *main_cluster = bundle->get_main_cluster();
  Opflash *flash = bundle->get_flash();

  // calculate the length ...
  
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  std::vector<int> range_v1 = main_cluster->get_uvwt_range();
  cluster_length = sqrt(2./3. * (pow(pitch_u*range_v1.at(0),2) + pow(pitch_v*range_v1.at(1),2) + pow(pitch_w*range_v1.at(2),2)) + pow(time_slice_width*range_v1.at(3),2))/units::cm;

  double total_pred_pe = 0;
  std::vector<double>& pred_pe = bundle->get_pred_pmt_light();
  for (size_t i=0;i!=pred_pe.size();i++){
    total_pred_pe += pred_pe.at(i);
  }
  double total_flash_pe = flash->get_total_PE();

  //std::cout << main_cluster->get_cluster_id() << " " << cluster_length << std::endl;
  
  if (total_pred_pe < 25 || cluster_length < 10)
    return 1;

  bool flag_anode = bundle->get_flag_close_to_PMT();
  bool flag_boundary = bundle->get_flag_at_x_boundary();

  double ks_dis = bundle->get_ks_dis();
  double chi2 = bundle->get_chi2();
  double ndf = bundle->get_ndf();

  double meas_pe[32];
  double max_meas_pe = 0;
  for (int i=0;i!=32;i++){
    meas_pe[i] = flash->get_PE(i);
    if (max_meas_pe < meas_pe[i]) max_meas_pe = meas_pe[i];
  }

  if (flash->get_type()==2){
    if (!flag_boundary){
      if (!( log10(total_pred_pe/total_flash_pe)>-0.55 &&
	     ks_dis<0.25 &&
	     ks_dis-0.15/1.4*log10(total_pred_pe/total_flash_pe)<0.32) )
	return 2;
      
    }else{
      if (!(ks_dis<0.6 && (log10(total_pred_pe/total_flash_pe)>-1.4 && flag_anode || log10(total_pred_pe/total_flash_pe)>-0.55) &&
	    log10(total_flash_pe)+1.8*log10(total_pred_pe/total_flash_pe)>1.25 &&
	    max_meas_pe/total_flash_pe+0.17/0.5*log10(total_pred_pe/total_flash_pe)>-0.05))
	return 2;
    }
  }
  
  
  return 0;
}

int WCP2dToy::ToyFiducial::check_LM_cuts(WCP::FlashTPCBundle *bundle, double& cluster_length){

  Opflash *flash = bundle->get_flash();
  /* comment out for now... to make direct comparison to Xin's*/
  /*
  if(flash->get_type()!=2){
    std::cout<<"NOT beam discriminator flash"<<std::endl;
    return 0;
  }
  */
  
  /* calculate the length */
  PR3DCluster *main_cluster = bundle->get_main_cluster();  
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  std::vector<int> range_v1 = main_cluster->get_uvwt_range();
  cluster_length = sqrt(2./3. * (pow(pitch_u*range_v1.at(0),2) + pow(pitch_v*range_v1.at(1),2) + pow(pitch_w*range_v1.at(2),2)) + pow(time_slice_width*range_v1.at(3),2))/units::cm;

  /* predicted PE */
  double total_pred_pe = 0;
  std::vector<double>& pred_pe = bundle->get_pred_pmt_light();
  for (size_t i=0;i!=pred_pe.size();i++){
    total_pred_pe += pred_pe.at(i);
  }

  if(total_pred_pe < 25 || cluster_length < 10){
    return 1; /* low energy event */
  }

  /* temporary... to make direct comparison to Xin's*/
  if(flash->get_type()!=2) return 0;
  
  double total_meas_pe = flash->get_total_PE();
  bool flag_anode = bundle->get_flag_close_to_PMT();
  bool flag_boundary = bundle->get_flag_at_x_boundary();
  double ks_dis = bundle->get_ks_dis();

  if(flag_boundary){ /* at anode or cathode */
    if(flag_anode){ /* at anode */
      if( !(log(total_pred_pe/total_meas_pe)>-1.8 && ks_dis<0.8) ){
	return 2; /* light mismatch */
      }
    }
    else{ /* at cathode */
      if( !(log(total_pred_pe/total_meas_pe)>-1.8 && ks_dis<0.45) ){
	return 2;
      }
    }
  }
  else{ /* NOT at anode or cathode */
    if( !(log(total_pred_pe/total_meas_pe)>-1.4 && ks_dis<0.25) ){
      return 2;
    }
  }

  return 0;
}

int WCP2dToy::ToyFiducial::check_LM_bdt(WCP::FlashTPCBundle *bundle, double& cluster_length){

  Opflash *flash = bundle->get_flash();
  /* comment out for now... to make direct comparison to Xin's*/
  /*
  if(flash->get_type()!=2){
    std::cout<<"NOT beam discriminator flash"<<std::endl;
    return 0;
  }
  */
  
  /* calculate the length */
  PR3DCluster *main_cluster = bundle->get_main_cluster();  
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  std::vector<int> range_v1 = main_cluster->get_uvwt_range();
  cluster_length = sqrt(2./3. * (pow(pitch_u*range_v1.at(0),2) + pow(pitch_v*range_v1.at(1),2) + pow(pitch_w*range_v1.at(2),2)) + pow(time_slice_width*range_v1.at(3),2))/units::cm;

  /* predicted PE */
  float total_pred_pe = 0;
  std::vector<double>& pred_pe = bundle->get_pred_pmt_light();
  for (size_t i=0;i!=pred_pe.size();i++){
    total_pred_pe += (float)pred_pe.at(i);
  }

  if(total_pred_pe < 25 || cluster_length < 10){
    return 1; /* low energy event */
  }

  /* temporary... to make direct comparison to Xin's*/
  if(flash->get_type()!=2) return 0;
    
  float total_meas_pe = (float)flash->get_total_PE();
  int flag_anode = (int)bundle->get_flag_close_to_PMT();
  int flag_boundary = (int)bundle->get_flag_at_x_boundary();
  float ks_dis = (float)bundle->get_ks_dis();
  float chi2 = (float)bundle->get_chi2();
  float ndf = (float)bundle->get_ndf();
  float cl = (float)cluster_length;
  
  float meas_pe[32];
  float max_meas_pe = 0;
  for (int i=0;i!=32;i++){
    meas_pe[i] = (float)flash->get_PE(i);
    if (max_meas_pe < meas_pe[i]) max_meas_pe = meas_pe[i];
  }

  /* TEMPORARY*/
  //TFile f("lm_bdt.root","READ");
  //TH1D *sig = (TH1D*)f.Get("MVA_BDT_effS");
  //TH1D *bgd = (TH1D*)f.Get("MVA_BDT_effB");
  std::ifstream inParams("input_data_files/lmHistParams.txt");
  std::ifstream inSigFile("input_data_files/lmSigEff.txt");
  std::ifstream inBgdFile("input_data_files/lmBgdEff.txt");
  int bins = 0;
  double binL = 0., binH = 0.;
  for(int i=1; i<=3; i++){
    if(i==1){ inParams >> bins; }
    if(i==2){ inParams >> binL; }
    if(i==3){ inParams >> binH; }
  }
  
  TH1D *sig = new TH1D("sig","",bins,binL,binH);
  TH1D *bgd = new TH1D("bgd","",bins,binL,binH);
  double eff = 0.;
  for(int i=1; i<=10000; i++){
    inSigFile >> eff;
    sig->SetBinContent(i, eff);
    inBgdFile >> eff;
    bgd->SetBinContent(i, eff);
  }
  
  int temp = -100;
  WCP::LMBDT lm(total_pred_pe,total_meas_pe,max_meas_pe, ks_dis,
		     chi2,ndf,cl,temp,flag_anode,flag_boundary); 
  
  //  if( !lm.isSignal(lm.get_BDT_score_max_significance(sig,bgd)) ){
  if( !lm.isSignal(lm.get_BDT_score_eff_background(0.005,bgd)) ){
    //f.Close();
    return 2; /* light mismatch event */
  }

  return 0;
}

bool WCP2dToy::ToyFiducial::check_fully_contained(WCP::FlashTPCBundle *bundle, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud,std::map<PR3DCluster*, PR3DCluster*>& old_new_cluster_map, unsigned int* fail_mode, int flag){
  PR3DCluster *main_cluster; 
  PR3DCluster *main_cluster1; 
  if (flag==1){ // check the current
    main_cluster = bundle->get_main_cluster();
    main_cluster1  = main_cluster;
    //replace it with the better one, which takes into account the dead channels ... 
    if (old_new_cluster_map.find(main_cluster)!=old_new_cluster_map.end())
      main_cluster1 = old_new_cluster_map[main_cluster];
  }else{
    main_cluster = bundle->get_orig_cluster();
    main_cluster1  = main_cluster;
    //replace it with the better one, which takes into account the dead channels ... 
    if (old_new_cluster_map.find(main_cluster)!=old_new_cluster_map.end())
      main_cluster1 = old_new_cluster_map[main_cluster];
    if (main_cluster==0) return true;
  }
  
  
  Opflash *flash = bundle->get_flash();
  std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = main_cluster1->get_extreme_wcps();

  TVector3 drift_dir(1,0,0);
  // hard coded for U and V plane ... 
  TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
  TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
  TVector3 W_dir(0,1,0);

  Vector main_dir = main_cluster->get_PCA_axis(0);
  TVector3 dir_main(main_dir.x,main_dir.y,main_dir.z);
  
  for (size_t i=0;i!=out_vec_wcps.size();i++){
    // check all the points ... 
    for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
      Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
      if (!inside_fiducial_volume(p1,offset_x)){
        if(fail_mode) *fail_mode |= 1U<<2;
        return false;
      }
    }


    Point p1(out_vec_wcps.at(i).at(0).x,out_vec_wcps.at(i).at(0).y,out_vec_wcps.at(i).at(0).z);
    TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
    dir *= (-1);

    // check U and V and W
    TVector3 dir_1(0,dir.Y(),dir.Z());
    double angle1 = dir_1.Angle(U_dir);
    TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
    double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;
    
    double angle2 = dir_1.Angle(V_dir);
    TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
    double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;
    
    double angle3 = dir_1.Angle(W_dir);
    TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
    double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;
    
    // not added for now, need to check to add this one in when more events are available ...
    // XQ, 7/11/2018
    // double angle4 = fabs(3.1415926/2.-dir.Angle(drift_dir))/3.1415926*180.;


    //std::cout << "A: " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << angle1_1 << " " << angle2_1 << " " << angle3_1 << std::endl;

    if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
      if (!check_signal_processing(p1,dir,ct_point_cloud,1*units::cm,offset_x)){
        if(fail_mode) *fail_mode |= 1U<<1;
        return false;
      }
    }
	    
    if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 ){
      if (!check_dead_volume(p1,dir,1*units::cm,offset_x)){
        if(fail_mode) *fail_mode |= 1U;
        return false;
      }
    }
    
  }



  
  return true;
}

bool WCP2dToy::ToyFiducial::check_tgm(WCP::FlashTPCBundle *bundle, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud,std::map<PR3DCluster*, PR3DCluster*>& old_new_cluster_map, int flag){

  PR3DCluster *main_cluster; 
  PR3DCluster *main_cluster1; 

  if (flag==1){ // check the current main cluster 
    main_cluster = bundle->get_main_cluster();
    main_cluster1 = main_cluster;
  //replace it with the better one, which takes into account the dead channels ... 
    if (old_new_cluster_map.find(main_cluster)!=old_new_cluster_map.end())
      main_cluster1 = old_new_cluster_map[main_cluster];
  }else if (flag==2){
    main_cluster = bundle->get_orig_cluster();
    main_cluster1 = main_cluster;
  //replace it with the better one, which takes into account the dead channels ... 
    if (old_new_cluster_map.find(main_cluster)!=old_new_cluster_map.end())
      main_cluster1 = old_new_cluster_map[main_cluster];
    if (main_cluster==0) return false;
  }

  
  Opflash *flash = bundle->get_flash();

  std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = main_cluster1->get_extreme_wcps();

  // if (main_cluster->get_cluster_id()==6)
  //std::cout << main_cluster << " " << main_cluster1 << std::endl;
  

  // int max_group = 0;
  // int max_count = out_vec_wcps.at(max_group).size();

  // for (size_t i=1; i!=out_vec_wcps.size();i++){
  //   if (out_vec_wcps.at(i).size() > max_count){
  //     max_group = i;
  //     max_count = out_vec_wcps.at(max_group).size();
  //   }
  // }

  TVector3 drift_dir(1,0,0);
  // hard coded for U and V plane ... 
  TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
  TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
  TVector3 W_dir(0,1,0);

  double length_limit = sqrt(pow(out_vec_wcps.at(0).at(0).x-out_vec_wcps.at(1).at(0).x,2)+
			     pow(out_vec_wcps.at(0).at(0).y-out_vec_wcps.at(1).at(0).y,2)+
			     pow(out_vec_wcps.at(0).at(0).z-out_vec_wcps.at(1).at(0).z,2));
  
  // std::cout << "Flash: " << flash->get_flash_id() << std::endl;

  // take a look at the first point ...
  for (size_t i=0;i!=out_vec_wcps.size();i++){
    bool flag_p1_inside = true;
    int p1_index = -1;
    for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
      Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
      flag_p1_inside = flag_p1_inside && inside_fiducial_volume(p1,offset_x);
      if (!flag_p1_inside){
	p1_index = j;
	break;
      }
    }
    
    
    // loop through the remaining groups and check ...
    for (size_t k=i+1;k!=out_vec_wcps.size();k++){
      bool flag_p2_inside = true;
      int p2_index = -1;
      for(size_t j=0;j!=out_vec_wcps.at(k).size();j++){
	Point p2(out_vec_wcps.at(k).at(j).x,out_vec_wcps.at(k).at(j).y,out_vec_wcps.at(k).at(j).z);
	flag_p2_inside = flag_p2_inside && inside_fiducial_volume(p2,offset_x);
	if (!flag_p2_inside){
	  p2_index = j;
	  break;
	}
      }
      

      // if (main_cluster->get_cluster_id()==7){
      //  	std::cout << main_cluster->get_cluster_id() << " " << i << " " <<
      //  	  out_vec_wcps.at(i).at(0).x/units::cm << " " << out_vec_wcps.at(i).at(0).y/units::cm << " " << out_vec_wcps.at(i).at(0).z/units::cm << " " << 
      //  	  k << " " << out_vec_wcps.at(k).at(0).x/units::cm << " " << out_vec_wcps.at(k).at(0).y/units::cm << " " << out_vec_wcps.at(k).at(0).z/units::cm <<
      //  	  " " << p1_index << " " << p2_index << " " << flag_p1_inside << " " << flag_p2_inside << std::endl;
      // 	// for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
      // 	//   Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
      // 	//   std::cout << j << " A " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << std::endl;
      // 	// }
      // 	// for(size_t j=0;j!=out_vec_wcps.at(k).size();j++){
      // 	//   Point p2(out_vec_wcps.at(k).at(j).x,out_vec_wcps.at(k).at(j).y,out_vec_wcps.at(k).at(j).z);
      // 	//   std::cout << j << " B " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << std::endl;
      // 	// }
      // }
	
      if ((!flag_p1_inside) && (!flag_p2_inside)){
	// if not a neutrino candidate ... to be worked out ...
	 // Point p1(out_vec_wcps.at(i).at(p1_index).x,out_vec_wcps.at(i).at(p1_index).y,out_vec_wcps.at(i).at(p1_index).z);
	 // Point p2(out_vec_wcps.at(k).at(p2_index).x,out_vec_wcps.at(k).at(p2_index).y,out_vec_wcps.at(k).at(p2_index).z);
	//	std::cout << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << inside_fiducial_volume(p1,offset_x) << " A " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << inside_fiducial_volume(p2,offset_x) << " " << offset_x/units::cm << std::endl;


	// std::cout << main_cluster->get_cluster_id() << " " << (out_vec_wcps.at(i).at(p1_index).x-offset_x)/units::cm <<3=r-t=-[\] " " << out_vec_wcps.at(i).at(p1_index).y/units::cm << " " << out_vec_wcps.at(i).at(p1_index).z/units::cm << " " ;
	// std::cout << (out_vec_wcps.at(k).at(p2_index).x-offset_x)/units::cm << " " << out_vec_wcps.at(k).at(p2_index).y/units::cm << " " << out_vec_wcps.at(k).at(p2_index).z/units::cm << " " <<  flag_p1_inside << " " << flag_p2_inside << " " << out_vec_wcps.size() << std::endl;

	//std::cout << p1_index << " " << p2_index << std::endl;
	
	// check two points in between
	bool flag_check = false;
	for (int kk=0;kk!=3;kk++){
	  Point p3(out_vec_wcps.at(i).at(p1_index).x+ (out_vec_wcps.at(k).at(p2_index).x - out_vec_wcps.at(i).at(p1_index).x)/4.*(kk+1),
		   out_vec_wcps.at(i).at(p1_index).y+ (out_vec_wcps.at(k).at(p2_index).y - out_vec_wcps.at(i).at(p1_index).y)/4.*(kk+1),
		   out_vec_wcps.at(i).at(p1_index).z+ (out_vec_wcps.at(k).at(p2_index).z - out_vec_wcps.at(i).at(p1_index).z)/4.*(kk+1));
	  flag_check = flag_check || inside_fiducial_volume(p3,offset_x);
	}
	
	// if (main_cluster->get_cluster_id()==10)
	//   std::cout << flag_check << " " << out_vec_wcps.size() << std::endl;

	// middle points are inside fiducial volume
	if (flag_check){
	  if (flash->get_type()==2){
	    
	    double temp_length = sqrt(pow(out_vec_wcps.at(i).at(p1_index).x-out_vec_wcps.at(k).at(p2_index).x,2)+
				      pow(out_vec_wcps.at(i).at(p1_index).y-out_vec_wcps.at(k).at(p2_index).y,2)+
				      pow(out_vec_wcps.at(i).at(p1_index).z-out_vec_wcps.at(k).at(p2_index).z,2));

	    //std::cout << temp_length/units::cm << " " << length_limit/units::cm << std::endl;

	    if (i==0&&k==1){
	      if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud,true)) &&
		   temp_length > 0.45*length_limit)
		return true;
	    }else{
	      if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud)) &&
		   temp_length > 0.45*length_limit)
		return true;
	    }
	    
	  }else{
	    return true; // through going muon ...
	  }
	  // no middle points are inside the fiducial volume
	}else{
	  
	  if (out_vec_wcps.size()==2){
	    return true;
	  }else{
	    // if (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x))
	    //   return true;

	    bool flag_check_again = false;
	    for (int kkk = 0;kkk!=out_vec_wcps.size(); kkk++){
	      if (kkk==i || kkk==k) continue;
	      for (int kk=0;kk!=4;kk++){
	    	Point p3(out_vec_wcps.at(i).at(p1_index).x+ (out_vec_wcps.at(kkk).at(0).x - out_vec_wcps.at(i).at(p1_index).x)/4.*(kk+1),
	    		 out_vec_wcps.at(i).at(p1_index).y+ (out_vec_wcps.at(kkk).at(0).y - out_vec_wcps.at(i).at(p1_index).y)/4.*(kk+1),
	    		 out_vec_wcps.at(i).at(p1_index).z+ (out_vec_wcps.at(kkk).at(0).z - out_vec_wcps.at(i).at(p1_index).z)/4.*(kk+1));
	    	flag_check_again = flag_check_again || inside_fiducial_volume(p3,offset_x);
	      }
	      
	      for (int kk=0;kk!=3;kk++){
	    	Point p3(out_vec_wcps.at(kkk).at(0).x+ (out_vec_wcps.at(k).at(p2_index).x - out_vec_wcps.at(i).at(0).x)/4.*(kk+1),
	    		 out_vec_wcps.at(kkk).at(0).y+ (out_vec_wcps.at(k).at(p2_index).y - out_vec_wcps.at(i).at(0).y)/4.*(kk+1),
	    		 out_vec_wcps.at(kkk).at(0).z+ (out_vec_wcps.at(k).at(p2_index).z - out_vec_wcps.at(i).at(0).z)/4.*(kk+1));
	    	flag_check_again = flag_check_again || inside_fiducial_volume(p3,offset_x);
	      }
	    }
	    if (!flag_check_again){
	      //find the longest one ...
	      double temp_length = sqrt(pow(out_vec_wcps.at(i).at(p1_index).x-out_vec_wcps.at(k).at(p2_index).x,2)+
					pow(out_vec_wcps.at(i).at(p1_index).y-out_vec_wcps.at(k).at(p2_index).y,2)+
					pow(out_vec_wcps.at(i).at(p1_index).z-out_vec_wcps.at(k).at(p2_index).z,2));
	      // check neutrino candidates ...
	      if (i==0&&k==1){
		if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud,true)) && temp_length > 0.45*length_limit)
		  return true;
	      }else{
		if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud)) && temp_length > 0.45*length_limit)
		  return true;
	      }
	    }	    
	  }
	}
      }else{
	Vector main_dir = main_cluster->get_PCA_axis(0);
	TVector3 dir_main(main_dir.x,main_dir.y,main_dir.z);
	TVector3 dir_test(out_vec_wcps.at(i).at(0).x-out_vec_wcps.at(k).at(0).x,
			  out_vec_wcps.at(i).at(0).y-out_vec_wcps.at(k).at(0).y,
			  out_vec_wcps.at(i).at(0).z-out_vec_wcps.at(k).at(0).z);

	//	std::cout << main_cluster->get_cluster_id() << " " << fabs((3.1415926/2.-dir_test.Angle(dir_main))/3.1415926*180.) << std::endl;

	// std::cout << main_cluster->get_cluster_id() << " " << (out_vec_wcps.at(i).at(0).x-offset_x)/units::cm << " " << out_vec_wcps.at(i).at(0).y/units::cm << " " << out_vec_wcps.at(i).at(0).z/units::cm << " " ;
	// std::cout << (out_vec_wcps.at(k).at(0).x-offset_x)/units::cm << " " << out_vec_wcps.at(k).at(0).y/units::cm << " " << out_vec_wcps.at(k).at(0).z/units::cm << " " <<  flag_p1_inside << " " << flag_p2_inside << " " << out_vec_wcps.size() << std::endl;
	
	if (fabs((3.1415926/2.-dir_test.Angle(dir_main))/3.1415926*180.)>75 || i==0 && k==1){
	  // check dead region ...
	  bool flag_p1_inside_p = flag_p1_inside;
	  if (flag_p1_inside_p){
	    Point p1(out_vec_wcps.at(i).at(0).x,out_vec_wcps.at(i).at(0).y,out_vec_wcps.at(i).at(0).z);
	    TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
	    dir *= (-1);

	    if (dir.Angle(dir_test) > 3.1415926*2./3.) continue;
	    
	    //	    std::cout << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << dir.X() << " " << dir.Y() << dir.Z() << std::endl;

	    // check U and V and W
	    TVector3 dir_1(0,dir.Y(),dir.Z());
	    double angle1 = dir_1.Angle(U_dir);
	    TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
	    double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;

	    double angle2 = dir_1.Angle(V_dir);
	    TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
	    double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;

	    double angle3 = dir_1.Angle(W_dir);
	    TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
	    double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;

	    // not added for now, need to check to add this one in when more events are available ...
	    // XQ, 7/11/2018
	    double angle4 = fabs(3.1415926/2.-dir.Angle(drift_dir))/3.1415926*180.;


	    //std::cout << "A: " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << angle1_1 << " " << angle2_1 << " " << angle3_1 << std::endl;

	    if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
	      flag_p1_inside_p = flag_p1_inside_p && check_signal_processing(p1,dir,ct_point_cloud,1*units::cm,offset_x);
	    }
	    
	    if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 )
	      flag_p1_inside_p= flag_p1_inside_p && check_dead_volume(p1,dir,1*units::cm,offset_x);
	  }
	  
	  bool flag_p2_inside_p = flag_p2_inside;
	  if (flag_p2_inside_p){
	    Point p2(out_vec_wcps.at(k).at(0).x,out_vec_wcps.at(k).at(0).y,out_vec_wcps.at(k).at(0).z);
	    TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
	    dir *= (-1);

	    if (dir.Angle(dir_test) < 3.1415926/3.) continue;
	    
	    //	    std::cout << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << dir.X() << " " << dir.Y() << dir.Z() << std::endl;
	    TVector3 dir_1(0,dir.Y(),dir.Z());
	    double angle1 = dir_1.Angle(U_dir);
	    TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
	    double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;

	    double angle2 = dir_1.Angle(V_dir);
	    TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
	    double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;

	    double angle3 = dir_1.Angle(W_dir);
	    TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
	    double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;

	    //	    std::cout << "B: " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " <<  angle1_1 << " " << angle2_1 << " " << angle3_1 << std::endl;

	    if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
	      flag_p2_inside_p = flag_p2_inside_p && check_signal_processing(p2,dir,ct_point_cloud,1*units::cm,offset_x);
	    }
	    
	    if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 )
	      flag_p2_inside_p=flag_p2_inside_p && check_dead_volume(p2,dir,1*units::cm,offset_x);
	  }
	  
	  if ((!flag_p1_inside_p) && (!flag_p2_inside_p)){
	    if (flash->get_type()==2){
	      double temp_length = sqrt(pow(out_vec_wcps.at(i).at(0).x-out_vec_wcps.at(k).at(0).x,2)+
					pow(out_vec_wcps.at(i).at(0).y-out_vec_wcps.at(k).at(0).y,2)+
					pow(out_vec_wcps.at(i).at(0).z-out_vec_wcps.at(k).at(0).z,2));
	      if (i==0&&k==1){
		if ((!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(0),out_vec_wcps.at(k).at(0),offset_x,ct_point_cloud,true)) &&
		    temp_length > 0.45*length_limit
		    )
		  return true;
	      }else{
		if ((!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(0),out_vec_wcps.at(k).at(0),offset_x,ct_point_cloud)) &&
		    temp_length > 0.45*length_limit
		    )
		  return true;
	      }
	    }else{
	      return true;
	    }
	  }
	  // check signal processing ...

	// {
	// 	if (flag_p1_inside)
	// 	  ;
	
	// 	if (flag_p2_inside)
	// 	  ;
	
	// }

	}
      }
    }
  }
  

  return false;

  // also check against the dead channel ...  
  // // check the fiducial volume ...
  // std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = main_cluster->get_main_axis_wcps();
  // Point p1(wcps.first.x,wcps.first.y,wcps.first.z);
  // Point p2(wcps.second.x,wcps.second.y,wcps.second.z);
  // //double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
  // bool flag_inside_p1 = inside_fiducial_volume(p1,offset_x);
  // bool flag_inside_p2 = inside_fiducial_volume(p2,offset_x);
  // //std::cout << main_cluster->get_cluster_id() << " " << (p1.x-offset_x)/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << (p2.x-offset_x)/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << fid->inside_fiducial_volume(p1,offset_x) << " " << fid->inside_fiducial_volume(p2,offset_x) << std::endl;

  
  
  // // check the dead region ... 
  // if (flag_inside_p1){
  //   // define a local direction ...
  //   TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
  //   dir *= (-1);
  //   flag_inside_p1=check_dead_volume(p1,dir,1*units::cm,offset_x);
  // }
  // if (flag_inside_p2){
  //   // define a  local direction ...
  //   TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
  //   dir *= (-1);
  //   flag_inside_p2=check_dead_volume(p2,dir,1*units::cm,offset_x);
  // }
  // return (!flag_inside_p1)&&(!flag_inside_p2);


  // bool flag_2nd = true;
       // {
	 
	 
       // 	 if ((!flag_inside_p1)&&(!flag_inside_p2)){
       // 	   event_type |= 1UL << 3; // through going muon ... 
       // 	   flag_2nd = false;
       // 	 }
	 
       // 	 if (flag_2nd && ((!flag_inside_p1)|| (!flag_inside_p2) )){
       // 	   // check the fiducial volume ...
       // 	   std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = main_cluster->get_extreme_wcps();
	   
       // 	   Point p1(wcps.first.x,wcps.first.y,wcps.first.z);
       // 	   Point p2(wcps.second.x,wcps.second.y,wcps.second.z);
	   
       // 	   flag_inside_p1 = fid->inside_fiducial_volume(p1,offset_x);
       // 	   flag_inside_p2 = fid->inside_fiducial_volume(p2,offset_x);

       // 	   std::cout << main_cluster->get_cluster_id() << " " << (p1.x-offset_x)/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << (p2.x-offset_x)/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << fid->inside_fiducial_volume(p1,offset_x) << " " << fid->inside_fiducial_volume(p2,offset_x) << std::endl;
	   
       // 	   // check the dead region ...
       // 	   if (flag_inside_p1){
       // 	     // define a local direction ...
       // 	     TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
       // 	     dir *= (-1);
       // 	     flag_inside_p1=fid->check_dead_volume(p1,dir,1*units::cm,offset_x);
       // 	   }
       // 	   if (flag_inside_p2){
       // 	     // define a  local direction ...
       // 	     TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
       // 	     dir *= (-1);
       // 	     flag_inside_p2=fid->check_dead_volume(p2,dir,1*units::cm,offset_x);
       // 	   }
	   
       // 	   if ((!flag_inside_p1)&&(!flag_inside_p2)){
       // 	     event_type |= 1UL << 3; // through going muon ... 
       // 	   }
       // 	 }

  
  return false;
}

bool WCP2dToy::ToyFiducial::check_neutrino_candidate(WCP::PR3DCluster *main_cluster,WCPointCloud<double>::WCPoint& wcp1 ,WCPointCloud<double>::WCPoint& wcp2, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud,bool flag_2view_check){
  main_cluster->Create_graph(ct_point_cloud);
  
  main_cluster->dijkstra_shortest_paths(wcp1);
  main_cluster->cal_shortest_path(wcp2);

  std::list<WCPointCloud<double>::WCPoint>& path_wcps = main_cluster->get_path_wcps();
  
  
  
  PointVector path_wcps_vec;
  PointVector path_wcps_vec1;  
  double low_dis_limit = 0.5*units::cm;
  for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
    if (path_wcps_vec.size()==0){
      Point p((*it).x,(*it).y,(*it).z);
      path_wcps_vec.push_back(p);
      path_wcps_vec1.push_back(p);
    }else{
      double dis = sqrt(pow((*it).x - path_wcps_vec.back().x,2)
			+pow((*it).y - path_wcps_vec.back().y,2)
			+pow((*it).z - path_wcps_vec.back().z,2));
      if (dis > low_dis_limit){
	Point p((*it).x,(*it).y,(*it).z);
	path_wcps_vec.push_back(p);
      }

      dis = sqrt(pow((*it).x - path_wcps_vec1.back().x,2)
		 +pow((*it).y - path_wcps_vec1.back().y,2)
		 +pow((*it).z - path_wcps_vec1.back().z,2));
      if (dis <= 2*low_dis_limit){
	Point p((*it).x,(*it).y,(*it).z);
	path_wcps_vec1.push_back(p);
      }else{
	int nseg = dis/2./low_dis_limit+1;
	for (int i=0;i!=nseg;i++){
	  Point temp_p;
	  temp_p.x = path_wcps_vec1.back().x + (i+1.)*((*it).x - path_wcps_vec1.back().x)/nseg;
	  temp_p.y = path_wcps_vec1.back().y + (i+1.)*((*it).y - path_wcps_vec1.back().y)/nseg;
	  temp_p.z = path_wcps_vec1.back().z + (i+1.)*((*it).z - path_wcps_vec1.back().z)/nseg;
	  path_wcps_vec1.push_back(temp_p);
	}
      }
    }
  }  

  
  //check whether path is good ... 
  {
    int num_nth = 0;
    double min_dis = 1e9;

    // bool flag_2view_check = true;
    
    // U and V induction view checks
    if (flag_2view_check){
      TVector3 drift_dir(1,0,0);
      // hard coded for U and V plane ... 
      TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
      TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
      TVector3 W_dir(0,1,0);

      TVector3 dir(wcp2.x-wcp1.x,wcp2.y-wcp1.y,wcp2.z-wcp1.z);
      
      TVector3 dir_1(0,dir.Y(),dir.Z());
      double angle1 = dir_1.Angle(U_dir);
      TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
      double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;

      double angle2 = dir_1.Angle(V_dir);
      TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
      double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;
      
      double angle3 = dir_1.Angle(W_dir);
      TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
      double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;

      double angle4 = fabs(3.1415926/2.-drift_dir.Angle(dir))/3.1415926*180.;
      
      
      if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5 || angle4 < 5.)){
	flag_2view_check = false;
      }
      
    }

    //int num_total_dead = 0;
    
    int num_bad = 0;
    
    for (int i=0;i!=path_wcps_vec1.size();i++){
      WCP::CTPointCloud<double> cloud_u = ct_point_cloud.get_closest_points(path_wcps_vec1.at(i),low_dis_limit*2,0);
      WCP::CTPointCloud<double> cloud_v = ct_point_cloud.get_closest_points(path_wcps_vec1.at(i),low_dis_limit*2,1);
      WCP::CTPointCloud<double> cloud_w = ct_point_cloud.get_closest_points(path_wcps_vec1.at(i),low_dis_limit*2,2);

      bool flag_reset = false;

      if (flag_2view_check){
	if (cloud_u.pts.size() >0 && cloud_v.pts.size() > 0 || // require two planes to be good ...
	    cloud_u.pts.size() >0 && cloud_w.pts.size() > 0 ||
	    cloud_v.pts.size() >0 && cloud_w.pts.size() > 0){
	  flag_reset = true;
	}else{
	  if (inside_dead_region(path_wcps_vec1.at(i)))
	    flag_reset = true;
	}
      }else{
	if (cloud_u.pts.size() >0  || // require one plane to be good ...
	    cloud_v.pts.size() >0  ||
	    cloud_w.pts.size() >0 ){
	  flag_reset = true;
	}else{
	  if (inside_dead_region(path_wcps_vec1.at(i)))
	    flag_reset = true;
	}
      }

      if (!ct_point_cloud.is_good_point(path_wcps_vec1.at(i)))
	num_bad ++;
      
      // std::cout << "O: " << path_wcps_vec1.at(i).x/units::cm << " " 
      // 		<< path_wcps_vec1.at(i).y/units::cm << " "
      // 		<< path_wcps_vec1.at(i).z/units::cm << " " << flag_reset << " " << cloud_u.pts.size() << " " << cloud_v.pts.size() << " "<< cloud_w.pts.size() << std::endl;
      
      if (flag_reset){
  	num_nth =0;
	min_dis = 1e9;
	num_bad = 0;
      }else{
	if (inside_fiducial_volume(path_wcps_vec1.at(i),offset_x)){
	  double dis1 = sqrt(pow(path_wcps_vec1.at(i).x-wcp1.x,2)+pow(path_wcps_vec1.at(i).y-wcp1.y,2)+pow(path_wcps_vec1.at(i).z-wcp1.z,2));
	  double dis2 = sqrt(pow(path_wcps_vec1.at(i).x-wcp2.x,2)+pow(path_wcps_vec1.at(i).y-wcp2.y,2)+pow(path_wcps_vec1.at(i).z-wcp2.z,2));
	  if (dis1 < min_dis) min_dis = dis1;
	  if (dis2 < min_dis) min_dis = dis2;
	  num_nth ++;
	  // num_total_dead ++;
	  
	  // if (main_cluster->get_cluster_id()==7)
	  //   std::cout << main_cluster->get_cluster_id() << " " << min_dis/units::cm << " " << flag_2view_check << " " << path_wcps_vec1.at(i).x/units::cm << " " 
	  // 	      << path_wcps_vec1.at(i).y/units::cm << " "
	  // 	      << path_wcps_vec1.at(i).z/units::cm << " "
	  // 	      << num_nth << " " << cloud_u.pts.size() << " " << cloud_v.pts.size() << " "<< cloud_w.pts.size() << " " << num_bad << std::endl;
	  
	}
	//std::cout << num_nth << std::endl;

	if (num_nth > 7 && min_dis < 25*units::cm && num_bad > 7) return true; // too big a gap ... 4 cm cut ...
	//if (num_nth > 7 && num_bad > 7) return true; // too big a gap ... 4 cm cut ...
      }
    }
  }
  
  
  // if (cluster_id == 13){
  //   std::cout << wcp1.x/units::cm << " " << wcp1.y/units::cm << " " << wcp1.z/units::cm << " " << wcp2.x/units::cm << " " << wcp2.y/units::cm << " " << wcp2.z/units::cm << std::endl;
  int count = 0;
  double max_angle=0;
  Point max_point(0,0,0);
  TVector3 drift_dir(1,0,0);
  for (size_t i=5;i+5<path_wcps_vec.size();i++){
    TVector3 dir1(path_wcps_vec.at(i).x - path_wcps_vec.at(i-5).x,
		  path_wcps_vec.at(i).y - path_wcps_vec.at(i-5).y,
		  path_wcps_vec.at(i).z - path_wcps_vec.at(i-5).z);
    TVector3 dir2(path_wcps_vec.at(i).x - path_wcps_vec.at(i+5).x,
		  path_wcps_vec.at(i).y - path_wcps_vec.at(i+5).y,
		  path_wcps_vec.at(i).z - path_wcps_vec.at(i+5).z);
    
    TVector3 dir3, dir4, dir5, dir6;
    {
      PointVector pts;
      double temp_x = 0;
      double temp_y = 0;
      double temp_z = 0;
      double temp_count = 0;
      for (size_t j=1;j!=15;j++){
       	if (i>=j){
	  Point pt(path_wcps_vec.at(i-j).x,path_wcps_vec.at(i-j).y,path_wcps_vec.at(i-j).z);
	  
	  if (j<=12&&j>2){
	    temp_x += pt.x;
	    temp_y += pt.y;
	    temp_z += pt.z;
	    temp_count ++;
	  }
	  pts.push_back(pt);
	}
	Point pt(path_wcps_vec.at(i).x,path_wcps_vec.at(i).y,path_wcps_vec.at(i).z);
	dir3 = main_cluster->calc_PCA_dir(pt,pts);
	dir5.SetXYZ(temp_x/temp_count - path_wcps_vec.at(i).x,
		    temp_y/temp_count - path_wcps_vec.at(i).y,
		    temp_z/temp_count - path_wcps_vec.at(i).z);
	if (dir3.Angle(dir1)>3.1415926/2.)
	  dir3 *= -1;
      }
    }
    {
      PointVector pts;
      double temp_x = 0;
      double temp_y = 0;
      double temp_z = 0;
      double temp_count = 0;
      for (size_t j=1;j!=15;j++){
	if (i+j<path_wcps_vec.size()){
	  Point pt(path_wcps_vec.at(i+j).x,path_wcps_vec.at(i+j).y,path_wcps_vec.at(i+j).z);
	  if (j<=12&&j>2){
	    temp_x += pt.x;
	    temp_y += pt.y;
	    temp_z += pt.z;
	    temp_count ++;
	  }
	  pts.push_back(pt);
	}
      }
      Point pt(path_wcps_vec.at(i).x,path_wcps_vec.at(i).y,path_wcps_vec.at(i).z);
      dir4 = main_cluster->calc_PCA_dir(pt,pts);
      dir6.SetXYZ(temp_x/temp_count - path_wcps_vec.at(i).x,
      		  temp_y/temp_count - path_wcps_vec.at(i).y,
      		  temp_z/temp_count - path_wcps_vec.at(i).z);
      if (dir4.Angle(dir2)>3.1415926/2.)
	dir4 *= -1;
    }

    int cut1 = 0;
    if ((3.1415926 - dir1.Angle(dir2))/3.1415926*180.>25) cut1++;
    if ((3.1415926 - dir3.Angle(dir4))/3.1415926*180.>25) cut1++;
    if ((3.1415926 - dir5.Angle(dir6))/3.1415926*180.>25) cut1++;
    int cut2 = 0;
    if (fabs(3.1415926/2.-drift_dir.Angle(dir1-dir2))/3.1415926*180. > 5) cut2++;
    if (fabs(3.1415926/2.-drift_dir.Angle(dir3-dir4))/3.1415926*180. > 5) cut2++;
    if (fabs(3.1415926/2.-drift_dir.Angle(dir5-dir6))/3.1415926*180. > 5) cut2++;

    
    // if (main_cluster->get_cluster_id()==7)
    //   std::cout << i << " " << path_wcps_vec.at(i).x/units::cm << " " << path_wcps_vec.at(i).y/units::cm << " " << path_wcps_vec.at(i).z/units::cm << " " << (3.1415926 - dir1.Angle(dir2))/3.1415926*180. << " " << (3.1415926 - dir3.Angle(dir4))/3.1415926*180. << " " << (3.1415926 - dir5.Angle(dir6))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(dir1-dir2))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(dir3-dir4))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(dir5-dir6))/3.1415926*180. << " " << cut1 << " " << cut2 << std::endl;
  
   
    
    
    if (cut1>=3 && cut2>=2){
      if ((3.1415926 - dir3.Angle(dir4))/3.1415926*180. > max_angle){
	max_angle = (3.1415926 - dir3.Angle(dir4))/3.1415926*180.;
	max_point = path_wcps_vec.at(i);
      }
      
      count ++;
      if (count >=3){
	TVector3 temp1(path_wcps_vec.at(i).x-wcp1.x,
		       path_wcps_vec.at(i).y-wcp1.y,
		       path_wcps_vec.at(i).z-wcp1.z);
	TVector3 temp2(path_wcps_vec.at(i).x-wcp2.x,
		       path_wcps_vec.at(i).y-wcp2.y,
		       path_wcps_vec.at(i).z-wcp2.z);

	// if (main_cluster->get_cluster_id()==7)
	//   std::cout << "A: " << (3.1415926-temp1.Angle(temp2))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180.<< " " << temp1.Mag()/units::cm << " " << temp2.Mag()/units::cm << std::endl;


	if (((3.1415926-temp1.Angle(temp2))/3.1415926*180. >35 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5|| (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 ) ||
	    ((3.1415926-temp1.Angle(temp2))/3.1415926*180. >32 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5|| (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 )&& temp1.Mag()>10*units::cm && temp2.Mag()>10*units::cm ||
	    ((3.1415926-temp1.Angle(temp2))/3.1415926*180. >25 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5 || (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 )
	    && temp1.Mag()>15*units::cm && temp2.Mag()>15*units::cm){

	  // if (main_cluster->get_cluster_id()==7)
	  //   std::cout << "B: " <<  (!inside_fiducial_volume(max_point,offset_x)) << " " << inside_dead_region(max_point) << std::endl;
	  
	  if ((!inside_fiducial_volume(max_point,offset_x)) || // must be in fiducial
	      inside_dead_region(max_point)&&(3.1415926-temp1.Angle(temp2))/3.1415926*180<45 // not in dead_volume
	      ){ // should not too close to anode 
	  }else{
	    return true;
	  }
	}
      } else if (count>=1){
      	TVector3 temp1(path_wcps_vec.at(i).x-wcp1.x,
      		       path_wcps_vec.at(i).y-wcp1.y,
      		       path_wcps_vec.at(i).z-wcp1.z);
      	TVector3 temp2(path_wcps_vec.at(i).x-wcp2.x,
      		       path_wcps_vec.at(i).y-wcp2.y,
      		       path_wcps_vec.at(i).z-wcp2.z);

	//	std::cout << "B: " << path_wcps_vec.at(i).x/units::cm << " " << path_wcps_vec.at(i).y/units::cm << " " << path_wcps_vec.at(i).z/units::cm << " " << (3.1415926-temp1.Angle(temp2))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. << " " << temp1.Mag()/units::cm << " " << temp2.Mag()/units::cm << std::endl;
	
      	if (((3.1415926-temp1.Angle(temp2))/3.1415926*180. >35 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5 ||
	     (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 )
	    && temp1.Mag()>5*units::cm && temp2.Mag()>5*units::cm 
	    ){
	  //	  std::cout << "AAA" << std::endl; 

      	  if ((!inside_fiducial_volume(max_point,offset_x)) || // must be in fiducial
      	      inside_dead_region(max_point)&&(3.1415926-temp1.Angle(temp2))/3.1415926*180<45 // not in dead_volume
      	      ){ // should not too close to anode 
      	  }else{
      	    return true;
      	  }
      	}
      }
    }else{
      count = 0 ;
      max_angle = 0;
      max_point.x = 0;
      max_point.y = 0;
      max_point.z = 0;
    }
  }
  
  return false;
}


WCP2dToy::ToyFiducial::~ToyFiducial(){
}

bool WCP2dToy::ToyFiducial::check_signal_processing(WCP::Point& p, TVector3& dir, WCP::ToyCTPointCloud& ct_point_cloud, double step, double offset_x){

  if (dir.Mag()==0){
    return true;
  }else{
    dir *= 1./dir.Mag();
    Point temp_p = p;

    int num_points = 0;
    int num_points_dead = 0;

    //  std::cerr << temp_p.x/units::cm << " " << temp_p.y/units::cm << " " << temp_p.z/units::cm << std::endl;
    
    while(inside_fiducial_volume(temp_p,offset_x)){
      num_points ++;
      //if (inside_dead_region(temp_p))
      //	num_points_dead ++;

      //      std::cerr << temp_p.x/units::cm << " " << temp_p.y/units::cm << " " << temp_p.z/units::cm << " ";
      
      WCP::CTPointCloud<double> cloud_u = ct_point_cloud.get_closest_points(temp_p,1.2*units::cm,0);
      WCP::CTPointCloud<double> cloud_v = ct_point_cloud.get_closest_points(temp_p,1.2*units::cm,1);
      WCP::CTPointCloud<double> cloud_w = ct_point_cloud.get_closest_points(temp_p,1.2*units::cm,2);
      
      //      std::cerr << cloud_u.pts.size() << " " << cloud_v.pts.size() << " " << cloud_w.pts.size() << std::endl;

      if (cloud_u.pts.size()>0 || cloud_v.pts.size()>0 || cloud_w.pts.size() > 0 || inside_dead_region(temp_p))
      	num_points_dead++;
      
      if (num_points - num_points_dead >=5) return true;
	
      temp_p.x += dir.X() * step;
      temp_p.y += dir.Y() * step;
      temp_p.z += dir.Z() * step;
    }

    //    std::cout << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << " " << num_points << " " << num_points_dead << std::endl;
    
    if (num_points_dead > 0.8*num_points){
    	return false;
    }else{
    	return true;
    }
  }
  
  return true;
}

bool WCP2dToy::ToyFiducial::check_dead_volume(WCP::Point& p, TVector3& dir, double step, double offset_x){
  if (!inside_fiducial_volume(p,offset_x)){
    return false;
  }else{
    if (dir.Mag()==0){
      return true;
    }else{
      dir *= 1./dir.Mag();
      Point temp_p = p;
      int num_points = 0;
      int num_points_dead = 0;
      while(inside_fiducial_volume(temp_p,offset_x)){

	num_points ++;
	if (inside_dead_region(temp_p))
	  num_points_dead ++;

	if (num_points - num_points_dead >=4) return true;
	
	temp_p.x += dir.X() * step;
	temp_p.y += dir.Y() * step;
	temp_p.z += dir.Z() * step;

	//	std::cout << temp_p.x/units::cm << " " << temp_p.y/units::cm << " " << temp_p.z/units::cm << std::endl;
      }

      // std::cout << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << " " << num_points << " " << num_points_dead << std::endl;
      
      if (num_points_dead > 0.81*num_points){
	return false;
      }else{
	return true;
      }
      
    }
  }
}

/*
bool WCP2dToy::ToyFiducial::inside_fiducial_volume(WCP::Point& p, double offset_x){

  int c1 = pnpoly(boundary_xy_x, boundary_xy_y, p.x-offset_x, p.y);
  int c2 = pnpoly(boundary_xz_x, boundary_xz_z, p.x-offset_x, p.z);

  //  std::cout << (p.x-offset_x)/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << std::endl;
  //std::cout << c1 << " " << c2 << std::endl;
  
  if (c1 && c2){
    return true;
  }else{
    return false;
  }
}
*/
bool WCP2dToy::ToyFiducial::inside_fiducial_volume(WCP::Point& p, double offset_x, std::vector<double>* tolerance_vec){

	int c1=0;
	int c2=0;
	int index_y = floor((p.y/units::cm+116)/24);
	int index_z = floor(p.z/units::m);
	if(index_y<0){index_y=0;} else if(index_y>9){index_y=9;}
	if(index_z<0){index_z=0;} else if(index_z>9){index_z=9;}
	if(tolerance_vec==NULL){
//		c1 = pnpoly(boundary_xy_x_array[index_z], boundary_xy_y_array[index_z], p.x-offset_x, p.y);
//		c2 = pnpoly(boundary_xz_x_array[index_y], boundary_xz_z_array[index_y], p.x-offset_x, p.z);
		c1 = pnpoly(boundary_xy_x, boundary_xy_y, p.x-offset_x, p.y);
		c2 = pnpoly(boundary_xz_x, boundary_xz_z, p.x-offset_x, p.z);
	} else{
		double tx =     tolerance_vec->at(0);
		double ty_bot = tolerance_vec->at(1);
		double ty_top = tolerance_vec->at(2);
		double tz =     tolerance_vec->at(3);
		//Adjust boundaries for tolerance, defined so that positive tolerance increases volume.
								boundary_SCB_xy_y_array[index_z][0] -= ty_bot;
								boundary_SCB_xy_y_array[index_z][1] -= ty_bot;
		boundary_SCB_xy_x_array[index_z][2] += tx;	boundary_SCB_xy_y_array[index_z][2] -= ty_bot;
		boundary_SCB_xy_x_array[index_z][3] += tx;	boundary_SCB_xy_y_array[index_z][3] += ty_top;
								boundary_SCB_xy_y_array[index_z][4] += ty_top;
								boundary_SCB_xy_y_array[index_z][5] += ty_top;
		
								boundary_SCB_xz_z_array[index_y][0] -= tz;
								boundary_SCB_xz_z_array[index_y][1] -= tz;
		boundary_SCB_xz_x_array[index_y][2] += tx;	boundary_SCB_xz_z_array[index_y][2] -= tz;
		boundary_SCB_xz_x_array[index_y][3] += tx;	boundary_SCB_xz_z_array[index_y][3] += tz;
								boundary_SCB_xz_z_array[index_y][4] += tz;
								boundary_SCB_xz_z_array[index_y][5] += tz;		

		c1 = pnpoly(boundary_SCB_xy_x_array[index_z], boundary_SCB_xy_y_array[index_z], p.x-offset_x, p.y);
		c2 = pnpoly(boundary_SCB_xz_x_array[index_y], boundary_SCB_xz_z_array[index_y], p.x-offset_x, p.z);

		//Revert tolerance shift
								boundary_SCB_xy_y_array[index_z][0] += ty_bot;
								boundary_SCB_xy_y_array[index_z][1] += ty_bot;
		boundary_SCB_xy_x_array[index_z][2] -= tx;	boundary_SCB_xy_y_array[index_z][2] += ty_bot;
		boundary_SCB_xy_x_array[index_z][3] -= tx;	boundary_SCB_xy_y_array[index_z][3] -= ty_top;
								boundary_SCB_xy_y_array[index_z][4] -= ty_top;
								boundary_SCB_xy_y_array[index_z][5] -= ty_top;
		
								boundary_SCB_xz_z_array[index_y][0] += tz;
								boundary_SCB_xz_z_array[index_y][1] += tz;
		boundary_SCB_xz_x_array[index_y][2] -= tx;	boundary_SCB_xz_z_array[index_y][2] += tz;
		boundary_SCB_xz_x_array[index_y][3] -= tx;	boundary_SCB_xz_z_array[index_y][3] -= tz;
								boundary_SCB_xz_z_array[index_y][4] -= tz;
								boundary_SCB_xz_z_array[index_y][5] -= tz;

	}

  //  std::cout << (p.x-offset_x)/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << std::endl;
  //std::cout << c1 << " " << c2 << std::endl;
  
  if (c1 && c2){
    return true;
  }else{
    return false;
  }
}

bool WCP2dToy::ToyFiducial::inside_dead_region(WCP::Point& p){
  // convert the position into U, V, W, and T number ...
  int time_slice = p.x * slope_t + offset_t;
  double pos_u = cos(angle_u) * p.z - sin(angle_u) *p.y;
  double pos_v = cos(angle_v) * p.z - sin(angle_v) *p.y;
  double pos_w = cos(angle_w) * p.z - sin(angle_w) *p.y;
  int ch_u = pos_u * slope_u + offset_u;
  int ch_v = pos_v * slope_v + offset_v + 2400;
  int ch_w = pos_w * slope_w + offset_w + 4800;

  //std::cout << ch_u << " " << ch_v << " " << ch_w << " " << time_slice << std::endl;
  //  std::cout << slope_w << " " << offset_w << " " << pos_w << std::endl;
  
  if (time_slice <0 || time_slice >=2398) return false;
  if (ch_u <0 || ch_u>=2400)  return false;
  if (ch_v <2400 || ch_v>=4800)  return false;
  if (ch_w <4800 || ch_w>=8256)  return false;

  std::set<SlimMergeGeomCell*> dead_u_mcells;
  std::set<SlimMergeGeomCell*> dead_v_mcells;
  std::set<SlimMergeGeomCell*> dead_w_mcells;

  if (ch_mcell_set_map.find(ch_u)!=ch_mcell_set_map.end())
    dead_u_mcells = ch_mcell_set_map[ch_u];
  if (ch_mcell_set_map.find(ch_v)!=ch_mcell_set_map.end())
    dead_v_mcells = ch_mcell_set_map[ch_v];
  if (ch_mcell_set_map.find(ch_w)!=ch_mcell_set_map.end())
    dead_w_mcells = ch_mcell_set_map[ch_w];
  
  // std::cout << ch_u << " " << ch_v << " " << ch578d-_w << " " << dead_u_mcells.size() << " " << dead_v_mcells.size() << " " << dead_w_mcells.size() << std::endl;
  
  // find the dead region given the U, V, and W number
  std::set<SlimMergeGeomCell*> results;
  for (auto it = dead_u_mcells.begin(); it!=dead_u_mcells.end(); it++){
    if (dead_v_mcells.find(*it)!=dead_v_mcells.end()){
      // compare UV sets
      results.insert(*it);
    }else if (dead_w_mcells.find(*it)!=dead_w_mcells.end()){
      // compare UW sets
      results.insert(*it);
    }
  }
  // compare VW sets ...
  for (auto it = dead_v_mcells.begin(); it!=dead_v_mcells.end(); it++){
    if (dead_w_mcells.find(*it)!=dead_w_mcells.end()){
      // compare UW sets
      results.insert(*it);
    }
  }
  
  
  // Check the T number for the remaining T ... 
  for (auto it = results.begin(); it!=results.end(); it++){
    if (mcell_time_map[*it].first <= time_slice &&
	time_slice <= mcell_time_map[*it].second)
      return true;
  }
  
  return false;
}


void WCP2dToy::ToyFiducial::AddDeadRegion(WCP::SlimMergeGeomCell* mcell, std::vector<int>& time_slices){

  mcells.push_back(mcell);
  int start_time = time_slices.front() - dead_region_ch_ext ;
  int end_time = time_slices.back() + dead_region_ch_ext;
  mcell_time_map[mcell] = std::make_pair(start_time, end_time);

  GeomWireSelection& uwires = mcell->get_uwires();
  GeomWireSelection& vwires = mcell->get_vwires();
  GeomWireSelection& wwires = mcell->get_wwires();

  // std::cout << uwires.size() << " " << vwires.size() << " " << wwires.size() << " " << start_time << " " << end_time << std::endl;
  
  std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();

  if (find(bad_planes.begin(), bad_planes.end(), WirePlaneType_t(0))!=bad_planes.end()){
    int start_ch = uwires.front()->channel() - dead_region_ch_ext;
    if (start_ch <0) start_ch = 0;
    int end_ch = uwires.back()->channel() + dead_region_ch_ext;
    if (end_ch>=2400) end_ch = 2399;
    for (int i = start_ch; i<=end_ch;i++){
      if (ch_mcell_set_map.find(i)==ch_mcell_set_map.end()){
	std::set<SlimMergeGeomCell*> mcells_set;
	mcells_set.insert(mcell);
	ch_mcell_set_map[i] = mcells_set;
      }else{
	ch_mcell_set_map[i].insert(mcell);
      }
    }
  }
  if (find(bad_planes.begin(), bad_planes.end(), WirePlaneType_t(1))!=bad_planes.end()){
    int start_ch = vwires.front()->channel() - dead_region_ch_ext;
    if (start_ch <2400) start_ch = 2400;
    int end_ch = vwires.back()->channel() + dead_region_ch_ext;
    if (end_ch>=4800) end_ch = 4799;
    for (int i = start_ch; i<=end_ch;i++){
      if (ch_mcell_set_map.find(i)==ch_mcell_set_map.end()){
	std::set<SlimMergeGeomCell*> mcells_set;
	mcells_set.insert(mcell);
	ch_mcell_set_map[i] = mcells_set;
      }else{
	ch_mcell_set_map[i].insert(mcell);
      }
    }
  }
  if (find(bad_planes.begin(), bad_planes.end(), WirePlaneType_t(2))!=bad_planes.end()){
    int start_ch = wwires.front()->channel() - dead_region_ch_ext;
    if (start_ch <4800) start_ch = 4800;
    int end_ch = wwires.back()->channel() + dead_region_ch_ext;
    if (end_ch>=8256) end_ch = 8255;
    for (int i = start_ch; i<=end_ch;i++){
      if (ch_mcell_set_map.find(i)==ch_mcell_set_map.end()){
	std::set<SlimMergeGeomCell*> mcells_set;
	mcells_set.insert(mcell);
	ch_mcell_set_map[i] = mcells_set;
      }else{
	ch_mcell_set_map[i].insert(mcell);
      }
    }
  }
  
  
}
