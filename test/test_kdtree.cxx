#include "WireCellNanoflann/nanoflann.h"
#include "WireCellData/WCPointCloud.h"

#include "TRandom.h"

using namespace WireCell;

int main()
{
  WCPointCloud<double> cloud;
  cloud.pts.resize(1000);
  for (size_t i = 0; i < 1000; i++)
    {
      cloud.pts[i].x = gRandom->Uniform(-100,100);
      cloud.pts[i].y = gRandom->Uniform(-100,100);
      cloud.pts[i].z = gRandom->Uniform(-100,100);
      SlimMergeGeomCell *cell = new SlimMergeGeomCell(i);
      cloud.pts[i].mcell = cell;
    }

  typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, WCPointCloud<double> > ,
    WCPointCloud<double>,
    3 /* dim */
    > my_kd_tree_t;

  my_kd_tree_t   index(3 /*dim*/, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index.buildIndex();

  const double query_pt[3] = { 0.5, 0.5, 0.5};

  {
    size_t num_results = 5;
    std::vector<size_t>   ret_index(num_results);
    std::vector<double>   out_dist_sqr(num_results);
    
    num_results = index.knnSearch(&query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);
    
    // In case of less points in the tree than requested:
    ret_index.resize(num_results);
    out_dist_sqr.resize(num_results);
    
    std::cout << "knnSearch(): num_results=" << num_results << "\n";
    for (size_t i = 0; i < num_results; i++){
      // int abc = cloud.pts[ret_index[i]].mcell;
      std::cout << "idx["<< i << "]=" << ret_index[i] << " dist["<< i << "]=" << out_dist_sqr[i] << " " << cloud.pts[ret_index[i]].mcell->GetIdent() << std::endl;
    }
    std::cout << "\n";
    
  }

  // ----------------------------------------------------------------
  // radiusSearch():  Perform a search for the N closest points
  // ----------------------------------------------------------------
  {
  	const double search_radius = 100;
  	std::vector<std::pair<size_t,double> >   ret_matches;
  
  	nanoflann::SearchParams params;
  	//params.sorted = false;
  
  	const size_t nMatches = index.radiusSearch(&query_pt[0], search_radius, ret_matches, params);
  
	std::cout << "radiusSearch(): radius=" << search_radius << " -> " << nMatches << " matches\n";
  	for (size_t i = 0; i < nMatches; i++)
	  std::cout << "idx["<< i << "]=" << ret_matches[i].first << " dist["<< i << "]=" << ret_matches[i].second << " " << cloud.pts[ret_matches[i].first].mcell->GetIdent()  << std::endl;
	std::cout << "\n";
  }

  
  return 0;
}
