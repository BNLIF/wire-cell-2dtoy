# ExamineBundles Function Documentation

## Overview
The ExamineBundles system consists of three main functions that work together to analyze and potentially separate clusters in TPC (Time Projection Chamber) data. These functions handle bundles of flash-TPC associations and perform graph-based cluster analysis.

## 1. ExamineBundles Function

### Purpose
This is the top-level function that processes multiple FlashTPCBundles.

### Input
- `bundles`: A FlashTPCBundleSelection (collection of FlashTPCBundle objects)
- `ct_point_cloud`: A ToyCTPointCloud object containing point cloud data

### Output
- Returns a new FlashTPCBundleSelection containing processed bundles

### Logic Flow
1. Creates a set to track used cluster IDs
2. Iterates through each bundle in the input collection
3. Processes each bundle using ExamineBundle
4. Cleans up original clusters and bundles
5. Returns the new collection of bundles

## 2. ExamineBundle Function

### Purpose
Examines a single FlashTPCBundle and potentially separates its clusters based on connectivity.

### Input
- `bundle`: A pointer to a FlashTPCBundle
- `used_cluster_ids`: Reference to a set of already used cluster IDs
- `ct_point_cloud`: Reference to a ToyCTPointCloud

### Output
- Returns a new FlashTPCBundle pointer

### Logic Flow
1. Determines the next available cluster ID
2. Collects all original clusters and their merged cells (MCells)
3. Creates a new temporary cluster containing all MCells
4. Uses Examine_graph to potentially separate the cluster
5. Creates new clusters based on the separation
6. Identifies the main cluster based on overlap with original main cluster
7. Creates a new bundle with the separated clusters
8. Transfers all properties from original bundle to new bundle

## 3. Examine_graph Function

### Purpose
Performs the actual graph-based analysis to determine if a cluster should be separated into multiple components.

### Input
- `ct_point_cloud`: Reference to a ToyCTPointCloud object

### Output
- Returns a vector of SMGCSelection (collections of SlimMergeGeomCells)

### Logic Flow
1. Cleanup and initialization:
   - Deletes existing graph if present
   - Creates point cloud if not already created
2. Graph construction:
   - Creates a new graph with N vertices (N = number of points)
   - Establishes connections between close points
   - Adds additional protection against over-clustering
3. Component analysis:
   - Performs connected components analysis on the graph
   - Identifies distinct connected components
4. Result organization:
   - Creates separate collections for each component
   - Associates MCells with their respective components
   - Returns the separated collections

## Key Concepts

### FlashTPCBundle
A data structure that associates:
- Flash data (light detection)
- TPC clusters (charge detection)
- Various matching quality indicators

### Clustering
The process uses graph theory to:
1. Represent points in 3D space as vertices
2. Connect nearby points with edges
3. Identify connected components
4. Separate clusters based on connectivity

### MCells (Merged Cells)
- Basic units that make up clusters
- Contain charge deposition information
- Are assigned to components based on point cloud connectivity

## Important Notes
1. Memory management is handled carefully with proper deletion of original clusters and bundles
2. The system maintains cluster ID uniqueness through the used_cluster_ids set
3. Properties from original bundles are preserved in new bundles
4. The main cluster is determined by maximum overlap with the original main cluster

## Code Examples

### Flow Diagram

The complete algorithm flow is visualized in [ExamineBundles_logic.md](ExamineBundles_logic.md):

### 1. Basic Usage Example
```cpp
// Initialize required objects
WCP::FlashTPCBundleSelection bundles;
WCP::ToyCTPointCloud ct_point_cloud;

// Add some bundles (simplified example)
FlashTPCBundle* bundle1 = create_sample_bundle();  // Implementation depends on your data
bundles.push_back(bundle1);

// Process the bundles
FlashTPCBundleSelection new_bundles = WCP2dToy::ExamineBundles(bundles, ct_point_cloud);

// Work with the results
for (auto bundle : new_bundles) {
    std::cout << "Main cluster ID: " << bundle->get_main_cluster()->get_cluster_id() << std::endl;
    std::cout << "Number of other clusters: " << bundle->get_other_clusters().size() << std::endl;
}
```

### 2. Examining a Single Bundle
```cpp
// Set up necessary objects
std::set<int> used_cluster_ids;
WCP::ToyCTPointCloud ct_point_cloud;
FlashTPCBundle* original_bundle = /* your bundle */;

// Examine single bundle
FlashTPCBundle* new_bundle = WCP2dToy::ExamineBundle(
    original_bundle, 
    used_cluster_ids, 
    ct_point_cloud
);

// Check the results
if (new_bundle->get_other_clusters().size() > 0) {
    std::cout << "Bundle was split into multiple clusters" << std::endl;
}
```

### 3. Graph Analysis Example
```cpp
// Assuming we have a PR3DCluster object
PR3DCluster* cluster = /* your cluster */;
WCP::ToyCTPointCloud ct_point_cloud;

// Examine the cluster's graph structure
std::vector<SMGCSelection> separated_clusters = cluster->Examine_graph(ct_point_cloud);

// Process the results
for (size_t i = 0; i < separated_clusters.size(); i++) {
    std::cout << "Component " << i << " contains " 
              << separated_clusters[i].size() << " cells" << std::endl;
    
    // Access individual cells in each component
    for (auto cell : separated_clusters[i]) {
        // Work with individual cells
        int time_slice = cell->GetTimeSlice();
        // ... other operations
    }
}
```

### 4. Complete Processing Pipeline
```cpp
#include "WCP2dToy/ExamineBundles.h"

// Initialize the system
WCP::FlashTPCBundleSelection input_bundles;
WCP::ToyCTPointCloud ct_point_cloud;
std::set<int> used_cluster_ids;

// Load or create your bundles
// ... bundle creation code ...

// Process all bundles
FlashTPCBundleSelection output_bundles = WCP2dToy::ExamineBundles(input_bundles, ct_point_cloud);

// Process results
for (auto bundle : output_bundles) {
    // Access the main cluster
    PR3DCluster* main_cluster = bundle->get_main_cluster();
    
    // Get associated flash information
    WCP::Flash* flash = bundle->get_flash();
    
    // Check matching quality
    double chi2 = bundle->get_chi2();
    double ks_dist = bundle->get_ks_dis();
    
    // Access other clusters if any were separated
    PR3DClusterSelection& other_clusters = bundle->get_other_clusters();
    
    // Process cluster information
    std::cout << "Bundle processing results:" << std::endl;
    std::cout << "  Main cluster cells: " << main_cluster->get_mcells().size() << std::endl;
    std::cout << "  Other clusters: " << other_clusters.size() << std::endl;
    std::cout << "  Chi2: " << chi2 << std::endl;
    std::cout << "  KS distance: " << ks_dist << std::endl;
}

// Cleanup is handled by ExamineBundles for input_bundles
// Make sure to clean up output_bundles when done
for (auto bundle : output_bundles) {
    delete bundle;
}
```

These examples demonstrate typical usage patterns of the ExamineBundles system, from basic bundle processing to complete pipeline implementation. They show how to:
- Process multiple bundles
- Examine individual bundles
- Work with the graph analysis results
- Handle the separated clusters
- Access bundle properties
- Properly manage memory

Remember that these examples assume the existence of certain helper functions and proper initialization of objects. In a real implementation, you would need to ensure proper data loading and error handling.