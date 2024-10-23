# Event Classification System in ToyFiducial for MicroBooNE
## 1. Overview

The event classification system in ToyFiducial implements multiple layers of event categorization, primarily focused on:
- Neutrino candidate identification
- Through-going muon (TGM) detection
- Light mismatch (LM) identification
- Signal processing quality assessment
- Cosmic ray tagging

## 2. Core Classification Functions

### 2.1 Neutrino Candidate Check
```cpp
bool check_neutrino_candidate(WCP::PR3DCluster *main_cluster, 
                            WCP::WCPointCloud<double>::WCPoint& wcp1,
                            WCP::WCPointCloud<double>::WCPoint& wcp2,
                            double offset_x,
                            WCP::ToyCTPointCloud& ct_point_cloud,
                            bool flag_2view_check = true)
```

#### Implementation Details:
1. **Path Analysis**
   ```cpp
   main_cluster->Create_graph(ct_point_cloud);
   main_cluster->dijkstra_shortest_paths(wcp1);
   main_cluster->cal_shortest_path(wcp2);
   ```

2. **Quality Checks**
   ```cpp
   // Point spacing control
   double low_dis_limit = 0.5*units::cm;
   
   // Track sampling
   for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++) {
       if (path_wcps_vec.size()==0) {
           Point p((*it).x,(*it).y,(*it).z);
           path_wcps_vec.push_back(p);
       } else {
           double dis = sqrt(pow((*it).x - path_wcps_vec.back().x,2)
                           +pow((*it).y - path_wcps_vec.back().y,2)
                           +pow((*it).z - path_wcps_vec.back().z,2));
           if (dis > low_dis_limit) {
               Point p((*it).x,(*it).y,(*it).z);
               path_wcps_vec.push_back(p);
           }
       }
   }
   ```

3. **Two-View Validation**
   ```cpp
   if (flag_2view_check) {
       // Angular checks for U and V planes
       TVector3 dir_1(0,dir.Y(),dir.Z());
       double angle1 = dir_1.Angle(U_dir);
       double angle2 = dir_1.Angle(V_dir);
       double angle3 = dir_1.Angle(W_dir);
       
       if ((angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)) {
           flag_2view_check = false;
       }
   }
   ```

### 2.2 Through-Going Muon Check
```cpp
bool check_tgm(WCP::FlashTPCBundle *bundle, 
               double offset_x,
               WCP::ToyCTPointCloud& ct_point_cloud,
               std::map<WCP::PR3DCluster*, WCP::PR3DCluster*>& old_new_cluster_map,
               int flag)
```

#### Key Features:
1. **Cluster Analysis**
   ```cpp
   std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = 
       main_cluster1->get_extreme_wcps();
   ```

2. **Track Length Validation**
   ```cpp
   double length_limit = sqrt(pow(out_vec_wcps.at(0).at(0).x-out_vec_wcps.at(1).at(0).x,2)+
                            pow(out_vec_wcps.at(0).at(0).y-out_vec_wcps.at(1).at(0).y,2)+
                            pow(out_vec_wcps.at(0).at(0).z-out_vec_wcps.at(1).at(0).z,2));
   ```

### 2.3 Light Mismatch Checks

#### Basic LM Check
```cpp
int check_LM(WCP::FlashTPCBundle *bundle, double& cluster_length)
```
- Evaluates basic light-charge matching quality
- Considers cluster length and PE predictions

#### Cut-based LM Check
```cpp
int check_LM_cuts(WCP::FlashTPCBundle *bundle, double& cluster_length)
```

Key thresholds:
```cpp
if(total_pred_pe < 25 || cluster_length < 10) {
    return 1; // low energy event
}

if(flag_boundary) { // at anode or cathode
    if(flag_anode) {
        if(!(log(total_pred_pe/total_meas_pe)>-1.8 && ks_dis<0.8)) {
            return 2; // light mismatch
        }
    }
}
```

#### BDT-based LM Check
```cpp
int check_LM_bdt(WCP::FlashTPCBundle *bundle, double& cluster_length)
```
Uses machine learning approach:
```cpp
WCP::LMBDT lm(total_pred_pe, total_meas_pe, max_meas_pe, ks_dis,
              chi2, ndf, cl, temp, flag_anode, flag_boundary); 

if(!lm.isSignal(lm.get_BDT_score_eff_background(0.005,bgd))) {
    return 2; // light mismatch event
}
```

### 2.4 Signal Processing Check
```cpp
bool check_signal_processing(WCP::Point& p, 
                           TVector3& dir,
                           WCP::ToyCTPointCloud& ct_point_cloud,
                           double step,
                           double offset_x)
```

Key features:
```cpp
// Point cloud analysis for each plane
WCP::CTPointCloud<double> cloud_u = ct_point_cloud.get_closest_points(temp_p,1.2*units::cm,0);
WCP::CTPointCloud<double> cloud_v = ct_point_cloud.get_closest_points(temp_p,1.2*units::cm,1);
WCP::CTPointCloud<double> cloud_w = ct_point_cloud.get_closest_points(temp_p,1.2*units::cm,2);

// Quality assessment
if (cloud_u.pts.size()>0 || cloud_v.pts.size()>0 || cloud_w.pts.size() > 0 || 
    inside_dead_region(temp_p))
    num_points_dead++;
```

## 3. Cosmic Tagging System

### 3.1 Main Cosmic Tagger
```cpp
void cosmic_tagger(double eventTime,
                  WCP::OpflashSelection& flashes,
                  FlashTPCBundleSelection *matched_bundles,
                  // ... other parameters ...)
```

#### Key Parameters:
```cpp
// Tolerances
double ks_main_stm_tol = 0.08;    // STM flash match tolerance
double ks_main_tgm_tol = 0.00;    // TGM flash match tolerance
double ks_main_tgm_nof_tol = 0.08;
double stm_flash_tol = 3.0*units::m;
double tgm_flash_tol = 3.0*units::m;
```

#### Classification Logic:
1. Flash-Cluster Matching
   ```cpp
   // PE calculations
   double total_pred_pe = 0;
   std::vector<double>& pred_pmt_light = new_bundle->get_pred_pmt_light();
   
   // Centroid calculations
   double pred_pe_z_centroid = 0;
   double flash_pe_z_centroid = 0;
   ```

2. Boundary Analysis
   ```cpp
   int boundary_num = check_boundary(extreme_points, offset_x, &tol_vec);
   ```

3. Tag Assignment
   ```cpp
   if(boundary_num==2) {tag_type = 22;}                          // TGM
   else if(boundary_num==3 && tag_type != 22) {tag_type = 24;}  // TGM no flash
   else if(boundary_num==1 && tag_type != 22 && 
           tag_type != 24) {tag_type = 23;}                      // STM
   ```

## 4. Helper Functions

### 4.1 Boundary Check
```cpp
int check_boundary(std::vector<std::vector<WCP::WCPointCloud<double>::WCPoint>> extreme_points,
                  double offset_x,
                  std::vector<double>* tol_vec)
```

#### Implementation:
```cpp
bool front_flag = false;
bool back_flag = false;

// Check extreme points
for(int i=0; i<int(extreme_points.size()); i++) {
    for(int ii=0; ii<int(extreme_points[i].size()); ii++) {
        Point p(extreme_points[i][ii].x,
               extreme_points[i][ii].y,
               extreme_points[i][ii].z);
               
        // Extended fiducial volume check
        std::vector<double> neg_tol_vec = {-1*tol_vec->at(0),
                                         -1*tol_vec->at(1),
                                         -1*tol_vec->at(2),
                                         -1*tol_vec->at(3)};
                                         
        if(!inside_fiducial_volume(p,offset_x,tol_vec)) {
            return -1;
        } else if(!inside_fiducial_volume(p,offset_x,&neg_tol_vec)) {
            if(i==0) {front_flag = true;}
            if(i==1) {back_flag = true;}
        }
    }
}

return front_flag + back_flag;
```

## 5. Classification Process Flow

1. **Initial Screening**
   - Track length check (>10cm)
   - PE threshold check (>25 PE)
   - Flash matching quality assessment

2. **Geometric Analysis**
   - Boundary intersection checks
   - Track topology analysis
   - Signal processing quality validation

3. **Light Matching**
   - PE ratio analysis
   - Spatial correlation checks
   - BDT-based classification

4. **Final Classification**
   - TGM identification
   - STM tagging
   - Light mismatch flagging

## 6. Performance Considerations

### 6.1 Optimization Techniques
1. Early exit conditions in checks
2. Efficient point cloud operations
3. Cached boundary calculations

### 6.2 Critical Parameters
1. Angular tolerances
2. PE matching thresholds
3. Boundary extension values

### 6.3 Error Prevention
1. Boundary condition validation
2. Flash matching quality checks
3. Signal processing validation