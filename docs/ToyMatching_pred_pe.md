# calculate_pred_pe Function Documentation

## Overview
The `calculate_pred_pe` function calculates predicted photoelectron (PE) signals for TPC clusters, considering their position and timing relative to PMTs. This is a crucial component in matching TPC clusters with light flashes.

## Function Signature
```cpp
void calculate_pred_pe(
    int run_no,                    // Run number
    double eventTime,              // Event timestamp
    int time_offset,               // Time offset in microseconds
    int nrebin,                    // Number of rebinning steps
    double time_slice_width,       // Width of time slice
    WCP::Photon_Library *pl,       // Photon library
    FlashTPCBundle* bundle,        // Flash-TPC bundle
    std::vector<double>* pred_pmt_light,  // Predicted PMT light
    std::vector<std::pair<WCP::PR3DCluster*,double>>* additional_clusters,  // Additional clusters
    WCP::PR3DClusterSelection* other_clusters,  // Other clusters
    WCP::PR3DClusterSelection* more_clusters,   // More clusters
    bool &flag_good_bundle,        // Output flag for good bundle
    bool flag_data,                // Data vs MC flag
    bool flag_timestamp            // Timestamp processing flag
);
```

## Input Parameters Explained

### Essential Parameters
- `run_no`: Run number for the experiment
- `eventTime`: Timestamp of the event being processed
- `time_offset`: Time offset in microseconds
- `nrebin`: Number of rebinning steps for time slices
- `time_slice_width`: Width of each time slice

### Configuration Objects
- `pl`: Pointer to Photon Library containing:
  - PMT mapping information
  - Light response data
  - Scaling factors and errors

### Bundle and Cluster Information
- `bundle`: Main Flash-TPC bundle being processed
- `additional_clusters`: Vector of additional clusters with their weights
- `other_clusters`: Collection of other related clusters
- `more_clusters`: Additional cluster collection for extended matching

### Control Flags
- `flag_data`: Boolean indicating if processing real data (true) or simulation (false)
- `flag_timestamp`: Flag for timestamp-based processing

## Output Parameters

- `pred_pmt_light`: Vector of predicted light signals for each PMT
- `flag_good_bundle`: Boolean flag indicating if bundle passes quality criteria

## Key Variables and Boundaries

```cpp
// Boundary definitions
double high_x_cut = 256 * units::cm;
double high_x_cut_ext1 = + 1.2*units::cm;
double high_x_cut_ext2 = - 2.0*units::cm;
double low_x_cut = 0*units::cm;
double low_x_cut_ext1 = - 2*units::cm;
double low_x_cut_ext2 = + 4.0*units::cm;
```

## Main Algorithm Flow

### Flow Diagram

The complete algorithm flow is visualized in [ToyMatching_pred_pe_logic.md](ToyMatching_pred_pe_logic.md):

### 1. Position Validation
First, the function checks if the cluster's position is within valid boundaries:

```cpp
double first_pos_x = (*((main_cluster->get_time_cells_set_map().begin())
                     ->second.begin()))->get_sampling_points().front().x;
double last_pos_x = (*((main_cluster->get_time_cells_set_map().rbegin())
                    ->second.begin()))->get_sampling_points().front().x;

if (first_pos_x-offset_x > low_x_cut + low_x_cut_ext1 -1.0*units::cm &&
    last_pos_x-offset_x > low_x_cut &&
    last_pos_x-offset_x < high_x_cut + high_x_cut_ext1 &&
    first_pos_x-offset_x < high_x_cut) {
    // Process valid position
}
```

### 2. Boundary Condition Checks
The function identifies clusters near boundaries:

```cpp
if (first_pos_x-offset_x <= low_x_cut + low_x_cut_ext2 && 
    first_pos_x-offset_x > low_x_cut + low_x_cut_ext1 - 1.0*units::cm) {
    bundle->set_flag_close_to_PMT(true);
    bundle->set_flag_at_x_boundary(true);
}
```

### 3. Light Prediction Calculation
For each valid cluster, calculates predicted light:

```cpp
for (auto it4 = mcells.begin(); it4!=mcells.end(); it4++){
    SlimMergeGeomCell *mcell = (*it4);
    if (mcell->get_q()>0){
        PointVector& pts = mcell->get_sampling_points();
        float charge = mcell->get_q()/pts.size();
        
        // Process each point
        for (size_t i=0; i!=pts.size(); i++){
            // Convert position
            p.x = pts.at(i).x - offset_x;
            p.y = pts.at(i).y;
            p.z = pts.at(i).z;
            
            // Get voxel ID and calculate light
            int voxel_id = convert_xyz_voxel_id(p);
            std::list<std::pair<int,float>>& pmt_list = 
                photon_library->at(voxel_id);
                
            // Sum up predictions
            for (auto it5 = pmt_list.begin(); it5!=pmt_list.end(); it5++){
                pred_pmt_light->at((*map_lib_pmt)[it5->first]) += 
                    charge * it5->second / elifetime_ratio;
            }
        }
    }
}
```

### 4. PMT Response Normalization
Applies normalization factors to PMT responses:

```cpp
double norm_factor[32];
for (int i=0; i!=32; i++){
    norm_factor[i] = 1;
}
if (flag_data){
    if ((run_no >= 12809 && (!flag_timestamp)) || 
        (flag_timestamp && eventTime >= 1505170407))
        norm_factor[17] = 0;
}
```

### 5. Final Scaling and Validation
Applies final scaling and checks bundle quality:

```cpp
for (size_t i=0; i!=32; i++){
    pred_pmt_light->at(i) *= scaling_light_mag * norm_factor[i];
    sum1 += flash->get_PE(i);
    sum2 += pred_pmt_light->at(i);
    if (pred_pmt_light->at(i) > max_pe)
        max_pe = pred_pmt_light->at(i);
}

if (valid_conditions) {
    flag_good_bundle = true;
}
```

## Special Considerations

1. **Electron Lifetime Correction**
   - Corrects for finite electron lifetime in the detector
   - Uses attenuation ratio based on drift distance

2. **PMT Mapping**
   - Handles mapping between library and PMT channel numbers
   - Accounts for different channel numbering schemes

3. **Boundary Effects**
   - Special handling for clusters near detector boundaries
   - Additional checks for clusters close to PMTs

4. **Data vs. Simulation**
   - Different processing for real data vs. simulation
   - Includes specific PMT normalizations for data

## Error Handling

- Validates input parameters and cluster positions
- Checks for boundary conditions
- Handles PMT response normalization
- Validates final light predictions

## Performance Impact

- Critical for flash matching accuracy
- Computationally intensive due to point-by-point calculations
- Memory efficient using reference parameters