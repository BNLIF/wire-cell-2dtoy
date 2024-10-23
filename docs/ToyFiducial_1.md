# Fiducial Volume Management System Analysis
## Overview
The fiducial volume management system in the Wire-Cell Prototype is a sophisticated component handling detector volume definitions, dead regions, and volume-based event classification. The system integrates boundary definitions, coordinate transformations, and multiple checking mechanisms.

## 1. Core Variables

### 1.1 Boundary Definition Variables
```cpp
// Basic detector boundaries
double m_top;        // Top distance
double m_bottom;     // Bottom distance
double m_upstream;   // Upstream boundary
double m_downstream; // Downstream boundary
double m_anode;      // Anode position
double m_cathode;    // Cathode position

// Space charge boundaries
double m_sc_bottom_1_x, m_sc_bottom_1_y; // Bottom space charge boundary point 1
double m_sc_bottom_2_x, m_sc_bottom_2_y; // Bottom space charge boundary point 2
double m_sc_top_1_x, m_sc_top_1_y;       // Top space charge boundary point 1
double m_sc_top_2_x, m_sc_top_2_y;       // Top space charge boundary point 2
double m_sc_upstream_1_x, m_sc_upstream_1_z;   // Upstream space charge point 1
double m_sc_upstream_2_x, m_sc_upstream_2_z;   // Upstream space charge point 2
double m_sc_downstream_1_x, m_sc_downstream_1_z; // Downstream space charge point 1
double m_sc_downstream_2_x, m_sc_downstream_2_z; // Downstream space charge point 2
```

### 1.2 Boundary Arrays
```cpp
// Boundary vectors for polygon checks
std::vector<double> boundary_xy_x, boundary_xy_y;  // XY plane boundaries
std::vector<double> boundary_xz_x, boundary_xz_z;  // XZ plane boundaries
std::vector<double> boundary_SCB_xy_x, boundary_SCB_xy_y;  // Space charge XY boundaries
std::vector<double> boundary_SCB_xz_x, boundary_SCB_xz_z;  // Space charge XZ boundaries

// Multi-segment boundary arrays
std::vector<std::vector<double>> boundary_xy_x_array, boundary_xy_y_array;
std::vector<std::vector<double>> boundary_xz_x_array, boundary_xz_z_array;
std::vector<std::vector<double>> boundary_SCB_xy_x_array, boundary_SCB_xy_y_array;
std::vector<std::vector<double>> boundary_SCB_xz_x_array, boundary_SCB_xz_z_array;
```

### 1.3 Coordinate Transformation Parameters
```cpp
// Time to position conversion
double offset_t, slope_t;  // (time_slice - offset_t) / slope_t = position_x

// Wire plane coordinate conversions
double offset_u, slope_u;  // U-plane conversion
double offset_v, slope_v;  // V-plane conversion
double offset_w, slope_w;  // W-plane conversion
double angle_u, angle_v, angle_w;  // Wire plane angles

int dead_region_ch_ext;  // Dead region channel extension
```

## 2. Core Functions

### 2.1 Fiducial Volume Check
```cpp
bool inside_fiducial_volume(WCP::Point& p, double offset_x, std::vector<double>* tolerance_vec)
```
#### Implementation Details:
1. Determines segment indices based on position:
```cpp
int index_y = floor((p.y/units::cm+116)/24);   // some of these are hard-coded numbers to match the input space charge boundary
int index_z = floor(p.z/units::m);
```

2. Boundary checking with tolerances:
```cpp
if(tolerance_vec==NULL) {
    c1 = pnpoly(boundary_xy_x, boundary_xy_y, p.x-offset_x, p.y);
    c2 = pnpoly(boundary_xz_x, boundary_xz_z, p.x-offset_x, p.z);
} else {
    // Applies tolerances to boundaries
    double tx = tolerance_vec->at(0);
    double ty_bot = tolerance_vec->at(1);
    double ty_top = tolerance_vec->at(2);
    double tz = tolerance_vec->at(3);
    // Adjusts boundaries and performs check
}
```

### 2.2 Dead Region Management
```cpp
bool inside_dead_region(WCP::Point& p)
```
#### Implementation Details:
1. Coordinate conversion:
```cpp
int time_slice = p.x * slope_t + offset_t;
double pos_u = cos(angle_u) * p.z - sin(angle_u) * p.y;
double pos_v = cos(angle_v) * p.z - sin(angle_v) * p.y;
double pos_w = cos(angle_w) * p.z - sin(angle_w) * p.y;
```

2. Channel mapping:
```cpp
int ch_u = pos_u * slope_u + offset_u;
int ch_v = pos_v * slope_v + offset_v + 2400;
int ch_w = pos_w * slope_w + offset_w + 4800;
```

3. Dead region checking:
```cpp
void AddDeadRegion(WCP::SlimMergeGeomCell* mcell, std::vector<int>& time_slices)
{
    // Time window extension
    int start_time = time_slices.front() - dead_region_ch_ext;
    int end_time = time_slices.back() + dead_region_ch_ext;
    
    // Channel mapping
    for (int i = start_ch; i <= end_ch; i++) {
        if (ch_mcell_set_map.find(i) == ch_mcell_set_map.end()) {
            std::set<SlimMergeGeomCell*> mcells_set;
            mcells_set.insert(mcell);
            ch_mcell_set_map[i] = mcells_set;
        } else {
            ch_mcell_set_map[i].insert(mcell);
        }
    }
}
```

### 2.3 Volume Checking Functions
```cpp
bool check_dead_volume(WCP::Point& p, TVector3& dir, double step, double offset_x)
```
#### Implementation:
1. Step-wise volume checking:
```cpp
while(inside_fiducial_volume(temp_p,offset_x)) {
    num_points++;
    if (inside_dead_region(temp_p))
        num_points_dead++;
        
    if (num_points - num_points_dead >= 4) 
        return true;
        
    temp_p.x += dir.X() * step;
    temp_p.y += dir.Y() * step;
    temp_p.z += dir.Z() * step;
}
```

### 2.4 Containment Check
```cpp
bool check_fully_contained(WCP::FlashTPCBundle *bundle, double offset_x, 
    WCP::ToyCTPointCloud& ct_point_cloud, 
    std::map<WCP::PR3DCluster*, WCP::PR3DCluster*>& old_new_cluster_map,
    unsigned int* fail_mode, int flag)
```
#### Key Features:
1. Cluster validation
2. Boundary intersection checks
3. Dead region consideration
4. Signal processing quality assessment

## 3. Helper Functions

### 3.1 Polygon Point Containment
```cpp
int pnpoly(std::vector<double>& vertx, std::vector<double>& verty, 
    double testx, double testy)
```
Implementation of point-in-polygon algorithm:
```cpp
int c = 0;
for (i = 0, j = vertx.size()-1; i < vertx.size(); j = i++) {
    if (((verty[i]>testy) != (verty[j]>testy)) &&
        (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / 
        (verty[j]-verty[i]) + vertx[i]))
        c = !c;
}
```

## 4. Integration Features

### 4.1 Data/MC Configuration
- Different boundary definitions for data and simulation
- Configurable tolerances and extensions
- Adjustable dead region parameters

### 4.2 Coordinate Systems
- Detector coordinates (x, y, z)
- Wire plane coordinates (u, v, w)
- Time slice coordinates

## 5. Performance Considerations

### 5.1 Optimization Techniques
1. Cached boundary arrays for quick lookup
2. Efficient point-in-polygon testing
3. Early exit conditions in volume checks

### 5.2 Memory Management
1. Use of STL containers for dynamic storage
2. Efficient mapping structures for dead regions
3. Optimized boundary representation

## 6. Usage Guidelines

### 6.1 Initialization
```cpp
ToyFiducial fid(dead_region_ch_ext, offset_t, offset_u, offset_v, offset_w,
    slope_t, slope_u, slope_v, slope_w, angle_u, angle_v, angle_w,
    boundary_dis_cut, top, bottom, upstream, downstream, anode, cathode,
    flag_data);
```

### 6.2 Best Practices
1. Initialize with appropriate boundary parameters
2. Consider detector conditions when setting tolerances
3. Regularly update dead region maps
4. Validate boundary definitions

### 6.3 Error Handling
1. Boundary condition validation
2. Coordinate range checking
3. Dead region overlap handling