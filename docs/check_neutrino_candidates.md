# check_neutrino_candidates() Function Documentation

## Purpose
The function determines whether a track segment has characteristics consistent with being a neutrino interaction candidate. It analyzes the trajectory and properties of a particle track to distinguish neutrino interactions from cosmic ray backgrounds.

## Function Signature
```cpp
bool check_neutrino_candidate(
    WCP::PR3DCluster *main_cluster,
    WCPointCloud<double>::WCPoint& wcp1,
    WCPointCloud<double>::WCPoint& wcp2, 
    double offset_x,
    WCP::ToyCTPointCloud& ct_point_cloud,
    bool flag_2view_check = true
)
```

## Input Parameters

- `main_cluster`: A 3D cluster object representing the particle track
- `wcp1`: Starting point of the track segment
- `wcp2`: End point of the track segment  
- `offset_x`: X-coordinate offset for drift time corrections
- `ct_point_cloud`: Point cloud object containing track points
- `flag_2view_check`: Boolean flag to control whether to check 2 wire plane views (default: true)

## Return Value
- Returns `true` if the track segment is likely from a neutrino interaction
- Returns `false` if the track segment is likely not from a neutrino interaction

## Key Analysis Steps

### Flow Diagram

The complete algorithm flow is visualized in [check_neutrino_candidates_logic.md](check_neutrino_candidates_logic.md):
![check_neutrino_candidates Logic Flow](check_neutrino_candidates_logic.png)

1. **Path Analysis**
   - Creates a graph representation of points in the track
   - Finds shortest path between start and end points
   - Samples points along path at regular intervals

2. **2-View Consistency Check** (if flag_2view_check is true)
   - Checks track visibility in different wire plane views (U, V, W planes)
   - Requires good signals in at least 2 planes for track segments
   - Analyzes track angles relative to wire plane orientations

3. **Track Gap Analysis** 
   - Looks for gaps in the track that could indicate reconstruction issues
   - Counts consecutive points without good signals
   - Flags tracks with gaps larger than thresholds

4. **Angle Analysis**
   - Calculates angles between track segments
   - Analyzes track direction changes
   - Looks for sharp kinks that could indicate neutrino interactions

5. **Dead Region Check**
   - Verifies if track segments pass through detector dead regions
   - Considers impact on track reconstruction quality

## Key Thresholds and Criteria

- Minimum track length: 10 cm
- Gap threshold: 7 points
- Minimum distance from gaps: 25 cm
- Angle thresholds:
  - Sharp turn detection: > 25 degrees
  - Multiple turn requirement: ≥ 3 sharp turns
  - Drift direction angle: > 5 degrees
  
## Neutrino Candidate Criteria

A track segment is considered a neutrino candidate if it:

1. Has consistent signals in required wire views
2. Does not have large unexplained gaps
3. Shows characteristic angle changes consistent with neutrino interactions
4. Is not fully contained within detector dead regions
5. Has sufficient length and track quality

## Usage Example
```cpp
bool is_neutrino = check_neutrino_candidate(
    cluster,
    start_point,
    end_point,
    drift_offset,
    point_cloud,
    true
);
```

## Code Examples

### 1. Path Analysis
```cpp
// Create graph and find shortest path
main_cluster->Create_graph(ct_point_cloud);
main_cluster->dijkstra_shortest_paths(wcp1);
main_cluster->cal_shortest_path(wcp2);

// Get path points and sample at regular intervals
std::list<WCPointCloud<double>::WCPoint>& path_wcps = main_cluster->get_path_wcps();
PointVector path_wcps_vec;
double low_dis_limit = 0.5*units::cm;

// Sample points with minimum distance between them
for (auto it = path_wcps.begin(); it != path_wcps.end(); it++) {
    if (path_wcps_vec.size() == 0) {
        Point p((*it).x, (*it).y, (*it).z);
        path_wcps_vec.push_back(p);
    } else {
        double dis = sqrt(pow((*it).x - path_wcps_vec.back().x, 2) +
                         pow((*it).y - path_wcps_vec.back().y, 2) +
                         pow((*it).z - path_wcps_vec.back().z, 2));
        if (dis > low_dis_limit) {
            Point p((*it).x, (*it).y, (*it).z);
            path_wcps_vec.push_back(p);
        }
    }
}
```

### 2. Wire Plane View Analysis
```cpp
// Define wire plane directions
TVector3 drift_dir(1, 0, 0);
TVector3 U_dir(0, cos(60./180.*3.1415926), sin(60./180.*3.1415926));
TVector3 V_dir(0, cos(60./180.*3.1415926), -sin(60./180.*3.1415926));
TVector3 W_dir(0, 1, 0);

// Check track angle relative to wire planes
TVector3 dir(wcp2.x-wcp1.x, wcp2.y-wcp1.y, wcp2.z-wcp1.z);
TVector3 dir_1(0, dir.Y(), dir.Z());

// Calculate angles with each plane
double angle1 = dir_1.Angle(U_dir);
double angle2 = dir_1.Angle(V_dir);
double angle3 = dir_1.Angle(W_dir);

// Project onto drift direction
TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1), 0);
double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;

// Check if track is too parallel to any plane
if (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5) {
    // Track might be difficult to reconstruct in this plane
}
```

### 3. Track Gap Detection
```cpp
// Check for gaps in track reconstruction
int num_nth = 0;  // Counter for points without good signals
double min_dis = 1e9;  // Distance to track endpoints

// Analyze each point along track
for (int i = 0; i < path_wcps_vec1.size(); i++) {
    // Get nearby points in each wire plane
    WCP::CTPointCloud<double> cloud_u = ct_point_cloud.get_closest_points(path_wcps_vec1.at(i), 1.2*units::cm, 0);
    WCP::CTPointCloud<double> cloud_v = ct_point_cloud.get_closest_points(path_wcps_vec1.at(i), 1.2*units::cm, 1);
    WCP::CTPointCloud<double> cloud_w = ct_point_cloud.get_closest_points(path_wcps_vec1.at(i), 1.2*units::cm, 2);

    bool flag_reset = false;
    
    // Check if point has good signals
    if (cloud_u.pts.size() > 0 || cloud_v.pts.size() > 0 || cloud_w.pts.size() > 0) {
        flag_reset = true;
    }
    
    if (flag_reset) {
        num_nth = 0;  // Reset gap counter
    } else {
        num_nth++;  // Increment gap counter
        
        // Calculate distance to endpoints
        double dis1 = sqrt(pow(path_wcps_vec1.at(i).x-wcp1.x,2) + 
                         pow(path_wcps_vec1.at(i).y-wcp1.y,2) + 
                         pow(path_wcps_vec1.at(i).z-wcp1.z,2));
        if (dis1 < min_dis) min_dis = dis1;
        
        // Check if gap is too large
        if (num_nth > 7 && min_dis < 25*units::cm) {
            return true;  // Gap indicates potential neutrino interaction
        }
    }
}
```

### 4. Angle Analysis for Interaction Vertex
```cpp
// Look for sharp angle changes that might indicate interaction vertex
for (size_t i = 5; i + 5 < path_wcps_vec.size(); i++) {
    // Calculate direction vectors for segments before and after point
    TVector3 dir1(path_wcps_vec.at(i).x - path_wcps_vec.at(i-5).x,
                  path_wcps_vec.at(i).y - path_wcps_vec.at(i-5).y,
                  path_wcps_vec.at(i).z - path_wcps_vec.at(i-5).z);
    TVector3 dir2(path_wcps_vec.at(i).x - path_wcps_vec.at(i+5).x,
                  path_wcps_vec.at(i).y - path_wcps_vec.at(i+5).y,
                  path_wcps_vec.at(i).z - path_wcps_vec.at(i+5).z);
    
    // Calculate angle between segments
    double angle = (3.1415926 - dir1.Angle(dir2))/3.1415926*180.;
    
    // Check for sharp turn
    if (angle > 25) {
        // Analyze angle relative to drift direction
        double drift_angle = fabs(3.1415926/2. - drift_dir.Angle(dir1-dir2))/3.1415926*180.;
        if (drift_angle > 5) {
            // This could be an interaction vertex
        }
    }
}
```

## Notes
- Function is part of neutrino event selection in the MicroBooNE detector
- Designed to work with liquid argon TPC detector data
- Uses both spatial and signal quality information
- Critical for rejecting cosmic ray backgrounds
- Code examples show key analysis techniques used in the algorithm