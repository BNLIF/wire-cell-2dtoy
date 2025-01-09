flowchart TD
    A[Start inside_fiducial_volume] --> B[Calculate index_y and index_z]
    
    B --> C{Bound Check indices}
    C -->|Out of bounds| D[Clamp indices to valid range]
    C -->|Within bounds| E{Check tolerance_vec}
    D --> E
    
    E -->|NULL| F1[Use standard boundaries]
    E -->|Not NULL| F2[Adjust boundaries with tolerances]
    
    F1 --> G[Calculate point in polygon for XY projection]
    F2 --> G
    
    G --> H[Calculate point in polygon for XZ projection]
    
    H --> I{Both checks pass?}
    I -->|Yes| J[Return true]
    I -->|No| K[Return false]
    
    F2 --> L[Restore original boundaries]
    L --> M[Continue to polygon checks]
    M --> G