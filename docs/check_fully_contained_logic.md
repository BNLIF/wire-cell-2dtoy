flowchart TD
    A[Start] --> B{Check flag}
    
    B -->|flag == 1| C1[Get main cluster from bundle]
    B -->|flag == 2| C2[Get original cluster from bundle]
    
    C1 --> D[Check old_new_cluster_map]
    C2 --> D
    
    D -->|Better cluster exists| E[Use new cluster]
    D -->|No better cluster| F[Use current cluster]
    
    E --> G[Get extreme points]
    F --> G
    
    G --> H[For each extreme point]
    
    H --> I{Check inside_fiducial_volume}
    I -->|false| J[Set fail_mode bit 2]
    J --> K[Return false]
    
    I -->|true| L[Calculate angles with U/V/W planes]
    
    L --> M{Angle < threshold?<br>U/V < 10° or W < 5°}
    
    M -->|yes| N{Check signal_processing}
    N -->|fail| O[Set fail_mode bit 1]
    O --> K

    M -->|no| P{Check angle with<br>main direction}
    N -->|pass| P
    
    P -->|> 60°| Q{Check dead_volume}
    Q -->|fail| R[Set fail_mode bit 0]
    R --> K
    
    P -->|≤ 60°| S[Continue to<br>next point]
    Q -->|pass| S
    
    S --> |More points| H
    S --> |No more points| T[Return true]
    
    style K fill:#ffcccc
    style T fill:#ccffcc