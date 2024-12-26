flowchart TD
    A[Start] --> B[Is point inside fiducial volume?]
    B --> |No| C[Return false]
    B --> |Yes| D[Is direction magnitude zero?]
    
    D --> |Yes| E[Return true]
    D --> |No| F[Normalize direction vector]
    
    F --> G[Initialize tracking:
    temp_p = p
    num_points = 0
    num_points_dead = 0]
    
    G --> H[Check position]
    
    H --> |Inside fiducial| I[num_points++]
    H --> |Outside fiducial| J[Check dead ratio]
    
    I --> K[Check dead region]
    
    K --> |In dead region| L[num_points_dead++]
    K --> |Not dead| M[Continue]
    
    L --> N[Check live points]
    M --> N
    
    N --> |≥ 4 live points| O[Return true]
    N --> |< 4 live points| P[Update position:
    Move step size in direction]
    
    P --> H
    
    J --> |> 81% dead| Q[Return false]
    J --> |≤ 81% dead| R[Return true]

    style B fill:#f9f
    style D fill:#f9f
    style H fill:#f9f
    style K fill:#f9f
    style N fill:#f9f
    style J fill:#f9f