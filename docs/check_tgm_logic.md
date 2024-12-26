``` mermaid
flowchart TD
    A[Start check_tgm] --> B{Check flag}
    B -->|flag=1| C1[Get current main cluster]
    B -->|flag=2| C2[Get original cluster]
    
    C1 --> D[Get extreme points]
    C2 --> D
    
    D --> E[Calculate length limit]
    
    E --> F[Loop through point groups i]
    F --> G[Check if points in group i are inside fiducial volume]
    
    G --> H[Loop through remaining groups k]
    H --> I[Check if points in group k are inside fiducial volume]
    
    I --> J{Both groups\noutside fiducial?}
    J -->|Yes| K[Check middle points]
    
    K --> L{Middle points\ninside fiducial?}
    L -->|Yes| M{Is flash type 2?}
    
    M -->|Yes| N[Calculate track length]
    N --> O{Length > 0.45 * limit?}
    
    O -->|Yes| P{Is neutrino\ncandidate?}
    P -->|No| Q[Return true]
    P -->|Yes| R[Continue checking]
    
    L -->|No| S{Track length meets\nminimum?}
    S -->|Yes| Q
    S -->|No| R
    
    J -->|No| T{Check track\nangles}
    T --> U{Angles meet\nTGM criteria?}
    U -->|Yes| V[Check dead regions]
    U -->|No| R
    
    V --> W{Dead region\ncriteria met?}
    W -->|Yes| Q
    W -->|No| R
    
    R --> F
    
    F -->|All checked| X[Return false]