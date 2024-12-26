flowchart TD
    Start([Start]) --> CheckDir{"Is dir.Mag()\n== 0?"}
    CheckDir -->|Yes| ReturnTrue1[Return true]
    CheckDir -->|No| NormDir[Normalize direction vector]
    NormDir --> InitCount1[Initialize point counter]
    NormDir --> InitCount2[Initialize dead point counter]
    InitCount1 --> SetPoint[Set temp_p = input point p]
    
    SetPoint --> CheckVol{"Is temp_p inside\nfiducial volume?"}
    
    CheckVol -->|Yes| IncrPoints[Increment num_points]
    CheckVol -->|No| CheckRatio[Check final ratio]
    
    IncrPoints --> GetU[Get closest points in U plane]
    GetU --> GetV[Get closest points in V plane]
    GetV --> GetW[Get closest points in W plane]
    
    GetW --> CheckPlanes{"Any plane has points\nOR point in dead region?"}
    
    CheckPlanes -->|Yes| IncrDead[Increment num_points_dead]
    CheckPlanes -->|No| Continue[Continue]
    
    IncrDead --> CheckPoints{"num_points -\nnum_points_dead >= 5?"}
    Continue --> CheckPoints
    
    CheckPoints -->|Yes| ReturnTrue2[Return true]
    CheckPoints -->|No| UpdatePos[Update temp_p position]
    
    UpdatePos --> CheckVol
    
    CheckRatio --> FinalCheck{"num_points_dead >\n0.8 * num_points?"}
    
    FinalCheck -->|Yes| ReturnFalse[Return false]
    FinalCheck -->|No| ReturnTrue3[Return true]

    classDef process fill:#000000,stroke:#a8d5ff,stroke-width:2px
    classDef decision fill:#000000,stroke:#ffcccc,stroke-width:2px
    classDef terminator fill:#000000,stroke:#d1e7dd,stroke-width:2px

    class Start,ReturnTrue1,ReturnTrue2,ReturnTrue3,ReturnFalse terminator
    class CheckDir,CheckVol,CheckPlanes,CheckPoints,FinalCheck decision
    class NormDir,InitCount1,InitCount2,SetPoint,GetU,GetV,GetW,IncrPoints,IncrDead,Continue,UpdatePos process