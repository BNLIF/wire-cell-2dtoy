::: mermaid
flowchart TB
    Start([Start]) --> GetPos[Get Cluster Position]
    GetPos --> ValidatePos{Validate\nPosition}
    
    subgraph PosCheck[Position Validation]
        ValidatePos -->|Invalid| SetFlagFalse[Set flag_good_bundle = false]
        ValidatePos -->|Valid| CheckBoundary{Check Boundary\nConditions}
        CheckBoundary -->|Near Boundary| SetFlags[Set Boundary Flags]
        CheckBoundary --> ProcessClusters[Process Clusters]
    end
    
    subgraph LightCalc[Light Calculation]
        ProcessClusters --> LoopCells[Loop Through Cells]
        LoopCells --> CalcCharge[Calculate Charge]
        CalcCharge --> LoopPoints[Loop Through Points]
        LoopPoints --> ConvertPos[Convert Position]
        ConvertPos --> GetVoxel[Get Voxel ID]
        GetVoxel --> CalcLight[Calculate Light Response]
        CalcLight --> UpdatePred[Update Predictions]
        UpdatePred -->|More Points| LoopPoints
        UpdatePred -->|Done| NormalizePMT[Normalize PMT Response]
    end
    
    subgraph FinalProcess[Final Processing]
        NormalizePMT --> ApplyScale[Apply Scaling Factors]
        ApplyScale --> ValidateResults{Validate\nResults}
        ValidateResults -->|Valid| SetFlagTrue[Set flag_good_bundle = true]
        ValidateResults -->|Invalid| SetFlagFalse
    end
    
    SetFlagTrue --> Return([Return])
    SetFlagFalse --> Return

    classDef process fill:#000000,stroke:#a8d5ff,stroke-width:2px;
    classDef decision fill:#000000,stroke:#ffcccc,stroke-width:2px;
    classDef terminator fill:#000000,stroke:#d1e7dd,stroke-width:2px;
    classDef subgraphBox fill:#000000,stroke:#f5f5f5,stroke-width:2px;

    class Start,Return terminator;
    class ValidatePos,CheckBoundary,ValidateResults decision;
    class ProcessClusters,CalcCharge,CalcLight process;
    class PosCheck,LightCalc,FinalProcess subgraphBox;

:::