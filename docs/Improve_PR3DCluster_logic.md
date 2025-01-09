``` mermaid
flowchart TD
    Start([Start]) --> Init[Initialize Time-Channel Maps]
    Init --> ProcessMCells[Process Existing Merged Cells]
    
    subgraph DataCollection[Initial Data Collection]
        ProcessMCells --> CollectU[Collect U Wire Data]
        CollectU --> CollectV[Collect V Wire Data]
        CollectV --> CollectW[Collect W Wire Data]
        CollectW --> StoreMaps[Store Charge Maps]
    end

    StoreMaps --> GetPath[Get Cluster Path Points]
    
    subgraph TrajectoryAnalysis[Trajectory Analysis & Gap Filling]
        GetPath --> CheckDist{Distance > 0.3cm?}
        CheckDist -->|Yes| Interpolate[Interpolate Points]
        CheckDist -->|No| AddPoint[Add Point]
        Interpolate --> ProcessPoint[Process Point]
        AddPoint --> ProcessPoint
        
        ProcessPoint --> Convert[Convert to Time-Channel]
        Convert --> CheckSurrounding[Check Surrounding Channels]
        CheckSurrounding --> AddMissing[Add Missing Associations]
    end

    AddMissing --> CreateHolder[Create WCPHolder]
    
    subgraph CellCreation[New Cell Creation & Filtering]
        CreateHolder --> CreateTiling[Create LowmemTiling]
        CreateTiling --> InitCells[Initialize Good Cells]
        InitCells --> FilterCells[Filter New Cells]
        
        FilterCells --> CheckOverlap{Overlaps with\nExisting Cells?}
        CheckOverlap -->|Yes| AddToMap[Add to New Map]
        CheckOverlap -->|No| MarkUnused[Mark as Unused]
        
        AddToMap --> CheckMore{More Cells?}
        MarkUnused --> CheckMore
        CheckMore -->|Yes| FilterCells
    end
    
    CheckMore -->|No| CreateCluster[Create New Cluster]
    
    subgraph FinalCreation[Final Cluster Creation]
        CreateCluster --> AddCells[Add Valid Cells]
        AddCells --> Cleanup[Cleanup Resources]
    end
    
    Cleanup --> Return([Return New Cluster])

    style DataCollection fill:#f9f,stroke:#333,stroke-width:2px
    style TrajectoryAnalysis fill:#bbf,stroke:#333,stroke-width:2px
    style CellCreation fill:#bfb,stroke:#333,stroke-width:2px
    style FinalCreation fill:#ffb,stroke:#333,stroke-width:2px