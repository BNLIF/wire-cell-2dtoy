flowchart TB
    Start([Start]) --> InitMaps[Initialize Bundle Maps]
    
    subgraph GroupingPhase["Initial Bundle Grouping"]
        InitMaps --> GroupBundles[Group Bundles by Flash]
        GroupBundles --> SeparateBundles[Separate Good/Remaining Bundles]
        SeparateBundles --> EvalConsistent{Has Consistent\nBundles?}
    end
    
    subgraph ConflictResolution["Conflict Resolution"]
        EvalConsistent -->|Yes| ProcessRemaining[Process Remaining Bundles]
        ProcessRemaining --> CheckCompat{Check Bundle\nCompatibility}
        CheckCompat -->|Compatible| KeepBundle[Keep Bundle]
        CheckCompat -->|Incompatible| MarkRemoval[Mark for Removal]
        EvalConsistent -->|No| SecondRound[To Second Round]
    end
    
    subgraph SecondPhase["Second Round Processing"]
        SecondRound --> EvalSecond[Evaluate Second Round Bundles]
        EvalSecond --> TryMatch[Try to Match\nwith Existing]
        TryMatch -->|Success| UpdateBundle[Update Bundle Maps]
        TryMatch -->|Fail| ThirdRound[To Third Round]
    end
    
    subgraph FinalPhase["Final Processing"]
        ThirdRound --> ProcessThird[Process Third Round Bundles]
        ProcessThird --> CleanupBundles[Cleanup Removed Bundles]
        CleanupBundles --> CreateFinal[Create Final Bundle Set]
    end
    
    CreateFinal --> Return([Return])
    
    classDef process fill:#333,stroke:#a8d5ff,stroke-width:2px
    classDef decision fill:#333,stroke:#ffcccc,stroke-width:2px
    classDef terminator fill:#333,stroke:#d1e7dd,stroke-width:2px
    classDef default fill:#333,stroke:#f5f5f5,stroke-width:2px
    
    class Start,Return terminator
    class EvalConsistent,CheckCompat,TryMatch decision
    class ProcessRemaining,EvalSecond,ProcessThird process
    class GroupingPhase,ConflictResolution,SecondPhase,FinalPhase default