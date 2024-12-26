::: mermaid
flowchart TB
    Start([Start]) --> Init[Initialize Bundle Collections]
    Init --> CreateBundles[Create Initial Flash-TPC Bundles]
    
    subgraph BundleCreation[Bundle Creation and Initial Filtering]
        CreateBundles --> CalcPE[Calculate Predicted PE]
        CalcPE --> CheckPos{Check Position\nConstraints}
        CheckPos -->|Valid| AddBundle[Add to Bundle Collection]
        CheckPos -->|Invalid| SkipBundle[Skip Bundle]
    end

    subgraph BundleEval[Bundle Evaluation]
        AddBundle --> ExamineBundles[Examine Bundles]
        ExamineBundles --> CalcMetrics[Calculate KS Distance\nand Chi-Square]
        CalcMetrics --> CheckMetrics{Check Quality\nMetrics}
        CheckMetrics -->|Good| MarkConsistent[Mark as Consistent]
        CheckMetrics -->|Poor| KeepForLater[Keep for Further Analysis]
    end

    subgraph Optimization[Match Optimization]
        MarkConsistent --> PrepLasso[Prepare LASSO Model]
        KeepForLater --> PrepLasso
        PrepLasso --> BuildMatrix[Build Optimization Matrices]
        BuildMatrix --> RunLasso[Run LASSO Optimization]
        RunLasso --> ProcessResults[Process Optimization Results]
    end

    subgraph FinalSelection[Final Selection]
        ProcessResults --> CreateMatches[Create Final Matches] 
        CreateMatches --> ResolveConflicts[Resolve Conflicts]
        ResolveConflicts --> OrganizeBundles[Organize Final Bundles]
    end

    OrganizeBundles --> Return([Return Bundle Selection])

    classDef process fill:#a8d5ff,stroke:#333,stroke-width:2px;
    classDef decision fill:#ffcccc,stroke:#333,stroke-width:2px;
    classDef terminator fill:#d1e7dd,stroke:#333,stroke-width:2px;
    classDef groupbox fill:#f5f5f5,stroke:#666,stroke-width:2px;

    class Start,Return terminator;
    class CheckPos,CheckMetrics decision;
    class CalcPE,ExamineBundles,BuildMatrix,RunLasso,CreateMatches process;
    class BundleCreation,BundleEval,Optimization,FinalSelection groupbox;
:::
