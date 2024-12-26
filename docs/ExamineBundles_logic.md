:::mermaid
flowchart TD
    subgraph ExamineBundles
        A[Start ExamineBundles] --> B[Initialize new_bundles]
        B --> C{Process each bundle}
        C --> |For each bundle| D[Call ExamineBundle]
        D --> E[Add new bundle to collection]
        E --> C
        C --> |Done| F[Cleanup original clusters]
        F --> G[Return new bundles]
    end

    subgraph ExamineBundle
        H[Start ExamineBundle] --> I[Get next cluster ID]
        I --> J[Collect original clusters & MCells]
        J --> K[Create temporary cluster]
        K --> L[Call Examine_graph]
        L --> M[Create new clusters]
        M --> N[Identify main cluster]
        N --> O[Create new bundle]
        O --> P[Transfer properties]
        P --> Q[Return new bundle]
    end

    subgraph Examine_graph
        R[Start Examine_graph] --> S[Initialize/cleanup graph]
        S --> T[Create point cloud]
        T --> U[Build graph]
        U --> V[Establish connections]
        V --> W[Protect against overclustering]
        W --> X[Find connected components]
        X --> Y[Organize MCells by component]
        Y --> Z[Return separated MCells]
    end

    C --> |Call| H
    L --> |Call| R
:::