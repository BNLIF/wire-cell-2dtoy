flowchart TD
    Start([Start]) --> Init[Create graph and calculate shortest path]
    Init --> SamplePoints[Sample points along path]
    SamplePoints --> Check2View{2-View check enabled?}
    
    Check2View -- Yes --> CalcWirePlaneAngles[Calculate angles relative to wire planes]
    Check2View -- No --> PathAnalysis
    
    CalcWirePlaneAngles --> WireAngleCheck{Angle < threshold?<br/>U: 10°<br/>V: 10°<br/>W: 5°}
    WireAngleCheck -- Yes --> SignalCheck[Check signal processing and dead volume]
    WireAngleCheck -- No --> PathAnalysis
    
    SignalCheck --> PathAnalysis[Begin path point analysis]
    
    PathAnalysis --> PointLoop[Analyze next point in path]
    PointLoop --> GetClosestPoints[Get closest points in U,V,W planes]
    
    GetClosestPoints --> CheckSignals{Has good signals<br/>in required planes?}
    CheckSignals -- Yes --> ResetCounters[Reset gap counters]
    CheckSignals -- No --> IncrementGap[Increment gap counter]
    
    ResetCounters --> ContinueCheck{More points<br/>to check?}
    IncrementGap --> GapCheck{Gap > 7 points AND<br/>distance < 25cm?}
    
    GapCheck -- Yes --> ReturnTrue([Return True])
    GapCheck -- No --> ContinueCheck
    
    ContinueCheck -- Yes --> PointLoop
    ContinueCheck -- No --> AngleAnalysis[Begin angle analysis]
    
    AngleAnalysis --> SegmentLoop[Check next track segment]
    SegmentLoop --> CalcAngles[Calculate segment angles]
    
    CalcAngles --> AngleCheck{Sharp turn detected?<br/>>25° AND drift angle >5°}
    AngleCheck -- Yes --> CountTurns[Increment turn counter]
    AngleCheck -- No --> NextSegment{More segments<br/>to check?}
    
    CountTurns --> TurnCheck{Enough sharp turns?<br/>count >= 3}
    TurnCheck -- Yes --> CheckPosition{Point in fiducial<br/>volume?}
    TurnCheck -- No --> NextSegment
    
    CheckPosition -- Yes --> ReturnTrue
    CheckPosition -- No --> NextSegment
    
    NextSegment -- Yes --> SegmentLoop
    NextSegment -- No --> ReturnFalse([Return False])

    style ReturnTrue fill:#9f9,stroke:#393
    style ReturnFalse fill:#f99,stroke:#933
    style Start fill:#99f,stroke:#339