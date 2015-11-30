namespace WireCell2dToy {
  class TPCParams {
    double m_pitch_u; // wire pitch u
    double m_pitch_v; // wire pitch v
    double m_pitch_w; // wire pitch w
    double m_ts_width; // time slice width 2 us * 1.6 mm/us ~ 3.2 mm
  public:
    // set defaults
  TPCParams() : m_pitch_u(3*units::mm)
      , m_pitch_v(3*units::mm)
      , m_pitch_w(3*units::mm)
      , m_ts_width(3.2*units::mm)
      {};
    
    // set/get u pitches
    void set_pitch_u(double p) { m_pitch_u = p; }
    double get_pitch_u() { return m_pitch_u; }

    // set/get v pitches
    void set_pitch_v(double p) { m_pitch_v = p; }
    double get_pitch_v() { return m_pitch_v; }
    
    // set/get w pitches
    void set_pitch_w(double p) { m_pitch_w = p; }
    double get_pitch_w() { return m_pitch_w; }
    
    double get_pitch(){return (m_pitch_u + m_pitch_v + m_pitch_w)/3.;};
    
    void set_ts_width(double p){ m_ts_width = p;}
    double get_ts_width(){return m_ts_width;}

    // etc for other parameters you need
  };
 }
