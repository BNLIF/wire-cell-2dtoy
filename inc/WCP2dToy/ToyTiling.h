#ifndef WIRECELL_TOYTILING_H
#define WIRECELL_TOYTILING_H

#include "WCPTiling/TilingBase.h"
#include "WCPData/GeomWCPMap.h"
#include "WCPNav/SliceDataSource.h"
#include "WCPSst/GeomDataSource.h"
#include "WCPNav/DetectorGDS.h"

#include "Rtypes.h"

namespace WCP2dToy {
    /**
     *  A bogus tiling class that doesn't do anything.
     */
    class ToyTiling : public WCP::TilingBase { 
    public:
      ToyTiling();
      // fixme: tiling should not know about slices data.
      
      // ToyTiling(const WCP::Slice& slice, WCP::GeomDataSource& gds);
      ToyTiling(const WCP::Slice& slice, WCP::GeomDataSource& gds, float rel_u = 0.05, float rel_v=0.05, float rel_w=0.05, float noise_u=14000*0.05, float noise_v=14000*0.03, float noise_w=14000*0.02, std::vector<float>* uplane_rms = 0, std::vector<float>* vplane_rms = 0, std::vector<float>* wplane_rms = 0);
      ToyTiling(const WCP::Slice& slice, WCP::DetectorGDS& gds, float rel_u = 0.05, float rel_v=0.05, float rel_w=0.05, float noise_u=14000*0.05, float noise_v=14000*0.03, float noise_w=14000*0.02, std::vector<float>* uplane_rms = 0, std::vector<float>* vplane_rms = 0, std::vector<float>* wplane_rms = 0);
      
      void twoplane_tiling(int time, int nrebin, WCP::GeomDataSource& gds, std::vector<float>& uplane_rms, std::vector<float>& vplane_rms, std::vector<float>& wplane_rms, WCP::ChirpMap& uplane_map, WCP::ChirpMap& vplane_map, WCP::ChirpMap& wplane_map);
      
      void AddCell(WCP::GeomDataSource& gds, WCP::GeomCell *cell, int u_index, int v_index, int w_index, 
		   float u_charge, float v_charge, float w_charge,
		   float u_charge_err, float v_charge_err, float w_charge_err);

      void AddCell(WCP::DetectorGDS& gds, int cryo, int apa, WCP::GeomCell *cell, int u_index, int v_index, int w_index, 
		   float u_charge, float v_charge, float w_charge,
		   float u_charge_err, float v_charge_err, float w_charge_err);


      void CreateCell(float tolerance, const WCP::GeomDataSource& gds, int face, int n_tpc, WCP::GeomWireSelection& temp_uwire, 
		      WCP::GeomWireSelection& temp_vwire, 
		      WCP::GeomWireSelection& temp_wwire);
      
      virtual ~ToyTiling();
      
      WCP::GeomWireSelection wires(const WCP::GeomCell& cell) const;
      WCP::GeomCellSelection cells(const WCP::GeomWire& wire) const;
      const WCP::GeomCell* cell(const WCP::GeomWireSelection& wires) const;
      
      WCP::GeomCellSelection get_allcell(){ return cell_all;}
      WCP::GeomWireSelection get_allwire(){ return wire_all;}
      int get_ncell(){return ncell;}

      WCP::GeomWireSelection get_wire_u(){ return wire_u;}
      WCP::GeomWireSelection get_wire_v(){ return wire_v;}
      WCP::GeomWireSelection get_wire_w(){return wire_w;}
      WCP::GeomWireSelection get_wire_all(){ return wire_all;}

      WCP::WireChargeMap& wcmap(){return wirechargemap;};
      WCP::WireChargeMap& wcemap(){return wirecharge_errmap;};

      std::map<int,float>& ccmap(){return channelchargemap;};
      std::map<int,float>& ccemap(){return channelcharge_errmap;};

      WCP::GeomCellMap& cmap(){return cellmap;};
      WCP::GeomWireMap& wmap(){return wiremap;};
      
      float get_ave_charge(){return ave_charge;};

    protected:
      WCP::GeomWireSelection wire_u;
      WCP::GeomWireSelection wire_v;
      WCP::GeomWireSelection wire_w;
      WCP::GeomWireSelection wire_all;
      int ncell;
      WCP::GeomCellSelection cell_all;
      



      WCP::GeomCellMap cellmap;
      WCP::GeomWireMap wiremap;
      
      WCP::WireChargeMap wirechargemap;
      WCP::WireChargeMap wirecharge_errmap;

      std::map<int,float> channelchargemap;
      std::map<int,float> channelcharge_errmap;
      
      //
      float ave_charge;
      
      
      ClassDef(ToyTiling,1);

    };

}



#endif
