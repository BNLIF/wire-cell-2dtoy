#ifndef WIRECELL_TOYMATRIX_H
#define WIRECELL_TOYMATRIX_H

#include "TMatrixD.h"
#include "TVectorD.h"
#include "WCPData/GeomCell.h"
#include "WCPData/GeomWire.h"
#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"
#include "WCP2dToy/MergeToyTiling.h"

namespace WCP2dToy{
  class ToyMatrix {
  public:
    ToyMatrix();
    ToyMatrix(WCP2dToy::ToyTiling& toytiling, WCP2dToy::MergeToyTiling& mergetiling, int svd_flag, int recon_t = 2000);
    
    ToyMatrix(WCP2dToy::ToyTiling& toytiling, WCP2dToy::MergeToyTiling& mergetiling, int recon_t=2000);
    ToyMatrix(const WCP::DetectorGDS& gds, WCP2dToy::ToyTiling& toytiling, WCP2dToy::MergeToyTiling& mergetiling, int recon_t=2000);
    

    

    virtual ~ToyMatrix();

    

    double Get_Chi2(){return chi2;};
    float Get_Solve_Flag(){return solve_flag;};
    int Get_mcindex(){return mcindex;};
    int Get_mwindex(){return mwindex;};
    int Get_swindex(){return swindex;};

    const TMatrixD* Get_MA(){return MA;};
    const TMatrixD* Get_MB(){return MB;};
    const TMatrixD* Get_MAT(){return MAT;};
    const TMatrixD* Get_MBT(){return MBT;};
    const TMatrixD* Get_Vy(){return Vy;};
    const TMatrixD* Get_VBy(){return VBy;};
    const TMatrixD* Get_Vx(){return Vx;};
    const TMatrixD* Get_VBy_inv(){return VBy_inv;};
    const TMatrixD* Get_Vx_inv(){return Vx_inv;};
    const TMatrixD* Get_MC(){return MC;};
    const TMatrixD* Get_MC_inv(){return MC_inv;};
    const TVectorD* Get_Wy(){return Wy;};
    
    double Get_Cell_Charge( const WCP::GeomCell *cell, int flag = 1) ;
    int Get_mcindex(const WCP::GeomCell *cell){ return mcimap[cell];};
    double Get_residual(const WCP::GeomCell *cell);

    void Set_chi2(double chi2t){chi2 = chi2t;};
    void Set_value(double value, int index){ (*Cx)[index] = value;};
    void Set_error(double value, int index){ (*dCx)[index] = value;};

    double Get_value(int index){return (*Cx)[index];};
    double  Get_error(int index){return (*dCx)[index];};

    void Set_ndf(double ndf1){ndf = ndf1;};
    void Set_Solve_Flag(float solve_flag1){solve_flag=solve_flag1;};
    void Print();
    int Get_ndf(){return ndf;};
    
    //void Set_blob(int num1){num_blob = num1;};
    int Get_blob(){return num_blob;};
    bool GetSimpleBlobReduction(){return simple_blob_reduction;};
    void SetSimpleBlobReduction(bool val){simple_blob_reduction = val;};
    void JudgeSimpleBlob(WCP2dToy::ToyTiling& toytiling, WCP2dToy::MergeToyTiling& mergetiling);
    

    void Update_pred();

    std::vector<int>& Get_svd_removed(){return svd_removed;};
    
    

  protected:
    int recon_threshold;
    std::vector<int> svd_removed;

    int num_blob;
    bool simple_blob_reduction;


    int Solve();
    int Solve_SVD();

    int svd_flag;

    TMatrixD *MA, *MB, *MAT, *MBT;
    TMatrixD *Vy, *VBy, *Vx, *VBy_inv, *Vx_inv;

    TMatrixD *MC, *MC_inv;

    TVectorD *Wy, *Cx, *dCx;
    
    TVectorD *MWy_pred, *MWy;

    double chi2;
    int ndf;
    int solve_flag;

    void Buildup_index(WCP2dToy::MergeToyTiling& mergetiling);
    void Buildup_index(const WCP::DetectorGDS& gds, WCP2dToy::MergeToyTiling& mergetiling);
   

    WCP::WireIndexMap mwimap, swimap;
    WCP::CellIndexMap mcimap;
    
    std::map<int,int> scimap;
    

    int mcindex;
    int mwindex;
    int swindex;
  };
}

#endif
