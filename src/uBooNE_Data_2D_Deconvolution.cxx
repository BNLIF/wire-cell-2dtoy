#include "WireCell2dToy/uBooNE_Data_2D_Deconvolution.h"

#include "WireCellData/GeomWire.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "TRandom.h"
#include "TF1.h"
#include "TSpectrum.h"

#include <set>
#include <algorithm>

using namespace WireCell;

WireCell2dToy::uBooNEData2DDeconvolutionFDS::uBooNEData2DDeconvolutionFDS(TH2I *hu_decon, TH2I *hv_decon, TH2I *hw_decon, TTree *T_bad, const WireCell::GeomDataSource& gds)
  : fds (0)
  , gds (gds)
  , max_frames(0)
  , nbin(0)
  , time_offset_uv(0)
  , time_offset_uw(0)
  , overall_time_offset(0)
  , nwire_u(0)
  , nwire_v(0)
  , nwire_w(0)
  , gu_2D_g(0)
  , gv_2D_g(0)
  , gw_2D_g(0)
  , hu_2D_g(0)
  , hv_2D_g(0)
  , hw_2D_g(0)
  , hu_2D_gg(0)
  , hv_2D_gg(0)
  , hw_2D_gg(0)
  , scale_u_2d(0)
  , scale_v_2d(0)
{
  load_results_from_file = true;

  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  Int_t chid=0, plane=0;
  Int_t start_time=0, end_time=0;
  T_bad->SetBranchAddress("chid",&chid);
  T_bad->SetBranchAddress("plane",&plane);
  T_bad->SetBranchAddress("start_time",&start_time);
  T_bad->SetBranchAddress("end_time",&end_time);
  for (Int_t i=0;i!=T_bad->GetEntries();i++){
    T_bad->GetEntry(i);
    
    std::pair<int,int> abc(start_time,end_time);
    if (plane == 0){
      umap[chid] = abc;
    }else if (plane == 1){
      vmap[chid-nwire_u] = abc;
    }else if (plane == 2){
      wmap[chid-nwire_u-nwire_v] = abc;
    }
  }
  
  frame.index =0;
  frame.clear();		// win or lose, we start anew
  
  bins_per_frame = hu_decon->GetNbinsY();

  // U plane
  for (size_t ind=0; ind < hu_decon->GetNbinsX(); ++ind) {
    WireCell::Trace trace;
    trace.chid = ind;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hu_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }
  
  // V plane
  for (size_t ind=0; ind < hv_decon->GetNbinsX(); ++ind) {
    WireCell::Trace trace;
    trace.chid = ind + nwire_u;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }

  // W plane
  for (size_t ind=0; ind < hw_decon->GetNbinsX(); ++ind) {
    WireCell::Trace trace;
    trace.chid = ind + nwire_u + nwire_v;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hw_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }

  // std::cout << umap.size() << std::endl;

// #include "calib_resp_v1.txt"
//   for (int i=0;i!=100;i++){
//     ref_resp_func[i] = ref_ele[i];
//   }
//   for (int j=0;j!=2400;j++){
//     for (int i=0;i!=100;i++){
//       u_resp_func[j][i] = calib_ele_chan[j][i];
//     }
//   }
//   for (int j=0;j!=2400;j++){
//     for (int i=0;i!=100;i++){
//       v_resp_func[j][i] = calib_ele_chan[j+2400][i];
//     }
//   }
//   for (int j=0;j!=3456;j++){
//     for (int i=0;i!=100;i++){
//       w_resp_func[j][i] = calib_ele_chan[j+4800][i];
//     }
//   }

}


WireCell2dToy::uBooNEData2DDeconvolutionFDS::uBooNEData2DDeconvolutionFDS(WireCell::FrameDataSource& fds1, const WireCell::GeomDataSource& gds,WireCell::ChirpMap& umap, WireCell::ChirpMap& vmap, WireCell::ChirpMap& wmap, int nframes_total, float time_offset_uv, float time_offset_uw, float overall_time_offset)
  : fds(&fds1)
  , gds(gds)
  , max_frames(nframes_total)
  , time_offset_uv(time_offset_uv)
  , time_offset_uw(time_offset_uw)
  , overall_time_offset(overall_time_offset)
  , umap(umap)
  , vmap(vmap)
  , wmap(wmap)
  , nbin(0)
  , nwire_u(0)
  , nwire_v(0)
  , nwire_w(0)
  , gu_2D_g(0)
  , gv_2D_g(0)
  , gw_2D_g(0)
  , hu_2D_g(0)
  , hv_2D_g(0)
  , hw_2D_g(0)
  , hu_2D_gg(0)
  , hv_2D_gg(0)
  , hw_2D_gg(0)
  , scale_u_2d(0)
  , scale_v_2d(0)
{  

  load_results_from_file = false;

  bins_per_frame = fds->Get_Bins_Per_Frame();
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  hu_2D_g = new TH2I("hu_2D_g","hu_2D_g",nwire_u,0,nwire_u,bins_per_frame,0,bins_per_frame);
  hv_2D_g = new TH2I("hv_2D_g","hv_2D_g",nwire_v,0,nwire_v,bins_per_frame,0,bins_per_frame);
  hw_2D_g = new TH2I("hw_2D_g","hw_2D_g",nwire_w,0,nwire_w,bins_per_frame,0,bins_per_frame);

  hu_2D_gg = new TH2I("hu_2D_gg","hu_2D_gg",nwire_u,0,nwire_u,bins_per_frame,0,bins_per_frame);
  hv_2D_gg = new TH2I("hv_2D_gg","hv_2D_gg",nwire_v,0,nwire_v,bins_per_frame,0,bins_per_frame);
  hw_2D_gg = new TH2I("hw_2D_gg","hw_2D_gg",nwire_w,0,nwire_w,bins_per_frame,0,bins_per_frame);

  #include "data_70_2D_11.txt"
  
  gu_2D_g = new TGraph*[11];
  gu_2D_g[0] = new TGraph(5000,u_2D_g_0_x,u_2D_g_0_y);
  gu_2D_g[1] = new TGraph(5000,u_2D_g_1_x,u_2D_g_1_y);
  gu_2D_g[2] = new TGraph(5000,u_2D_g_2_x,u_2D_g_2_y);
  gu_2D_g[3] = new TGraph(5000,u_2D_g_3_x,u_2D_g_3_y);
  gu_2D_g[4] = new TGraph(5000,u_2D_g_4_x,u_2D_g_4_y);
  gu_2D_g[5] = new TGraph(5000,u_2D_g_5_x,u_2D_g_5_y);
  gu_2D_g[6] = new TGraph(5000,u_2D_g_6_x,u_2D_g_6_y);
  gu_2D_g[7] = new TGraph(5000,u_2D_g_7_x,u_2D_g_7_y);
  gu_2D_g[8] = new TGraph(5000,u_2D_g_8_x,u_2D_g_8_y);
  gu_2D_g[9] = new TGraph(5000,u_2D_g_9_x,u_2D_g_9_y);
  gu_2D_g[10] = new TGraph(5000,u_2D_g_10_x,u_2D_g_10_y);

  gv_2D_g = new TGraph*[11];
  gv_2D_g[0] = new TGraph(5000,v_2D_g_0_x,v_2D_g_0_y);
  gv_2D_g[1] = new TGraph(5000,v_2D_g_1_x,v_2D_g_1_y);
  gv_2D_g[2] = new TGraph(5000,v_2D_g_2_x,v_2D_g_2_y);
  gv_2D_g[3] = new TGraph(5000,v_2D_g_3_x,v_2D_g_3_y);
  gv_2D_g[4] = new TGraph(5000,v_2D_g_4_x,v_2D_g_4_y);
  gv_2D_g[5] = new TGraph(5000,v_2D_g_5_x,v_2D_g_5_y);
  gv_2D_g[6] = new TGraph(5000,v_2D_g_6_x,v_2D_g_6_y);
  gv_2D_g[7] = new TGraph(5000,v_2D_g_7_x,v_2D_g_7_y);
  gv_2D_g[8] = new TGraph(5000,v_2D_g_8_x,v_2D_g_8_y);
  gv_2D_g[9] = new TGraph(5000,v_2D_g_9_x,v_2D_g_9_y);
  gv_2D_g[10] = new TGraph(5000,v_2D_g_10_x,v_2D_g_10_y);

  
  gw_2D_g = new TGraph*[11];
  gw_2D_g[0] = new TGraph(5000,w_2D_g_0_x,w_2D_g_0_y);
  gw_2D_g[1] = new TGraph(5000,w_2D_g_1_x,w_2D_g_1_y);
  gw_2D_g[2] = new TGraph(5000,w_2D_g_2_x,w_2D_g_2_y);
  gw_2D_g[3] = new TGraph(5000,w_2D_g_3_x,w_2D_g_3_y);
  gw_2D_g[4] = new TGraph(5000,w_2D_g_4_x,w_2D_g_4_y);
  gw_2D_g[5] = new TGraph(5000,w_2D_g_5_x,w_2D_g_5_y);
  gw_2D_g[6] = new TGraph(5000,w_2D_g_6_x,w_2D_g_6_y);
  gw_2D_g[7] = new TGraph(5000,w_2D_g_7_x,w_2D_g_7_y);
  gw_2D_g[8] = new TGraph(5000,w_2D_g_8_x,w_2D_g_8_y);
  gw_2D_g[9] = new TGraph(5000,w_2D_g_9_x,w_2D_g_9_y);
  gw_2D_g[10] = new TGraph(5000,w_2D_g_10_x,w_2D_g_10_y);

  scale_u_2d = 1.0;
  scale_v_2d = 1.0;


  // #include "calib_resp_v1.txt"
  // for (int i=0;i!=100;i++){
  //   ref_resp_func[i] = ref_ele[i];
  // }
  // for (int j=0;j!=2400;j++){
  //   for (int i=0;i!=100;i++){
  //     u_resp_func[j][i] = calib_ele_chan[j][i];
  //   }
  // }
  // for (int j=0;j!=2400;j++){
  //   for (int i=0;i!=100;i++){
  //     v_resp_func[j][i] = calib_ele_chan[j+2400][i];
  //   }
  // }
  // for (int j=0;j!=3456;j++){
  //   for (int i=0;i!=100;i++){
  //     w_resp_func[j][i] = calib_ele_chan[j+4800][i];
  //   }
  // }
}

WireCell2dToy::uBooNEData2DDeconvolutionFDS::~uBooNEData2DDeconvolutionFDS()
{

  if (!load_results_from_file){
    delete hu_2D_g;
    delete hv_2D_g;
    delete hw_2D_g;
    
    delete hu_2D_gg;
    delete hv_2D_gg;
    delete hw_2D_gg;
    
    for (int i=0;i!=11;i++){
      delete gu_2D_g[i];
      delete gv_2D_g[i];
      delete gw_2D_g[i];
    }
    
    delete [] gu_2D_g;
    delete [] gv_2D_g;
    delete [] gw_2D_g;
  }
}

int WireCell2dToy::uBooNEData2DDeconvolutionFDS::jump(int frame_number){
  
  if (load_results_from_file) return frame_number;

  if (frame.index == frame_number) {
    return frame_number;
  }
  
  
  frame.clear();
  fds->jump(frame_number);

  TVirtualFFT::SetTransform(0);


  //#include "calib_resp_v1.txt"
  // for (int i=0;i!=100;i++){
  //   ref_resp_func[i] = ref_ele[i];
  // }
  // for (int j=0;j!=2400;j++){
  //   for (int i=0;i!=100;i++){
  //     u_resp_func[j][i] = calib_ele_chan[j][i];
  //   }
  // }
  // for (int j=0;j!=2400;j++){
  //   for (int i=0;i!=100;i++){
  //     v_resp_func[j][i] = calib_ele_chan[j+2400][i];
  //   }
  // }
  // for (int j=0;j!=3456;j++){
  //   for (int i=0;i!=100;i++){
  //     w_resp_func[j][i] = calib_ele_chan[j+4800][i];
  //   }
  // }


  std::cout << "Deconvolution with garfield field response for 2-D W Plane" << std::endl;
  Deconvolute_2D(2);
 
  std::cout << "Deconvolution with garfield field response for 2-D V Plane" << std::endl;
  Deconvolute_2D(1);

  std::cout << "Deconvolution with garfield field response for 2-D U Plane" << std::endl;
  Deconvolute_2D(0);


  // load data to be added ... 
  for (int i=0;i!=nwire_u;i++){
    Trace t;
    t.chid = i;
    t.tbin = 0;
    t.charge.resize(bins_per_frame, 0.0);
    for (int j=0;j!= bins_per_frame; j++){
      t.charge.at(j) = hu_2D_g->GetBinContent(i+1,j+1);
    }
    frame.traces.push_back(t);
  }

  for (int i=0;i!=nwire_v;i++){
    Trace t;
    t.chid = i+nwire_u;
    t.tbin = 0;
    t.charge.resize(bins_per_frame, 0.0);
    for (int j=0;j!= bins_per_frame; j++){
      t.charge.at(j) = hv_2D_g->GetBinContent(i+1,j+1);
    }
    frame.traces.push_back(t);
  }

  for (int i=0;i!=nwire_w;i++){
    Trace t;
    t.chid = i + nwire_u + nwire_v;
    t.tbin = 0;
    t.charge.resize(bins_per_frame, 0.0);
    for (int j=0;j!= bins_per_frame; j++){
      t.charge.at(j) = hw_2D_g->GetBinContent(i+1,j+1);
    }
    frame.traces.push_back(t);
  }
  

  frame.index = frame_number;
  return frame.index;
}

int WireCell2dToy::uBooNEData2DDeconvolutionFDS::size() const{
  return max_frames;
}

void WireCell2dToy::uBooNEData2DDeconvolutionFDS::Deconvolute_2D(int plane){
  const Frame& frame1 = fds->get();
  size_t ntraces = frame1.traces.size();

  const int nticks = bins_per_frame;
  nbin = nticks;
  double value_re[nticks];
  double value_im[nticks];
  for (int i=0;i!=nticks;++i){
    value_re[i] = 0;
    value_im[i] = 0;
  }
  
  TF1 *filter_time = new TF1("filter_time","(x>0.0)*exp(-0.5*pow(x/[0],[1]))");
  int nwires = 0;
  float scale = 1;
  TGraph **gfield = 0;
  double time_offset = 0;
  if (plane == 0){
    //double par[2]={1.43555e+01/200.*2.,4.95096e+00}; // note filter
    //double par[2]={1.78695e+01/200.*2.,5.33129e+00}; //  2.2 us + MIP + 200 ticks
    //double par[2]={1.86783e+01/200.*2.,5.49510e+00}; //  2.0 us + MIP + 200 ticks
    //double par[2] = {5.77121e+01/800.*2,3.94872e+00}; //2.0 us  + MIP + 800 ticks
    double par[2] = {5.75416e+01/800.*2,4.10358e+00}; // 2.2 us  + MIP + 800 ticks
    // double par[2] = {1.11408e-01,2}; // Gaussian note


    filter_time->SetParameters(par);
    nwires = nwire_u;
    scale = scale_u_2d;
    gfield = gu_2D_g;
  }else if (plane == 1){
    //double par[2]={1.47404e+01/200.*2.,4.97667e+00}; // note filter
    //double par[2]={1.84666e+01/200.*2.,5.60489e+00}; //  2.2 us + MIP + 200 ticks
    //double par[2]={1.79987e+01/200.*2.,5.15181e+00}; //  2.0 us + MIP + 200 ticks
    //double par[2] = {6.05906e+01/800.*2,4.11203e+00}; //2.0 us  + MIP + 800 ticks
    double par[2] = {5.99306e+01/800.*2,4.20820e+00}; // 2.2 us  + MIP + 800 ticks
    // double par[2] = {1.11408e-01,2}; // Gaussian note

    filter_time->SetParameters(par);
    nwires = nwire_v;
    scale = scale_v_2d;
    gfield = gv_2D_g;
    time_offset = - time_offset_uv;// -0.113;
  }else if (plane == 2){
    //double par[2]={1.45874e+01/200.*2.,5.02219e+00}; // note filter
    //double par[2]={1.83044e+01/200.*2.,5.44945e+00}; //  2.2 us + MIP + 200 ticks
    //double par[2]={1.84494e+01/200.*2.,5.26641e+00}; //  2.0 us + MIP + 200 ticks
    //double par[2] = {5.89055e+01/800.*2,3.99122e+00}; //2.0 us  + MIP + 800 ticks
    double par[2] = {5.88802e+01/800.*2,4.17455e+00}; // 2.2 us  + MIP + 800 ticks
    //double par[2] = {1.11408e-01,2}; // Gaussian note
    
    filter_time->SetParameters(par);
    nwires = nwire_w;
    gfield = gw_2D_g;
    time_offset = -time_offset_uw;// + 0.803;
  }
  
  TF1 *filter_g = new TF1("filter_g","(x>0.0)*exp(-0.5*pow(x/[0],2))");
  double par3[1]={1.11408e-01};
  filter_g->SetParameters(par3);

  TF1 *filter_wire = new TF1("filter_wire","exp(-0.5*pow(x/[0],2))");
  double par4[1] = {1.0/sqrt(3.1415926)*1.4};
  filter_wire->SetParameters(par4);

  TF1 *filter_wire1 = new TF1("filter_wire1","exp(-0.5*pow(x/[0],2))");
  double par5[1] = {1.0/sqrt(3.1415926)*3.0};
  filter_wire1->SetParameters(par5);
  
  TF1 *filter_low = new TF1("filter_low","1-exp(-pow(x/0.02,2))");

  TH1F *hr = new TH1F("hr1","hr1",nbin,0,nbin); // half us tick
  const int nchannels = nwires;
  
  double rho_res[21][nticks], phi_res[21][nticks];
  for (int i=0;i!=21;++i){
    for (int j=0;j!=nticks;++j){
      rho_res[i][j] = 0.0;
      phi_res[i][j] = 0.0;
    }
  }

  TH1F hmr_u("hmr_u","hmr_u",nbin,0,nbin); 
  TH1F hpr_u("hpr_u","hpr_u",nbin,0,nbin); 

  for (int j=0;j!=21;j++){ 
    TGraph *gtemp=0;
    if (j==0 || j==20){
      gtemp = gfield[10];
    }else if (j==1 || j==19){
      gtemp = gfield[9];
    }else if (j==2 || j==18){
      gtemp = gfield[8];
    }else if (j==3 || j==17){
      gtemp = gfield[7];
    }else if (j==4 || j==16){
      gtemp = gfield[6];
    }else if (j==5 || j==15){
      gtemp = gfield[5];
    }else if (j==6 || j==14){
      gtemp = gfield[4];
    }else if (j==7 || j==13){
      gtemp = gfield[3];
    }else if (j==8 || j==12){
      gtemp = gfield[2];
    }else if (j==9 || j==11){
      gtemp = gfield[1];
    }else if (j==10){
      gtemp = gfield[0];
    }

    for (int i=0; i!=nbin; i++){  
      double time = hr->GetBinCenter(i+1)/2.-90 ;
      //*** scale factors for 70kV ***//
      double x = time + time_offset;
      if (x > -84 && x  < 15.8){
	hr->SetBinContent(i+1,gtemp->Eval(x)/scale); //70kV
      }else{
	hr->SetBinContent(i+1,0);
      }
    }
    hmr_u.Reset();
    hpr_u.Reset();
    hr->FFT(&hmr_u,"MAG");
    hr->FFT(&hpr_u,"PH");
    
    for (Int_t i=0;i!=nticks;i++){
      rho_res[j][i] = hmr_u.GetBinContent(i+1);
      phi_res[j][i] = hpr_u.GetBinContent(i+1);
    }

    // std::cout << "abc " << j << std::endl;
    //delete hmr_u;
    //delete hpr_u;
  }

    // deconvolution ... 
   std::vector<std::vector<double>> rho_u, phi_u,result_re,result_im;
   for (int i=0;i!=nchannels;i++){
     std::vector<double> temp,temp1,temp2,temp3;
     temp.resize(nticks,0);
     rho_u.push_back(temp);
     temp1.resize(nticks,0);
     phi_u.push_back(temp1);
     temp2.resize(nticks,0);
     result_re.push_back(temp2);
     temp3.resize(nticks,0);
     result_im.push_back(temp3);
   }
   int tbin_save[nchannels];
   for (int i=0;i!=nchannels;++i){
     tbin_save[i] = 0;
   }
   // do the nominal electronics FFT
   #include "calib_resp_v1.txt"

   TH1F hm1("hm1","hm1",nbin,0,nbin);
   TH1F hp1("hp1","hp1",nbin,0,nbin);
   TH1F *htemp1 = new TH1F("htemp1","htemp1",nbin,0,nbin);
   for (int j=0;j!=100;j++){
     htemp1->SetBinContent(j+1,ref_ele[j]);
   }
   htemp1->FFT(&hm1,"MAG");
   htemp1->FFT(&hp1,"PH");
   

   TH1F hm("hm","hm",nbin,0,nbin);
   TH1F hp("hp","hp",nbin,0,nbin);

   TH1F hm2("hm2","hm2",nbin,0,nbin);
   TH1F hp2("hp2","hp2",nbin,0,nbin);

   for (size_t ind=0; ind<ntraces; ++ind) {
     const Trace& trace = frame1.traces[ind];
     int tbin = trace.tbin;
     int chid = trace.chid;
     int nbins = trace.charge.size();
     if (chid >=nwire_u && plane == 0) continue;
     if (plane == 1 && (chid >=nwire_u+nwire_v || chid < nwire_u)) continue;
     if (plane == 2 && chid < nwire_u+nwire_v) continue;

     TH1F *htemp = new TH1F("htemp","htemp",nbin,0,nbin);
     for (int i = tbin;i!=tbin+nbins;i++){
       int tt = i+1;
       htemp->SetBinContent(tt,trace.charge.at(i));
     }
     htemp->FFT(&hm,"MAG");
     htemp->FFT(&hp,"PH");
     
     //do FFT for individual electronics response
     TH1F *htemp2 = new TH1F("htemp2","htemp2",nbin,0,nbin);
     for (int j=0;j!=100;j++){
       if (plane == 0 || plane == 1){
     	 htemp2->SetBinContent(j+1,ref_ele1_ind[j]); // 3.5% lower induction wire plane @ 2.2 us
     	 //htemp2->SetBinContent(j+1,ref_ele[j]*1.1*0.965); // 3.5% lower induction wire plane @ 2 us
	 // htemp2->SetBinContent(j+1,calib_ele_chan[chid][j]); // calibrated resp
       }else{
     	 htemp2->SetBinContent(j+1,ref_ele1_ind[j]); // 3.5% lower induction wire plane @ 2.2 us
     	 //htemp2->SetBinContent(j+1,ref_ele[j]*1.1); // for test purpose @ 2.0 us
     	 //htemp2->SetBinContent(j+1,calib_ele_chan[chid][j]); // calibrated resp
       }
     }
     htemp2->FFT(&hm2,"MAG");
     htemp2->FFT(&hp2,"PH");


     if (plane == 0){
       for (Int_t j=0;j!=nticks;j++){ // Note this 1.1 is trying to take care the difference between 2.0 us vs. 2.2 us ... 
	 // rho_u[chid][j] = hm.GetBinContent(j+1)/1.1/0.965;// * hm1.GetBinContent(j+1) / hm2.GetBinContent(j+1);
	 // phi_u[chid][j] = hp.GetBinContent(j+1);// + hp1.GetBinContent(j+1) - hp2.GetBinContent(j+1);

	 rho_u[chid][j] = hm.GetBinContent(j+1) * hm1.GetBinContent(j+1) / hm2.GetBinContent(j+1);
	 phi_u[chid][j] = hp.GetBinContent(j+1) + hp1.GetBinContent(j+1) - hp2.GetBinContent(j+1);
       }
       tbin_save[chid] = tbin;
     }else if (plane == 1){
       for (Int_t j=0;j!=nticks;j++){
	 // rho_u[chid-nwire_u][j] = hm.GetBinContent(j+1)/1.1/0.965;// * hm1.GetBinContent(j+1) / hm2.GetBinContent(j+1);
	 // phi_u[chid-nwire_u][j] = hp.GetBinContent(j+1);// + hp1.GetBinContent(j+1) - hp2.GetBinContent(j+1);

	 rho_u[chid-nwire_u][j] = hm.GetBinContent(j+1) * hm1.GetBinContent(j+1) / hm2.GetBinContent(j+1);
	 phi_u[chid-nwire_u][j] = hp.GetBinContent(j+1) + hp1.GetBinContent(j+1) - hp2.GetBinContent(j+1);
       }
       tbin_save[chid-nwire_u] = tbin;
     }else if (plane == 2){
       for (Int_t j=0;j!=nticks;j++){
	 // rho_u[chid-nwire_u-nwire_v][j] = hm.GetBinContent(j+1)/1.1;// * hm1.GetBinContent(j+1) / hm2.GetBinContent(j+1);
	 // phi_u[chid-nwire_u-nwire_v][j] = hp.GetBinContent(j+1);// + hp1.GetBinContent(j+1) - hp2.GetBinContent(j+1);

	 rho_u[chid-nwire_u-nwire_v][j] = hm.GetBinContent(j+1) * hm1.GetBinContent(j+1) / hm2.GetBinContent(j+1);
	 phi_u[chid-nwire_u-nwire_v][j] = hp.GetBinContent(j+1) + hp1.GetBinContent(j+1) - hp2.GetBinContent(j+1);
       }
       tbin_save[chid-nwire_u-nwire_v] = tbin;
     }


    
     //delete htemp2;
     
     delete htemp;
   }
   delete htemp1;

   double resp_re[nchannels], resp_im[nchannels];
   for (int i=0;i!=nchannels;i++){
     resp_re[i] = 0.0;
     resp_im[i] = 0.0;
   }
   for (Int_t i=0;i!=nticks;i++){
     Double_t freq=0;
     if (i < nticks/2.){
       freq = i/(1.*nticks)*2.;
     }else{
       freq = (nticks - i)/(1.*nticks)*2.;
     }
     
     for (Int_t j=0;j!=nchannels;j++){
       value_re[j] = rho_u[j][i]*cos(phi_u[j][i]);
       value_im[j] = rho_u[j][i]*sin(phi_u[j][i]);
       if (j<21){
   	 resp_re[j] = rho_res[j][i]*cos(phi_res[j][i]);
   	 resp_im[j] = rho_res[j][i]*sin(phi_res[j][i]);
       }else{
   	 resp_re[j] = 0.;
   	 resp_im[j] = 0.;
       }
     }
     
     Int_t m=nchannels;


   TVirtualFFT *ifft = TVirtualFFT::FFT(1,&m,"C2CFORWARD M K");
     ifft->SetPointsComplex(value_re,value_im);
     ifft->Transform();
     Double_t temp_re[nchannels],temp_im[nchannels];
     ifft->GetPointsComplex(temp_re,temp_im);
     
     ifft->SetPointsComplex(resp_re,resp_im);
     ifft->Transform();
     Double_t temp1_re[nchannels],temp1_im[nchannels];
     ifft->GetPointsComplex(temp1_re,temp1_im);
     
     Double_t temp2_re[nchannels],temp2_im[nchannels];
     for (Int_t j=0;j!=nchannels;j++){
       Double_t freq_wire = 0;
       if (j < nchannels/2.){
   	 freq_wire = j/(1.*nchannels)*2.;
       }else{
   	 freq_wire = (nchannels -j)/(1.*nchannels)*2.;
       }
       
       
       if (temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j]>0){
	 if (plane!=2){
	   temp2_re[j] = (temp_re[j]*temp1_re[j]+temp_im[j]*temp1_im[j])/m/
	     (temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j])*filter_wire->Eval(freq_wire);
	   temp2_im[j] = (temp_im[j]*temp1_re[j]-temp_re[j]*temp1_im[j])
	     /m/(temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j])*filter_wire->Eval(freq_wire);
	 }else{
	   temp2_re[j] = (temp_re[j]*temp1_re[j]+temp_im[j]*temp1_im[j])/m/
	     (temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j])*filter_wire1->Eval(freq_wire);
	   temp2_im[j] = (temp_im[j]*temp1_re[j]-temp_re[j]*temp1_im[j])
	     /m/(temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j])*filter_wire1->Eval(freq_wire);
	 }
       }else{
   	 temp2_re[j] = 0;
   	 temp2_im[j] = 0;
       }
     }
     
     TVirtualFFT *ifft3 = TVirtualFFT::FFT(1,&m,"C2CBACKWARD M K");
     ifft3->SetPointsComplex(temp2_re,temp2_im);
     ifft3->Transform();
     Double_t temp3_re[nchannels],temp3_im[nchannels];
     ifft3->GetPointsComplex(temp3_re,temp3_im);
     
     for (Int_t j=0;j!=nchannels;j++){
       Int_t shift = j - 10;
       if (shift <0) shift += nchannels;
       result_re[j][i] = temp3_re[shift]/nticks;//*filter_u->Eval(freq);
       result_im[j][i] = temp3_im[shift]/nticks;//*filter_u->Eval(freq);
     }
     
     delete ifft;
     delete ifft3;
   }

   for (Int_t chid=0;chid!=nchannels;chid++){
     int n = nticks;
     TVirtualFFT *ifft2 = TVirtualFFT::FFT(1,&n,"C2R M K");
     TVirtualFFT *ifft2_g = TVirtualFFT::FFT(1,&n,"C2R M K");
     double temp_re[nticks],temp_im[nticks];
     double tempg_re[nticks],tempg_im[nticks];
     for (int j=0;j!=nticks;++j){
       temp_re[j] = 0;
       temp_im[j] = 0;

       tempg_re[j] = 0;
       tempg_im[j] = 0;
     }
     for (int j=0;j!=nticks;j++){
       Double_t freq=0;
       if (j < nticks/2.){
   	 freq = j/(1.*nticks)*2.;
       }else{
   	 freq = (nticks - j)/(1.*nticks)*2.;
       }
       temp_re[j] = result_re[chid][j]*filter_time->Eval(freq);// *filter_low->Eval(freq);
       temp_im[j] = result_im[chid][j]*filter_time->Eval(freq);// *filter_low->Eval(freq);
       
       tempg_re[j] = result_re[chid][j]*filter_g->Eval(freq);// *filter_low->Eval(freq);
       tempg_im[j] = result_im[chid][j]*filter_g->Eval(freq);// *filter_low->Eval(freq);
     }
     ifft2->SetPointsComplex(temp_re,temp_im);
     ifft2->Transform();
     ifft2_g->SetPointsComplex(tempg_re,tempg_im);
     ifft2_g->Transform();
     
     TH1 *fb = 0;
     fb = TH1::TransformHisto(ifft2,fb,"Re");
     
     TH1 *fb_g = 0;
     fb_g = TH1::TransformHisto(ifft2_g,fb_g,"Re");
     

     TH1F *htemp = new TH1F("htemp","htemp",nbin,0,nbin);
     for (int i=0;i!=nbin;i++){
       htemp->SetBinContent(i+1,fb->GetBinContent(i+1)/( 14.*1.2*4096./2000.));
     }

     TH1F *htemp_g = new TH1F("htemp_g","htemp_g",nbin,0,nbin);
     for (int i=0;i!=nbin;i++){
       htemp_g->SetBinContent(i+1,fb_g->GetBinContent(i+1)/( 14.*1.2*4096./2000.));
     }

     delete ifft2;
     delete fb;
     delete ifft2_g;
     delete fb_g;

     int start = -1, end = -1;
     if (plane == 0){
       if (umap.find(chid)!=umap.end()){
	 start = umap[chid].first;
	 end = umap[chid].second;
       }
     }else if (plane == 1){
       if (vmap.find(chid)!=vmap.end()){
	 start = vmap[chid].first;
	 end = vmap[chid].second;
       }
     }else if (plane == 2){
       if (wmap.find(chid)!=wmap.end()){
	 start = wmap[chid].first;
	 end = wmap[chid].second;
       }
     }

     if (start !=-1){
       for (int j=start;j<=end;j++){
	 int bin = j-180;
	 if (bin <0 ) bin += bins_per_frame;
	 htemp->SetBinContent(bin+1,0);
	 htemp_g->SetBinContent(bin+1,0);
       }
     }

     // correct baseline 
     restore_baseline(htemp);
     restore_baseline(htemp_g);
     
    
     
     // put results back into the 2-D histogram
     for (int j=0;j!=bins_per_frame;j++){
       int bin = j+180;
       if (bin >= bins_per_frame) bin -= bins_per_frame;
       if (bin >=start && bin <=end){
       }else{
	 if (plane == 0){
	   hu_2D_g->SetBinContent(chid+1,bin+1,int(htemp->GetBinContent(j+1)));
	   hu_2D_gg->SetBinContent(chid+1,bin+1,int(htemp_g->GetBinContent(j+1)));
	 }else if (plane == 1){
	   hv_2D_g->SetBinContent(chid+1,bin+1,int(htemp->GetBinContent(j+1)));
	   hv_2D_gg->SetBinContent(chid+1,bin+1,int(htemp_g->GetBinContent(j+1)));
	 }else if (plane == 2){
	   hw_2D_g->SetBinContent(chid+1,bin+1,int(htemp->GetBinContent(j+1)));
	   hw_2D_gg->SetBinContent(chid+1,bin+1,int(htemp_g->GetBinContent(j+1)));
	 }
       }
     }
     delete htemp;
     delete htemp_g;
   }
   
  delete filter_low;
  delete filter_wire;
  delete filter_wire1;
  delete filter_time;
  delete filter_g;
  delete hr;
  
}

void WireCell2dToy::uBooNEData2DDeconvolutionFDS::Clear(){
  hu_2D_g->Reset();
  hv_2D_g->Reset();
  hw_2D_g->Reset();

  hu_2D_gg->Reset();
  hv_2D_gg->Reset();
  hw_2D_gg->Reset();

  frame.clear();
}


void WireCell2dToy::uBooNEData2DDeconvolutionFDS::restore_baseline(TH1F *htemp){
  //correct baseline 
  double max = htemp->GetMaximum();
  double min = htemp->GetMinimum();
  int nbin_b = max - min;
  if (nbin_b ==0) nbin_b = 1;
  TH1F *h1 = new TH1F("h1","h1",nbin_b,min,max);
  for (int j=0;j!=nbin;j++){
    //    if (j%100==0) std::cout << j << " " << nbin << std::endl;
    if (htemp->GetBinContent(j+1)!=0)
      h1->Fill(htemp->GetBinContent(j+1));
  }
  float ped = h1->GetMaximumBin()*(max-min)/(nbin_b*1.) + min;
  float ave=0,ncount = 0;
  for (int j=0;j!=nbin;j++){
    if (fabs(htemp->GetBinContent(j+1)-ped)<400 && htemp->GetBinContent(j+1)!=0){
      ave +=htemp->GetBinContent(j+1);
      ncount ++;
      }
    }
    if (ncount==0) ncount=1;
    ave = ave/ncount;
    
    for (int j=0;j!=nbin;j++){
      double content = htemp->GetBinContent(j+1);
      content -= ave;
      if (htemp->GetBinContent(j+1)!=0)
	htemp->SetBinContent(j+1,content);
    }
    delete h1;
}
