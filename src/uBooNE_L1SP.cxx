#include "WireCell2dToy/uBooNE_L1SP.h"

#include "WireCellRess/LassoModel.h"
#include "WireCellRess/ElasticNetModel.h"
#include <Eigen/Dense>


using namespace Eigen;
using namespace WireCell;

WireCell2dToy::uBooNE_L1SP::uBooNE_L1SP(TH2F *hv_raw, TH2F *hv_decon, TH2F *hv_decon_g, int nrebin)
  : hv_raw(hv_raw)
  , hv_decon(hv_decon)
  , hv_decon_g(hv_decon_g)
  , nrebin(nrebin)
{
  #include "data_70_2D_11.txt"
  TGraph **gv_2D_g = new TGraph*[11];
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

  TGraph **gw_2D_g = new TGraph*[11];
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

  gv = new TGraph();
  gw = new TGraph();
  
  int scaling = 14 * 1.2 * 4096/2000.;
  
  for (Int_t i=0;i!=gv_2D_g[0]->GetN();i++){
    Double_t x,y;
    gv_2D_g[0]->GetPoint(i,x,y);
    for (Int_t j=1;j!=11;j++){
      y += gv_2D_g[j]->Eval(x) * 2;
    }
    y = y * scaling;
    gv->SetPoint(i,x,y);
    
    y = gw_2D_g[0]->Eval(x+3.0);
    for (Int_t j=1;j!=11;j++){
      y += gw_2D_g[j]->Eval(x+3.0) * 2;
    }
    y = y * scaling;
    gw->SetPoint(i,x,y);
  }

  for (int i=0;i!=11;i++){
    delete gv_2D_g[i];
    delete gw_2D_g[i];
  }
  delete [] gv_2D_g;
  delete [] gw_2D_g;
}

WireCell2dToy::uBooNE_L1SP::~uBooNE_L1SP(){
  delete gv;
  delete gw;
}

void WireCell2dToy::uBooNE_L1SP::AddWireTime_Raw(){
  for (int wire_index= 1166; wire_index<=1905;wire_index++){
    for (int time_slice = 0; time_slice != hv_raw->GetNbinsY();time_slice++){
      float content = hv_raw->GetBinContent(wire_index+1,time_slice+1);
      if (content>10){
	if (init_map.find(wire_index)!=init_map.end()){
	  init_map[wire_index].push_back(time_slice);
	}else{
	  std::vector<int> times;
	  times.push_back(time_slice);
	  init_map[wire_index] = times;
	}
      }
    }
  }
}

void WireCell2dToy::uBooNE_L1SP::AddWires(int time_slice, GeomWireSelection& wires){
  if (wires.size() >0){
    
    // save +- 2 time slices ... 
    init_time_slice_set.insert(time_slice);
    init_time_slice_set.insert(time_slice+1);
    init_time_slice_set.insert(time_slice-1);
    init_time_slice_set.insert(time_slice+2);
    init_time_slice_set.insert(time_slice-2);
      
    
    for (auto it = wires.begin(); it!=wires.end(); it++){
      int wire_index = (*it)->index();
      if (init_map.find(wire_index)!=init_map.end()){
	init_map[wire_index].push_back(time_slice);
      }else{
	std::vector<int> times;
	times.push_back(time_slice);
	init_map[wire_index] = times;
      }
    }
  }
}

void WireCell2dToy::uBooNE_L1SP::Form_rois(int pad){
  
  for (auto it = init_map.begin(); it!=init_map.end(); it++){
    int wire_index = it->first;
    std::vector<int> time_slices = it->second;
    std::sort(time_slices.begin(), time_slices.end());

    std::vector<std::pair<int,int>> rois;
    std::vector<std::pair<int,int>> rois_save;
    
    rois.push_back(std::make_pair(time_slices.front(),time_slices.front()));
    for (size_t i=1; i<time_slices.size();i++){
      if (time_slices.at(i) - rois.back().second <= pad*2){
	rois.back().second = time_slices.at(i);
      }else{
	rois.push_back(std::make_pair(time_slices.at(i),time_slices.at(i)));
      }
    }
    
    // extend the rois to both side according to the bin content
    for (auto it = rois.begin(); it!= rois.end();  it++){
      int start_bin = it->first;
      int end_bin = it->second;
      //std::cout << start_bin << " " << end_bin << " " ;
      while(hv_decon_g->GetBinContent(wire_index+1,start_bin)>0){
	start_bin --;
	if (start_bin<=0) break;
      }
      while(hv_decon_g->GetBinContent(wire_index+1,end_bin)>0){
	end_bin ++;
	if (end_bin >= hv_decon_g->GetNbinsX()-1) break;
      }

      // add padding ...
      start_bin = start_bin - pad;
      while(hv_decon_g->GetBinContent(wire_index+1,start_bin)>0){
	start_bin -=pad;
	if (start_bin<=0) break;
      }
      end_bin = end_bin + pad;
      while(hv_decon_g->GetBinContent(wire_index+1,end_bin)>0){
	end_bin +=pad;
	if (end_bin >= hv_decon_g->GetNbinsX()-1) break;
      }
      
      if (start_bin <0) start_bin = 0;
      if (end_bin>hv_decon_g->GetNbinsX()-1) end_bin = hv_decon_g->GetNbinsX()-1;
      it->first = start_bin;
      it->second = end_bin;
      //std::cout << start_bin << " " << end_bin << std::endl;
    }
    
    // std::cout << wire_index << " " << rois.size() << " ";    
    // merge them ...
    for (auto it = rois.begin(); it!= rois.end();  it++){
      if (rois_save.size()==0){
	rois_save.push_back(*it);
      }else if (it->first <= rois_save.back().second){
	rois_save.back().second = it->second;
      }else{
	rois_save.push_back(*it);
      }
    }

    for (auto it = rois_save.begin(); it!=rois_save.end(); it++){
      L1_fit(wire_index, it->first * nrebin, it->second * nrebin + nrebin);
      //   std::cout << wire_index << " " << it->first * nrebin << " " << it->second * nrebin + nrebin <<  " " << (it->second - it->first + 1) * nrebin << std::endl;
    }
    
    // std::cout << wire_index << " " << rois_save.size() << std::endl;
    // std::cout << wire_index << " " << time_slices.size() << " " << time_slices.front() << " " << time_slices.back() << " " << hv_decon_g->GetBinContent(wire_index+1,time_slices.front()+1)<< std::endl;
  }
}


void WireCell2dToy::uBooNE_L1SP::L1_fit(int wire_index, int start_tick, int end_tick){
  const int nbin_fit = end_tick-start_tick;
  // fill the data ... 
  VectorXd W = VectorXd::Zero(nbin_fit);

  double temp_sum = 0;
  double temp1_sum = 0;
  for (int i=0; i!= nbin_fit; i++){
    W(i) = hv_raw->GetBinContent(wire_index+1,start_tick+i+1);
    if (fabs(W(i))>6) {
      temp_sum += W(i);
      temp1_sum += fabs(W(i));
    }
  }

  std::cout << nbin_fit << " " << wire_index+2400 << " " << start_tick/4. << " " << (end_tick-start_tick)/4. << " " << temp_sum << " " << temp1_sum << std::endl;

  int flag_l1 = 0; // do nothing
  // 1 do L1 

  if (temp_sum/(temp1_sum*90./nbin_fit)>0.5&& temp1_sum>160){
    flag_l1 = 1; // do L1 ...
  }else if (0){
    flag_l1 = 2; //remove signal ... 
  }

  
  
  if (flag_l1==1){
    //for matrix G
    MatrixXd G = MatrixXd::Zero(nbin_fit, nbin_fit*2);
    for (int i=0;i!=nbin_fit;i++){ 
      // X 
      Double_t t1 = i/2.; // us, measured time  
      for (int j=0;j!=nbin_fit;j++){
	// Y ... 
	Double_t t2 = j/2.; // us, real signal time
	double delta_t = t1 - t2;
	if (delta_t >-15 && delta_t < 10){
	  G(i,j) = gw->Eval(delta_t) *500;
	  G(i,nbin_fit+j) = gv->Eval(delta_t)*500; 
	}
      }
    }
    
    double lambda = 5;//1/2.;
    WireCell::LassoModel m2(lambda, 100000, 0.05);
    m2.SetData(G, W);
    m2.Fit();
    
    VectorXd beta = m2.Getbeta();
    
    double sum1 = 0;
    double sum2 = 0;
    
    for (int i=0;i!=nbin_fit;i++){
      sum1 += beta(i);
      sum2 += beta(nbin_fit+i);
    }
    
    //std::cout << sum1 << " " << sum2 << std::endl;
    
    if (sum1 >6 ){
      // replace it in the decon_v ...
      for (int i=0;i<nbin_fit/nrebin;i++){
	double content = 0;
	for (int j=0;j!=nrebin;j++){
	  content += beta(nrebin*i+j) + beta(nbin_fit + nrebin*i + j);
	}
	content *= 500;
	hv_decon->SetBinContent(wire_index+1,start_tick/nrebin+1+i,content);
	hv_decon_g->SetBinContent(wire_index+1,start_tick/nrebin+1+i,content);
	
	if (content > 0 && init_time_slice_set.find(start_tick/nrebin+i)!=init_time_slice_set.end()){
	  time_slice_set.insert(start_tick/nrebin+i);
	}
	//std::cout << hv_decon->GetBinContent(wire_index+1,start_tick/nrebin+1+i) << " " << content << std::endl;
      }
      //  std::cout << sum1 << " " << sum2 << std::endl;
    }
  }else if (flag_l1==2){
    for (int i=0;i<nbin_fit/nrebin;i++){
      double content = 0;
      content *= 500;
      hv_decon->SetBinContent(wire_index+1,start_tick/nrebin+1+i,content);
      hv_decon_g->SetBinContent(wire_index+1,start_tick/nrebin+1+i,content);
    }
  }
}
