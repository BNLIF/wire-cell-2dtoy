#include "WireCell2dToy/uBooNE_L1SP.h"

#include "WireCellRess/LassoModel.h"
#include "WireCellRess/ElasticNetModel.h"
#include <Eigen/Dense>

#include "TVirtualFFT.h"

using namespace Eigen;
using namespace WireCell;

WireCell2dToy::uBooNE_L1SP::uBooNE_L1SP(TH2F *hv_raw, TH2F *hv_decon, TH2F *hv_decon_g, int nrebin, double time_offset)
  : hv_raw(hv_raw)
  , hv_decon(hv_decon)
  , hv_decon_g(hv_decon_g)
  , nrebin(nrebin)
  , time_offset(time_offset)
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
    y = y * scaling/2.;
    gv->SetPoint(i,x-time_offset,y);
    
    y = gw_2D_g[0]->Eval(x+3.0);
    for (Int_t j=1;j!=11;j++){
      y += gw_2D_g[j]->Eval(x+3.0) * 2;
    }
    y = y * scaling/2.;
    gw->SetPoint(i,x-time_offset,y);
  }

  for (int i=0;i!=11;i++){
    delete gv_2D_g[i];
    delete gw_2D_g[i];
  }
  delete [] gv_2D_g;
  delete [] gw_2D_g;

  filter_g = new TF1("filter_g","exp(-0.5*pow(x/[0],2))");
  double par3[1]={1.11408e-01};
  filter_g->SetParameters(par3);
  
}

WireCell2dToy::uBooNE_L1SP::~uBooNE_L1SP(){
  delete gv;
  delete gw;
  delete filter_g;
}

void WireCell2dToy::uBooNE_L1SP::AddWireTime_Raw(){
  TH1F *h1 = new TH1F("h1","h1",200,-50,50);
  for (int wire_index= 1166; wire_index<=1905;wire_index++){
    h1->Reset();
    for (int time_tick = 0; time_tick != hv_raw->GetNbinsY();time_tick++){
      h1->Fill(hv_raw->GetBinContent(wire_index+1,time_tick+1));
    }
    double par[3];
    double xq = 0.5;
    h1->GetQuantiles(1,&par[1],&xq);
    xq = 0.5 + 0.34;
    h1->GetQuantiles(1,&par[0],&xq);
    xq = 0.5 - 0.34;
    h1->GetQuantiles(1,&par[2],&xq);

    double cut = 4.2 *  sqrt((pow(par[0]-par[1],2)+pow(par[2]-par[1],2))/2.); // 4.2 sigma cut ...
    if (cut < 9) cut = 9;
    
    for (int time_tick = 0; time_tick != hv_raw->GetNbinsY();time_tick++){
      float content = hv_raw->GetBinContent(wire_index+1,time_tick+1);

      if (content > cut && hv_decon->GetBinContent(wire_index+1,int(time_tick/nrebin+1))<=0){
	double content1 = hv_raw->GetBinContent(wire_index+1,int(time_tick/nrebin)*nrebin+1);
	if (content1 < hv_raw->GetBinContent(wire_index+1,int(time_tick/nrebin)*nrebin+2))
	  content1 = hv_raw->GetBinContent(wire_index+1,int(time_tick/nrebin)*nrebin+2);
	if (content1 < hv_raw->GetBinContent(wire_index+1,int(time_tick/nrebin)*nrebin+3))
	  content1 = hv_raw->GetBinContent(wire_index+1,int(time_tick/nrebin)*nrebin+3);
	if (content1 < hv_raw->GetBinContent(wire_index+1,int(time_tick/nrebin)*nrebin+4))
	  content1 = hv_raw->GetBinContent(wire_index+1,int(time_tick/nrebin)*nrebin+4);
       	content1 *= 500/4.*nrebin;
	hv_decon->SetBinContent(wire_index+1,int(time_tick/nrebin+1),content1);
	hv_decon_g->SetBinContent(wire_index+1,int(time_tick/nrebin+1),content1);
      }
      
      if (content>cut){
	if (init_map.find(wire_index)!=init_map.end()){
	  for (int time_slice = time_tick/nrebin; time_slice < time_tick/nrebin+5; time_slice ++){
	    //if (find(init_map[wire_index].begin(),init_map[wire_index].end(),time_slice)==init_map[wire_index].end())
	    init_map[wire_index].insert(time_slice);
	  }
	}else{
	  std::set<int> times;
	  for (int time_slice = time_tick/nrebin; time_slice < time_tick/nrebin +5; time_slice ++){
	    //	    times.push_back(time_slice);
	    times.insert(time_slice);
	  }
	  init_map[wire_index] = times;
	}
      }
    }
  }
  delete h1;
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
	init_map[wire_index].insert(time_slice);
      }else{
	std::set<int> times;
	times.insert(time_slice);
	init_map[wire_index] = times;
      }
    }
  }
}

void WireCell2dToy::uBooNE_L1SP::Form_rois(int pad){
  
  for (auto it = init_map.begin(); it!=init_map.end(); it++){
    int wire_index = it->first;
    std::set<int>& time_slices_set = it->second;
    std::vector<int> time_slices;
    std::copy(time_slices_set.begin(), time_slices_set.end(), std::back_inserter(time_slices));
    //  std::sort(time_slices.begin(), time_slices.end());

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
      //if (it->second * nrebin + nrebin - it->first * nrebin < 150)
      L1_fit(wire_index, it->first * nrebin, it->second * nrebin + nrebin);
      //std::cout << wire_index << " " << it->first * nrebin << " " << it->second * nrebin + nrebin <<  " " << (it->second - it->first + 1) * nrebin << std::endl;
    }
    
    // std::cout << wire_index << " " << rois_save.size() << std::endl;
    // std::cout << wire_index << " " << time_slices.size() << " " << time_slices.front() << " " << time_slices.back() << " " << hv_decon_g->GetBinContent(wire_index+1,time_slices.front()+1)<< std::endl;
  }
}


void WireCell2dToy::uBooNE_L1SP::L1_fit(int wire_index, int start_tick, int end_tick){
  const int nbin_fit = end_tick-start_tick;
  // std::cout << start_tick << " " << end_tick << std::endl;
  
  // fill the data ... 
  VectorXd init_W = VectorXd::Zero(nbin_fit);
  VectorXd final_beta = VectorXd::Zero(nbin_fit*2);
  
  double temp_sum = 0;
  double temp1_sum = 0;
  for (int i=0; i!= nbin_fit; i++){
    init_W(i) = hv_raw->GetBinContent(wire_index+1,start_tick+i+1);
    if (fabs(init_W(i))>6) {
      temp_sum += init_W(i);
      temp1_sum += fabs(init_W(i));
    }
  }
  
  int flag_l1 = 0; // do nothing
  // 1 do L1 
  if (temp_sum/(temp1_sum*90./nbin_fit)>0.2&& temp1_sum>160){
    flag_l1 = 1; // do L1 ...
  }else if (temp1_sum*90./nbin_fit < 50.){
    flag_l1 = 2; //remove signal ... 
  }

  //std::cout << nbin_fit << " " << wire_index+2400 << " " << start_tick/4. << " " << (end_tick-start_tick)/4. << " " << temp_sum << " " << temp1_sum << " " << flag_l1 << std::endl;
  //  std::cout << flag_l1 << std::endl;
  //flag_l1 = 0;
  
  if (flag_l1==1){
    int n_section = std::round(nbin_fit/120.);
    if (n_section ==0) n_section =1;
    std::vector<int> boundaries;
    for (int i=0;i!=n_section;i++){
      boundaries.push_back(int(i * nbin_fit /n_section));
    }
    boundaries.push_back(nbin_fit);
    
    for (int i=0;i!=n_section;i++){
      // std::cout << i << " " << boundaries.at(i+1)  << " " <<  boundaries.at(i) << " " << boundaries.at(i+1)  - boundaries.at(i) << std::endl;
      
      int temp_nbin_fit = boundaries.at(i+1)-boundaries.at(i);
      VectorXd W = VectorXd::Zero(temp_nbin_fit);
      for (int j=0;j!=temp_nbin_fit;j++){
	W(j) = init_W(boundaries.at(i)+j);
      }
            
      //for matrix G
      MatrixXd G = MatrixXd::Zero(temp_nbin_fit, temp_nbin_fit*2);
      for (int i=0;i!=temp_nbin_fit;i++){ 
	// X 
	Double_t t1 = i/2.; // us, measured time  
	for (int j=0;j!=temp_nbin_fit;j++){
	  // Y ... 
	  Double_t t2 = j/2.; // us, real signal time
	  double delta_t = t1 - t2;
	  if (delta_t >-15-time_offset && delta_t < 10-time_offset){
	    G(i,j) = gw->Eval(delta_t) *500;
	    G(i,temp_nbin_fit+j) = gv->Eval(delta_t)*500; 
	  }
	}
      }
      
      double lambda = 5;//1/2.;
      WireCell::LassoModel m2(lambda, 100000, 0.05);
      m2.SetData(G, W);
      m2.Fit();
      VectorXd beta = m2.Getbeta();
      for (int j=0;j!=temp_nbin_fit;j++){
	final_beta(boundaries.at(i)+j) = beta(j);
	final_beta(nbin_fit + boundaries.at(i)+j) = beta(temp_nbin_fit+j);
      }
    }


    double sum1 = 0;
    double sum2 = 0;
    for (int i=0;i!=nbin_fit;i++){
      sum1 += final_beta(i);
      sum2 += final_beta(nbin_fit+i);
    }
    
    
    
    //std::cout << sum1 << " " << sum2 << std::endl;

   
    
    if (sum1 >6 ){
       
      // add the software filter ... 
      TH1F hl1_signal("hl1_signal","hl1_signal",nbin_fit,0,nbin_fit);
      // the factor of 4 for the induction is based on Brian's new study ...
      // data somehow favors 1 as this factor ... 
      for (int j=0;j!=nbin_fit;j++){
	hl1_signal.SetBinContent(j+1,final_beta(j)*1.15+final_beta(nbin_fit+j)*0.5);
      }
      
      TH1F hl2_signal("hl2_signal","hl2_signal",nbin_fit,0,nbin_fit);
      // faster ... 
      Double_t smear_matrix[21]={0.000305453, 0.000978027, 0.00277049, 0.00694322, 0.0153945,
				 0.0301973, 0.0524048, 0.0804588, 0.109289, 0.131334, 
				 0.139629, 0.131334, 0.109289, 0.0804588, 0.0524048, 
				 0.0301973, 0.0153945, 0.00694322, 0.00277049, 0.000978027,
				 0.000305453};
      for (int j=0;j!=nbin_fit;j++){
	Double_t content = hl1_signal.GetBinContent(j+1);
	if (content>0){
	  for (int k=0;k!=21;k++){
	    int bin = j+k-10;
	    if (bin>=0&&bin<nbin_fit)
	      hl2_signal.SetBinContent(bin+1,hl2_signal.GetBinContent(bin+1)+content*smear_matrix[k]);
	  }
	}
      }
      for (int j=0;j!=nbin_fit;j++){
      	if (hl2_signal.GetBinContent(j+1)<50/500.){ // 50 electrons
      	  hl1_signal.SetBinContent(j+1,0);
      	}else{
      	  hl1_signal.SetBinContent(j+1,hl2_signal.GetBinContent(j+1));
      	}
      }
      
      // // convolute a gaussian filter
      // TH1 *hm = hl1_signal.FFT(0,"MAG");
      // TH1 *hp = hl1_signal.FFT(0,"PH");
      // Double_t value_re[1000],value_im[1000];
      // Int_t nticks = nbin_fit;
      // for (int j=0;j!=nticks;j++){
      // 	Double_t freq=0;
      // 	if (j < nticks/2.){
      // 	  freq = j/(1.*nticks)*2.;
      // 	}else{
      // 	  freq = (nticks - j)/(1.*nticks)*2.;
      // 	}
	
      // 	value_re[j] = hm->GetBinContent(j+1) * cos(hp->GetBinContent(j+1)) * filter_g->Eval(freq)/nticks;
      // 	value_im[j] = hm->GetBinContent(j+1) * sin(hp->GetBinContent(j+1)) * filter_g->Eval(freq)/nticks;
      // }
      // TVirtualFFT *ifft = TVirtualFFT::FFT(1,&nticks,"C2R M K");
      // ifft->SetPointsComplex(value_re,value_im);
      // ifft->Transform();
      // TH1 *fb = TH1::TransformHisto(ifft,0,"Re");
      // for (int j=0;j!=nticks;j++){
      // 	if (fb->GetBinContent(j+1)<50/500.){ // 50 electrons
      // 	  hl1_signal.SetBinContent(j+1,0);
      // 	}else{
      // 	  hl1_signal.SetBinContent(j+1,fb->GetBinContent(j+1));
      // 	}
      // }
      // delete hm;
      // delete hp;
      // delete fb;
      // delete ifft;
      // //

      
      // replace it in the decon_v ...
      for (int i=0;i<nbin_fit/nrebin;i++){
	double content = 0;
	for (int j=0;j!=nrebin;j++){
	  content += hl1_signal.GetBinContent(nrebin*i+j+1);
	  //content += final_beta(nrebin*i+j) + final_beta(nbin_fit + nrebin*i + j) * 2.0;
	}
	content *= 500;

	double content1 = hv_decon->GetBinContent(wire_index+1,start_tick/nrebin+1+i);
	//double content2 = hv_decon_g->GetBinContent(wire_index+1,start_tick/nrebin+1+i);
	
	hv_decon->SetBinContent(wire_index+1,start_tick/nrebin+1+i,content);
	hv_decon_g->SetBinContent(wire_index+1,start_tick/nrebin+1+i,content);
	
	if ((content > 0 || content1>0 ) && init_time_slice_set.find(start_tick/nrebin+i)!=init_time_slice_set.end()){
	  time_slice_set.insert(start_tick/nrebin+i);
	}
	//std::cout << hv_decon->GetBinContent(wire_index+1,start_tick/nrebin+1+i) << " " << content << std::endl;
      }
      //  std::cout << sum1 << " " << sum2 << std::endl;
    }
  }else if (flag_l1==2){
    for (int i=0;i<nbin_fit/nrebin;i++){
      if (hv_decon->GetBinContent(wire_index+1,start_tick/nrebin+1+i)!=0 && init_time_slice_set.find(start_tick/nrebin+i)!=init_time_slice_set.end()){
	  time_slice_set.insert(start_tick/nrebin+i);
      }
      double content = 0;
      hv_decon->SetBinContent(wire_index+1,start_tick/nrebin+1+i,content);
      hv_decon_g->SetBinContent(wire_index+1,start_tick/nrebin+1+i,content);
      
    }
  }
}
