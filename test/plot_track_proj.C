#include <vector>

void plot_track_proj(TString runno="0_0_0", int cluster_id=1, int min_time = 0, int max_time = 9600, int min_u1_ch = 0, int max_u1_ch = 2399, int min_v1_ch = 2400, int max_v1_ch = 4799, int min_w1_ch = 4800, int max_w1_ch = 8256){
  TString filename = "tracking_" + runno + ".root";
  TChain *T = new TChain("T_proj","T_proj");
  T->AddFile(filename);

  // figure out the boundaries for U, V, W planes ...  
  std::vector<int> *proj_cluster_id = new std::vector<int>;
  std::vector<std::vector<int>> *proj_cluster_channel = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_cluster_timeslice= new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_cluster_charge= new std::vector<std::vector<int>>;
  T->SetBranchAddress("cluster_id",&proj_cluster_id);
  T->SetBranchAddress("channel",&proj_cluster_channel);
  T->SetBranchAddress("time_slice",&proj_cluster_timeslice);
  T->SetBranchAddress("charge",&proj_cluster_charge);

  T->GetEntry(0);
  int saved_id = -1;
  for (size_t i=0;i!=proj_cluster_id->size();i++){
    if (proj_cluster_id->at(i)==cluster_id){
      saved_id = i;
      break;
    }
  }
  if (saved_id != -1){
    //std::cout << saved_id << std::endl;
    // figure out the range of U, V, W ...
    int min_ch_u=2400, max_ch_u=0;
    int min_ch_v=4800, max_ch_v=2400;
    int min_ch_w=8256, max_ch_w=4800;
    int min_time_u=9600, max_time_u=0;
    int min_time_v=9600, max_time_v=0;
    int min_time_w=9600, max_time_w=0;

    for (size_t i=0;i!=proj_cluster_channel->at(saved_id).size();i++){
      int ch = proj_cluster_channel->at(saved_id).at(i);
      int time = proj_cluster_timeslice->at(saved_id).at(i);

      if (time >=min_time && time <= max_time ){
	if (ch >= min_u1_ch && ch <= max_u1_ch){
	  if (time < min_time_u) min_time_u = time;
	  if (time > max_time_u) max_time_u = time;
	  
	  if (ch < min_ch_u) min_ch_u = ch;
	  if (ch > max_ch_u) max_ch_u = ch;
	}else if (ch >= min_v1_ch && ch <= max_v1_ch){
	  if (time < min_time_v) min_time_v = time;
	  if (time > max_time_v) max_time_v = time;
	  
	  if (ch < min_ch_v) min_ch_v = ch;
	  if (ch > max_ch_v) max_ch_v = ch;
	}else if (ch >= min_w1_ch && ch <= max_w1_ch){
	  if (time < min_time_w) min_time_w = time;
	  if (time > max_time_w) max_time_w = time;

	  if (ch < min_ch_w) min_ch_w = ch;
	  if (ch > max_ch_w) max_ch_w = ch;
	}
      }
    }

    TCanvas *c1 = new TCanvas("c1","c1",1200,400);
    c1->Divide(3,1);
    c1->cd(1);

    std::cout << min_time_u << " " << max_time_u << " " << min_ch_u << " " << max_ch_u << std::endl;
    
    if (max_time_u > min_time_u && max_ch_u > min_ch_u){
      TH2F *hu = new TH2F("hu","hu",max_ch_u-min_ch_u+4, min_ch_u-2, max_ch_u+2, max_time_u - min_time_u + 4, min_time_u -2, max_time_u + 2);
      T->Project("hu","time_slice:channel",Form("charge*(cluster_id==%d)",cluster_id));
      hu->Draw("COLZ");
      hu->SetTitle("U plane (time x channel)");
    }

    c1->cd(2);
    if (max_time_v > min_time_v && max_ch_v > min_ch_v){
      TH2F *hv = new TH2F("hv","hv",max_ch_v-min_ch_v+4, min_ch_v-2, max_ch_v+2, max_time_v - min_time_v + 4, min_time_v -2, max_time_v + 2);
      T->Project("hv","time_slice:channel",Form("charge*(cluster_id==%d)",cluster_id));
      hv->Draw("COLZ");
      hv->SetTitle("V plane (time x channel)");
    }
    

    c1->cd(3);
    
    if (max_time_w > min_time_w && max_ch_w > min_ch_w){
      TH2F *hw = new TH2F("hw","hw",max_ch_w-min_ch_w+4, min_ch_w-2, max_ch_w+2, max_time_w - min_time_w + 4, min_time_w -2, max_time_w + 2);
      T->Project("hw","time_slice:channel",Form("charge*(cluster_id==%d)",cluster_id));
      hw->Draw("COLZ");
      hw->SetTitle("W plane (time x channel)");
    }
  }

  
  //  std::cout << proj_cluster_id->size() << std::endl;
}
