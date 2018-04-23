#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TVector3.h"

#include <iostream>
#include <vector>

using namespace std;

bool boundary(TVector3 v)
{
    bool isBound1 = true;
    bool isBound2 = true;
    Double_t x = v.X();
    Double_t y = v.Y();
    Double_t z = v.Z();
    // X-Y plane: bottom - (80, -117) --- (256, -98), top - (100, 118) --- (256, 103)
    // 6 distances (sign) to boundary
    Double_t d1 = x;
    Double_t d2 = y-118;
    Double_t d3 = y+117;
    Double_t d4 = x-256;
    Double_t d5 = ((103-118)*x-(256-100)*y+256*118-103*100)/TMath::Sqrt((103-118)*(103-118)+(256-100)*(256-100));
    Double_t d6 = ((-98+117)*x-(256-80)*y+256*(-117)-(-98)*80)/TMath::Sqrt((-98+117)*(-98+117)+(256-80)*(256-80));
   
    //cout<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<" "<<d5<<" "<<d6<<endl;
    Double_t fv = 5; // FC cut, cm
    if(d1>fv && d2<-fv && d3>fv && d4<-fv && d5>fv && d6<-fv) isBound1 = false; // contained 

    // X-Z plane: left - (0, 120) --- (11, 256), right - (1026, 256) --- (1037, 120)
    Double_t dd1 = z;
    Double_t dd2 = x-256;
    Double_t dd3 = x;
    Double_t dd4 = z - 1037;
    Double_t dd5 = ((120-256)*x-(1037-1026)*y+1037*256-120*1026)/TMath::Sqrt((120-256)*(120-256)+(1037-1026)*(1037-1026));
    Double_t dd6 = ((256-120)*x-(11-0)*y+11*120-256*0)/TMath::Sqrt((256-120)*(256-120)+(11-0)*(11-0));
    //cout<<dd1<<" "<<dd2<<" "<<dd3<<" "<<dd4<<" "<<dd5<<" "<<dd6<<endl;
    if(dd1>fv && dd2<-fv && dd3>fv && dd4<-fv && dd5>fv && dd6>fv) isBound2 = false; // contained 

    //cout<<"Bound1: "<<isBound1<<endl;
    //cout<<"Bound2: "<<isBound2<<endl;
    return isBound1||isBound2;
}


int main(int argc, char* argv[])
{
    if(argc < 2){
        std::cerr << "usage: apps /path/to/match.root" << std::endl;
        return 1;
    }
    TH1::AddDirectory(kFALSE);
    
    TString fileinput = argv[1];
    TFile *file = new TFile(fileinput, "READ");

    TTree *proj = (TTree*)file->Get("T_proj");
    std::vector<int> *cluster_id = new std::vector<int>;
    std::vector<std::vector<int>> *channel_vec = new std::vector<std::vector<int>>;
    std::vector<std::vector<int>> *time_slice_vec = new std::vector<std::vector<int>>;
    std::vector<std::vector<int>> *charge_vec = new std::vector<std::vector<int>>;
    proj->SetBranchAddress("cluster_id", &cluster_id);
    proj->SetBranchAddress("channel", &channel_vec);
    proj->SetBranchAddress("time_slice", &time_slice_vec);
    proj->SetBranchAddress("charge", &charge_vec);

    // Trun info
    TTree* Trun = (TTree*)file->Get("Trun");
    int eventNo;
    int runNo;
    int subRunNo;
    unsigned int triggerBits;
    Trun->SetBranchAddress("eventNo", &eventNo);
    Trun->SetBranchAddress("runNo", &runNo);
    Trun->SetBranchAddress("subRunNo", &subRunNo);
    Trun->SetBranchAddress("triggerBits", &triggerBits);
    Trun->GetEntry(0);

    double lowerwindow = 0;
    double upperwindow = 0;
    if(triggerBits==2048) { lowerwindow = 3; upperwindow = 6; }// bnb
    if(triggerBits==512) { lowerwindow = 3.45; upperwindow = 5.45; } // extbnb

    // flash id and time
    TTree *flash = (TTree*)file->Get("T_flash");
    int flashes_id;
    double flash_time;
    flash->SetBranchAddress("flash_id", &flashes_id);
    flash->SetBranchAddress("time", &flash_time);
    //flash id - time
    std::map<int, double> flash_time_map;
    for(int i=0; i<flash->GetEntries(); i++)
    {
        flash->GetEntry(i);
        if(flash_time<upperwindow && flash_time>lowerwindow) flash_time_map[flashes_id] = flash_time;
    }


    // tpc-flash matching
    TTree *match = (TTree*)file->Get("T_match");
    int tpc_cluster_id;
    int flash_id;
    match->SetBranchAddress("tpc_cluster_id", &tpc_cluster_id);
    match->SetBranchAddress("flash_id", &flash_id);
    
    
    // filtering in-time flash
    std::map<int, int> intime_cluster_id; // cluster_id, flash_id

    // flash time filter
    for(int i=0; i<match->GetEntries(); i++)
    {
        match->GetEntry(i);
        if(flash_id != -1){
            auto it = flash_time_map.find(flash_id);
            if(it != flash_time_map.end()) // within (2, 6) us
            {
                intime_cluster_id[tpc_cluster_id] = flash_id;
            }
        }
    }


    // matched cluster id
    for(auto &it : intime_cluster_id)
    {
        //cout<<"Cluster: "<<it.first<<endl;
    }


    // clusters
    TTree *clusters = (TTree*)file->Get("T_cluster");
    TTree *intime_cluster;
    bool through = true; //any part of cluster touch the boundary within 10 cm
    bool thiscluster = false;
    bool matched = false;

    // channel range
    int uchannel_min = 10000;
    int uchannel_max = 0;
    int vchannel_min = 10000;
    int vchannel_max = 0;
    int wchannel_min = 10000;
    int wchannel_max = 0;
    
    int uchannel0_min = 10000;
    int uchannel0_max = 0;
    int vchannel0_min = 10000;
    int vchannel0_max = 0;
    int wchannel0_min = 10000;
    int wchannel0_max = 0;

    // time range
    int utime_min = 3000;
    int utime_max = 0;
    int vtime_min = 3000;
    int vtime_max = 0;
    int wtime_min = 3000;
    int wtime_max = 0;
    
    int utime0_min = 3000;
    int utime0_max = 0;
    int vtime0_min = 3000;
    int vtime0_max = 0;
    int wtime0_min = 3000;
    int wtime0_max = 0;

    //cout<<"READ: "<<proj->GetEntries()<<endl; // should be 1
    for(int entry=0; entry<proj->GetEntries(); entry++)
    {
        proj->GetEntry(entry);

        // # of charge hits
        //std::cout<<"Cluster: "<<cluster_id->size()<<" Channel: "<<channel_vec->size()<<" Time slices: "<<time_slice_vec->size()<<" Charge: "<<charge_vec->size()<<std::endl;
        
        for(int i=0; i<cluster_id->size(); i++)
        {
            //cout<<i<<endl;
            if(intime_cluster_id.find(cluster_id->at(i)) != intime_cluster_id.end())
            {
            std::vector<int> charge = charge_vec->at(i);
            std::vector<int> time = time_slice_vec->at(i);
            std::vector<int> channel = channel_vec->at(i);          

            cout<<"Cluster: "<<cluster_id->at(i)<<" Time: "<<flash_time_map[intime_cluster_id[cluster_id->at(i)]]<<endl;
            //cout<<"# of Channels: "<<channel.size()<<endl;
            //cout<<"# of time slices: "<<time.size()<<endl;
            //cout<<"# of charges : "<<charge.size()<<endl;

            if(charge.size()>1000){
		    matched = true;
            thiscluster=false;
            for(int p=0; p<charge.size(); p++)
            {
                if(channel.at(p)<2400) {
                    //cout<<"U plane: "<<channel.at(p)<<endl;
                    uchannel_min = uchannel_min<channel.at(p)?uchannel_min:channel.at(p);
                    uchannel_max = uchannel_max>channel.at(p)?uchannel_max:channel.at(p);
                    utime_min = utime_min<time.at(p)?utime_min:time.at(p);
                    utime_max = utime_max>time.at(p)?utime_max:time.at(p);
                }
                if(channel.at(p)>=2400 && channel.at(p)<4800) {
                    //cout<<"V plane: "<<channel.at(p)<<endl;
                    vchannel_min = vchannel_min<channel.at(p)?vchannel_min:channel.at(p);
                    vchannel_max = vchannel_max>channel.at(p)?vchannel_max:channel.at(p);
                    vtime_min = vtime_min<time.at(p)?vtime_min:time.at(p);
                    vtime_max = vtime_max>time.at(p)?vtime_max:time.at(p);
                }
                if(channel.at(p)>=4800) {
                    //cout<<"W plane: "<<channel.at(p)<<endl;
                    wchannel_min = wchannel_min<channel.at(p)?wchannel_min:channel.at(p);
                    wchannel_max = wchannel_max>channel.at(p)?wchannel_max:channel.at(p);
                    wtime_min = wtime_min<time.at(p)?wtime_min:time.at(p);
                    wtime_max = wtime_max>time.at(p)?wtime_max:time.at(p);
                }
            }
 
            // on boundary
            TString String_clusterid;
            String_clusterid.Form("cluster_id==%d", cluster_id->at(i)); 
            gROOT->cd();
            intime_cluster = clusters->CopyTree(String_clusterid.Data());
            //cout<<intime_cluster->GetEntries()<<endl;
            double cx;
            double cy;
            double cz;
            intime_cluster->SetBranchAddress("x", &cx);
            intime_cluster->SetBranchAddress("y", &cy);
            intime_cluster->SetBranchAddress("z", &cz);

            for(int i=0; i<intime_cluster->GetEntries(); i++)
            {
                intime_cluster->GetEntry(i);
                TVector3 v0(cx, cy, cz);
                if(boundary(v0)){
                    thiscluster=true;
                    cout<<"This cluter is on boundary!"<<endl;
                    break;
                }
            }
            // update cluster range if not on boundary
            if(!thiscluster){
                uchannel0_min = uchannel0_min<uchannel_min?uchannel0_min:uchannel_min;
                vchannel0_min = vchannel0_min<vchannel_min?vchannel0_min:vchannel_min;
                wchannel0_min = wchannel0_min<wchannel_min?wchannel0_min:wchannel_min;
                utime0_min = utime0_min<utime_min?utime0_min:utime_min;
                vtime0_min = vtime0_min<vtime_min?vtime0_min:vtime_min;
                wtime0_min = wtime0_min<wtime_min?wtime0_min:wtime_min; 
                uchannel0_max = uchannel0_max>uchannel_max?uchannel0_max:uchannel_max;
                vchannel0_max = vchannel0_max>vchannel_max?vchannel0_max:vchannel_max;
                wchannel0_max = wchannel0_max>wchannel_max?wchannel0_max:wchannel_max;
                utime0_max = utime0_max>utime_max?utime0_max:utime_max;
                vtime0_max = vtime0_max>vtime_max?vtime0_max:vtime_max;
                wtime0_max = wtime0_max>wtime_max?wtime0_max:wtime_max; 
            }

            // reset individual cluster range
            uchannel_min = 10000;
            uchannel_max = 0;
            vchannel_min = 10000;
            vchannel_max = 0;
            wchannel_min = 10000;
            wchannel_max = 0;
            utime_min = 3000;
            utime_max = 0;
            vtime_min = 3000;
            vtime_max = 0;
            wtime_min = 3000;
            wtime_max = 0;

            through = through&&thiscluster; 
            } // cluster size
            else{
                cout<<"Extremely low visible energy event!"<<endl;
            }
            } // intime cluster
        }
    }
    if(through && matched){ 
        cout<<"All clusters at boundary!"<<endl;
    }

    if(!through && matched){
        cout<<"One candidate: "<<endl;
        cout<<"U channel range: "<<uchannel0_min<<" - "<<uchannel0_max<<endl;
        cout<<"V channel range: "<<vchannel0_min<<" -"<<vchannel0_max<<endl;
        cout<<"W channel range: "<<wchannel0_min<<" - "<<wchannel0_max<<endl;
        cout<<"U time tick (2 us) range: "<<utime0_min<<" - "<<utime0_max<<endl;
        cout<<"V time tick (2 us) range: "<<vtime0_min<<" - "<<vtime0_max<<endl;
        cout<<"W time tick (2 us) range: "<<wtime0_min<<" - "<<wtime0_max<<endl;
    }

    file->Close();
    
    return 0;
}
