#include "WireCellQuickhull/QuickHull.h"
#include "WireCellQuickhull/MathUtils.h"

#include <vector>
#include <set>
#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>

int main(){
  quickhull::QuickHull<float> qh;
  
  std::vector<quickhull::Vector3<float>> pc;
  for (int h=0;h<1000;h++) {
    pc.emplace_back(gRandom->Uniform(-100,100),gRandom->Uniform(-100,100),gRandom->Uniform(-100,100));
  }

  quickhull::ConvexHull<float> hull = qh.getConvexHull(pc,false,false);
  std::cout << hull.getIndexBuffer().size() << " " << hull.getVertexBuffer().size() << " " << pc.size() << std::endl;

  std::set<int> indices;
  
  for (size_t i=0;i!=hull.getIndexBuffer().size();i++){
    indices.insert(hull.getIndexBuffer().at(i));
    //    std::cout << hull.getIndexBuffer().at(i) << std::endl;
  }
  std::cout << indices.size() << std::endl;
  for (size_t i=0;i!=hull.getVertexBuffer().size();i++){
    //  std::cout << hull.getVertexBuffer()[i].x << std::endl;
    //    std::cout << hull.getIndexBuffer().at(i) << std::endl;
  }

  TFile *file = new TFile("temp.root","RECREATE");
  TTree *T = new TTree("T_cluster","T_cluster");
  Double_t x,y,z;
  Int_t cluster_id;
  T->Branch("cluster_id",&cluster_id,"cluster_id/I");
  T->Branch("x",&x,"x/D");
  T->Branch("y",&y,"y/D");
  T->Branch("z",&z,"z/D");
  T->SetDirectory(file);
  cluster_id = 0;
  for (int i=0;i!=1000;i++){
    x = pc.at(i).x;
    y = pc.at(i).y;
    z = pc.at(i).z;
    T->Fill();
  }
  cluster_id = 1;
  for (auto it = indices.begin(); it!=indices.end(); it++){
    x = pc.at(*it).x;
    y = pc.at(*it).y;
    z = pc.at(*it).z;
    T->Fill();
  }

   TTree *Trun = new TTree("Trun","Trun");
  Trun->SetDirectory(file);

  int detector=0, event_no=0, subrun_no=0,run_no=0;
  
  
  Trun->Branch("detector",&detector,"detector/I");

  Trun->Branch("eventNo",&event_no,"eventNo/I");
  Trun->Branch("runNo",&run_no,"runNo/I");
  Trun->Branch("subRunNo",&subrun_no,"runRunNo/I");
  
  Trun->Fill();


  
  
  file->Write();
  file->Close();
  
  
  
  return 0;
}
