#include "WireCell2dToy/WCCosmic.h"
#include "TVector3.h"

using namespace WireCell;

WireCell2dToy::WCCosmic::WCCosmic(WireCell2dToy::ToyTrackingSelection& toytrackings)
  : toytrackings(toytrackings)
{
  center.x = 0;
  center.y = 0;
  center.z = 0;
  
  int sum = 0;
  for (int i=0;i!=toytrackings.size();i++){
    ToyTracking *tracking = toytrackings.at(i);
    for (int j=0;j!=tracking->get_good_tracks().size();j++){
      WCTrack *track = tracking->get_good_tracks().at(j);
      for (int k=0;k!=track->get_centerVP_cells().size();k++){
	MergeSpaceCell *mcell = track->get_centerVP_cells().at(k);
	auto it = find(mcells.begin(),mcells.end(),mcell);
	if (it == mcells.end()){
	  mcells.push_back(mcell);
	  center.x += mcell->Get_all_spacecell().size() * mcell->Get_Center().x;
	  center.y += mcell->Get_all_spacecell().size() * mcell->Get_Center().y;
	  center.z += mcell->Get_all_spacecell().size() * mcell->Get_Center().z;
	  sum += mcell->Get_all_spacecell().size();
	}
      }
    }
  }

  //calculate the center
  center.x = center.x/sum;
  center.y = center.y/sum;
  center.z = center.z/sum;

  
  //calculate the angle
  ct = new ClusterTrack(mcells.at(0));
  for (int i=1;i< mcells.size();i++){
    ct->AddMSCell_anyway(mcells.at(i));
  }
  ct->SC_Hough(center);
  theta = ct->Get_Theta();
  phi = ct->Get_Phi();

  Sort();
}



void WireCell2dToy::WCCosmic::Sort(){
  //Sorting 
  std::vector<WireCell2dToy::MSC_Struct> msc_vector;
  for (int i=0;i!=mcells.size();i++){
    msc_vector.push_back(WireCell2dToy::MSC_Struct(cal_pos(mcells.at(i)),mcells.at(i)));
  }
  
  std::sort(msc_vector.begin(),msc_vector.end(),WireCell2dToy::less_than_key());
  
  mcells.clear();
  for (int i=0;i!=msc_vector.size();i++){
    mcells.push_back(msc_vector.at(i).mcell);
  }
}

float WireCell2dToy::WCCosmic::cal_pos(MergeSpaceCell *mcell1){
  TVector3 dir(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  TVector3 dir1(mcell1->Get_Center().x - center.x,
		mcell1->Get_Center().y - center.y,
		mcell1->Get_Center().z - center.z);
  float dis1 = dir.Dot(dir1);
  return dis1;
}


float WireCell2dToy::WCCosmic::cal_dist(WireCell2dToy::WCCosmic *cosmic){
  TVector3 dir1(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  TVector3 dir2(sin(cosmic->get_theta())*cos(cosmic->get_phi()),
		sin(cosmic->get_theta())*sin(cosmic->get_phi()),
		cos(cosmic->get_theta()));
  TVector3 pos1(center.x, center.y, center.z);
  TVector3 pos2(cosmic->get_center().x, cosmic->get_center().y, cosmic->get_center().z);
  
  TVector3 perp1 = dir1.Cross(dir2);
  TVector3 perp2 = pos1 - pos2;
  float abc = perp1.Dot(perp2);

  float dis=1e9;
  if (perp1.Mag() == 0){
    TVector3 perp3 = perp2.Cross(dir1);
    dis = perp3.Mag()/dir1.Mag();
  }else{
    dis = fabs(abc/perp1.Mag());
  }

  return dis;
}

float WireCell2dToy::WCCosmic::cal_costh(WireCell2dToy::WCCosmic *cosmic){
  TVector3 dir1(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  TVector3 dir2(sin(cosmic->get_theta())*cos(cosmic->get_phi()),
		sin(cosmic->get_theta())*sin(cosmic->get_phi()),
		cos(cosmic->get_theta()));
  return dir1.Dot(dir2)/dir1.Mag()/dir2.Mag();
}


WireCell2dToy::WCCosmic::~WCCosmic(){
  delete ct;
}
