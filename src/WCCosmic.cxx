#include "WireCell2dToy/WCCosmic.h"
#include "TVector3.h"

using namespace WireCell;

WireCell2dToy::WCCosmic::WCCosmic(WireCell2dToy::ToyTrackingSelection& toytrackings)
  : toytrackings(toytrackings)
{
  center.x = 0;
  center.y = 0;
  center.z = 0;

  cosmic_flag = false;
  
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

  ct = 0;

  //std::cout << mcells.size() << std::endl;
  if (sum >1){
    //calculate the center
    center.x = center.x/sum;
    center.y = center.y/sum;
    center.z = center.z/sum;
    
    //std::cout << "Hough " << std::endl;
    //calculate the angle
    ct = new ClusterTrack(mcells.at(0));
    for (int i=1;i< mcells.size();i++){
      ct->AddMSCell_anyway(mcells.at(i));
    }
    ct->SC_Hough(center);
    theta = ct->Get_Theta();
    phi = ct->Get_Phi();
    
    //std::cout << "Sort " << std::endl; 
    Sort();
    //std::cout << "Fill Points" << std::endl;
    fill_points();
    //std::cout << mcells.size() << " " << cal_pos(mcells.front()) << " " << cal_pos(mcells.back()) << std::endl;
    judge_cosmic();
  }
}

void WireCell2dToy::WCCosmic::judge_cosmic(){
  cosmic_flag = false;
  
  // if anythng is outside x
  for (int i=0;i!=points.size();i++){
    if (points.at(i).x < 1*units::cm || points.at(i).x > (256 - 1) * units::cm){
      // std::cout << points.at(i) << std::endl;
      cosmic_flag = true;
      break;
    }
  }

  if (!cosmic_flag){
    bool front_boundary = false;
    bool back_boundary = false;

    if (fabs(points.front().y) > 233 * units::cm/2. - 10*units::cm) front_boundary = true;
    if (fabs(points.back().y) > 233 * units::cm/2. - 10*units::cm) back_boundary = true;
    if (points.front().z < 10 * units::cm || points.front().z > 10.3692 *100 * units::cm - 10 * units::cm) front_boundary = true;
    if (points.back().z < 10 * units::cm || points.back().z > 10.3692 *100 * units::cm - 10 * units::cm) back_boundary = true;

    if (front_boundary && back_boundary) cosmic_flag = true;

    if ((front_boundary || back_boundary) && (!cosmic_flag)){
      // look at the angle
      TVector3 dir1(0,1,0);
      TVector3 dir3(0,0,1);
      TVector3 dir2(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
      //     std::cout << dir1.Dot(dir2)/dir2.Mag() << std::endl;
      if (fabs(dir1.Dot(dir2)/dir2.Mag()) > 0.9) cosmic_flag = true;
      TVector3 dir4;
      Point vertex_p;
      if (front_boundary){
	dir4.SetXYZ(points.front().x - points.back().x, points.front().y - points.back().y, points.front().z - points.back().z);

	vertex_p.x = points.back().x;
	vertex_p.y = points.back().y;
	vertex_p.z = points.back().z;

	if (dir4.Dot(dir2) < 0 ) 
	  dir2.SetXYZ(-dir2.x(), -dir2.y(), -dir2.z());
      }else{
	dir4.SetXYZ(points.back().x - points.front().x, points.back().y - points.front().y, points.back().z - points.front().z);

	vertex_p.x = points.front().x;
	vertex_p.y = points.front().y;
	vertex_p.z = points.front().z;

	if (dir4.Dot(dir2) < 0)
	  dir2.SetXYZ(-dir2.x(), -dir2.y(), -dir2.z());
      }
      
      if (dir3.Dot(dir2)/dir2.Mag() < -0.5) cosmic_flag = true;
      //std::cout << dir3.Dot(dir2)/dir2.Mag() << " " << front_boundary << " " << back_boundary << std::endl;
      
      if (!cosmic_flag){
	// count how many good tracks were in ... 
	int ntrack = 0;
	for (int i=0;i!=toytrackings.size();i++){
	  WireCell2dToy::ToyTracking *toytracking = toytrackings.at(i);
	  for (int j=0;j!=toytracking->get_vertices().size();j++){
	    WCVertex *vertex = toytracking->get_vertices().at(j);
	    if (vertex->get_ntracks() >=2){
	      for (int k=0;k!=vertex->get_ntracks();k++){
		WCTrack *track = vertex->get_tracks().at(k);
		std::vector<float> pp = track->get_direction();
		TVector3 dir5(pp.at(0),pp.at(1),pp.at(2));
		float dis = sqrt(pow(vertex_p.x - vertex->Center().x,2) + pow(vertex_p.y - vertex->Center().y,2) + pow(vertex_p.z - vertex->Center().z,2));
		
		if (dis < 5*units::cm){
		  if (fabs(dir5.Dot(dir2)/dir5.Mag()/dir2.Mag()) < 0.9)
		    ntrack ++;
		}else{
		  if (track->get_range()/units::cm > 5 && fabs(dir5.Dot(dir2)/dir5.Mag()/dir2.Mag()) < 0.9)
		    ntrack++;
		}
		
		//		if (fabs(dir5.Dot(dir2)/dir5.Mag()/dir2.Mag()) < 0.9)
		//ntrack ++;
		//		std::cout << i << " " << j << " " << dir5.Dot(dir2)/dir5.Mag()/dir2.Mag() << " " << track->get_range()/units::cm<< std::endl;
	      }
	    }
	    // std::cout << vertex->get_ntracks() << std::endl;
	    // ntrack += toytrackings.at(i)->get_good_tracks().size();
	    //	for (int j=0;j!=toytrackings.at(i)->get_good_tracks().size()
	  }
	}

	if (ntrack == 0) cosmic_flag = true;
	
      }
      //std::cout << theta << " " << phi << " " << ntrack << std::endl;
    }
    
  }
  
}


void WireCell2dToy::WCCosmic::Add(WCCosmic *cosmic){
  toytrackings.insert(toytrackings.end(), cosmic->get_trackings().begin(),
		      cosmic->get_trackings().end());
  mcells.insert(mcells.end(),cosmic->get_mcells().begin(),cosmic->get_mcells().end());
  center.x = 0;
  center.y = 0;
  center.z = 0;
  
  int sum = 0;
  for (int i=0;i!=mcells.size();i++){
    MergeSpaceCell *mcell = mcells.at(i);
    center.x += mcell->Get_all_spacecell().size() * mcell->Get_Center().x;
    center.y += mcell->Get_all_spacecell().size() * mcell->Get_Center().y;
    center.z += mcell->Get_all_spacecell().size() * mcell->Get_Center().z;
    sum += mcell->Get_all_spacecell().size();
  }
  
  if (sum >1){
    //calculate the center
    center.x = center.x/sum;
    center.y = center.y/sum;
    center.z = center.z/sum;

    delete ct;
    ct = new ClusterTrack(mcells.at(0));
    for (int i=1;i< mcells.size();i++){
      ct->AddMSCell_anyway(mcells.at(i));
    }
    ct->SC_Hough(center);
    theta = ct->Get_Theta();
    phi = ct->Get_Phi();
    
    //std::cout << "Sort " << std::endl; 
    Sort();
    //std::cout << "Fill Points" << std::endl;
    fill_points();
    judge_cosmic();
  }

}


void WireCell2dToy::WCCosmic::fill_points(){
  points.clear();
  
  //find the smallest point and largest point
  float min_dis = 20*units::m;
  Point min_point;
  float max_dis = -20*units::m;
  Point max_point;

  // std::cout << "qx: q "  << mcells.size() << std::endl;
  for (int i=0;i!=mcells.size();i++){
    MergeSpaceCell *mcell = mcells.at(i);
    if (i < 5){
      for (int j=0;j!=mcell->Get_all_spacecell().size();j++){
	SpaceCell *cell = mcell->Get_all_spacecell().at(j);
	float dis = cal_pos(cell);
	if (dis < min_dis){
	  min_dis = dis;
	  min_point.x = cell->x();
	  min_point.y = cell->y();
	  min_point.z = cell->z();
	}
      }
    }
    
    if (i +5 >= mcells.size()){
      //std::cout << "qx: q " << i << " " << mcells.size() << std::endl;
      for (int j=0;j!=mcell->Get_all_spacecell().size();j++){
	SpaceCell *cell = mcell->Get_all_spacecell().at(j);
	float dis = cal_pos(cell);
	//std::cout << "qx: " << dis/units::cm <<std::endl;
	if (dis > max_dis){
	  max_dis = dis;
	  max_point.x = cell->x();
	  max_point.y = cell->y();
	  max_point.z = cell->z();
	}
      }
    }
  }
  

  
  
  //at the beginning
  points.push_back(min_point);
  
  int nbin = (max_dis - min_dis)/(5*units::cm);
  
  if (nbin >0){
    // calculate middle points
    TH1F *hx = new TH1F("hx","hx",nbin,min_dis,max_dis);
    TH1F *hy = new TH1F("hy","hy",nbin,min_dis,max_dis);
    TH1F *hz = new TH1F("hz","hz",nbin,min_dis,max_dis);
    TH1F *hc = new TH1F("hc","hc",nbin,min_dis,max_dis);
    
    for (int i=0;i!=mcells.size();i++){
      MergeSpaceCell *mcell = mcells.at(i);
      for (int j=0;j!=mcell->Get_all_spacecell().size();j++){
	SpaceCell *cell = mcell->Get_all_spacecell().at(j);
	float dis = cal_pos(cell);
	hx->Fill(dis,cell->x());
	hy->Fill(dis,cell->y());
	hz->Fill(dis,cell->z());
	hc->Fill(dis,1);
      }
    }
    
    Point p;
    for (int i=0;i!=nbin;i++){
      if (hc->GetBinContent(i+1)!=0){
	p.x = hx->GetBinContent(i+1)/hc->GetBinContent(i+1);
	p.y = hy->GetBinContent(i+1)/hc->GetBinContent(i+1);
	p.z = hz->GetBinContent(i+1)/hc->GetBinContent(i+1);
	points.push_back(p);
      }
    }
    delete hx;
    delete hy;
    delete hz;
    delete hc;
  }

  //at the end
  points.push_back(max_point);

  // std::cout << "qx: " << min_dis << " " << min_point.x/units::cm << " " << min_point.y/units::cm << " " << min_point.z/units::cm << " " 
  // 	    << max_dis << " " << max_point.x/units::cm << " " << max_point.y/units::cm << " " << max_point.z/units::cm << std::endl;
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

float WireCell2dToy::WCCosmic::cal_pos(SpaceCell *cell){
  TVector3 dir(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  TVector3 dir1(cell->x() - center.x,
		cell->y() - center.y,
		cell->z() - center.z);
  float dis1 = dir.Dot(dir1);
  return dis1;
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
  if (ct != 0)
    delete ct;
}
