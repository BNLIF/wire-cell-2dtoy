
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLegend.h"

#include "WireCellRess/GaussProcess.h"



#include <iostream>
using namespace std;
using namespace WireCell;

int main(int argc, char* argv[])
{
  GaussProcess gp(1);
  
   std::vector<double> pars;
  pars.push_back(1);
  pars.push_back(1);
  gp.set_parameters(pars);

  std::vector<std::tuple<double,double,double> > data;
  data.push_back(std::make_tuple(1,  0.3,0.1));
  data.push_back(std::make_tuple(1.5,0.2,0.05));
  data.push_back(std::make_tuple(2,  0,0.2));
  data.push_back(std::make_tuple(3,  0.2, 0.05));
  data.push_back(std::make_tuple(4,  -0.3,0.005));
  data.push_back(std::make_tuple(4.5,-0.4,0.1));
  data.push_back(std::make_tuple(5,  0.2,0.2));
  data.push_back(std::make_tuple(5.5,0.3,0.05));
  data.push_back(std::make_tuple(6.5,-0.5,0.01));
  data.push_back(std::make_tuple(8,  -0.2,0.01));
  data.push_back(std::make_tuple(9,  -0.3,0.1));
  data.push_back(std::make_tuple(10, -0.1,0.01));

  gp.set_measurements(data);

  std::vector<double> vec_x;
  for (size_t i=0;i!=91;i++){
    vec_x.push_back(1+0.1*i);
  }

  std::vector<double> vec_y = gp.cal_conditional_mean(vec_x);
  std::vector<double> vec_var_y = gp.cal_conditional_variance(vec_x);
  //for (size_t i=0;i!=vec_y.size();i++){
    //   std::cout << vec_x.at(i) << " "<< vec_y.at(i) << " " << vec_var_y.at(i) << std::endl;
  //}
  
  TApplication theApp("theApp",&argc,argv);
  theApp.SetReturnFromRun(true);
    
  TCanvas c1("ToyMC","ToyMC",800,600);
  c1.Draw();

  TGraphErrors *g_meas = new TGraphErrors();
  for (size_t i=0;i!=data.size();i++){
    g_meas->SetPoint(i,std::get<0>(data.at(i)), std::get<1>(data.at(i)));
    g_meas->SetPointError(i,0,std::get<2>(data.at(i)));
  }
  g_meas->Draw("A*L");
  
  TGraphErrors *g_pred = new TGraphErrors();
  for (size_t i=0;i!=vec_x.size();i++){
    g_pred->SetPoint(i,vec_x.at(i), vec_y.at(i));
    g_pred->SetPointError(i,0,sqrt(vec_var_y.at(i)));
  }
  g_pred->Draw("AL*");
  g_pred->SetMarkerStyle(21);
  g_pred->SetMarkerColor(2);
  g_pred->SetLineColor(2);
  g_meas->Draw("L*same");
  g_meas->SetMarkerStyle(20);
  
  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->SetFillColor(10);
  le1->AddEntry(g_meas,"Measurements","pl");
  le1->AddEntry(g_pred,"Predictions","pl");
  le1->Draw();
  
  theApp.Run();
}
