#ifdef HAVE_XDATA               // only build if user has Xdata package

/// Convert "shower3D" file to a Wire Cell Exchange Data ROOT format.

#include "WireCellData/GeomCell.h"

#include "WireCellXdataRoot/Writer.h"
#include "WireCellXdataRoot/Reader.h"
#include "WireCellXdataRoot/CloneHelper.h"
#include "WireCellXdataRoot/WireDB.h"

#include "TFile.h"
#include "TTree.h"

#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <chrono>

using namespace std;
using namespace WireCell;
namespace Xdata = WireCellXdataRoot;


int main(int argc, const char* argv[])
{
  string filename = argv[1];
  Xdata::Reader xreader(filename);
 
  TString filename1 = argv[2];
  TFile *file1 = new TFile(filename1,"RECREATE");
  TTree *T = new TTree("T","T");
  T->SetDirectory(file1);
  Int_t length = 0;
  T->Branch("length",&length,"length/I");

  std::map<int,std::set<int>> signals;

  while (xreader.frame_reader.Next()) {
    Xdata::Frame& f = *xreader.frame;
    Xdata::CloneHelper<Xdata::Deco> decoca(*f.decos);
    for (Int_t i=0;i!=decoca.size();i++){
      Xdata::Deco *deco = decoca.get(i);
      int chid = deco->chanid;
      int time = deco->slice;
      
      if (signals.find(chid)==signals.end()){
	std::set<int> times;
	times.insert(time);
	signals[chid] = times;
      }else{
	signals[chid].insert(time);
      }
    }
    
    for (auto it = signals.begin(); it!= signals.end(); it++){
      // start a new channel
      int channel = (*it).first;
      length = 0;
      int time, prev_time=-1;
      for (auto it1 = (*it).second.begin(); it1 != (*it).second.end(); it1++){
	time = (*it1);

	//std::cout << channel << " " << time << std::endl;

	if (prev_time == -1){
	  length ++;
	}else{
	  if (time - prev_time == 1){
	    length ++;
	  }else{
	    T->Fill();
	    //std::cout << "Xin: " << length << std::endl;
	    length = 1;
	  }
	}
	prev_time = time;
      }
      
      if (length >0)
	T->Fill();
	//std::cout << "Xin: " << length << std::endl;

    }
    // cerr << "Frame: " << f.ident << " trig=" << f.trigid << " geom=" << f.geomid
    // 	 << " #decos=" << decoca.size() << endl;

  }
  file1->Write();
  file1->Close();

  return 0;

}
    
#else  // if user doesn't have Xdata package
int main() { return 0; }
#endif

