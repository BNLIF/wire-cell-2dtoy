#include <vector>
void convert3(){
  #include "2dtoy/src/calib_resp_v1.txt"
  ofstream outfile("run_calib_resp_v1.txt");
  
  for (int i=0;i!=310;i++){ outfile << ref_ele[i] << " ";} outfile << endl;
  for (int i=0;i!=310;i++){ outfile << ref_ele1_ind[i] << " ";} outfile << endl;
  std::cout << calib_ele_chan.size() << " " << calib_ele_chan.at(0).size() << endl;
  for (size_t i=0;i!=calib_ele_chan.size();i++){
    for (size_t j=0;j!=calib_ele_chan.at(i).size();j++){
      outfile << calib_ele_chan.at(i).at(j) << " ";
    }
    outfile << endl;
  }
}
