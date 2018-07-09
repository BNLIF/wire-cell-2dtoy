void convert1(){
  #include "./2dtoy/src/data_70_2D_11.txt"

  ofstream outfile("run_data_70_2D_11.txt");
  // Double_t u_1D_c_x[5000]
  //Double_t u_1D_c_y[5000]
  //Double_t v_1D_c_x[5000]
  //Double_t v_1D_c_y[5000]
  //Double_t w_1D_c_x[5000]
  //Double_t w_1D_c_y[5000]
  //Double_t u_2D_g_0_x[5000]
  //Double_t u_2D_g_0_y[5000]
  //Double_t v_2D_g_0_x[5000]
  //Double_t v_2D_g_0_y[5000]
  //Double_t w_2D_g_0_x[5000]
  // Double_t w_2D_g_0_y[5000]
  // al the way to 10 ...

  for (int i=0;i!=5000;i++){ outfile << u_1D_c_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << u_1D_c_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_1D_c_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_1D_c_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_1D_c_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_1D_c_y[i] << " "; } outfile << endl;

  for (int i=0;i!=5000;i++){ outfile << u_2D_g_0_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << u_2D_g_0_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_0_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_0_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_0_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_0_y[i] << " "; } outfile << endl;

  for (int i=0;i!=5000;i++){ outfile << u_2D_g_1_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << u_2D_g_1_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_1_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_1_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_1_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_1_y[i] << " "; } outfile << endl;

  for (int i=0;i!=5000;i++){ outfile << u_2D_g_2_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << u_2D_g_2_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_2_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_2_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_2_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_2_y[i] << " "; } outfile << endl;

  for (int i=0;i!=5000;i++){ outfile << u_2D_g_3_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << u_2D_g_3_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_3_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_3_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_3_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_3_y[i] << " "; } outfile << endl;

  for (int i=0;i!=5000;i++){ outfile << u_2D_g_4_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << u_2D_g_4_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_4_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_4_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_4_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_4_y[i] << " "; } outfile << endl;

  for (int i=0;i!=5000;i++){ outfile << u_2D_g_5_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << u_2D_g_5_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_5_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_5_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_5_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_5_y[i] << " "; } outfile << endl;

  for (int i=0;i!=5000;i++){ outfile << u_2D_g_6_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << u_2D_g_6_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_6_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_6_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_6_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_6_y[i] << " "; } outfile << endl;

  for (int i=0;i!=5000;i++){ outfile << u_2D_g_7_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << u_2D_g_7_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_7_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_7_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_7_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_7_y[i] << " "; } outfile << endl;

  for (int i=0;i!=5000;i++){ outfile << u_2D_g_8_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << u_2D_g_8_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_8_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_8_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_8_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_8_y[i] << " "; } outfile << endl;

  for (int i=0;i!=5000;i++){ outfile << u_2D_g_9_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << u_2D_g_9_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_9_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_9_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_9_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_9_y[i] << " "; } outfile << endl;

  for (int i=0;i!=5000;i++){ outfile << u_2D_g_10_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << u_2D_g_10_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_10_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << v_2D_g_10_y[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_10_x[i] << " "; } outfile << endl;
  for (int i=0;i!=5000;i++){ outfile << w_2D_g_10_y[i] << " "; } outfile << endl;
}
