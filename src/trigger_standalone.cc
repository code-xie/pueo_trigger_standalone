#include <iostream>
#include "TGraph.h"
#include "TGraph2D.h"
#include "TTree.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TRootCanvas.h"
#include "TStyle.h"
#include "TApplication.h"
#include <vector>
#include <math.h>
#include <complex.h>
#include <iostream>
#include <tuple>
#include <string>
#include <fftw3.h>
#include <algorithm>
#include <random>
#include "TAxis.h" 
#include "TH1D.h" 
#include "TH1.h"
#include "FTPair.h"
#include "TFitResult.h"
#include "TF1.h"
#include "trigger.h"


double e_field(double theta_offcone, double freq, double normalisation ) {
  if (freq < 1 || freq > 2000) {
    return 0;
  } else {
    double theta_cerenkov = 41; //degrees
    double cone_width = 2.2 * 1000 / freq; //degrees
    double theta_view =  theta_cerenkov - theta_offcone;
    double ice_attenu = exp(-500/(440 + (-120)/(1000-200)*(freq-200)));
    return ice_attenu * normalisation * sin(theta_view/180*M_PI) / sin(theta_cerenkov/180*M_PI)   * exp(-pow(theta_offcone/cone_width, 2)) * 2.53E-7 * freq / 1150 / (1+pow(freq / 1150,1.44));
  }
}


void visualiseFields() {
  std::vector<double> off_cone_angles = {0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 6.0};
  std::vector<double> freq_list;
  for(int i=0; i<2001; i+=10){
    freq_list.push_back(i);
  }
  std::vector<std::vector<double>> e_field_array_all;
  for (auto const &theta_offcone : off_cone_angles) {
    std::vector<double> e_field_array;
    for (auto const &freq : freq_list) {
      e_field_array.push_back(e_field(theta_offcone, freq, 10));
    }
    double max_ele = *std::max_element(e_field_array.begin(), e_field_array.end());
    e_field_array_all.push_back(e_field_array);
  }
 
  
  TCanvas *c_fields = new TCanvas("c_fields","electric field with attenuation",500,500,600,400);
  TMultiGraph *mg = new TMultiGraph();
  //std::cout << "Digitised data:" << "\n";

  int off_cone_angle_count = 0;
  for (auto const &theta_offcone : off_cone_angles) {

    TGraph *gr3 = new TGraph (freq_list.size(), &freq_list[0], &e_field_array_all.at(off_cone_angle_count)[0]);
    gr3->SetTitle(std::to_string(theta_offcone).c_str());
    mg->Add(gr3);
    off_cone_angle_count++;
  }
  c_fields->cd();
  mg->Draw("A pmc plc");
  c_fields->BuildLegend();
}


double * getImpulse(double theta_offcone_deg, int expected_size) {

  //get PUEO impulse response from csv (based on A3 IR with filter limiting to 0.3-1.2Ghz)
   

  TTree *t = new TTree("t", "tree from PUEO_signal_3ghz.csv");
  t->ReadFile("./PUEO_signal_3ghz.csv", "t_ns/D:signal/D");
  //t->Draw("signal : t_ns");

  Double_t IR_y;
  Double_t IR_time;
  t->SetBranchAddress("signal",&IR_y);
  t->SetBranchAddress("t_ns",&IR_time);

  int IR_size = t->GetEntries();
  int IR_fft_size = IR_size/2 + 1; 

  if (IR_size != expected_size) {
    throw std::invalid_argument( "Impulse response file not the expected number of values");
  }

  double IR_in [IR_size]; 
  std::complex<double> IR_fft [IR_fft_size];
  double  fft_freq [IR_fft_size];

  //proccess in fourier domain
  //load IR into input for fftw
  for(int i=0; i<IR_size; i++){
    //std::cout<<"time (ns): "<<IR_time<< "\t" <<"signal: "<<signal_y<<std::endl;
    t->GetEntry(i);
    IR_in[i] = (IR_y);
  }

  //frequencies for plotting in fourier domain
  for(int i=0; i<IR_fft_size; i++){
    fft_freq[i] = (i * 2.94912E9 / IR_size);
  }
 
  //do fft
  fftw_plan p = fftw_plan_dft_r2c_1d(IR_size, &IR_in[0] , reinterpret_cast<fftw_complex*>(&IR_fft[0]) , FFTW_ESTIMATE);
  fftw_execute(p);

  fftw_destroy_plan(p);
  
  //plot norm for impulse response FFT
  double out_norm[IR_fft_size];
  for(int i = 0; i < IR_fft_size; i++) {
    out_norm[i] = std::norm(IR_fft[i]);
    //std::cout << ((double)out_norm[i]) << std::endl;
  }

  //TCanvas *c_sig_ft = new TCanvas("c_sig_ft","fourier transformed signal",500,500,600,400);
  //TGraph *gr_sig_ft = new TGraph (IR_fft_size, &fft_freq[0], &out_norm[0]);
  //c_sig_ft->cd();
  //gr_sig_ft->Draw();

  //combine with electric field
  std::complex<double> signal_before_norm_fft[IR_fft_size];
  for (int i=0; i < IR_fft_size; i++) {
   signal_before_norm_fft[i] = (exp(M_PI/2) * e_field(theta_offcone_deg, fft_freq[i] / 1E6, 1E7) * IR_fft[i]);
   //std::cout <<  signal_before_norm_fft[i] << std::endl;
  }

  //irfft
  double signal_before_norm [IR_size]; 
  p = fftw_plan_dft_c2r_1d(IR_size, reinterpret_cast<fftw_complex*>(&signal_before_norm_fft[0]) , &signal_before_norm[0] , FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  
  //normalise to peak-to-peak/2 = 1
  double peak_to_peak =  *std::max_element(signal_before_norm, signal_before_norm + IR_size) - *std::min_element(signal_before_norm, signal_before_norm + IR_size);
  double * signal_normalised = new double [expected_size];  
  for(int i=0; i<IR_size; i++){
    signal_normalised[i] = (signal_before_norm[i] / peak_to_peak * 2);
    //std::cout <<  signal_normalised[i] << std::endl;
  }


  //plot signal before norm 
  //times for plotting 
  double sig_times [IR_size];
  for(int i=0; i<IR_size; i++){
    t->GetEntry(i);
    sig_times[i] = IR_time;
  }

  //TCanvas *c_sig = new TCanvas("c_sig","signal",500,500,600,400);
  //TGraph *gr_sig = new TGraph (IR_size, &sig_times[0], &signal_normalised[0]);
  //c_sig->cd();
  //gr_sig->Draw();


  
  //fftw_free(IR_in); 
  //fftw_free(out); 

  delete t;
  delete[] signal_normalised;

  return signal_normalised;
}


double * generate_noise_flat(double norm, int size, double min_freq_Mhz, double max_freq_Mhz, std::mt19937& gen) {
  int freq_bins = ceil(size / 2.0) + 1;
  double sample_freq_Mhz = 2949.12;
  double nyq_freq_Mhz = sample_freq_Mhz /2;
  double scaling = norm / sqrt(size);

  std::complex<double> * freq_output = new std::complex<double>[freq_bins]; 

  //std::random_device rd;  // Will be used to obtain a seed for the random number engine
  //std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0.0, 1.0);

  for (int f = 0; f < freq_bins; f++) {
    if ((double)f / (double)freq_bins > min_freq_Mhz / nyq_freq_Mhz && (double)f / (double)freq_bins < max_freq_Mhz / nyq_freq_Mhz) {
      freq_output[f] = scaling * sqrt(-2 * log(dis(gen))) * exp(std::complex<double>(0, 1) * std::complex<double>(dis(gen) * 2 * M_PI, 0)  );
    } else {
       freq_output[f] = std::complex<double>(0, 0);
    }
  }

  int time_bins = (freq_bins - 1) * 2;

  double * noise_output = new double[time_bins]; 
  fftw_plan p = fftw_plan_dft_c2r_1d(time_bins, reinterpret_cast<fftw_complex*>(&freq_output[0]) , &noise_output[0] , FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);

  double *  noise_output_xaxis = new double[time_bins];
  for (int i = 0; i < time_bins; i++) {
    noise_output_xaxis[i] = i;
  }


  //TCanvas *c_noise = new TCanvas("c_noise","c_noise",500,500,600,400);
  //TGraph *gr_noise = new TGraph (time_bins,&noise_output_xaxis[0], &noise_output[0]);
  //c_noise->cd();
  //gr_noise->Draw();

  delete[] noise_output_xaxis;
  delete[] freq_output;
  return noise_output;
}

TGraph * combine_signal_noise(double signal_short[], double noise[], int signal_size, int signal_loc, int size, double SNR, double multiplier, double freq_Ghz) {

  int size_fft = ceil(size / 2) + 1;
  if (signal_size + signal_loc > size) {
    throw std::invalid_argument( "signal and location too large to fit within total profile");
  }


  double * signal = new double[size];
  std::complex<double> * sig_fft = new std::complex<double>[size_fft];
  std::complex<double> * noise_fft = new std::complex<double>[size_fft];
  std::complex<double> * sig_and_noise_fft = new std::complex<double>[size_fft];
  double * sig_and_noise = new double[size];
  double * time_bins = new double[size];

  for (int i = 0; i < size; i++) {
    if (i < signal_loc || i >= signal_loc + signal_size) {
      signal[i] = 0;
    } else {
      signal[i] = signal_short[i-signal_loc];
    }
  }

  for (int i = 0; i < size; i++) {
    time_bins[i] = i * 1.0 / freq_Ghz / 1E9;
  }

  double * freq_bins = new double[size_fft];
  for (int i = 0; i < size_fft; i++) {
    freq_bins[i] = i;
  }

  //check RMS of noise
  //double sum_of_noise_sq = .0;
  //for(int i = 0; i < size; ++i) {
  //  sum_of_noise_sq += pow(noise[i],2);
  //}  
  //double noise_rms = pow(sum_of_noise_sq/size, 1/2) * 1000000;
  //std::cout << "RMS of noise " << noise_rms;

  //check peak to peak of signal
  //double max_of_sig = .0;
  //double min_of_sig = .0;
  //for(int i = 0; i < size; ++i) {
  //  if (signal[i] > max_of_sig) {
  //    max_of_sig = signal[i];
  //  }
  //  if (signal[i] < min_of_sig) {
  //    min_of_sig = signal[i];
  //  }
  //}  
  //double ptp_2 = (max_of_sig - min_of_sig)/2;
  //std::cout << "signal peak to peak /2 = " << ptp_2;


  //do fft on signal and noise to add in momentum space
  fftw_plan p = fftw_plan_dft_r2c_1d(size, &signal[0] , reinterpret_cast<fftw_complex*>(&sig_fft[0]) , FFTW_ESTIMATE);
  fftw_execute(p);

  fftw_destroy_plan(p);

  p = fftw_plan_dft_r2c_1d(size, &noise[0] , reinterpret_cast<fftw_complex*>(&noise_fft[0]) , FFTW_ESTIMATE);
  fftw_execute(p);

  fftw_destroy_plan(p);

  //plot norm
  //double * noise_fft_norm = new double[size_fft];
  //for (int i = 0; i < size_fft; i++) {
  //  noise_fft_norm[i] = std::norm(noise_fft[i]);
  //}
  //TCanvas *c_sig_test = new TCanvas("c_sig_test","c_sig_test",500,500,600,400);
  //c_sig_test->cd();
  //TGraph *sig_tgraph = new TGraph (size_fft,&freq_bins[0], &noise_fft_norm[0]);
  //sig_tgraph -> Draw();
  
  for (int i = 0; i < size_fft; i++) {
    sig_and_noise_fft[i]  = (multiplier / size) * (SNR * sig_fft[i] + noise_fft[i]) ; //FFTW does unnormalised DFT, introducing scaling by n 
  }

  p = fftw_plan_dft_c2r_1d(size,reinterpret_cast<fftw_complex*>(&sig_and_noise_fft[0]), &sig_and_noise[0] ,FFTW_ESTIMATE);
  fftw_execute(p);

  fftw_destroy_plan(p);
  
  TGraph *sig_noise_tgraph = new TGraph (size,&time_bins[0], &sig_and_noise[0]);
  
  //TGraph *noise_tgraph = new TGraph (size,&time_bins[0], &noise[0]);
  
  //TCanvas *c_signoise_test = new TCanvas("c_signoise_test","c_signoise_test",500,500,600,400);
  //c_signoise_test->cd();
  //sig_noise_tgraph -> Draw();

  delete[] signal;
  delete[] sig_fft;
  delete[] noise_fft;
  delete[] sig_and_noise_fft;
  delete[] sig_and_noise;
  delete[] time_bins;
  delete[] freq_bins;
  
  return sig_noise_tgraph;
  
  //return noise_tgraph;
}


TGraph * combine_noise_zero_signal(double noise[], int size, double multiplier, double freq_Ghz) {

  int size_fft = ceil(size / 2) + 1;
  std::complex<double> * noise_fft = new std::complex<double>[size_fft];
  std::complex<double> * sig_and_noise_fft = new std::complex<double>[size_fft];
  double * sig_and_noise = new double[size];
  double * time_bins = new double[size];

  for (int i = 0; i < size; i++) {
    time_bins[i] = i * 1.0 / freq_Ghz / 1E9;
  }

  double * freq_bins = new double[size_fft];
  for (int i = 0; i < size_fft; i++) {
    freq_bins[i] = i;
  }


  fftw_plan p = fftw_plan_dft_r2c_1d(size, &noise[0] , reinterpret_cast<fftw_complex*>(&noise_fft[0]) , FFTW_ESTIMATE);
  fftw_execute(p);

  fftw_destroy_plan(p);

  
  for (int i = 0; i < size_fft; i++) {
    sig_and_noise_fft[i]  = (multiplier / size) * (noise_fft[i]) ; //FFTW does unnormalised DFT, introducing scaling by n 
  }

  p = fftw_plan_dft_c2r_1d(size,reinterpret_cast<fftw_complex*>(&sig_and_noise_fft[0]), &sig_and_noise[0] ,FFTW_ESTIMATE);
  fftw_execute(p);

  fftw_destroy_plan(p);
  
  TGraph *sig_noise_tgraph = new TGraph (size,&time_bins[0], &sig_and_noise[0]);
  
  //TGraph *noise_tgraph = new TGraph (size,&time_bins[0], &noise[0]);
  
  //TCanvas *c_signoise_test = new TCanvas("c_signoise_test","c_signoise_test",500,500,600,400);
  //c_signoise_test->cd();
  //sig_noise_tgraph -> Draw();

  delete[] noise_fft;
  delete[] sig_and_noise_fft;
  delete[] sig_and_noise;
  delete[] time_bins;
  delete[] freq_bins;
  
  return sig_noise_tgraph;
  
  //return noise_tgraph;
}

std::vector<nicemc::FTPair> signal_gen(std::mt19937& gen, double snr, int n_samples, bool noise_only) {

  //choosing source
  double theta = 10.0 / 180.0 * M_PI;
  double phi   = 10.0 / 180.0 * M_PI;


  //initialising antenna data
  double x[n_samples], y[pueoSim::pueoTrigger::n_ant_L2][n_samples];
  double std = 6;
  double pos = n_samples /2;
  double peak = 16;


  double w1 = 0.6; //horizontal separation in m between 2 adjacent azimuthal sectors
  double w = 3 * 0.6; //horizontal separation in m between azimuthal sectors furthest apart in L2 sector
  double h1 = 0.74; //vertical separation between adjacent antennas for lower 3 antennas in a //phi sector
  double h2 = 3.0; //vertical separation between adjacent antennas between top two antennas in //a given phi sector
  double h = 2 * h1 + h2;
  double c = 3e8;
  double delta_t = 1/(2.6e9); //time between samples, for 2.6Ghz sample frequency

 /****
   * label antennas in L2 sector as below
   * top
   * 0  4   8   12
   * 1  5   9   13
   * 2  6   10  14
   * 3  7   11  15
   * bottom
   **/

  TGraph * received_signal_pre_interp = new TGraph[16];

  double sample_padding = 64;
  double n_sample_generated = n_samples+sample_padding * 2; //size of noise generated must be smaller than size of combined signal combined for trigger

  //std::random_device rd;  // Will be used to obtain a seed for the random number engine
  //std::mt19937 gen1(rd()); // Standard mersenne_twister_engine seeded with rd()


  for (int i = 0 ; i < 16; i++) {
    double * noise = generate_noise_flat(0.8343, n_sample_generated, 240, 1300, gen); //generate noise with RMS=1 and with frequencies in this range
    //double sum = 0;
    //for(int j = 0; j < n_sample_generated; j++) {
    //  sum += pow(noise[j], 2);
    //}
    //std::cout << "RMS: " << sqrt(sum/n_sample_generated) << std::endl;
    
    //std::cout << "peak-to-peak/2: " << (*std::max_element(signal, signal+128) - *std::min_element(signal, signal+128)) / 2 << std::endl;
    
    TGraph * temp_tgraph;
    
    if (noise_only) {
      temp_tgraph = combine_noise_zero_signal(noise, n_sample_generated, 1., 2.94912);
    } else {
      double * signal = getImpulse(1, 128);
      temp_tgraph = combine_signal_noise(signal, noise, 128, 128, n_sample_generated, snr, 1., 2.94912); //double signal_short[], double noise[], int signal_size, int signal_loc, int size, double SNR, double multiplier, double freq_Ghz
    
    } 
    received_signal_pre_interp[i] = *temp_tgraph;
    delete temp_tgraph;
    
    delete [] noise;
  }

  for (int i=0; i<n_samples; i++) {
    
    x[i] = i;
    y[0][i] =  received_signal_pre_interp[0]  .Eval((i+ sample_padding + 0    * sin (theta) / (c * delta_t) - 0    * sin (phi)  / (c * delta_t) )/2.94912E9  );
    y[1][i] =  received_signal_pre_interp[1]  .Eval((i+ sample_padding + h2   * sin (theta) / (c * delta_t) - 0    * sin (phi)  / (c * delta_t) )/2.94912E9  );
    y[2][i] =  received_signal_pre_interp[2]  .Eval((i+ sample_padding + h1+h2* sin (theta) / (c * delta_t) - 0    * sin (phi)  / (c * delta_t) )/2.94912E9  );
    y[3][i] =  received_signal_pre_interp[3]  .Eval((i+ sample_padding + h    * sin (theta) / (c * delta_t) - 0    * sin (phi)  / (c * delta_t) )/2.94912E9  );

    y[4][i] =  received_signal_pre_interp[4]  .Eval((i+ sample_padding + 0    * sin (theta) / (c * delta_t) - w1   * sin (phi)  / (c * delta_t) )/2.94912E9  );
    y[5][i] =  received_signal_pre_interp[5]  .Eval((i+ sample_padding + h2   * sin (theta) / (c * delta_t) - w1   * sin (phi)  / (c * delta_t) )/2.94912E9  );
    y[6][i] =  received_signal_pre_interp[6]  .Eval((i+ sample_padding + h1+h2* sin (theta) / (c * delta_t) - w1   * sin (phi)  / (c * delta_t) )/2.94912E9  );
    y[7][i] =  received_signal_pre_interp[7]  .Eval((i+ sample_padding + h    * sin (theta) / (c * delta_t) - w1   * sin (phi)  / (c * delta_t) )/2.94912E9  );

    y[8][i] =  received_signal_pre_interp[8]  .Eval((i+ sample_padding + 0    * sin (theta) / (c * delta_t) - 2*w1 * sin (phi)  / (c * delta_t) )/2.94912E9  );
    y[9][i] =  received_signal_pre_interp[9]  .Eval((i+ sample_padding + h2   * sin (theta) / (c * delta_t) - 2*w1 * sin (phi)  / (c * delta_t) )/2.94912E9  );
    y[10][i] = received_signal_pre_interp[10] .Eval((i+ sample_padding + h1+h2* sin (theta) / (c * delta_t) - 2*w1 * sin (phi)  / (c * delta_t) )/2.94912E9  );
    y[11][i] = received_signal_pre_interp[11] .Eval((i+ sample_padding + h    * sin (theta) / (c * delta_t) - 2*w1 * sin (phi)  / (c * delta_t) )/2.94912E9  );

    y[12][i] = received_signal_pre_interp[12] .Eval((i+ sample_padding + 0    * sin (theta) / (c * delta_t) - 3*w1 * sin (phi)  / (c * delta_t) )/2.94912E9  );
    y[13][i] = received_signal_pre_interp[13] .Eval((i+ sample_padding + h2   * sin (theta) / (c * delta_t) - 3*w1 * sin (phi)  / (c * delta_t) )/2.94912E9  );
    y[14][i] = received_signal_pre_interp[14] .Eval((i+ sample_padding + h1+h2* sin (theta) / (c * delta_t) - 3*w1 * sin (phi)  / (c * delta_t) )/2.94912E9  );
    y[15][i] = received_signal_pre_interp[15] .Eval((i+ sample_padding + h    * sin (theta) / (c * delta_t) - 3*w1 * sin (phi)  / (c * delta_t) )/2.94912E9  );
  }

  delete[] received_signal_pre_interp;


  std::vector<nicemc::FTPair> generated_signals ;

  for (int i_ant=0; i_ant<pueoSim::pueoTrigger::n_ant_L2; i_ant++) {
    TGraph gr = TGraph (n_samples, x, y[i_ant]);
    nicemc::FTPair ftp(gr);
    generated_signals.push_back(ftp);
  }


  return generated_signals;
}

void visualiseTrigger(int argc, char **argv, int L1_threshold, int L2_threshold, double snr) {
  TApplication app("app", &argc, argv); //this allows interactive plots for cmake application

  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

  std::vector<nicemc::FTPair> internally_generated_signal = signal_gen(gen, snr, 512, false);
  pueoSim::pueoTrigger * ptrigger = new pueoSim::pueoTrigger(2.56E9);
  ptrigger->newSignal(internally_generated_signal);
  
  std::cout.precision(5);


  //evaluate_threshold(ptrigger, gen, true); //ptrigger object, random number generator, is_L2
  //evaluate_threshold(ptrigger, gen, false); //ptrigger object, random number generator, is_L2

  //For testing - different off axis signals
  //visualiseFields();
  
  //For testing - getting signal
  //double * test_signal = getImpulse(1, 128); //get signal at given off cone angle (in degrees), and expected number of values in impulse response file

  //For testing - getting noise
  //double * test_noise = generate_noise_flat(0.8343, ptrigger.n_samples, 240, 1300);
  //TGraph signal_noise_test = *combine_signal_noise(test_signal, test_noise, 128, 64, n_samples, 10, 1, 2.94912);



  //create digitised data
  ptrigger->digitize(4);

  //Plot digitised signals data
  TCanvas *c_digi = new TCanvas("c_digi","Discrete antenna data",500,500,600,400);
  TMultiGraph *mg = new TMultiGraph();
  //std::cout << "Digitised data:" << "\n";
  for (std::vector<TGraph>::iterator it=ptrigger->signals_discrete.begin();it!=ptrigger->signals_discrete.end(); ++it) {
    TGraph * gr3 = new TGraph(*it);
    gr3->SetTitle(std::to_string(it - ptrigger->signals_discrete.begin()).c_str());
    gr3->SetLineColor(it - ptrigger->signals_discrete.begin() + 1);
    mg->Add(gr3);
  }
  c_digi->cd();
  mg->Draw("A pmc plc");
  c_digi->BuildLegend();


  //do L1 trigger
  ptrigger->l1Trigger(8, 16, L1_threshold, 64); //args: ptrigger, step, window, threshold, edge size

  //do L2 trigger
  ptrigger->l2Trigger(8, 16, L2_threshold, 64); //args: ptrigger, step, window, threshold, edge size

  //print L1 triggers
  //std::cout << "L1 Triggers" <<"\n";
  //std::cout << "beam\t" <<"triggered\t" << "Max coherent sum\t" << "Max sum location" <<"\n";
  //for (int i=0; i<ptrigger.L1_triggers.size(); i++) {
  //  std::cout << i << "\t\t" << ptrigger.L1_triggers[i] << "\t\t" << ptrigger.L1_max_value[i]  <<"\n";
  //}

  ////print L2 triggers
  //std::cout << "L2 Triggers" <<"\n";
  //std::cout << "beam\t" <<"triggered\t" << "Max coherent sum\t" << "Max sum location" <<"\n";
  //for (int i=0; i<ptrigger.L2_max_value.size(); i++) {
  //  std::cout << i << "\t\t" << ptrigger.L2_max_value[i]  <<"\n";
  //}


  //visual represent triggers - L1

  double w1 = 0.6; //horizontal separation in m between 2 adjacent azimuthal sectors
  double w = 3 * 0.6; //horizontal separation in m between azimuthal sectors furthest apart in L2 sector
  double h1 = 0.74; //vertical separation between adjacent antennas for lower 3 antennas in a //phi sector
  double h2 = 3.0; //vertical separation between adjacent antennas between top two antennas in //a given phi sector
  double h = 2 * h1 + h2;
  double c = 3e8;
  double delta_t = 1/(2.6e9); //time between samples, for 2.6Ghz sample frequency

  TCanvas *c_L1 = new TCanvas("L1","L1",0,0,600,400);
  TGraph2D *dt_L1 = new TGraph2D();
  dt_L1->SetTitle("Beam heat map L1; Theta; Phi; Coherent Sum Value");

  //vertical separation
  int v_max = floor((h * sin(50.0 / 180.0 * M_PI))/ (delta_t * c));
  int v_min = ceil((h * sin(-50.0 / 180.0 * M_PI))/ (delta_t * c));
  //horizontal separation
  int h_max = floor((w1 * sin(50.0 / 180.0 * M_PI))/ (delta_t * c));
  int h_min = ceil((w1 * sin(-50.0 / 180.0 * M_PI))/ (delta_t * c));

  int beam_L1_count = 0;
  //std::cout  << "L1 heat map\n";
  //std::cout  << "Theta/Phi\t";
  //for (int j = h_min; j < h_max+1; j++) {
  //     std::cout << asin(j * (delta_t * c) / w1) / M_PI * 180 << "\t";
  //}
  //std::cout << "\n";
   for (int i = v_min; i < v_max+1; i++) {
     //std::cout << asin(i * (delta_t * c) / h) / M_PI * 180 << "\t\t";
      for (int j = h_min; j < h_max+1; j++) {
            //std::cout << ptrigger.L1_max_value[beam_L1_count] << "\t";
            dt_L1->SetPoint(beam_L1_count,asin(i * (delta_t * c) / h) / M_PI * 180, asin(j * (delta_t * c) / w1) / M_PI * 180,  ptrigger->L1_max_value[beam_L1_count]);
            beam_L1_count++;
      }
      //std::cout << "\n";
   }
  c_L1->cd();
  gStyle->SetPalette(1);
  dt_L1->Draw("surf1");
  gPad->Update();
  dt_L1->GetXaxis()-> SetTitleOffset(2); 
  dt_L1->GetYaxis()-> SetTitleOffset(2); 
  dt_L1->GetZaxis()-> SetTitleOffset(-0.5); 

  dt_L1->GetXaxis()->SetLabelSize(0.05);
  dt_L1->GetYaxis()->SetLabelSize(0.05);
  dt_L1->GetZaxis()->SetLabelSize(0.05);

  dt_L1->GetXaxis()->SetTitleSize(.05);
  dt_L1->GetYaxis()->SetTitleSize(.05);


  //visual represent triggers - L2

  TCanvas *c_L2 = new TCanvas("L2","L2",0,0,600,400);
  TGraph2D *dt_L2 = new TGraph2D();
  dt_L2->SetTitle("Beam heat map L2; Theta; Phi; Coherent Sum Value");

  //vertical separation
  v_max = floor((h * sin(50.0 / 180.0 * M_PI))/ (delta_t * c));
  v_min = ceil((h * sin(-50.0 / 180.0 * M_PI))/ (delta_t * c));
  //horizontal separation
  h_max = floor((w * sin(50.0 / 180.0 * M_PI))/ (delta_t * c));
  h_min = ceil((w * sin(-50.0 / 180.0 * M_PI))/ (delta_t * c));

  int beam_L2_count = 0;
  //std::cout  << "L2 heat map\n";
  //std::cout  << "Theta/Phi\t";
  //for (int j = h_min; j < h_max+1; j++) {
  //     std::cout << asin(j * (delta_t * c) / w) / M_PI * 180 << "\t";
  //}
  //std::cout << "\n";
   for (int i = v_min; i < v_max+1; i++) {
     //std::cout << asin(i * (delta_t * c) / h) / M_PI * 180 << "\t\t";
      for (int j = h_min; j < h_max+1; j++) {
            //std::cout << ptrigger.L2_max_value[beam_L2_count] << "\t";
            dt_L2->SetPoint(beam_L2_count,asin(i * (delta_t * c) / h) / M_PI * 180, asin(j * (delta_t * c) / w) / M_PI * 180,  ptrigger->L2_max_value[beam_L2_count]);
            beam_L2_count++;
      }
     //std::cout << "\n";
   }

  gStyle->SetPalette(1);
  dt_L2->Draw("surf1");
  gPad->Update();
  dt_L2->GetXaxis()-> SetTitleOffset(2); 
  dt_L2->GetYaxis()-> SetTitleOffset(2); 
  dt_L2->GetZaxis()-> SetTitleOffset(-0.5); 

  dt_L2->GetXaxis()->SetLabelSize(0.05);
  dt_L2->GetYaxis()->SetLabelSize(0.05);
  dt_L2->GetZaxis()->SetLabelSize(0.05);

  dt_L2->GetXaxis()->SetTitleSize(.05);
  dt_L2->GetYaxis()->SetTitleSize(.05);


  //reenable cout
  std::cout.clear();
  std::cout.precision(5);


  //connect TApplication to canvas - seems like one is enough
  c_L2->Modified(); c_L2->Update();
  TRootCanvas *rc_L2 = (TRootCanvas *)c_L2->GetCanvasImp();
  rc_L2->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

  app.Run();

}  
int main(int argc, char **argv) {
  
  //setting up random numbers for closed loop noise generation
  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

  
  
  int l1threshold = 368;
  int l2threshold = 1306;
  int repeats;
  float samplingFreqHz = 2.56E9;

  TApplication app("app", &argc, argv); //interactive plots for reviewing threshold eval. Also needs app.Run() later
  //L1 threshold evaluation - at least 1E4 iterations; fast

  std::vector<nicemc::FTPair> noise_export = signal_gen(gen, .0, 512, true);
  TGraph gr = noise_export.at(1).getTimeDomain();
  gr.SaveAs("exportNoise.csv",".csv");


  std::cout<< "\n" << "--L1 threshold evaluation--" << "\n";
  repeats = 1E3;
  pueoSim::triggerThreshold * tThresholdL1 = new pueoSim::triggerThreshold(samplingFreqHz);
  tThresholdL1->setTriggerScaling(1.0);
  for (int r=0; r< repeats ; r++) {
    //replace signal_gen with pueoSim noise, no signal, vector of 16 FTPairs, each 512 samples
    tThresholdL1->L1Threshold_addData(signal_gen(gen, .0, 512, true));
    if (r % (repeats/10) == 0) {
      std::cout << "L1 iteration " << r << " done" << "\n";
    }
  }
  l1threshold = tThresholdL1->L1Threshold_eval();
  std::cout<< "\n\n" << "L1 threshold set to " << l1threshold << "\n";
  delete tThresholdL1;

  app.Run(); //interactive plots for reviewing threshold eval
  //L2 threshold evaluation - should be at least 1E4; slow 
  std::cout<< "\n" << "--L2 threshold evaluation--" << "\n";
  repeats = 1E2;
  pueoSim::triggerThreshold * tThresholdL2 = new pueoSim::triggerThreshold(samplingFreqHz);
  tThresholdL2->setTriggerScaling(1.0);
  for (int r=0; r< repeats ; r++) {
    //replace signal_gen with pueoSim noise, no signal, vector of 16 FTPairs, each 512 samples
    tThresholdL2->L2Threshold_addData(signal_gen(gen, .0, 512, true), l1threshold);
    if (r % (repeats/10) == 0) {
      std::cout << "L2 iteration " << r << " done" << "\n";
    }
  }
  l2threshold = tThresholdL2->L2Threshold_eval();
  std::cout<< "\n" << "L2 threshold set to " << l2threshold << "\n";
  delete tThresholdL2;

  app.Run(); //interactive plots for reviewing threshold eval

  //Do runs of trigger with the evaluated thresholds 
  std::cout<< "\n" << "--Trigger on signals with evaluated threshold--" << "\n";
  int total_pueo_runs = 100;
  double snr = 1.0;
  int L2_triggered_count = 0;
  pueoSim::pueoTrigger * ptrigger = new pueoSim::pueoTrigger(samplingFreqHz);
  ptrigger->setScaling(1.1);
  for(int run = 0; run < total_pueo_runs; run++ ) {
    //replace signal_gen with pueoSim signal+noise,vector of 16 FTPairs, each 512 samples     
    ptrigger->newSignal(signal_gen(gen, snr, 512, false));
    ptrigger->digitize(4);
    ptrigger->l1Trigger(8, 16, l1threshold, 64); //args: ptrigger, step, window, threshold, edge size
    ptrigger->l2Trigger(8, 16, l2threshold, 64); //args: ptrigger, step, window, threshold, edge size
    if (ptrigger->L2_triggered == true) {
      L2_triggered_count++;
    }
    
  }

  delete ptrigger;

  std::cout << "\n" <<"Loop runs: " <<  total_pueo_runs << "\n";
  std::cout << "\n" <<"Loop triggers: " <<  L2_triggered_count << "\n";

  //visualiseTrigger(argc, argv, L1_threshold, L2_threshold, snr);
  
  //app.Run(); //interactive plots for reviewing threshold eval

  return 0;
}















