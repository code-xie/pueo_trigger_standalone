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




class PueoSimple {
  public:
  static const int n_samples=512;
  static const int n_ant_L1=8;
  static const int n_ant_L2=16;

  std::vector<bool> L1_triggers;
  std::vector<bool> L2_triggers;

  std::vector<int> L1_max_value;
  std::vector<int> L2_max_value;
  std::vector<int> L1_max_loc;
  std::vector<int> L2_max_loc;


  std::vector<TGraph> signals;
  std::vector<TGraph> signals_discrete;

  PueoSimple(std::mt19937& gen, double snr);
};


void Digitize(PueoSimple* pueo1, int bits, double multiplier) {
  int digitise_max = pow(2,bits-1)-1;
  int digitise_min = -1 * digitise_max;

  std::vector<TGraph>::iterator it_ant;
  it_ant=pueo1->signals.begin();

  for (it_ant=pueo1->signals.begin();it_ant!=pueo1->signals.end(); ++it_ant) {

    pueo1->signals_discrete.emplace_back(*it_ant);
    TGraph *gr_digi = & pueo1->signals_discrete.back();

    const int n = it_ant->GetN();
    int yout[n];

    for (int i=0; i<gr_digi->GetN(); i++) {
      int digitised_y;
      double scaled_y = gr_digi->GetPointY(i) * multiplier;
      if (scaled_y > digitise_max) {
        digitised_y = digitise_max;
      } else if (scaled_y < digitise_min) {
        digitised_y = digitise_min;
      } else {
        digitised_y = round(scaled_y);
      }

      gr_digi->SetPointY(i,digitised_y );
    }
  }

}


static void generate_beams_L1(std::vector<std::vector<int>> &L1_beams) {
  double w1 = 0.6; //horizontal separation in m between 2 phi sectors
  double h1 = 0.6; //vertical separation between adjacent antennas for lower 3 antennas in a phi sector
  double h2 = 2.8; //vertical separation between adjacent antennas between top two antennas in a given phi sector
  double h = 2 * h1 + h2;
  double c = 3e8;
  double delta_t = 1/(2.6e9); //time between samples, for 2.6Ghz sample frequency

  /****
   * label antennas in L1 sector as below
   * top
   * 0  4
   * 1  5
   * 2  6
   * 3  7
   * bottom
   **/

  //vertical separation
  int v_max = floor((h * sin(50.0 / 180.0 * M_PI))/ (delta_t * c));
  int v_min = ceil((h * sin(-50.0 / 180.0 * M_PI))/ (delta_t * c));

  std::vector<std::vector<int>> vertical;
  for (int i = v_min; i < v_max+1; i++) {
    std::vector<int> beam;
    beam.push_back(0);
    beam.push_back((i / h * h2));
    beam.push_back((i / h * (h1+h2)));
    beam.push_back(i);
    //std::cout << beam.at(0)  << "\t" << beam.at(1) << "\t"<< beam.at(2) << "\t"<< beam.at(3) << "\t" <<"\n";
    vertical.push_back(beam);
  }

  //horizontal separation
  int h_max = floor((w1 * sin(50.0 / 180.0 * M_PI))/ (delta_t * c));
  int h_min = ceil((w1 * sin(-50.0 / 180.0 * M_PI))/ (delta_t * c));
  std::vector<std::vector<int>> horizontal;
  for (int i = h_min; i < h_max+1; i++) {
    std::vector<int> beam;
    beam.push_back(0);
    beam.push_back(i);
    //std::cout << beam.at(0)  << "\t" << beam.at(1) << "\t" <<"\n";
    horizontal.push_back(beam);
  }

  //combine vertical and horizontal for beams

  std::cout << "**********L1 Beams**********" <<"\n";

  for (std::vector<std::vector<int>>::iterator it_v = vertical.begin(); it_v != vertical.end(); ++it_v) {
    for (std::vector<std::vector<int>>::iterator it_h = horizontal.begin(); it_h != horizontal.end(); ++it_h) {
      std::vector<int> beam;
      beam.push_back(round(it_v->at(0)-it_h->at(0)));
      beam.push_back(round(it_v->at(1)-it_h->at(0)));
      beam.push_back(round(it_v->at(2)-it_h->at(0)));
      beam.push_back(round(it_v->at(3)-it_h->at(0)));
      beam.push_back(round(it_v->at(0)-it_h->at(1)));
      beam.push_back(round(it_v->at(1)-it_h->at(1)));
      beam.push_back(round(it_v->at(2)-it_h->at(1)));
      beam.push_back(round(it_v->at(3)-it_h->at(1)));

      std::cout << "L1 Beam:\t" << L1_beams.size()<< "\t" << "Offsets:" << "\t"  << beam.at(0)  << "\t" << beam.at(1) << "\t"<< beam.at(2) << "\t"<< beam.at(3) << "\t"<< beam.at(4)  << "\t" << beam.at(5) << "\t"<< beam.at(6) << "\t"<< beam.at(7) << "\t" <<"\n";
      L1_beams.push_back(beam);
    }
  }


}


////first variable is vector, each element a vector representing an L2 beam, with 16 delay values in that vector
////second variable is vector, each element a vector representing an L1 beam, with corresponding L2 beams in that vector
static void generate_beams_L2(std::vector<std::vector<int>> &L2_beams, std::vector<std::vector<int>> &L1_L2_map) {
  double w1 = 0.6; //horizontal separation in m between 2 phi sectors
  double w = 3 * 0.6; //horizontal separation in m between 2 phi sectors
  double h1 = 0.6; //vertical separation between adjacent antennas for lower 3 antennas in a //phi sector
  double h2 = 2.8; //vertical separation between adjacent antennas between top two antennas in //a given phi sector
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

  //vertical separation
  int v_max = floor((h * sin(50.0 / 180.0 * M_PI))/ (delta_t * c));
  int v_min = ceil((h * sin(-50.0 / 180.0 * M_PI))/ (delta_t * c));

  std::vector<std::vector<int>> vertical;
  for (int i = v_min; i < v_max+1; i++) {
    std::vector<int> beam;
    beam.push_back(0);
    beam.push_back((i / h * h2));
    beam.push_back((i / h * (h1+h2)));
    beam.push_back(i);
    //std::cout << beam.at(0)  << "\t" << beam.at(1) << "\t"<< beam.at(2) << "\t"<< beam.at(3) << "\t" << "\n";
    vertical.push_back(beam);
  }

  //horizontal separation
  int h_max = floor((w * sin(50.0 / 180.0 * M_PI))/ (delta_t * c));
  int h_min = ceil((w * sin(-50.0 / 180.0 * M_PI))/ (delta_t * c));
  std::vector<std::vector<int>> horizontal;
  for (int i = h_min; i < h_max+1; i++) {
    std::vector<int> beam;
    beam.push_back(0);
    beam.push_back((i / 3));
    beam.push_back((i / 3 * 2));
    beam.push_back(i);
    //std::cout << beam.at(0)  << "\t" << beam.at(1) << "\t" << beam.at(2)  << "\t" << beam.at(3) << "\t" <<"\n";
    horizontal.push_back(beam);
  }

  //horizontal separation L1
  int h_max_L1 = floor((w1 * sin(50.0 / 180.0 * M_PI))/ (delta_t * c));
  int h_min_L1 = ceil((w1 * sin(-50.0 / 180.0 * M_PI))/ (delta_t * c));
  std::vector<int> horizontal_L1;
  for (int i = h_min_L1; i < h_max_L1+1; i++) {
    horizontal_L1.push_back(i);
    //std::cout << horizontal_L1.at(0)   <<"\n";
  }

  std::vector<int> empty_map;
  L1_L2_map = std::vector<std::vector<int>> (vertical.size() * horizontal_L1.size(), empty_map);
  //std::cout << "L1 L2 map intialised with size " << vertical.size() * horizontal_L1.size()<< "\n";


  //combine vertical and horizontal for beams

  int v_counter = 0;
  std::cout << "**********L2 Beams**********" <<"\n";

  for (std::vector<std::vector<int>>::iterator it_v = vertical.begin(); it_v != vertical.end(); ++it_v) {
    for (std::vector<std::vector<int>>::iterator it_h = horizontal.begin(); it_h != horizontal.end(); ++it_h) {
      std::vector<int> beam;
      beam.push_back(round(it_v->at(0)-it_h->at(0)));
      beam.push_back(round(it_v->at(1)-it_h->at(0)));
      beam.push_back(round(it_v->at(2)-it_h->at(0)));
      beam.push_back(round(it_v->at(3)-it_h->at(0)));
      beam.push_back(round(it_v->at(0)-it_h->at(1)));
      beam.push_back(round(it_v->at(1)-it_h->at(1)));
      beam.push_back(round(it_v->at(2)-it_h->at(1)));
      beam.push_back(round(it_v->at(3)-it_h->at(1)));
      beam.push_back(round(it_v->at(0)-it_h->at(2)));
      beam.push_back(round(it_v->at(1)-it_h->at(2)));
      beam.push_back(round(it_v->at(2)-it_h->at(2)));
      beam.push_back(round(it_v->at(3)-it_h->at(2)));
      beam.push_back(round(it_v->at(0)-it_h->at(3)));
      beam.push_back(round(it_v->at(1)-it_h->at(3)));
      beam.push_back(round(it_v->at(2)-it_h->at(3)));
      beam.push_back(round(it_v->at(3)-it_h->at(3)));

      std::cout << "L2 Beam:\t" << L2_beams.size() << "\t" << "Offsets:" << "\t" ;
      for (int i = 0; i < beam.size(); i++) {
        std::cout << beam.at(i)  << "\t" ;
      }
      std::cout <<"\n";
      L2_beams.push_back(beam);

      //std::cout << "Adding to L1 L2 map" <<"\n";
      //map L1 beams to L2 beams where the delay between horizontal neighbours are +-1
      int element_min = std::max(*min_element(horizontal_L1.begin(), horizontal_L1.end()),it_h->at(1)-1);
      int pos_min = std::distance(horizontal_L1.begin(), std::find(horizontal_L1.begin(), horizontal_L1.end(), element_min));

      int element_max = std::min(*max_element(horizontal_L1.begin(), horizontal_L1.end()),it_h->at(1)+1);
      int pos_max = std::distance(horizontal_L1.begin(), std::find(horizontal_L1.begin(), horizontal_L1.end(), element_max));

      //std::cout << "matching element ranges: " << element_min << "\t" << element_max  << "\t" << "\n";
      //std::cout << "matching positions ranges: " << pos_min << "\t" << pos_max  << "\t" << "\n";
      //std::cout << "matching L1 beams: ";
      for(int i = v_counter * horizontal_L1.size() + pos_min; i < v_counter * horizontal_L1.size() + pos_max + 1; i++){
        //std::cout << i << "\t";
        L1_L2_map.at(i).push_back(L2_beams.size()-1);
      }
    }
    //std::cout << "v_counter" << "\t" << v_counter << "\n";
    v_counter++;

  }

  std::cout << "**********L1 L2 Map**********" <<"\n";
  for (int i=0; i < L1_L2_map.size(); ++i) {
    std::cout << "L1: " << i << "\t"<< "Mapped L2: " << "\t"  ;
    for (int j=0; j < L1_L2_map.at(i).size(); ++j) {
      std::cout << L1_L2_map.at(i).at(j) << "\t" ;
    }
    std::cout << "\n";
  }

}

//writes to vector of L1 triggers; requires defn of make up and relative delay between antennas for L1 beams
void L1Trigger(PueoSimple* pueo1, std::vector<std::vector<int>> & L1_ants,std::vector<std::vector<int>> & L1_beams,  int step, int window, int threshold, int max_shift) {


  for(int i_beam=0; i_beam <  L1_beams.size(); i_beam += 1) {

    bool beam_triggered = false;
    int max_coherent = 0;
    int max_loc = 0;
    for(int wind_pos=max_shift; wind_pos < pueo1->n_samples - max_shift-window -16; wind_pos+=step) {

      int coherent_sum = 0;
      for (int samp_pos=0; samp_pos < window ; samp_pos +=1){
        int coherent_sum_sample = 0;
        for (int i_ant=0; i_ant<pueo1->n_ant_L1; i_ant+=1) {
          coherent_sum_sample += pueo1->signals_discrete.at(L1_ants.at(i_beam).at(i_ant)).GetPointY(wind_pos+samp_pos-L1_beams.at(i_beam).at(i_ant));
        }
        int coherent_sum_sqr_sample = coherent_sum_sample * coherent_sum_sample;
        coherent_sum += coherent_sum_sqr_sample;
      }
      if (coherent_sum > threshold) {
        beam_triggered = true;
      }
      if (coherent_sum > max_coherent) { 
        max_coherent = coherent_sum;
        max_loc = wind_pos;
      }
      

    }

  pueo1->L1_triggers.push_back(beam_triggered);
  pueo1->L1_max_value.push_back(max_coherent);
  pueo1->L1_max_loc.push_back(max_loc);
  }
}


//outputs boolean of L2 trigger; requires dfn of make up and relative delay between antennas for L2 beams
void L2Trigger(PueoSimple* pueo1, std::vector<std::vector<int>> & L2_ants,std::vector<std::vector<int>> & L2_beams, std::vector<std::vector<int>> & L1_L2_map, int step, int window, int threshold, int max_shift) {


  //get list of L2 beams that need to be triggered
  std::vector<bool> L2_from_L1 (L2_beams.size(), 0);
  for (int i = 0; i < pueo1->L1_triggers.size(); i++) {
    if (pueo1->L1_triggers.at(i) == true) {
      for (std::vector<int>::iterator j = L1_L2_map.at(i).begin(); j != L1_L2_map.at(i).end(); ++j) {
        L2_from_L1.at(*j) = true;
      }
    }
  }

  for(int i_beam=0; i_beam <  L2_beams.size(); i_beam += 1) {
    if (L2_from_L1.at(i_beam) == true)  {
      bool beam_triggered = false;
      int max_coherent = 0;
      int max_loc = 0;
      for(int wind_pos=max_shift; wind_pos < pueo1->n_samples - max_shift-window -16; wind_pos+=step) {

        int coherent_sum = 0;
        for (int samp_pos=0; samp_pos < window ; samp_pos +=1){
          int coherent_sum_sample = 0;
          for (int i_ant=0; i_ant<pueo1->n_ant_L2; i_ant+=1) {
            coherent_sum_sample += pueo1->signals_discrete.at(L2_ants.at(i_beam).at(i_ant)).GetPointY(wind_pos+samp_pos-L2_beams.at(i_beam).at(i_ant));
          }
          int coherent_sum_sqr_sample = coherent_sum_sample * coherent_sum_sample;
          coherent_sum += coherent_sum_sqr_sample;
        }
        if (coherent_sum > threshold) {
          beam_triggered = true;
        }
        if (coherent_sum > max_coherent) {
          max_coherent = coherent_sum;
          max_loc = wind_pos;
        }
      }

      pueo1->L2_triggers.push_back(beam_triggered);
      pueo1->L2_max_value.push_back(max_coherent);
      pueo1->L2_max_loc.push_back(max_loc);
    } else {
      pueo1->L2_triggers.push_back(false);
      pueo1->L2_max_value.push_back(-1);
      pueo1->L2_max_loc.push_back(-1);
    }

  }
}


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


//std::array<double, 128> getSignal(double theta_offcone_deg, double expected_size) {
double * getSignal(double theta_offcone_deg, int expected_size) {

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


  fftw_destroy_plan(p);
  //fftw_free(IR_in); 
  //fftw_free(out); 

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


double gaussian_pulse(double x, double c, double pos, double peak) { //x is sample count, pos is which sample is centre of signal, c is width in samples 
    return peak * exp(-pow((x - pos -c /2),2) / 2 / c /c) - peak * exp(-pow((x - pos + c /2),2) / 2 / c /c) ;
}

void evaluate_threshold(PueoSimple& pueo1, std::mt19937& gen, bool is_L2) {
  

  int bits = 4;
  double multiplier = 2;
  //std::vector<double> samples;
  int n_samples = 512;
  int repeats = 40000;

  int number_ants = 8;
  if (is_L2 == true) {
    number_ants = 16;
  }

  TH1D *h1 = new TH1D("Histogram for coherent values","Histogram for coherent values",100,0.,5000.);

  int digitise_max = pow(2,bits-1)-1;
  int digitise_min = -1 * digitise_max;

  double x[n_samples];
  for (int i=0; i<n_samples; i++) {
    x[i] = i;
  }

  double sample_padding = 64;
  double n_sample_generated =  n_samples+sample_padding * 2;

  for (int r=0; r< repeats ; r++) {
    std::vector<TGraph> noise_tgraphs;
    for (int i = 0 ; i < number_ants; i++) {
    double * noise = generate_noise_flat(0.8343, n_sample_generated, 240, 1300, gen);
    TGraph gr = TGraph (n_samples, x, noise);
    
    //digitise
    for (int j=0; j<gr.GetN(); j++) {
      int digitised_y;
      double scaled_y = gr.GetPointY(j) * multiplier;

      if (scaled_y > digitise_max) {
        digitised_y = digitise_max;
      } else if (scaled_y < digitise_min) {
        digitised_y = digitise_min;
      } else {
        digitised_y = round(scaled_y);
      }
      
      gr.SetPointY(j,digitised_y );
    }
    //add to vector of tgraphs
    noise_tgraphs.push_back(gr); 
    }

    int step = 8;
    int window = 16;
    int max_shift = 64;


    int max_coherent = 0;
    for(int wind_pos=max_shift; wind_pos < n_samples - max_shift-window -16; wind_pos+=step) {
      int coherent_sum = 0;
      for (int samp_pos=0; samp_pos < window ; samp_pos +=1){
        int coherent_sum_sample = 0;
        for (int i_ant=0; i_ant<number_ants; i_ant+=1) {
          coherent_sum_sample += noise_tgraphs.at(i_ant).GetPointY(wind_pos+samp_pos);
        }
        int coherent_sum_sqr_sample = coherent_sum_sample * coherent_sum_sample;
        coherent_sum += coherent_sum_sqr_sample;
      }
      //samples.push_back(coherent_sum);
      h1->Fill(coherent_sum);

    }
  }
  

  h1->Draw();
  gPad->SetLogy();
  
} 




TGraph * combine_signal_noise(double signal_short[], double noise[], int signal_size, int signal_loc, int size, double SNR, double multiplier, double freq_Ghz) {

  int size_fft = ceil(size / 2) + 1;
  if (signal_size + signal_loc > size) {
    throw std::invalid_argument( "signal and location to large to fit within total profile");
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


  fftw_plan p = fftw_plan_dft_r2c_1d(size, &signal[0] , reinterpret_cast<fftw_complex*>(&sig_fft[0]) , FFTW_ESTIMATE);
  fftw_execute(p);

  

  p = fftw_plan_dft_r2c_1d(size, &noise[0] , reinterpret_cast<fftw_complex*>(&noise_fft[0]) , FFTW_ESTIMATE);
  fftw_execute(p);

  double * noise_fft_norm = new double[size_fft];
  for (int i = 0; i < size_fft; i++) {
    noise_fft_norm[i] = std::norm(noise_fft[i]);
  }
  //TCanvas *c_sig_test = new TCanvas("c_sig_test","c_sig_test",500,500,600,400);
  //c_sig_test->cd();
  //TGraph *sig_tgraph = new TGraph (size_fft,&freq_bins[0], &noise_fft_norm[0]);
  //sig_tgraph -> Draw();
  
  for (int i = 0; i < size_fft; i++) {
    sig_and_noise_fft[i]  = (multiplier / size) * (SNR * sig_fft[i] + noise_fft[i]) ; //FFTW does unnormalised DFT, introducing scaling by n 
  }

  p = fftw_plan_dft_c2r_1d(size,reinterpret_cast<fftw_complex*>(&sig_and_noise_fft[0]), &sig_and_noise[0] ,FFTW_ESTIMATE);
  fftw_execute(p);

  TGraph *sig_noise_tgraph = new TGraph (size,&time_bins[0], &sig_and_noise[0]);

  //TCanvas *c_signoise_test = new TCanvas("c_signoise_test","c_signoise_test",500,500,600,400);
  //c_signoise_test->cd();
  //sig_noise_tgraph -> Draw();


  delete[] signal;
  delete[] sig_fft;
  delete[] noise_fft;
  delete[] sig_and_noise_fft;
  return sig_noise_tgraph;
}


PueoSimple::PueoSimple(std::mt19937& gen, double snr) {
  

  //choosing source
  double theta = 10.0 / 180.0 * M_PI;
  double phi   = 10.0 / 180.0 * M_PI;


  //initialising antenna data
  double x[n_samples], y[n_ant_L2][n_samples];
  double std = 6;
  double pos = n_samples /2;
  double peak = 16;

  double w1 = 0.6; //horizontal separation in m between 2 phi sectors
  double w = 3 * 0.6; //horizontal separation in m between 2 phi sectors
  double h1 = 0.6; //vertical separation between adjacent antennas for lower 3 antennas in a //phi sector
  double h2 = 2.8; //vertical separation between adjacent antennas between top two antennas in //a given phi sector
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
    double * noise = generate_noise_flat(0.8343, n_sample_generated, 240, 1300, gen); 
    //double sum = 0;
    //for(int j = 0; j < n_sample_generated; j++) {
    //  sum += pow(noise[j], 2);
    //}
    //std::cout << "RMS: " << sqrt(sum/n_sample_generated) << std::endl;
    double * signal = getSignal(1, 128);
    //std::cout << "peak-to-peak/2: " << (*std::max_element(signal, signal+128) - *std::min_element(signal, signal+128)) / 2 << std::endl;
    received_signal_pre_interp[i] = *combine_signal_noise(signal, noise, 128, 128, n_sample_generated, snr, 2., 2.94912); //double signal_short[], double noise[], int signal_size, int signal_loc, int size, double SNR, double multiplier, double freq_Ghz
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





  for (int i_ant=0; i_ant<n_ant_L2; i_ant++) {
    TGraph gr = TGraph (n_samples, x, y[i_ant]);
    signals.push_back(gr);
  }

}


int main(int argc, char **argv) {
  TApplication app("app", &argc, argv); //this allows interactive plots for cmake application
  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

  std::ios::sync_with_stdio(false);
  std::cout.setstate(std::ios_base::failbit);

  //init beams
  std::vector<std::vector<int>> L1_ants; //beam number, antenna within beam
  std::vector<std::vector<int>> L1_beams;
  std::vector<std::vector<int>> L2_ants;
  std::vector<std::vector<int>> L2_beams;
  std::vector<std::vector<int>> L1_L2_map;

  generate_beams_L1(L1_beams);
  std::vector<int> antennas_L1 = {0,1,2,3,4,5,6,7};
  for (int i=0; i<L1_beams.size(); i++) {
    L1_ants.push_back(antennas_L1);
  }
  generate_beams_L2(L2_beams, L1_L2_map);
  std::vector<int> antennas_L2 = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  for (int i=0; i<L2_beams.size(); i++) {
    L2_ants.push_back(antennas_L2);
  }

  //reenable cout
  std::cout.clear();
  std::cout.precision(3);

  std::cout << "\n" <<"L1 Beam count: " <<  L1_beams.size() << "\n";
  std::cout << "\n" <<"L2 Beam count: " <<  L2_beams.size() << "\n";

  //suppress cout
  std::ios::sync_with_stdio(false);
  std::cout.setstate(std::ios_base::failbit);


  //Do loop of runs
  int total_pueo_runs = 0;
  double snr = 1.1;

  //int L1_threshold = 0;
  int L1_threshold = 1050;
  int L2_threshold = 5223;
  
  int L2_triggered = 0;
  double total_L1_beams;
  double L1_beams_triggered;
  for(int run = 0; run < total_pueo_runs; run++ ) {
    PueoSimple pueo(gen, snr);
    //suppress cout
    
    //create digitised data
    Digitize(&pueo, 4, pow(2,0));
  
    //do L1 trigger
    L1Trigger(&pueo,L1_ants,L1_beams, 8, 16, L1_threshold, 64); //args: pueo, step, window, threshold, edge size

    for (bool b: pueo.L1_triggers) {
      total_L1_beams++;
      if (b == true) {
        L1_beams_triggered++;
      }
    }

    //do L2 trigger
    L2Trigger(&pueo,L2_ants,L2_beams,L1_L2_map, 8, 16, L2_threshold, 64); //args: pueo, step, window, threshold, edge size

    //Compare L2 max values to threshold
    if ( *max_element(pueo.L2_max_value.begin(), pueo.L2_max_value.end()) > L2_threshold) {
      L2_triggered++;
    }

    

  }
  //reenable cout
  std::cout.clear();
  std::cout.precision(3);

  std::cout << "\n" <<"Loop runs: " <<  total_pueo_runs << "\n";
  std::cout << "\n" <<"Loop triggers: " <<  L2_triggered << "\n";
  std::cout << "\n" <<"Total L1 beams: " <<  total_L1_beams << "\n";
  std::cout << "\n" <<"L1 beams triggered: " <<  L1_beams_triggered << "\n";




  //suppress cout
  std::ios::sync_with_stdio(false);
  std::cout.setstate(std::ios_base::failbit);




  //initialise pueo object - sets up signal and noise for each channels
  PueoSimple pueo(gen, snr);

  std::cout.precision(3);


  //evaluate_threshold(pueo, gen, true); //pueo object, random number generator, is_L2
  //evaluate_threshold(pueo, gen, false); //pueo object, random number generator, is_L2

  //For testing - different off axis signals
  //visualiseFields();
  
  //For testing - getting signal
  //double * test_signal = getSignal(1, 128); //get signal at given off cone angle (in degrees), and expected number of values in impulse response file

  //For testing - getting noise
  //double * test_noise = generate_noise_flat(0.8343, pueo.n_samples, 240, 1300);
  //TGraph signal_noise_test = *combine_signal_noise(test_signal, test_noise, 128, 64, n_samples, 10, 1, 2.94912);



  //create digitised data
  Digitize(&pueo, 4, pow(2,0));

  //Plot digitised signals data
  TCanvas *c_digi = new TCanvas("c_digi","Discrete antenna data",500,500,600,400);
  TMultiGraph *mg = new TMultiGraph();
  //std::cout << "Digitised data:" << "\n";
  for (std::vector<TGraph>::iterator it=pueo.signals_discrete.begin();it!=pueo.signals_discrete.end(); ++it) {
    TGraph * gr3 = new TGraph(*it);
    gr3->SetTitle(std::to_string(it - pueo.signals_discrete.begin()).c_str());
    gr3->SetLineColor(it - pueo.signals_discrete.begin() + 1);
    mg->Add(gr3);
  }
  c_digi->cd();
  mg->Draw("A pmc plc");
  c_digi->BuildLegend();


  


  //do L1 trigger
  L1Trigger(&pueo,L1_ants,L1_beams, 8, 16, L1_threshold, 64); //args: pueo, step, window, threshold, edge size

  //do L2 trigger
  L2Trigger(&pueo,L2_ants,L2_beams,L1_L2_map, 8, 16, L2_threshold, 64); //args: pueo, step, window, threshold, edge size

  //print L1 triggers
  std::cout << "L1 Triggers" <<"\n";
  std::cout << "beam\t" <<"triggered\t" << "Max coherent sum\t" << "Max sum location" <<"\n";
  for (int i=0; i<pueo.L1_triggers.size(); i++) {
    std::cout << i << "\t\t" << pueo.L1_triggers[i] << "\t\t" << pueo.L1_max_value[i]<< "\t\t" << pueo.L1_max_loc[i]  <<"\n";
  }

  //print L2 triggers
  std::cout << "L2 Triggers" <<"\n";
  std::cout << "beam\t" <<"triggered\t" << "Max coherent sum\t" << "Max sum location" <<"\n";
  for (int i=0; i<pueo.L2_triggers.size(); i++) {
    std::cout << i << "\t\t" << pueo.L2_triggers[i] << "\t\t" << pueo.L2_max_value[i]<< "\t\t" << pueo.L2_max_loc[i]  <<"\n";
  }




  //visual represent triggers - L1
  double w1 = 0.6; //horizontal separation in m between 2 phi sectors
  double w = 3 * 0.6; //horizontal separation in m between 2 phi sectors
  double h1 = 0.6; //vertical separation between adjacent antennas for lower 3 antennas in a //phi sector
  double h2 = 2.8; //vertical separation between adjacent antennas between top two antennas in //a given phi sector
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
  std::cout  << "L1 heat map\n";
  std::cout  << "Theta/Phi\t";
  for (int j = h_min; j < h_max+1; j++) {
       std::cout << asin(j * (delta_t * c) / w1) / M_PI * 180 << "\t";
  }
  std::cout << "\n";
   for (int i = v_min; i < v_max+1; i++) {
     std::cout << asin(i * (delta_t * c) / h) / M_PI * 180 << "\t\t";
      for (int j = h_min; j < h_max+1; j++) {
            std::cout << pueo.L1_max_value[beam_L1_count] << "\t";
            dt_L1->SetPoint(beam_L1_count,asin(i * (delta_t * c) / h) / M_PI * 180, asin(j * (delta_t * c) / w1) / M_PI * 180,  pueo.L1_max_value[beam_L1_count]);
            beam_L1_count++;
      }
      std::cout << "\n";
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
  std::cout  << "L2 heat map\n";
  std::cout  << "Theta/Phi\t";
  for (int j = h_min; j < h_max+1; j++) {
       std::cout << asin(j * (delta_t * c) / w) / M_PI * 180 << "\t";
  }
  std::cout << "\n";
   for (int i = v_min; i < v_max+1; i++) {
     std::cout << asin(i * (delta_t * c) / h) / M_PI * 180 << "\t\t";
      for (int j = h_min; j < h_max+1; j++) {
            std::cout << pueo.L2_max_value[beam_L2_count] << "\t";
            dt_L2->SetPoint(beam_L2_count,asin(i * (delta_t * c) / h) / M_PI * 180, asin(j * (delta_t * c) / w) / M_PI * 180,  pueo.L2_max_value[beam_L2_count]);
            beam_L2_count++;
      }
      std::cout << "\n";
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
  std::cout.precision(3);


  //connect TApplication to canvas - seems like one is enough
  c_L2->Modified(); c_L2->Update();
  TRootCanvas *rc_L2 = (TRootCanvas *)c_L2->GetCanvasImp();
  rc_L2->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

  app.Run();


  return 0;
}
















