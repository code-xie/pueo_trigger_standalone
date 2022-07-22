#include <random>
#include "TAxis.h"
#include "TH1D.h"
#include "TH1.h"
#include "FTPair.h"
#include "TFitResult.h"
#include "TF1.h"
#include "trigger.h"

pueoSim::pueoTrigger::pueoTrigger(float samplingFreqHz_input, int antenna_start){

  samplingFreqHz = samplingFreqHz_input;
  first_antenna = (antenna_start+96)%96;

  get_beamsL1_simpleSeparation(L1_beams);
  n_beams_L1 = L1_beams.size();

  get_beamsL2_simpleSeparation(L2_beams);
  n_beams_L2 = L2_beams.size();

  

  //reenable cout
  std::cout.clear();
  std::cout.precision(5);

  std::cout << "\n" <<"pueoTrigger initialised with " <<  n_beams_L1 << " L1 beams, " <<  n_beams_L2 <<" L2 beams" << "\n";
}


void pueoSim::pueoTrigger::setScaling(float multiplier) {
  scaling = multiplier;
}


//Define beams at given intervals of deltaPhi and deltaTheta
void pueoSim::pueoTrigger::get_beamsL1_simpleSeparation(std::vector<std::vector<int>> &L1_beams) {
  gErrorIgnoreLevel = kError; // Suppress warnings as top two lines not read from photogrammetry file
  TTree *tree = new TTree("ntuple","data from csv file");
  tree->ReadFile("../data/pueoPhotogrammetry_220617.csv","An:X(in):Y(in):Z(in):HorizDist(in):AzCenter(deg):AperAz(deg):AperElev(deg):AntSize(in):description/C",',');
  gErrorIgnoreLevel = kPrint;

  double c = 3e8;
  float inch_to_m = 0.0254;
  
  Float_t An, X, Y, Z, Az, r;
  //tree->Print();
  tree->SetBranchAddress("An",&An);
  tree->SetBranchAddress("X(in)",&X);
  tree->SetBranchAddress("Y(in)",&Y);
  tree->SetBranchAddress("Z(in)",&Z);
  tree->SetBranchAddress("AzCenter(deg)",&Az);
  tree->SetBranchAddress("HorizDist(in)",&r);

  //for (int i = 0, N = tree->GetEntries(); i < N; ++i) {
  //  tree->GetEntry(i) ;
  //  std::cout<< r;
  //  std::cout<< "\n";
  //}

  //Choose beams from at given multiples of azimuthal and elevation angles from the centre of the triggering sector
  //define centre of sector as elevation perpendicular to the payload, and azimuth the average of the two bottom antennas (for L1) and average of the left/right most two bottom antennas (for L2)

  //1. Find centre of sector by taking bottom antennas and average their azimuth (evaluated from cartesian coordinates)
  tree->GetEntry(first_antenna+3);
  float azimuth_1 = Az;
  tree->GetEntry((first_antenna +7)%96);
  float azimuth_2 = Az;
  float centre_azimuth = (azimuth_1 + azimuth_2)/2;
  //std::cout << "Centre of sector azimuth = " << centre_azimuth << "\n";
  float centre_elevation = -10; //positive is defined as above the horizon. The antennas all point down 10 degrees below horizon


  //2. Loop through azimuths and elevation angles, then for each antenna,  
  for (float azimuth = centre_azimuth - 40.; azimuth <= centre_azimuth + 40.01; azimuth +=10.) {
      //std::cout << "Azimuth = " << azimuth << "\n";
      
      for (float elevation = centre_elevation - 40.; elevation <= centre_elevation + 40.01; elevation +=2.) {
      //std::cout << "Elevation = " << elevation << "\n";
      
      
      tree->GetEntry(first_antenna );
      
      double x_0 = cos(elevation/180.*M_PI) * (inch_to_m*Z * tan(elevation/180.*M_PI) - inch_to_m*r * cos((azimuth - Az)/180.*M_PI)) ;
      
      std::vector<int> beam;
      for (int antenna = first_antenna; antenna < (first_antenna + 8) ; antenna++) {
        tree->GetEntry((antenna+96) % 96);
        double delta_x = x_0 - cos(elevation/180.*M_PI) * (inch_to_m*Z * tan(elevation/180.*M_PI) - inch_to_m*r * cos((azimuth - Az)/180.*M_PI)) ;
        double delta_samples = delta_x * (samplingFreqHz / c);
        //std::cout << delta_samples << "     \t";
        beam.push_back(int(round(delta_samples)));
      }
      L1_beams.push_back(beam);
      //std::cout <<  "\n";
    }
    //std::cout << "\n";
  }

}

//Define beams at given intervals of deltaPhi and deltaTheta
void pueoSim::pueoTrigger::get_beamsL2_simpleSeparation(std::vector<std::vector<int>> &L2_beams) {
  gErrorIgnoreLevel = kError; // Suppress warnings as top two lines not read from photogrammetry file
  TTree *tree = new TTree("ntuple","data from csv file");
  tree->ReadFile("../data/pueoPhotogrammetry_220617.csv","An:X(in):Y(in):Z(in):HorizDist(in):AzCenter(deg):AperAz(deg):AperElev(deg):AntSize(in):description/C",',');
  gErrorIgnoreLevel = kPrint;

  //TTreeReader myReader("ntuple", tree);
  //TTreeReaderValue<Float_t> myPx(myReader, "An");

  double c = 3e8;
  float inch_to_m = 0.0254;
  
  Float_t An, X, Y, Z, Az, r;
  //tree->Print();
  tree->SetBranchAddress("An",&An);
  tree->SetBranchAddress("X(in)",&X);
  tree->SetBranchAddress("Y(in)",&Y);
  tree->SetBranchAddress("Z(in)",&Z);
  tree->SetBranchAddress("AzCenter(deg)",&Az);
  tree->SetBranchAddress("HorizDist(in)",&r);

  //for (int i = 0, N = tree->GetEntries(); i < N; ++i) {
  //  tree->GetEntry(i) ;
  //  std::cout<< r;
  //  std::cout<< "\n";
  //}

  //Choose beams from at given multiples of azimuthal and elevation angles from the centre of the triggering sector
  //define centre of sector as elevation perpendicular to the payload, and azimuth the average of the two bottom antennas (for L1) and average of the left/right most two bottom antennas (for L2)

  //1. Find centre of sector by taking bottom antennas and average their azimuth (evaluated from cartesian coordinates)
  tree->GetEntry(first_antenna+3);
  float azimuth_1 = Az;
  tree->GetEntry((first_antenna +15)%96);
  float azimuth_2 = Az;
  float centre_azimuth = (azimuth_1 + azimuth_2)/2;
  //std::cout << "Centre of sector azimuth = " << centre_azimuth << "\n";
  float centre_elevation = -10; //positive is defined as above the horizon. The antennas all point down 10 degrees below horizon


  //2. Loop through azimuths and elevation angles, then for each antenna,  
  for (float azimuth = centre_azimuth - 40.; azimuth <= centre_azimuth + 40.01; azimuth +=5.) {
      //std::cout << "Azimuth = " << azimuth << "\n";
      
      for (float elevation = centre_elevation - 40.; elevation <= centre_elevation + 40.01; elevation +=2.) {
      //std::cout << "Elevation = " << elevation << "\n";
      
      
      tree->GetEntry(first_antenna );
      
      double x_0 = cos(elevation/180.*M_PI) * (inch_to_m*Z * tan(elevation/180.*M_PI) - inch_to_m*r * cos((azimuth - Az)/180.*M_PI)) ;
      
      std::vector<int> beam;
      for (int antenna = first_antenna; antenna < first_antenna + 16 ; antenna++) {
        tree->GetEntry((antenna+96) % 96);
        double delta_x = x_0 - cos(elevation/180.*M_PI) * (inch_to_m*Z * tan(elevation/180.*M_PI) - inch_to_m*r * cos((azimuth - Az)/180.*M_PI)) ;
        double delta_samples = delta_x * (samplingFreqHz / c);
        //std::cout << delta_samples << "     \t";
        beam.push_back(int(round(delta_samples)));
      }
      L2_beams.push_back(beam);
      //std::cout <<  "\n";
    }
    //std::cout << "\n";
  }

}


void pueoSim::pueoTrigger::newSignal(std::vector<nicemc::FTPair> input_signals) {
  signals.clear();
  signals_discrete.clear();
  L2_triggered = false;
  L1_max_value.clear();
  L2_max_value.clear();
  L1_triggered_windows.clear();
  L1_max_value = std::vector<int>(n_beams_L1);
  L2_max_value = std::vector<int>(n_beams_L2);
  
  for (int i_ant=0; i_ant<n_ant_L2; i_ant++) {
    //std::cout << "Looping over antenna" << i_ant << std::endl;
    TGraph gr = input_signals.at(i_ant).getTimeDomain();
    signals.push_back(gr);
  }

  n_samples = signals.at(0).GetN();
}

void pueoSim::pueoTrigger::digitize(int bits) {
  int digitise_max = pow(2,bits-1)-1;
  int digitise_min = -1 * digitise_max;

  std::vector<TGraph>::iterator it_ant;
  //it_ant=signals.begin();
  signals_discrete.clear();

  for (it_ant=signals.begin();it_ant!=signals.end(); ++it_ant) {

    signals_discrete.emplace_back(*it_ant);
    TGraph *gr_digi = & signals_discrete.back();

    const int n = it_ant->GetN();

    for (int i=0; i<gr_digi->GetN(); i++) {
      int digitised_y;
      double scaled_y = gr_digi->GetPointY(i) * scaling;
      if (scaled_y > digitise_max) {
        digitised_y = digitise_max;
      } else if (scaled_y < digitise_min) {
        digitised_y = digitise_min;
      } else {
        digitised_y = round(scaled_y);
        //digitised_y = scaled_y;
      }

      gr_digi->SetPointY(i,digitised_y );
    }
  }
}

void pueoSim::pueoTrigger::digitize_afterFilter(int bits) {
  int digitise_max = pow(2,bits-1)-1;
  int digitise_min = -1 * digitise_max;

  std::vector<TGraph>::iterator it_ant;
  //it_ant=signals.begin();
  signals_discrete.clear();

  for (it_ant=signals_filtered.begin();it_ant!=signals_filtered.end(); ++it_ant) {

    signals_discrete.emplace_back(*it_ant);
    TGraph *gr_digi = & signals_discrete.back();

    const int n = it_ant->GetN();

    for (int i=0; i<gr_digi->GetN(); i++) {
      int digitised_y;
      double scaled_y = gr_digi->GetPointY(i) * scaling;
      if (scaled_y > digitise_max) {
        digitised_y = digitise_max;
      } else if (scaled_y < digitise_min) {
        digitised_y = digitise_min;
      } else {
        digitised_y = round(scaled_y);
        //digitised_y = scaled_y;
        //std::cout << scaled_y << "\t\t" << round(scaled_y) << "\n";
      }

      gr_digi->SetPointY(i,digitised_y );
    }
  }
}

void pueoSim::pueoTrigger::firFilter() {

  double filter[] = {0, 0.0475, 0, -0.0938, 0, 0.3046, 0.4832, 0.3046, 0, -0.0938, 0, 0.0475, 0};
  int filter_size = 13;

  std::vector<TGraph>::iterator it_ant;
  signals_filtered.clear();

  for (it_ant=signals_discrete.begin();it_ant!=signals_discrete.end(); ++it_ant) {
    signals_filtered.emplace_back(*it_ant);
    TGraph *gr_fir = & signals_filtered.back();

    int signal_size = gr_fir->GetN();
    for (int i=0; i < signal_size; i++) {
      double  acc = 0;
      for (int j=0; j < filter_size; j++) {
        if ((i-j) > -1) {
            //std::cout << it_ant->GetPointY(i-j) << " ";
            acc += it_ant->GetPointY(i-j) * filter[j];
        }
      }
      //std::cout << std::endl;
      gr_fir->SetPointY(i, acc);
    //std::cout << acc << " ";
    }
    //std::cout << std::endl;

  }

  

}

void pueoSim::pueoTrigger::firFilter_signal_to_fir() {

  double filter[] = {0, 0.0475, 0, -0.0938, 0, 0.3046, 0.4832, 0.3046, 0, -0.0938, 0, 0.0475, 0};
  //double filter[] =  {0.0467, 0.0849, -0.0647, -0.0821, 0.0441, 0.3157, 0.4479, 0.3157, 0.0441, -0.0821, -0.0647, 0.0849, 0.0467};
  int filter_size = 13;

  std::vector<TGraph>::iterator it_ant;
  signals_filtered.clear();

  //testing variable
  bool testDone = false;
  double filter_output[512]; 
  double filter_input[512]; 

  for (it_ant=signals.begin();it_ant!=signals.end(); ++it_ant) {
    signals_filtered.emplace_back(*it_ant);
    TGraph *gr_fir = & signals_filtered.back();

    int signal_size = gr_fir->GetN();
    for (int i=0; i < signal_size; i++) {
      double  acc = 0;
      for (int j=0; j < filter_size; j++) {
        if ((i-j) > -1) {
            //std::cout << it_ant->GetPointY(i-j) << " ";
            acc += it_ant->GetPointY(i-j) * filter[j];
        }
      }
      //std::cout << std::endl;
      gr_fir->SetPointY(i, acc);
      if (!testDone) {
        filter_output[i] = acc;
      }

    //std::cout << acc << " ";
    }
    //std::cout << std::endl;
    //
    if (!testDone) {
      for (int i=0; i < signal_size; i++) {
        filter_input[i] = it_ant->GetPointY(i);
      }
    }
    testDone = true;
  }

  /*
  //plot tests
  int size = 512;
  int size_fft = ceil(size / 2) + 1;

  double * freq_bins = new double[size_fft];
  for (int i = 0; i < size_fft; i++) {
    freq_bins[i] = i;
  }

  double * time_bins = new double[size];
  double freq_hz = 2.56E9;
  for (int i = 0; i < size; i++) {
    time_bins[i] = i * 1.0 / freq_hz;
  }

  //input fft norm plot
  std::complex<double> * input_fft = new std::complex<double>[size_fft];
  double * input_fft_norm = new double[size_fft];

  fftw_plan p = fftw_plan_dft_r2c_1d(size, &filter_input[0] , reinterpret_cast<fftw_complex*>(&input_fft[0]) , FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);

  for (int i = 0; i < size_fft; i++) {
    input_fft_norm[i] = std::norm(input_fft[i]);
    //std::cout << input_fft_norm[i] << " ";
  }
  //std::cout << std::endl;

  TGraph *input_fft_norm_tgraph = new TGraph (size_fft,&freq_bins[0], &input_fft_norm[0]);
  TCanvas *d3 = new TCanvas("d3","d3",500,500,600,400);
  d3->cd();
  input_fft_norm_tgraph -> Draw();

  //output fft norm plot
  std::complex<double> * output_fft = new std::complex<double>[size_fft];
  double * output_fft_norm = new double[size_fft];

  p = fftw_plan_dft_r2c_1d(size, &filter_output[0] , reinterpret_cast<fftw_complex*>(&output_fft[0]) , FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);

  for (int i = 0; i < size_fft; i++) {
    output_fft_norm[i] = std::norm(output_fft[i]);
    //std::cout << input_fft_norm[i] << " ";
  }
  //std::cout << std::endl;

  TGraph *output_fft_norm_tgraph = new TGraph (size_fft,&freq_bins[0], &output_fft_norm[0]);
  TCanvas *d4 = new TCanvas("d4","d4",500,500,600,400);
  d4->cd();
  output_fft_norm_tgraph -> Draw();

  //input timedomain plot
  TGraph *input_tgraph = new TGraph (size,&time_bins[0], &filter_input[0]);
  TCanvas *d6 = new TCanvas("d6","d6",500,500,600,400);
  d6->cd();
  input_tgraph -> Draw();

  //output timedomain plot
  TGraph *output_tgraph = new TGraph (size,&time_bins[0], &filter_output[0]);
  TCanvas *d5 = new TCanvas("d5","d5",500,500,600,400);
  d5->cd();
  output_tgraph -> Draw();

  //get RMS
  double squared_sum = 0;
  for (int i = 0; i < size; i++) {
    squared_sum += filter_input[i] *filter_input[i] ;
  }
  std::cout << "Input RMS: " << sqrt(squared_sum/size);

  squared_sum = 0;
  for (int i = 0; i < size; i++) {
    squared_sum += filter_output[i] * filter_output[i];
  }
  std::cout << "Output RMS: " << sqrt(squared_sum/size);
  std::cout << std::endl;

  */


}



void pueoSim::pueoTrigger::l1Trigger(int step, int window, int threshold, int max_shift) {
//window is number of samples in a given coherent sum. Step is how many samples between subsequent windows.
//threshold is value compared with coherent sum for trigger
//max_shift accounts for edge effects where samples should not be formed around the start and end of the signal as samples have been shifted.

  //initialise max values to all zeros. This is used for visualisation purposes, not required for the simulation itself
  for(int i_beam=0; i_beam <  n_beams_L1; i_beam += 1) {
    L1_max_value.at(i_beam) = 0;
  }

  //initalize all windows to have NOT triggered (set to false):
  for(int wind_pos=max_shift; wind_pos < n_samples - max_shift-window -16; wind_pos+=step) {
    L1_triggered_windows.push_back(false);
  }



  //Do first L1 sector
  for(int i_beam=0; i_beam <  n_beams_L1; i_beam += 1) {
  int waveform_length = signals_discrete.at(0).GetN();
  int total_shifted[waveform_length]={0};

  //sum signals across antennas after shifting each based on beam definition
  for (int i_ant=0;i_ant<n_ant_L1;i_ant++){
    for (int samp_pos=0;samp_pos<waveform_length;samp_pos++){
      total_shifted[samp_pos]+= signals_discrete.at(i_ant).GetPointY(samp_pos-L1_beams.at(i_beam).at(i_ant));
      //std::cout<<"here! "<<samp_pos<<", "<<total_shifted[samp_pos]<<std::endl;

    }
  }

  //square so that total_shifted has power:
  for (int samp_pos=0;samp_pos<waveform_length;samp_pos++){
    total_shifted[samp_pos]=pow(total_shifted[samp_pos],2);
    //std::cout<<"shifted val is "<<total_shifted[samp_pos]<<std::endl;
  }
  //Now cycle through windows and record if a winow triggers. Note that as long as any beam triggers, the trigger for that window is set to //true:
  int window_count = 0;

  for(int wind_pos=max_shift; wind_pos < n_samples - max_shift-window -16; wind_pos+=step) {
    int coherent_sum = 0;

    for (int samp_pos=0;samp_pos<window; samp_pos++){
    coherent_sum+=total_shifted[wind_pos+samp_pos];
    }
    //std::cout<<"coherent sum: "<<coherent_sum <<", threshold: "<<threshold << std::endl;
    if(coherent_sum>threshold){
    //std::cout<<"triggered l1!"<<std::endl;
    L1_triggered_windows.at(window_count)=true;
    }

    //This is used for visualisation purposes, not required for the simulation itself
    if (coherent_sum > L1_max_value.at(i_beam)){
      L1_max_value.at(i_beam) = coherent_sum;
    }

    window_count++;
  }

  }

  //Do second L1 sector
  for(int i_beam=0; i_beam <  n_beams_L1; i_beam += 1) {
    int waveform_length = signals_discrete.at(0).GetN();
    int total_shifted[waveform_length]={0};

    for (int i_ant=0;i_ant<n_ant_L1;i_ant++){
      for (int samp_pos=0;samp_pos<waveform_length;samp_pos++){
        total_shifted[samp_pos]+= signals_discrete.at(i_ant+8).GetPointY(samp_pos-L1_beams.at(i_beam).at(i_ant));
        //std::cout<<"here! "<<samp_pos<<", "<<total_shifted[samp_pos]<<std::endl;

      }
    }

    //square so that total_shifted has power:
    for (int samp_pos=0;samp_pos<waveform_length;samp_pos++){
      total_shifted[samp_pos]=pow(total_shifted[samp_pos],2);
      //std::cout<<"shifted val is "<<total_shifted[samp_pos]<<std::endl;
    }
    //Now cycle through windows and record if a winow triggers. Note that as long as any beam triggers, the trigger for that window is set to //true:
    int window_count = 0;

    for(int wind_pos=max_shift; wind_pos < n_samples - max_shift-window -16; wind_pos+=step) {
      int coherent_sum = 0;

      for (int samp_pos=0;samp_pos<window; samp_pos++){
        coherent_sum+=total_shifted[wind_pos+samp_pos];
      }
      //std::cout<<"coherent sum: "<<coherent_sum <<", threshold: "<<threshold << std::endl;
      if(coherent_sum>threshold){
        //std::cout<<"triggered l1!"<<std::endl;
        L1_triggered_windows.at(window_count)=true;
      }
      window_count++;
    }

  }

//old way - less efficient
// for(int wind_pos=max_shift; wind_pos < n_samples - max_shift-window -16; wind_pos+=step) {
//      bool window_triggered = false;

//      //Do first L1 sector
//      for(int i_beam=0; i_beam <  n_beams_L1; i_beam += 1) {
//        int coherent_sum = 0;

//        for (int samp_pos=0; samp_pos < window ; samp_pos +=1){
//          int coherent_sum_sample = 0;
//          for (int i_ant=0; i_ant<n_ant_L1; i_ant+=1) {
//            coherent_sum_sample += signals_discrete.at(L1_ants.at(i_beam).at(i_ant)).GetPointY(wind_pos+samp_pos-L1_beams.at(i_beam).at(i_ant));
//          }
//          int coherent_sum_sqr_sample = coherent_sum_sample * coherent_sum_sample;
//          coherent_sum += coherent_sum_sqr_sample;
//        }
//      
//        if (coherent_sum > threshold) {
//          window_triggered = true;

//        }
//        if (coherent_sum > L1_max_value.at(i_beam)){
//          L1_max_value.at(i_beam) = coherent_sum;
//        }
//      
//      }

//      //Do second L1 sector
//      for(int i_beam=0; i_beam <  n_beams_L1; i_beam += 1) {
//        int coherent_sum = 0;

//        for (int samp_pos=0; samp_pos < window ; samp_pos +=1){
//          int coherent_sum_sample = 0;
//          for (int i_ant=0; i_ant<n_ant_L1; i_ant+=1) {
//            coherent_sum_sample += signals_discrete.at(L1_ants.at(i_beam).at(i_ant)+8).GetPointY(wind_pos+samp_pos-L1_beams.at(i_beam).at(i_ant));
//          }
//          int coherent_sum_sqr_sample = coherent_sum_sample * coherent_sum_sample;
//          coherent_sum += coherent_sum_sqr_sample;
//        }
//      
//        if (coherent_sum > threshold) {
//          //std::cout << "sum=" << coherent_sum;
//          window_triggered = true;
//        }
//        if (coherent_sum > L1_max_value.at(i_beam)){
//          L1_max_value.at(i_beam) = coherent_sum;
//        }
//      
//      }

//      //std::cout << "tr=" << window_triggered << "/n";
//      L1_triggered_windows.push_back(window_triggered);
//  }

  
}

void pueoSim::pueoTrigger::l2Trigger(int step, int window, int threshold, int max_shift) {

  //initialise max values to all zeros. This is used for visualisation purposes, not required for the simulation itself
  for(int i_beam=0; i_beam <  n_beams_L2; i_beam += 1) {
    L2_max_value.at(i_beam) = 0;
  }

  L2_triggered = false;

  for(int i_beam=0; i_beam <  n_beams_L2; i_beam += 1) {
    int waveform_length = signals_discrete.at(0).GetN();
    int total_shifted[waveform_length]={0};

    //sum signals across antennas after shifting each based on beam definition
    for (int i_ant=0;i_ant<n_ant_L2;i_ant++){
      for (int samp_pos=0;samp_pos<waveform_length;samp_pos++){
        total_shifted[samp_pos]+= signals_discrete.at(i_ant).GetPointY(samp_pos-L2_beams.at(i_beam).at(i_ant));
        //std::cout<<"here! "<<samp_pos<<", "<<total_shifted[samp_pos]<<std::endl;
        //std::cout << total_shifted[samp_pos] << " ";
      }
    }
    //std::cout<< "\n";

    //square so that total_shifted has power:
    for (int samp_pos=0;samp_pos<waveform_length;samp_pos++){

      total_shifted[samp_pos]=total_shifted[samp_pos]*total_shifted[samp_pos];

    }

    //Now cycle through windows and record if a winow triggers. Note that as long as any beam triggers, the trigger for that window is set to true:
    int window_count = 0;
    for(int wind_pos=max_shift; wind_pos < n_samples - max_shift-window -16; wind_pos+=step) {
      if (L1_triggered_windows.at(window_count) == true) {
        int coherent_sum = 0;

        for (int samp_pos=0;samp_pos<window; samp_pos+=1 ){
        coherent_sum+=total_shifted[wind_pos+samp_pos];
        }

        if(coherent_sum>threshold){
          L2_triggered = true;
        }

        if (coherent_sum > L2_max_value.at(i_beam)){
          L2_max_value.at(i_beam) = coherent_sum;
        }

        
      }
      window_count++;
    }
  }

  //old way - less efficient
//  int window_count = 0;
//  L2_triggered = false;
//  
//  for(int wind_pos=max_shift; wind_pos < n_samples - max_shift-window -16; wind_pos+=step) {
//    if (L1_triggered_windows.at(window_count) == true) {
//      for(int i_beam=0; i_beam <  n_beams_L2; i_beam += 1) {
//        int coherent_sum = 0;
//        for (int samp_pos=0; samp_pos < window ; samp_pos +=1){
//          int coherent_sum_sample = 0;
//          for (int i_ant=0; i_ant<n_ant_L2; i_ant+=1) {
//            coherent_sum_sample += signals_discrete.at(i_ant).GetPointY(wind_pos+samp_pos-L2_beams.at(i_beam).at(i_ant));
//          }
//          int coherent_sum_sqr_sample = coherent_sum_sample * coherent_sum_sample;
//          coherent_sum += coherent_sum_sqr_sample;
//        }
//        if (coherent_sum > threshold) {
//          L2_triggered = true;
//
//        }
//        if (coherent_sum > L2_max_value.at(i_beam)){
//          L2_max_value.at(i_beam) = coherent_sum;
//        }
//      }
//    }
//    window_count++;
//  }

}

pueoSim::triggerThreshold::triggerThreshold(float samplingFreqHz_input, int first_antenna) {
  window_count = 0;
  ptrigger = new pueoTrigger(samplingFreqHz_input, (first_antenna+96)%96);
  
  //h1 = new TH1D("Histogram for coherent values","Histogram for coherent values",100,0.,5000.);
}

void pueoSim::triggerThreshold::setTriggerScaling(float multiplier) {
  ptrigger->setScaling(multiplier);
}

//note: only 1 beam is generated, as otherwise the data points may not be independent. 
//note: it's also not a real beam, but this shouldn't matter for our purposes
void pueoSim::triggerThreshold::L1Threshold_addData(std::vector<nicemc::FTPair> input_signals) {
  ptrigger->newSignal(input_signals);
  ptrigger->digitize(4);
  ptrigger->firFilter_signal_to_fir();
  ptrigger->digitize_afterFilter(4);

  
  int n_samples = ptrigger->signals_discrete.at(0).GetN();
  int number_ants = 8;  
  int step = 8;
  int window = 16;
  int max_shift = 64;

  int waveform_length = ptrigger->signals_discrete.at(0).GetN();
  int total_shifted[waveform_length]={0};

  for (int i_ant=0;i_ant<number_ants;i_ant++){
    for (int samp_pos=0;samp_pos<waveform_length;samp_pos++){
      total_shifted[samp_pos]+= ptrigger->signals_discrete.at(i_ant).GetPointY(samp_pos);
      //std::cout<<"here! "<<samp_pos<<", "<<total_shifted[samp_pos]<<std::endl;
    }
  }

  //square so that total_shifted has power:
  for (int samp_pos=0;samp_pos<waveform_length;samp_pos++){
    total_shifted[samp_pos]=total_shifted[samp_pos]*total_shifted[samp_pos];
  }

  for(int wind_pos=max_shift; wind_pos < n_samples - max_shift-window -16; wind_pos+=step) {
    int coherent_sum = 0;
    for (int samp_pos=0;samp_pos<window; samp_pos++){
    coherent_sum+=total_shifted[wind_pos+samp_pos];
    }
    coherent_values.push_back(coherent_sum);
    window_count++;
  }


  //old way
  //for(int wind_pos=max_shift; wind_pos < n_samples - max_shift-window -16; wind_pos+=step) {
  //  int coherent_sum = 0;
  //  for (int samp_pos=0; samp_pos < window ; samp_pos +=1){
  //    int coherent_sum_sample = 0;
  //    for (int i_ant=0; i_ant<number_ants; i_ant+=1) {
  //      coherent_sum_sample += ptrigger->signals_discrete.at(i_ant).GetPointY(wind_pos+samp_pos);
  //    }
  //    int coherent_sum_sqr_sample = coherent_sum_sample * coherent_sum_sample;
  //    coherent_sum += coherent_sum_sqr_sample;
  //  }
  //  
  //  
  //  coherent_values.push_back(coherent_sum);
  //  window_count++;
  //}

}

int pueoSim::triggerThreshold::L1Threshold_eval(double sample_rate) {  

  double stdev =  TMath::StdDev(coherent_values.begin(), coherent_values.end() );
  double mean =  TMath::Mean(coherent_values.begin(), coherent_values.end() );
  double size =  coherent_values.size();
  std::cout << "\nCount = " << size << ", STDeV = " << stdev << ", Mean = " << mean << "\n\n";


  TCanvas *g1 = new TCanvas("g1","g1",500,500,600,400);
  g1->cd();

  TH1D * h1 = new TH1D("Histogram for L1","Histogram for L1",100,0.,8 * stdev + mean);
  for (int i = 0; i < coherent_values.size(); i++) {
    h1->Fill(coherent_values.at(i));
  }



  TF1 * f1 = new TF1("f1","expo");
  f1->SetParameters(1,1);
  h1->Fit("f1","","", mean + 3 * stdev, mean + 7 * stdev);
  TF1 * fitFunc = h1->GetFunction("f1");
  //std::cout << "p0=" << fitFunc->GetParameter(0) << " p1=" << fitFunc->GetParameter(1) ;
  
  h1->Draw();
  gPad->SetLogy();

  double sector_rate = 1E6;
  double step = 8;


  
  return round((std::log(sector_rate * step / sample_rate / ptrigger->n_beams_L1 *  window_count * (1 - exp((8 * stdev + mean)/100 *  fitFunc->GetParameter(1)))) - fitFunc->GetParameter(0) ) / fitFunc->GetParameter(1));

}


void  pueoSim::triggerThreshold::L2Threshold_addData(std::vector<nicemc::FTPair> input_signals, int L1Threshold) {
  ptrigger->newSignal(input_signals);
  ptrigger->digitize(4);
  ptrigger->firFilter_signal_to_fir();
  ptrigger->digitize_afterFilter(4);

  ptrigger->l1Trigger(8, 16, L1Threshold, 64);
 
  

  int n_samples = ptrigger->signals_discrete.at(0).GetN();
  int step = 8;
  int window = 16;
  int max_shift = 64;

// pre-shift approach has issue since we want to keep track of largest coherent value at the window, across beams. Needs array of coherent sum maxes
  int num_windows = 0;
  for(int wind_pos=max_shift; wind_pos < n_samples - max_shift-window -16; wind_pos+=step) {
    num_windows++;
  }
  int window_count_this_signal = 0;

/*
  //new way - about half the speed of old way due to nature of calculation
  int coherent_sum_max[num_windows] = {0};
  for(int i_beam=0; i_beam <  ptrigger->n_beams_L2; i_beam += 1) {
    int waveform_length = ptrigger->signals_discrete.at(0).GetN();
    int total_shifted[waveform_length]={0};
    

    for (int i_ant=0;i_ant<number_ants;i_ant++){
      for (int samp_pos=0;samp_pos<waveform_length;samp_pos++){
        total_shifted[samp_pos]+= ptrigger->signals_discrete.at(i_ant).GetPointY(samp_pos-ptrigger->L2_beams.at(i_beam).at(i_ant));
        //std::cout<<"here! "<<samp_pos<<", "<<total_shifted[samp_pos]<<std::endl;
      }
    }

    //square so that total_shifted has power:
    for (int samp_pos=0;samp_pos<waveform_length;samp_pos++){
      total_shifted[samp_pos]=total_shifted[samp_pos]*total_shifted[samp_pos];
    }

    //
    window_count_this_signal = 0;
    for(int wind_pos=max_shift; wind_pos < n_samples - max_shift-window -16; wind_pos+=step) {
      if (ptrigger->L1_triggered_windows.at(window_count_this_signal) == true) {
        int coherent_sum = 0;
        for (int samp_pos=0;samp_pos<window; samp_pos++){
          coherent_sum+=total_shifted[wind_pos+samp_pos];
        }
        if (coherent_sum > coherent_sum_max[window_count_this_signal]) {
          coherent_sum_max[window_count_this_signal] = coherent_sum;
        }
      }
      window_count_this_signal++;
    }
  }

  window_count_this_signal = 0;
  for(int wind_pos=max_shift; wind_pos < n_samples - max_shift-window -16; wind_pos+=step) {
    if (ptrigger->L1_triggered_windows.at(window_count_this_signal) == true) {
      coherent_values.push_back(coherent_sum_max[window_count_this_signal]);
      //std::cout << "coherent sum pushed = " << coherent_sum_max[window_count_this_signal] << std::endl;
     
    }
    window_count++; //every window counts for threshold evaluation regardless if triggered
    window_count_this_signal++;
    
  }

*/
  
///*
  for(int wind_pos=max_shift; wind_pos < n_samples - max_shift-window -16; wind_pos+=step) {
    if (ptrigger->L1_triggered_windows.at(window_count_this_signal) == true) { //the threshold value is only relevant if L1 triggered 
      int coherent_sum_max = 0;
      for(int i_beam=0; i_beam <  ptrigger->n_beams_L2; i_beam += 1) {
        int coherent_sum = 0;
        for (int samp_pos=0; samp_pos < window ; samp_pos +=1){
          int coherent_sum_sample = 0;
          for (int i_ant=0; i_ant<ptrigger->n_ant_L2; i_ant+=1) {
            coherent_sum_sample += ptrigger->signals_discrete.at(i_ant).GetPointY(wind_pos+samp_pos-ptrigger->L2_beams.at(i_beam).at(i_ant));
          }
          int coherent_sum_sqr_sample = coherent_sum_sample * coherent_sum_sample;
          coherent_sum += coherent_sum_sqr_sample;
        }
        if (coherent_sum > coherent_sum_max) {
          coherent_sum_max = coherent_sum;
        }
      }
      //h1->Fill(coherent_sum_max);
      coherent_values.push_back(coherent_sum_max);
    }
    window_count_this_signal++; //used to keep track of window vs. L1 triggers 
    window_count++; //every window counts for threshold evaluation regardless if triggered
  }
//*/
}

int  pueoSim::triggerThreshold::L2Threshold_eval(double sample_rate) {  

  double stdev =  TMath::StdDev(coherent_values.begin(), coherent_values.end() );
  double mean =  TMath::Mean(coherent_values.begin(), coherent_values.end() );
  double size =  coherent_values.size();
  std::cout << "\nCount = " << size << ", STDeV = " << stdev << ", Mean = " << mean << "\n\n";

  TCanvas *g2 = new TCanvas("g2","g2",500,500,600,400);
  g2->cd();


  TH1D * h1 = new TH1D("Histogram for L2","Histogram for L2",100,0.,8 * stdev + mean);
  for (int i = 0; i < coherent_values.size(); i++) {
    h1->Fill(coherent_values.at(i));
  }

  TF1 * f1 = new TF1("f1","expo");
  f1->SetParameters(1,1);
  h1->Fit("f1","","", mean + 1 * stdev, mean + 6 * stdev);
  TF1 * fitFunc = h1->GetFunction("f1");
  //std::cout << "p0=" << fitFunc->GetParameter(0) << " p1=" << fitFunc->GetParameter(1) ;

  double sector_rate = 12;
  double step = 8;
  int l2threshold = round((std::log(sector_rate * step / sample_rate *  window_count * (1 - exp((8 * stdev + mean)/100 *  fitFunc->GetParameter(1))) ) - fitFunc->GetParameter(0) ) / fitFunc->GetParameter(1));
  //std::cout<< "\n" << "L2 threshold " << l2threshold << "\n";

  h1->Draw();
  gPad->SetLogy();


  return l2threshold;
}