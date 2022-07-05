#include <random>
#include "TAxis.h"
#include "TH1D.h"
#include "TH1.h"
#include "FTPair.h"
#include "TFitResult.h"
#include "TF1.h"
#include "trigger.h"

pueoSim::pueoTrigger::pueoTrigger(float samplingFreqHz_input){
  get_beamsL1_simpleSeparation(L1_beams);

  generate_beams_L1(L1_beams);
  n_beams_L1 = L1_beams.size();
  std::vector<int> antennas_L1 = {0,1,2,3,4,5,6,7};
  for (int i=0; i<n_beams_L1; i++) {
	L1_ants.push_back(antennas_L1);
  }

  generate_beams_L2(L2_beams, L1_L2_map);
  n_beams_L2 = L2_beams.size();
  std::vector<int> antennas_L2 = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  for (int i=0; i<n_beams_L2; i++) {
	L2_ants.push_back(antennas_L2);
  }

  samplingFreqHz = samplingFreqHz_input;

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

  //TTreeReader myReader("ntuple", tree);
  //TTreeReaderValue<Float_t> myPx(myReader, "An");

  Float_t An, X, Y;
  tree->SetBranchAddress("An",&An);

  //for (int i = 0, N = tree->GetEntries(); i < N; ++i) {
  //  tree->GetEntry(i) ;
  //  std::cout<< An;
  //  std::cout<< "\n";
  //}
}

void pueoSim::pueoTrigger::generate_beams_L1(std::vector<std::vector<int>> &L1_beams) {
  double w1 = 0.6; //horizontal separation in m between 2 phi sectors
  double h1 = 0.74; //vertical separation between adjacent antennas for lower 3 antennas in a phi sector
  double h2 = 3.0; //vertical separation between adjacent antennas between top two antennas in a given phi sector
  double h = 2 * h1 + h2;
  double c = 3e8;
  double delta_t = 1/(2.56e9); //time between samples, for 2.6Ghz sample frequency

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
  for (int i = v_min; i < v_max+1; i+=2) {
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
  for (int i = h_min; i < h_max+1; i+=2) {
    std::vector<int> beam;
    beam.push_back(0);
    beam.push_back(i);
    //std::cout << beam.at(0)  << "\t" << beam.at(1) << "\t" <<"\n";
    horizontal.push_back(beam);
  }

  //combine vertical and horizontal for beams

  //std::cout << "**********L1 Beams**********" <<"\n";

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

      //std::cout << "L1 Beam:\t" << L1_beams.size()<< "\t" << "Offsets:" << "\t"  << beam.at(0)  << "\t" << beam.at(1) << "\t"<< beam.at(2) << "\t"<< beam.at(3) << "\t"<< beam.at(4)  << "\t" << beam.at(5) << "\t"<< beam.at(6) << "\t"<< beam.at(7) << "\t" <<"\n";
      L1_beams.push_back(beam);
    }
  }


}

  ////first variable is vector, each element a vector representing an L2 beam, with 16 delay values in that vector
  ////second variable is vector, each element a vector representing an L1 beam, with corresponding L2 beams in that vector
void pueoSim::pueoTrigger::generate_beams_L2(std::vector<std::vector<int>> &L2_beams, std::vector<std::vector<int>> &L1_L2_map) {
  double w1 = 0.6; //horizontal separation in m between 2 adjacent azimuthal sectors
  double w = 3 * 0.6; //horizontal separation in m between azimuthal sectors furthest apart in L2 sector
  double h1 = 0.74; //vertical separation between adjacent antennas for lower 3 antennas in a //phi sector
  double h2 = 3.0; //vertical separation between adjacent antennas between top two antennas in //a given phi sector
  double h = 2 * h1 + h2;
  double c = 3e8;
  double delta_t = 1/(2.56e9); //time between samples, for 2.6Ghz sample frequency


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
  for (int i = v_min; i < v_max+1; i+=2) {
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
  for (int i = h_min; i < h_max+1; i+=2) {
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
  //std::cout << "**********L2 Beams**********" <<"\n";

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

      //std::cout << "L2 Beam:\t" << L2_beams.size() << "\t" << "Offsets:" << "\t" ;
      //for (int i = 0; i < beam.size(); i++) {
      //  std::cout << beam.at(i)  << "\t" ;
      //}
      //std::cout <<"\n";
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

  //std::cout << "**********L1 L2 Map**********" <<"\n";
  //for (int i=0; i < L1_L2_map.size(); ++i) {
  //  std::cout << "L1: " << i << "\t"<< "Mapped L2: " << "\t"  ;
  //  for (int j=0; j < L1_L2_map.at(i).size(); ++j) {
  //    std::cout << L1_L2_map.at(i).at(j) << "\t" ;
  //  }
  //  std::cout << "\n";
  //}
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
  it_ant=signals.begin();

  for (it_ant=signals.begin();it_ant!=signals.end(); ++it_ant) {

    signals_discrete.emplace_back(*it_ant);
    TGraph *gr_digi = & signals_discrete.back();

    const int n = it_ant->GetN();
    int yout[n];

    for (int i=0; i<gr_digi->GetN(); i++) {
      int digitised_y;
      double scaled_y = gr_digi->GetPointY(i) * scaling;
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
      total_shifted[samp_pos]=total_shifted[samp_pos]*total_shifted[samp_pos];
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
        total_shifted[samp_pos]=total_shifted[samp_pos]*total_shifted[samp_pos];
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

    //Now cycle through windows and record if a winow triggers. Note that as long as any beam triggers, the trigger for that window is set to //true:
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

pueoSim::triggerThreshold::triggerThreshold(float samplingFreqHz_input) {
  window_count = 0;
  ptrigger = new pueoTrigger(samplingFreqHz_input);
  
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

int pueoSim::triggerThreshold::L1Threshold_eval() {  

  double stdev =  TMath::StdDev(coherent_values.begin(), coherent_values.end() );
  double mean =  TMath::Mean(coherent_values.begin(), coherent_values.end() );
  double size =  coherent_values.size();
  std::cout << "\nCount = " << size << ", STDeV = " << stdev << ", Mean = " << mean << "\n\n";

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
  double sample_rate = 2.56E9;
  double step = 8;

  
  return round((std::log(sector_rate * step / sample_rate / ptrigger->n_beams_L1 *  window_count ) - fitFunc->GetParameter(0) ) / fitFunc->GetParameter(1));

}


void  pueoSim::triggerThreshold::L2Threshold_addData(std::vector<nicemc::FTPair> input_signals, int L1Threshold) {
  ptrigger->newSignal(input_signals);
  ptrigger->digitize(4);
  ptrigger->l1Trigger(8, 16, L1Threshold, 64);
  
  

  int n_samples = ptrigger->signals_discrete.at(0).GetN();
  int number_ants = 16;  
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
          for (int i_ant=0; i_ant<number_ants; i_ant+=1) {
            coherent_sum_sample += ptrigger->signals_discrete.at(ptrigger->L2_ants.at(i_beam).at(i_ant)).GetPointY(wind_pos+samp_pos-ptrigger->L2_beams.at(i_beam).at(i_ant));
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

int  pueoSim::triggerThreshold::L2Threshold_eval() {  

  double stdev =  TMath::StdDev(coherent_values.begin(), coherent_values.end() );
  double mean =  TMath::Mean(coherent_values.begin(), coherent_values.end() );
  double size =  coherent_values.size();
  std::cout << "\nCount = " << size << ", STDeV = " << stdev << ", Mean = " << mean << "\n\n";

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
  double sample_rate = 2.56E9;
  double step = 8;
  int l2threshold = round((std::log(sector_rate * step / sample_rate *  window_count ) - fitFunc->GetParameter(0) ) / fitFunc->GetParameter(1));
  //std::cout<< "\n" << "L2 threshold " << l2threshold << "\n";

  h1->Draw();
  gPad->SetLogy();


  return l2threshold;
}