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
//#include <fftw3.h>
#include "fftw3.h"
#include <algorithm>
#include <random>
#include "TAxis.h" 
#include "TH1D.h" 
#include "TH1.h"
#include "FTPair.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TTree.h"

namespace pueoSim {

  /**
   * @class pueoTrigger
   * @brief Implements the PUEO trigger.
   */

  class pueoTrigger {
  public:
    pueoTrigger(float samplingFreqHz_input, int antenna_start, float theta_width_L1, float phi_width_L1, float theta_width_L2, float phi_width_L2, bool firOn);

    static const int n_ant_L1=8;
    static const int n_ant_L2=16;

    bool visualisationOn = false;
    bool firFilterOn = false;
    
    int first_antenna;
    int n_samples;
    int n_beams_L1;
    int n_beams_L2;
    std::vector<bool> L1_triggered_windows;
    bool L2_triggered = false;
    float scaling = 1.0;
    float samplingFreqHz;

    std::vector<int> L1_max_value;
    std::vector<int> L2_max_value;


    std::vector<TGraph> signals;
    std::vector<TGraph> signals_discrete;
    std::vector<TGraph> signals_filtered;

    std::vector<std::vector<int>> L1_beams;
    std::vector<std::vector<int>> L2_beams;
    std::vector<std::vector<int>> L1_L2_map;



    void setScaling(float multiplier);
    void get_beamsL1_simpleSeparation(std::vector<std::vector<int>> &L1_beams, float theta_width_L1, float phi_width_L1);
    void get_beamsL2_simpleSeparation(std::vector<std::vector<int>> &L2_beams, float theta_width_L2, float phi_width_L2);
    void newSignal(std::vector<nicemc::FTPair> input_signals);
    void digitize(int bits);
    void digitize_afterFilter(int bits);
    void firFilter();
    void firFilter_signal_to_fir();
    void l1Trigger(int step, int window, int threshold, int max_shift);
    void l2Trigger(int step, int window, int threshold, int max_shift);




  };

    /**
     * @class triggerThreshold
     * @brief Evaluate threshold for trigger.
     */

  class triggerThreshold {
    public:
        triggerThreshold(float samplingFreqHz_input, int antenna_start, float theta_width_L1, float phi_width_L1, float theta_width_L2, float phi_width_L2, bool firOn);


        pueoTrigger * ptrigger;
        std::vector<int> coherent_values;
        int window_count;

        void setTriggerScaling(float multiplier) ;
        void setFir(bool firFilterYes);

        void L1Threshold_addData(int step, int window, int max_shift, int digitize_bits, std::vector<nicemc::FTPair> input_signals);
        int L1Threshold_eval(double samplingFreqHz, int step, bool diagnosticOn);
        void L2Threshold_addData(int step, int window, int max_shift, int digitize_bits, std::vector<nicemc::FTPair> input_signals, int L1Threshold);
        int L2Threshold_eval(double samplingFreqHz, int step, bool diagnosticOn);
  };

  

}