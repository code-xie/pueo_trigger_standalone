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
    pueoTrigger(float samplingFreqHz_input);

    static const int n_ant_L1=8;
    static const int n_ant_L2=16;
    
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

    std::vector<std::vector<int>> L1_ants; //beam number, antenna within beam
    std::vector<std::vector<int>> L1_beams;
    std::vector<std::vector<int>> L2_ants;
    std::vector<std::vector<int>> L2_beams;
    std::vector<std::vector<int>> L1_L2_map;



    void setScaling(float multiplier);
    void generate_beams_L1(std::vector<std::vector<int>> &L1_beams);
    void get_beamsL1_simpleSeparation(std::vector<std::vector<int>> &L1_beams);
    void generate_beams_L2(std::vector<std::vector<int>> &L2_beams, std::vector<std::vector<int>> &L1_L2_map);
    void get_beamsL2_simpleSeparation(std::vector<std::vector<int>> &L2_beams);
    void newSignal(std::vector<nicemc::FTPair> input_signals);
    void digitize(int bits);
    void l1Trigger(int step, int window, int threshold, int max_shift);
    void l2Trigger(int step, int window, int threshold, int max_shift);




  };

    /**
     * @class triggerThreshold
     * @brief Evaluate threshold for trigger.
     */

  class triggerThreshold {
    public:
        triggerThreshold(float samplingFreqHz_input);


        pueoTrigger * ptrigger;
        std::vector<int> coherent_values;
        int window_count;

        void setTriggerScaling(float multiplier) ;
        void L1Threshold_addData(std::vector<nicemc::FTPair> input_signals);
        int L1Threshold_eval();
        void L2Threshold_addData(std::vector<nicemc::FTPair> input_signals, int L1Threshold);
        int L2Threshold_eval();
  };

  

}