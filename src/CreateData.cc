

#include <iostream>
#include "TGraph.h"
#include <vector>
#include <math.h>
#include <iostream>
#include <tuple>
#include <string>  
#include <istream>
#include <iterator>
#include <fstream>
#include <sstream>




class CSVRow
{
    public:
        std::string_view operator[](std::size_t index) const
        {
            return std::string_view(&m_line[m_data[index] + 1], m_data[index + 1] -  (m_data[index] + 1));
        }
        std::size_t size() const
        {
            return m_data.size() - 1;
        }
        void readNextRow(std::istream& str)
        {
            std::getline(str, m_line);

            m_data.clear();
            m_data.emplace_back(-1);
            std::string::size_type pos = 0;
            while((pos = m_line.find(',', pos)) != std::string::npos)
            {
                m_data.emplace_back(pos);
                ++pos;
            }
            // This checks for a trailing comma with no data after it.
            pos   = m_line.size();
            m_data.emplace_back(pos);
        }
    private:
        std::string         m_line;
        std::vector<int>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}   




int main() {
    
    //load ANITA3 Impulse Response

    float PUEO_IR_x[400], PUEO_IR_y[400];

    std::ifstream       file("C:/Users/Humma/Dropbox/PhD/Adhoc python scripts/PUEO_signal_10ghz.csv");
    CSVRow              row;
    int counter = 0;
    while(file >> row)
    {
        PUEO_IR_x[counter] = (std::stod(std::string((row[0]))));
        PUEO_IR_y[counter] = (std::stod(std::string((row[1]))));
        std::cout << std::to_string(counter) << "\t" << row[0] << "\t" << row[1] << "\n";
        counter++;
    }
        
   // TGraph * gr3 = new TGraph(400, PUEO_IR_x,PUEO_IR_y);
   // gr3 ->Draw();

    //TCanvas *c1 = new TCanvas("c1","A3 IR",500,500,600,400);
    //TMultiGraph *mg = new TMultiGraph();
    //gr3->SetTitle("Anita3 Impulse Response");
    //mg->Add(gr3);   
    //c1->cd();
    //mg->Draw("A pmc plc");  
    //c1->BuildLegend();




}


/****
//generate testing data - electric field at given angle and freq; frequency in Mhz
double e_field(double theta_offcone, double freq_mhz, double normalisation) {
  if (freq_mhz <1 || freq_mhz > 2000) {
    return 0;
  } else{
    double theta_cerenkov =  41; //degrees
    double cone_width = 2.2 * 1000 / freq_mhz; //degrees
    double theta_view = theta_cerenkov - theta_offcone;
    double ice_atten = exp(-500/(440 + (-120)/(1000-200)*(f-200)));
    return  normalisation * sin(theta_view/180*M_PI) / sin(theta_cerenkov/180*M_PI)   * exp(-(theta_offcone/cone_width)**2) * 2.53E-7 * freq_mhz / 1150 * pow((1+pow(freq_mhz / 1150, 1.44)),-1)
  }
        
}
  

void e_field_test() {
  
  std::vector<double> freq_list;
  for (double i=1, i<2e3,i+=10) {
    freq_list.push_back(i);
  }
  std::vector<double> off_cone_angles {0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 6.0}

  std::vector<vector<double>> 



}
****/
