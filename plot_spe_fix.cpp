#include <iostream>
#include <fstream> 
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TAxis.h>

#include <filesystem>
namespace fs = std::filesystem;

#include <utils.hh>
#include <readConfigFit.hh>
#include <my_data.hh>

//channel default to calculate spe
std::vector<std::string> list_channel = {"1050","1051","1060","1061","1070","1071","1080","1081","2030","2031","2040","2041","2070","2071","2080","2081"};
//separation default 

#define bins 1024

int main(int argc, char** argv)
{
    //home folder
    std::string base_dir = "/home/gabriel/Documents/protodune/protodune_vd/light_analysis/waveforms/";

    std::string run_number;
    for(int i=0;i<list_channel.size();i++)
    {
        run_number = get_map_spe2()[list_channel[i]];
        std::string file_name1,file_name2;
        file_name1=base_dir+run_number+"/waveform_"+list_channel[i]+".txt";
        file_name2=base_dir+run_number+"/waveform_"+list_channel[i]+"_filtered.txt";

        std::vector<double> x1, signal1;
        std::vector<double> x2, signal2;
        double t, s;
        //ler arquivos
        std::ifstream f1(file_name1);
        if (!f1.is_open())
        {
            std::cerr << "Erro ao abrir: " << file_name1 << std::endl;
            continue;
        }
       
        while (f1 >> t >> s)
        {
            x1.push_back(t);
            signal1.push_back(s);
        }
        f1.close();

          //ler arquivos
        std::ifstream f2(file_name2);
        if (!f2.is_open())
        {
            std::cerr << "Erro ao abrir: " << file_name2 << std::endl;
            continue;
        }
  
        while (f2 >> t >> s)
        {
            x2.push_back(t);
            signal2.push_back(s);
        }
        f2.close();


        TCanvas* c = new TCanvas("c", "Waveforms", 800, 600);
        TGraph* g1 = new TGraph(bins, x1.data(), signal1.data());
        g1->SetTitle(Form("Waveform - Channel %d", std::stoi(list_channel[i])));
        g1->GetXaxis()->SetTitle("Sample");
        g1->GetYaxis()->SetTitle("Amplitude (ADC)");
        g1->SetLineColor(kRed);
        g1->Draw("AL");

        TGraph* g2 = new TGraph(bins, x2.data(), signal2.data());
        g2->SetTitle(Form("Waveform - Channel %d", std::stoi(list_channel[i])));
        g2->GetXaxis()->SetTitle("Sample");
        g2->GetYaxis()->SetTitle("Amplitude (ADC)");
        g2->SetLineColor(kBlue);
        g2->Draw("SAME");

        c->SaveAs(std::string(base_dir+run_number+"/waveform_"+list_channel[i]+"_FIXED.png").c_str());

        delete g1;
        delete g2;
        delete c;
    }
 
    
  
   
    
    return 0;
}