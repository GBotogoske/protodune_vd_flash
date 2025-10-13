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
#include <my_data.hh>

std::vector<int> list_channel = {1050};
std::vector<double> separation = {-5.25901,415.99,825.73,1256.75};

int main(int argc, char** argv)
{
    std::string file_name;
    std::string base_dir = "/home/gabriel/Documents/protodune/data/VD/";
    std::string run_number = "null";
    std::string this_ch="39215";

    if (argc > 1) 
    {
        bool run_found = false;
        for (int i = 1; i < argc; ++i)
        {
            std::string arg = argv[i];
            if (arg == "-ch" && i + 1 < argc)
            {
                this_ch = argv[i + 1];
                run_number = get_map_spe2()[this_ch];
                run_found = true;
                break;
            }
        }
        if(run_found)
        {
            file_name=search_file(base_dir,run_number,std::string("_spe.root"));
        }
        else
        {
            file_name = std::string(argv[1]);
        }
    } 
    else 
    {
        file_name = std::string("/home/gabriel/Documents/protodune/data/VD/np02vd_raw_run039357_0000_df-s04-d0_dw_0_20250915T151645_myAnalyser.root");
    }

    std::cout << "FILE NAME: " << file_name << std::endl;
    std::cout << "RUN: " << run_number << std::endl;
    std::cout << "CH: " << this_ch << std::endl;

    TFile *file = TFile::Open(file_name.c_str(), "READ");
    if (!file || file->IsZombie()) 
    {
        std::cerr << "Erro ao abrir o arquivo!" << std::endl;
        return 1;
    }
    TTree* tree = (TTree*)file->Get("T1");
    if (!tree) 
    {
        std::cerr << "Não achei a TTree no arquivo!" << std::endl;
        return 1;
    }
    int n = tree->GetEntries();
    std::cout << "O arquivo contem " << n << " entradas\n";

    int norm_ch=0; 
    
    std::vector<double> signal = std::vector<double>(1024, 0.0);
    my_data* data = nullptr;
    tree->SetBranchAddress("Data", &data);
    int my_channel = std::stoi(this_ch);

    for(int i = 0; i < tree->GetEntries(); i++) 
    {
        tree->GetEntry(i);
        int ch=data->Channel;
       
        if (data->Channel==my_channel)
        {
            if(data->baseline>4000 && data->baseline<5000 && data->noise<10) // && data->amplitude<100)
            {
                double integral=data->integral;
                if(!(integral<separation[0] || integral>separation[separation.size()-1]))
                {
                    int norm;
                    for(norm=0;norm<separation.size()-1;norm++)
                    {
                        if(integral>separation[norm] && integral<separation[norm+1])
                        {
                            break;
                        }
                    }
                    norm=norm+1;
                    norm_ch+=norm;
                    for (int k = 0; k < data->adcs.size(); k++)
                    {
                        signal[k] += (data->adcs[k] - data->baseline);
                    }
                }
            }
                
        }
            
               
    }

    file->Close();
    delete data;

    std::vector<double> x(1024);
    for(int i=0;i<1024;i++)
    {
        x[i]=i;
    }

    for (int k = 0; k < 1024; k++)
    {
        signal[k] = signal[k]/norm_ch;
    }
  
    TCanvas* c = new TCanvas("c", "Waveforms", 800, 600);

    int ch = my_channel;
    TGraph* g = new TGraph(1024, x.data(), signal.data());
    g->SetTitle(Form("Waveform - Channel %d", ch));
    g->GetXaxis()->SetTitle("Sample");
    g->GetYaxis()->SetTitle("Amplitude (ADC)");
    g->Draw("AL");


    std::string out_base_dir;
    if (run_number == "null")
    {
        out_base_dir = "../waveforms/" + fs::path(file_name).filename().string();
    }
    else
    {
        out_base_dir = "../waveforms/" + run_number;
    }
    
    fs::create_directory(out_base_dir);


    std::string fig_name = out_base_dir+ Form("/waveform_%d.png", ch);
    c->SaveAs(fig_name.c_str());

    std::ofstream fout(out_base_dir+Form("/waveform_%d.txt", ch));
    for (int k = 0; k < 1024; k++)
    {
        fout << x[k] << "\t" << signal[k] << "\n";
    }
    fout.close();

    delete g;
    delete c;
    

    return 0;
}