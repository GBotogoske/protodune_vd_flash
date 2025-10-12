#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <TCanvas.h>
#include <TApplication.h>

#include <utils.hh>
#include <my_data.hh>

#define PRINT false

std::vector<int> ch_list={1050};

int main(int argc, char** argv)
{
    
    std::string file_name;
    if(argc > 1)
    {
        file_name = std::string(argv[1]);
        std::cout << file_name << std::endl;
    }
    else
    {
        file_name = std::string("/home/gabriel/Documents/protodune/data/VD/np02vd_raw_run039357_0000_df-s04-d0_dw_0_20250915T151645_myAnalyser.root");
    }
    TApplication app("app", &argc, argv);
    TFile *file = TFile::Open(file_name.c_str(), "READ");
    if (!file || file->IsZombie()) 
    {
        std::cerr << "Erro ao abrir o arquivo!" << std::endl;
        return 1;
    }

    TTree *tree = (TTree*)file->Get("waveAna/waveform_tree");
    if (!tree) 
    {
        std::cerr << "Não achei a TTree no arquivo!" << std::endl;
        return 1;
    }

    int n = tree->GetEntries();
    std::cout << "O arquivo contem " << n << " entradas\n";

    std::vector<short>* signal = nullptr; 
    int channel;
    long time;
    tree->SetBranchAddress("adc",&signal);
    tree->SetBranchAddress("offline_channel",&channel);
    tree->SetBranchAddress("timestamp",&time);

    auto data = new my_data();

    TFile *newfile = TFile::Open("data_analysed.root","RECREATE");
    TTree *tree_write = new TTree("T1","data");
    tree_write->Branch("Data", "my_data", &data); 

    auto c = new TCanvas("my Canvas", "Waveforms",800,600);

    for(int i=0;i<n;i++)
    {
        tree->GetEntry(i);
        if (std::find(ch_list.begin(), ch_list.end(), channel) != ch_list.end()) 
        {
            data->set_parameters(channel,time,*signal);
            data->calc_baseline(200);
            data->calc_noise(200);
            data->calc_t0(45,200);
            data->calc_tend(data->t0+15,1000);
            data->calc_amplitude(255,280);
            data->calc_integral(255,280);
            data->calc_prompt(100);
            //std::cout << signal->size() << std::endl;
            tree_write->Fill();

       
            if(PRINT)
            {
                
                print_waveform(c,data,i);
                std::cout << "Press ENTER para continuar, ou 'c' + ENTER para sair..." << std::endl;
                std::string input;
                std::getline(std::cin, input);
                if (input == "c" || input == "C") 
                {
                    break; 
                }
            }
        }
    }

    newfile->Write();
    newfile->Close();
 
    return 0;
}