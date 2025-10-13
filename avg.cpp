#include <iostream>
#include <fstream> 
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TAxis.h>

#include <utils.hh>
#include <my_data.hh>
#include <filesystem>
namespace fs = std::filesystem;

//lista de canais para calcular o sinal medio
std::vector<int> list_channel = {1050,1051,1060,1061,1070,1071,1080,1081,2030,2031,2040,2041,2050,2051,2080,2081};

int main(int argc, char** argv)
{
    std::string file_name;
    //base folder
    std::string base_dir = "/home/gabriel/Documents/protodune/data/VD/";
    std::string run_number = "null";

    if (argc > 1) 
    {
        bool run_found = false;
        for (int i = 1; i < argc; ++i)
        {
            std::string arg = argv[i];
            if (arg == "-run" && i + 1 < argc) // se rodou com -run procura o arquivo correspondente
            {
                run_number = argv[i + 1];
                run_found = true;
                break;
            }
        }
        if(run_found) // procura o arquivo
        {
            file_name=search_file(base_dir,run_number,std::string("_flash_filtered.root"));
        }
        else // se nao mandou flag, le o arquivo mandado
        {
            file_name = std::string(argv[1]);
        }
    } 
    else //arquivo default
    {
        file_name = std::string("/home/gabriel/Documents/protodune/data/VD/np02vd_raw_run039510_0000_df-s04-d0_dw_0_20250919T123428_myAnalyser.root");
    }
    std::cout << "FILE: " << file_name << std::endl;
    std::cout << "RUN: " << run_number << std::endl;

    //abre arquivo para leitura
    TFile *file = TFile::Open(file_name.c_str(), "READ");
    if (!file || file->IsZombie()) 
    {
        std::cerr << "Erro ao abrir o arquivo!" << std::endl;
        return 1;
    }
    //abre arvore para leitura
    TTree* tree = (TTree*)file->Get("T1");
    if (!tree) 
    {
        std::cerr << "Não achei a TTree no arquivo!" << std::endl;
        return 1;
    }
    //numero de entradas
    int n = tree->GetEntries();
    std::cout << "O arquivo contem " << n << " entradas\n";

    //cria um mapa èara vincular o canal do larsoft do canal, para o indice do canal do list_channel
    std::map<int,int> ch_map;
    int n_ch=list_channel.size();
    std::vector<double> norm_ch(n_ch, 0.0); 
    
    for(int i=0;i<n_ch;i++)
    {
        ch_map[list_channel[i]]=i;
    }

    //vetores nulo para cada canal
    std::vector<std::vector<double>> signals(n_ch, std::vector<double>(1024));
    my_data* data = nullptr;
    tree->SetBranchAddress("Data", &data);

    //varre a arvore
    for(int i = 0; i < tree->GetEntries(); i++) 
    {
        tree->GetEntry(i);
        int ch=data->Channel;
        bool found=false;
        for(int j=0;j<n_ch;j++)
        {
            //procura o canal na lista de canai
            if(ch==list_channel[j])
            {
                found=true;
                break;
            }
        }
        if(found) // se encontrou o canal
        {
            // adiciona na media
            for (int k = 0; k < data->adcs.size(); k++)
            {
                norm_ch[ch_map[ch]]++;
                signals[ch_map[ch]][k] += (data->adcs[k] - data->baseline);
            }
        }        
    }
    //fecha o arquivo
    file->Close();
    delete data;

    //vetor de tempo
    std::vector<double> x(1024);
    for(int i=0;i<1024;i++)
    {
        x[i]=i;
    }

    //normaliza as waveforms
    for(int i=0; i<n_ch; i++)
    {
        for (int k = 0; k < 1024; k++)
        {
            signals[i][k] = signals[i][k]/norm_ch[i];
        }
    }

    //dir base para salvar
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


    //cria o vetor medio para cada canal e salva png e txt
    for(int i=0; i<n_ch; i++)
    {
        TCanvas* c = new TCanvas("c", "Waveforms", 800, 600);

        int ch = list_channel[i];
        TGraph* g = new TGraph(1024, x.data(), signals[i].data());
        g->SetTitle(Form("Waveform - Channel %d", ch));
        g->GetXaxis()->SetTitle("Sample");
        g->GetYaxis()->SetTitle("Amplitude (ADC)");
        g->Draw("AL");

        std::string fig_name = out_base_dir+ Form("/waveform_%d.png", ch);
        c->SaveAs(fig_name.c_str());

        std::ofstream fout(out_base_dir+Form("/waveform_%d.txt", ch));
        for (int k = 0; k < 1024; k++)
        {
            fout << x[k] << "\t" << signals[i][k] << "\n";
        }
        fout.close();

        delete g;
        delete c;
    }

    return 0;
}