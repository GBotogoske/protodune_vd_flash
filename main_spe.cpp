#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <TCanvas.h>
#include <TApplication.h>

#include <utils.hh>
#include <my_data.hh>
#include <filesystem>
namespace fs = std::filesystem;

#define PRINT false


//channel list para calcular os paremetos, com oum valor default setado
std::vector<int> ch_list={1054};

int main(int argc, char** argv)
{
    std::string file_name="null";
    //diretorio base para procurar o arquivo
    std::string base_dir = "/home/gabriel/Documents/protodune/data/VD/";
    std::string run_number = "null";

    if (argc > 1) 
    {
        bool run_found = false;
        for (int i = 1; i < argc; ++i)
        {
            std::string arg = argv[i];
            if (arg == "-run" && i + 1 < argc)//procurar o arquivo caso -run seja dada
            {
                run_number = argv[i + 1];
                run_found = true;
                break;
            }
        }
        if(run_found)
        {
            file_name=search_file(base_dir,run_number,std::string("myAnalyser.root")); //procurar o arquivo caso -run seja dada
        }
        else
        {
            file_name = std::string(argv[1]); //caso contrario usar o arquivo dado
        }
    }  
    else //e caso contrario ainda, usar o nome padrao abaixo
    {
        file_name = std::string("/home/gabriel/Documents/protodune/data/VD/np02vd_raw_run039358_0000_df-s04-d0_dw_0_20250915T151645_myAnalyser.root");
    }

    std::cout << file_name << std::endl;
    TApplication app("app", &argc, argv);
    //abrir arquivo para leitura
    TFile *file = TFile::Open(file_name.c_str(), "READ");
    if (!file || file->IsZombie()) 
    {
        std::cerr << "Erro ao abrir o arquivo!" << std::endl;
        return 1;
    }

    //abrir ttree
    TTree *tree = (TTree*)file->Get("waveAna/waveform_tree");
    if (!tree) 
    {
        std::cerr << "Não achei a TTree no arquivo!" << std::endl;
        return 1;
    }
    //pegar o numero de entradas
    int n = tree->GetEntries();
    std::cout << "O arquivo contem " << n << " entradas\n";

    //variaveis de leitura
    std::vector<short>* signal = nullptr; 
    int channel;
    long time;
    //faz as variaveis apontarem no lugar correto da arvore
    tree->SetBranchAddress("adc",&signal);
    tree->SetBranchAddress("offline_channel",&channel);
    tree->SetBranchAddress("timestamp",&time);

    //cria uma class my data para salvar os paremetros
    auto data = new my_data();

    //cria arquivo de saida
    std::string file_output;
    if(run_number == "null")
    {
        file_output = fs::path(file_name).replace_extension("").string() + "_spe.root";
    }
    else
    {   
        file_output=fs::path(file_name).parent_path() / ("data_analysed_" + run_number + "_spe.root");
        //pega o ch_list adequado dessa run
        ch_list = get_map_spe()[run_number];
    }
    std::cout << file_output << std::endl;
    TFile* newfile = TFile::Open(file_output.c_str(), "RECREATE");
    //cria arvore de saida
    TTree *tree_write = new TTree("T1","data");
    //aponta corretamente o my_data 
    tree_write->Branch("Data", "my_data", &data); 

    //cria canvas que vai ser usado se a opcao print estiver setada
    auto c = new TCanvas("my Canvas", "Waveforms",800,600);

    //varre todas as entradas e calcula os paremetros
    for(int i=0;i<n;i++)
    {
        tree->GetEntry(i);
        if (std::find(ch_list.begin(), ch_list.end(), channel) != ch_list.end()) 
        {
            data->set_parameters(channel,time,*signal);
            data->calc_baseline(0,180,480,1024);
            data->calc_noise(200);
            data->calc_t0(45,200);
            data->calc_tend(data->t0+15,1000);
            data->calc_amplitude(255,1024);
            data->calc_integral(255,280);
            data->calc_prompt(100);
            //std::cout << signal->size() << std::endl;
            tree_write->Fill();

            //aqui printa a waveform se PRINT = true
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

    tree_write->Write("", TObject::kOverwrite);
    newfile->Write();
    newfile->Close();
 
    return 0;
}