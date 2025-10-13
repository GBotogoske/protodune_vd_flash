#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <string>

#include <utils.hh>
#include <my_data.hh>
#include <filesystem>
namespace fs = std::filesystem;

#define deltaT_max 10  // unidades do seu timestamp

// estrutura par organizar os hits
struct Hit
{
    int index;
    int channel;
    long time;
    std::vector<short> adc; 
};

int main(int argc, char** argv)
{
    std::string file_name;
    std::string base_dir = "/home/gabriel/Documents/protodune/data/VD/"; // base dir dos arquivos
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
            file_name=search_file(base_dir,run_number,std::string("myAnalyser.root"));
        }
        else
        {
            file_name = std::string(argv[1]); // caso tenha passado o nome do arquivo
        }
    } 
    else // arquivo default
    {
        file_name = std::string("/home/gabriel/Documents/protodune/data/VD/np02vd_raw_run039510_0000_df-s04-d0_dw_0_20250919T123428_myAnalyser.root");
    }
    std::cout << "FILE: " << file_name << std::endl;
    std::cout << "RUN: " << run_number << std::endl;

    //abre arquivo de leitura
    TFile* file = TFile::Open(file_name.c_str(), "READ");
    if (!file || file->IsZombie()) 
    {
        std::cerr << "Erro ao abrir o arquivo!\n";
        return 1;
    }

    //abre arvore de leitura
    TTree* tree = (TTree*)file->Get("waveAna/waveform_tree");
    if (!tree) 
    {
        std::cerr << "Não achei a TTree!\n";
        return 1;
    }

    //variasveis de leituras
    std::vector<short>* signal = nullptr;
    int channel = 0;
    long timestamp = 0;
    //faz as varivaeis apontarem para as variaveis da arvore
    tree->SetBranchAddress("adc", &signal);
    tree->SetBranchAddress("offline_channel", &channel);
    tree->SetBranchAddress("timestamp", &timestamp);
    //numero de entradas
    const int n = tree->GetEntries();
    std::cout << "O arquivo contém " << n << " entradas\n";

    // Lê headers para vetor "hits" e detecta se já está ordenado
    std::vector<Hit> hits;
    hits.reserve(n);

    bool already_sorted = true;
    Long64_t last_t = std::numeric_limits<Long64_t>::min();

    //le todas as entrradas e coloca no vetor de hits
    for (int i = 0; i < n; ++i) 
    {
        tree->GetEntry(i);
        hits.push_back({i, channel, timestamp, *signal});
        if (timestamp < last_t) 
        {
            already_sorted = false;
        }
        last_t = timestamp;
    }

    //ordenas os hits por tempo
    if (!already_sorted) 
    {
        std::sort(hits.begin(), hits.end(),[](const Hit& a, const Hit& b) { return a.time < b.time;});
    }

    // --- Cria arquivo de saida  ---
    std::string file_output;
    if(run_number == "null")
    {
        file_output = fs::path(file_name).replace_extension("").string() + "_flash.root";
    }
    else
    {   
        file_output=fs::path(file_name).parent_path() / ("data_analysed_" + run_number + "_spe.root");
    }
    TFile* newfile = TFile::Open(file_output.c_str(), "RECREATE");
    std::cout << newfile << std::endl;

    //cria arvore de saida
    auto* tree_write = new TTree("T1", "data");
    //aponta a variavel my_data para a arvore
    auto* data = new my_data();
    tree_write->Branch("Data", "my_data", &data);

    // --- Varre os hits, agrupando em flashs ---
    int flash_ID = 0;
    size_t start = 0;
    while (start < hits.size()) 
    {
        size_t end = start;
        //mapa para nao repetir canal
        std::unordered_set<int> this_flash_channels;
        //min_t do flash
        const Long64_t min_t = hits[start].time;

        //varre todos ou ate a diferenca for maior que o maximo aceitado
        while (end < hits.size() && (hits[end].time - min_t) <= deltaT_max) 
        {
            //checa se nao colocou antes
            if(this_flash_channels.insert(hits[end].channel).second)          
            {
                //calcula parametros
                data->set_parameters(hits[end].channel, hits[end].time, hits[end].adc);             
                data->calc_baseline(50);
                data->calc_noise(50);
                data->calc_t0(45,200);
                data->calc_tend(data->t0+15,1000);
                data->calc_amplitude();
                data->calc_integral();
                data->calc_prompt(100);
                data->flash = flash_ID;
                tree_write->Fill();
            }
            end++;
        }

        if ((flash_ID & 0xFF) == 0) 
        { 
            tree_write->AutoSave("SaveSelf");
        }
        flash_ID++;
        //comeco o proximo loop, onde parou o interno
        start = end;
    }

    tree_write->Write("", TObject::kOverwrite);
    newfile->Close();
    delete data;
    delete newfile;

    file->Close();
    delete file;

    return 0;
}
