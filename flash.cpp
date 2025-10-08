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

#define PRINT false
#define deltaT_max 10  // unidades do seu timestamp

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
    if (argc > 1) 
    {
        file_name = std::string(argv[1]);
        std::cout << file_name << std::endl;
    } 
    else 
    {
        file_name = std::string("/home/gabriel/Documents/protodune/data/VD/np02vd_raw_run039510_0000_df-s04-d0_dw_0_20250919T123428_myAnalyser.root");
    }

    TFile* file = TFile::Open(file_name.c_str(), "READ");
    if (!file || file->IsZombie()) 
    {
        std::cerr << "Erro ao abrir o arquivo!\n";
        return 1;
    }

    TTree* tree = (TTree*)file->Get("waveAna/waveform_tree");
    if (!tree) 
    {
        std::cerr << "Não achei a TTree!\n";
        return 1;
    }

    std::vector<short>* signal = nullptr;
    int channel = 0;
    long timestamp = 0;

    tree->SetBranchAddress("adc", &signal);
    tree->SetBranchAddress("offline_channel", &channel);
    tree->SetBranchAddress("timestamp", &timestamp);

    const int n = tree->GetEntries();
    std::cout << "O arquivo contém " << n << " entradas\n";

    // Lê headers para vetor "hits" e detecta se já está ordenado
    std::vector<Hit> hits;
    hits.reserve(n);

    bool already_sorted = true;
    Long64_t last_t = std::numeric_limits<Long64_t>::min();

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

    if (!already_sorted) 
    {
        std::sort(hits.begin(), hits.end(),[](const Hit& a, const Hit& b) { return a.time < b.time;});
    }

    // --- Saída ROOT otimizada ---
    TFile* newfile = TFile::Open("data_analysed_flash.root", "RECREATE");
    auto* tree_write = new TTree("T1", "data");
    auto* data = new my_data();
    tree_write->Branch("Data", "my_data", &data);

    // --- Janela deslizante para flash ---
    int flash_ID = 0;
    size_t start = 0;
    while (start < hits.size()) 
    {
        size_t end = start;
        std::unordered_set<int> this_flash_channels;
        const Long64_t min_t = hits[start].time;

        while (end < hits.size() && (hits[end].time - min_t) <= deltaT_max) 
        {
            if (this_flash_channels.insert(hits[end].channel).second)          
            {
                data->set_parameters(hits[end].channel, hits[end].time, hits[end].adc);             
                data->calc_baseline(50);
                data->calc_noise(50);
                data->calc_t0(45, 120);
                data->calc_tend(400, 1000);
                data->calc_amplitude();
                data->calc_integral();
                data->calc_prompt(100);
                data->flash = flash_ID;
                tree_write->Fill();
            }
            ++end;
        }

        if ((flash_ID & 0xFF) == 0) 
        { 
            tree_write->AutoSave("SaveSelf");
        }
        ++flash_ID;
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
