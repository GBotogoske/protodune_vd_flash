#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <cmath>
#include <filesystem>
#include <vector>
#include <numeric>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TApplication.h>
#include <TAxis.h>

#include <utils.hh>
#include <readConfigFit.hh>
#include <my_data.hh>

namespace fs = std::filesystem;

// Canal padrão
std::vector<int> list_channel = {1050};
// Separações padrão (picos SPE)
std::vector<double> separation = {-5.25901, 415.99, 825.73, 1256.75};

// =========================================================
// Função auxiliar: moda
// =========================================================
double calc_moda(const std::vector<double>& values, double bin_width = 0.5)
{
    if (values.empty()) return 0.0;

    std::unordered_map<int, int> freq;
    freq.reserve(values.size());

    for (double v : values)
    {
        int bin = static_cast<int>(std::floor(v / bin_width));
        freq[bin]++;
    }

    int best_bin = 0;
    int best_count = 0;
    for (auto& [bin, count] : freq)
    {
        if (count > best_count)
        {
            best_count = count;
            best_bin = bin;
        }
    }

    // retorna o centro do bin de maior frequência
    return (best_bin * bin_width) + (bin_width / 2.0);
}

// =========================================================
// Programa principal
// =========================================================
int main(int argc, char** argv)
{
    std::string file_name;
    std::string base_dir = "/home/gabriel/Documents/protodune/data/VD/";
    std::string run_number = "null";
    std::string this_ch = "-1";

    readConfigFit* fitConfig = nullptr;

    // ---------- Processa argumentos ----------
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

        if (run_found)
        {
            file_name = search_file(base_dir, run_number, "_spe.root");
            fitConfig = new readConfigFit(std::stoi(this_ch), 1);
        }
        else
        {
            file_name = std::string(argv[1]);
            fitConfig = new readConfigFit(-1);
        }
    }
    else
    {
        file_name = "/home/gabriel/Documents/protodune/data/VD/np02vd_raw_run039357_0000_df-s04-d0_dw_0_20250915T151645_myAnalyser.root";
        fitConfig = new readConfigFit(-1);
    }

    std::cout << "FILE NAME: " << file_name << std::endl;
    std::cout << "RUN: " << run_number << std::endl;
    std::cout << "CH: " << this_ch << std::endl;

    // ---------- Abre arquivo ROOT ----------
    TFile* file = TFile::Open(file_name.c_str(), "READ");
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
    std::cout << "O arquivo contém " << n << " entradas\n";

    // ---------- Inicializa variáveis ----------
    int norm_ch = 0;
    std::vector<std::vector<double>> all_signals(1024);

    my_data* data = nullptr;
    tree->SetBranchAddress("Data", &data);

    int my_channel = std::stoi(this_ch);
    float baseline_min = fitConfig->getParam("baseline_min");
    float baseline_max = fitConfig->getParam("baseline_max");
    float noise_max = fitConfig->getParam("noise_max");
    float amplitude_min = fitConfig->getParam("amplitude_min");
    float amplitude_max = fitConfig->getParam("amplitude_max");
    float pre_amplitude_min = fitConfig->getParam("pre_amplitude_min");
    float pre_amplitude_max = fitConfig->getParam("pre_amplitude_max");
    float post_amplitude_min = fitConfig->getParam("post_amplitude_min");
    float post_amplitude_max = fitConfig->getParam("post_amplitude_max");


    std::filesystem::path p(file_name);
    std::filesystem::path file_output = p.parent_path() / (p.stem().string() + "_" + this_ch +"_n.root");

    std::cout << "Creating new file: " << file_output << std::endl;
    // Cria novo arquivo
    TFile* newfile = TFile::Open(file_output.c_str(), "RECREATE");
    //cria arvore de saida
    TTree *tree_write = new TTree("T1","data");
    //aponta corretamente o my_data 
    tree_write->Branch("Data", "my_data", &data); 


    // ---------- Recarrega separações ----------
    if (my_channel != -1)
    {
        separation.clear();
        std::string file_sep = Form("/home/gabriel/Documents/protodune/protodune_vd/light_analysis/data/FIT_SPE/%d.txt", my_channel);
        double s;
        std::ifstream data_sp(file_sep);
        while (data_sp >> s)
        {
            std::cout << "Separation: " << s << std::endl;
            separation.push_back(s);
        }
        data_sp.close();
    }


    bool found=false;
    // ---------- Loop principal ----------
    for (int i = 0; i < tree->GetEntries(); i++)
    {
        tree->GetEntry(i);
        if (data->Channel != my_channel) continue;

        if (data->baseline > baseline_min && data->baseline < baseline_max &&
            data->noise < noise_max &&
            data->amplitude < amplitude_max && data->amplitudemin > amplitude_min &&
            data->preamplitude < pre_amplitude_max && data->preamplitudemin > pre_amplitude_min &&
            data->postamplitude < post_amplitude_max && data->postamplitudemin > post_amplitude_min)
        {
            double integral = data->integral;
            found = false;
            if (!(integral < separation.front() || integral > separation.back()))
            {
                int norm = 0;
                for (norm = 0; norm < (int)separation.size() - 1; norm=norm+2)
                {
                    if (integral > separation[norm] && integral < separation[norm + 1])
                    {
                        found = true;
                        break;
                    }
                        
                }
                if(found)
                {
                    norm = 1 + norm / 2;
                    norm_ch++;

                    for (int k = 0; k < (int)data->adcs.size(); k++)
                    {
                        all_signals[k].push_back((data->adcs[k] - data->baseline) / norm);

                    }
                    data->n_spe=norm;
                    tree_write->Fill();
                }
                
            }
        }
    }

    file->Close();
    tree_write->Write("", TObject::kOverwrite);
    newfile->Write();
    newfile->Close();
    delete data;

    // ---------- Calcula moda e média ao redor ----------
    std::vector<double> signal(1024, 0.0);
    for (int k = 0; k < 1024; k++)
    {
        const auto& vals = all_signals[k];
        if (vals.empty()) continue;


        double moda = calc_moda(vals, 0.1);
 
        double var = 0.0;
        for (double v : vals)
            var += (v - moda) * (v - moda);
        double sigma = std::sqrt(var / vals.size());

        sigma=0.3*sigma;
      
        std::vector<double> window;
        window.reserve(vals.size());
        for (double v : vals)
            if (std::fabs(v - moda) <= sigma)
                window.push_back(v);

        if (!window.empty())
            signal[k] = std::accumulate(window.begin(), window.end(), 0.0) / window.size();
        else
            signal[k] = moda; // fallback
    }

    // ---------- Gera eixo X ----------
    std::vector<double> x(1024);
    for (int i = 0; i < 1024; i++) x[i] = i;

    // ---------- Desenha e salva ----------
    TCanvas* c = new TCanvas("c", "Waveform (Moda+Média Local)", 800, 600);
    TGraph* g = new TGraph(1024, x.data(), signal.data());
    g->SetTitle(Form("Waveform (Moda + média local) - Channel %d", my_channel));
    g->GetXaxis()->SetTitle("Sample");
    g->GetYaxis()->SetTitle("Amplitude (ADC)");
    g->Draw("AL");

    std::string out_base_dir;
    if (run_number == "null")
        out_base_dir = "../waveforms/" + fs::path(file_name).filename().string();
    else
        out_base_dir = "../waveforms/" + run_number;

    fs::create_directories(out_base_dir);

    std::string fig_name = out_base_dir + Form("/waveform_modamedia_%d.png", my_channel);
    c->SaveAs(fig_name.c_str());

    std::ofstream fout(out_base_dir + Form("/waveform_modamedia_%d.txt", my_channel));
    for (int k = 0; k < 1024; k++)
        fout << x[k] << "\t" << signal[k] << "\n";
    fout.close();

    delete g;
    delete c;

    std::cout << "✅ Forma de onda (moda + média local, σ centrado na moda) salva em: " << fig_name << std::endl;

    return 0;
}
