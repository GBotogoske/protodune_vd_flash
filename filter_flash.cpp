#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <string>
#include <TSystem.h>
#include <TStyle.h>
#include <TColor.h>
#include <TH2F.h>
#include <TBox.h>
#include <utils.hh>
#include <my_data.hh>
#include <filesystem>

#define print_flash true

namespace fs = std::filesystem;

//estrutura para os flashs
typedef struct 
{
    std::vector<my_data> data;
    int n_ch=0;
}my_Flash;

//lista de canais para procurar os flashs ao mesmo tempo
std::vector<int> list_channel = {1050,1051,1060,1061,1070,1071,1080,1081,2030,2031,2040,2041,2080,2081};

//lista de spe_charge
std::vector<double> list_spe_q(list_channel.size(),1.0);

bool print_charge(std::vector<my_data> this_data,bool print = true)
{
    std::map<int,int> map_ch;
    for(int i=0;i<list_channel.size();i++)
    {
        map_ch[list_channel[i]]=i;
    }

    int n_ch = this_data.size();
    std::vector<double> values_q(list_channel.size(),0.0);

    for(int i=0;i<n_ch;i++)
    {
        int ch = this_data[i].Channel;
        auto it = map_ch.find(ch);
        if (it != map_ch.end())
        {
            int idx = it->second;
            values_q[idx]=this_data[i].integral!=-1000?this_data[i].integral:0.0;
        } 
    }

    std::vector<double> values_q_total(list_channel.size()/2,0.0);
    double first_max = -2e6;
    double second_max = -2e6;
    for(int i=0;i<values_q_total.size();i++)
    {
        values_q_total[i]=values_q[2*i]/list_spe_q[2*i]+values_q[2*i+1]/list_spe_q[2*i+1];
        if(i>=0 && i<4)
        {
            if(values_q_total[i]>first_max)
            {
                first_max=values_q_total[i];
            }
        }
        else
        {
            if(values_q_total[i]>second_max)
            {
                second_max=values_q_total[i];
            }
        }
    }
    std::vector<double> values_q_cathode(0);
    std::vector<double> values_q_membrane(0);
    for(int i=0;i<values_q_total.size();i++)
    {
        if(i>=0 && i<4)
        {
            values_q_total[i]/=first_max;
            values_q_cathode.push_back(values_q_total[i]);
        }
        else
        {
           values_q_total[i]/=second_max;
           values_q_membrane.push_back(values_q_total[i]);
        }
    }

    int n1=0;
    for(int i=0;i<values_q_membrane.size();i++)
    {
        if(values_q_membrane[i]>0.7)
        {
            n1++;
        }
    }

    if(print)
    {
       static TCanvas* c = nullptr;

        if (!c)
        {
            c = new TCanvas("event_display", "Event Display", 700, 800);
            c->Divide(1,2);   // apenas uma vez
        }

        gStyle->SetOptStat(0);

        // limpar
        c->cd(1); gPad->Clear();
        c->cd(2); gPad->Clear();

        // ----- CATHODE -----
        c->cd(1);

        TH2F* hcat = new TH2F(Form("hcat_%d", rand()), "Cathode;X;Y",
                            4, 0, 4,
                            4, 0, 4);

        // preto = -1
        for(int ix=1; ix<=4; ix++)
            for(int iy=1; iy<=4; iy++)
                hcat->SetBinContent(ix, iy, 0);

        // posições
        hcat->SetBinContent(2, 4, values_q_total[3]);
        hcat->SetBinContent(3, 3, values_q_total[2]);
        hcat->SetBinContent(1, 2, values_q_total[1]);
        hcat->SetBinContent(3, 1, values_q_total[0]);

        // paleta
        static bool palette_init = false;
        if (!palette_init)
        {
            // 3 cores:
            // stop=0.0 → PRETO (valor -1)
            // stop=0.5 → BRANCO (valor 0)
            // stop=1.0 → VERDE (valor 1)

            const int N = 3;
            double stops[N]   = {0.0, 0.5, 1.0};
            double red[N]     = {0.0, 1.0, 0.0};
            double green[N]   = {0.0, 1.0, 1.0};
            double blue[N]    = {0.0, 1.0, 0.0};

            // ROOT 6.32 cria e aplica automaticamente a paleta
            TColor::CreateGradientColorTable(N, stops, red, green, blue, 255);

            palette_init = true;
        }

        hcat->SetMaximum(1.0);
        hcat->SetMinimum(-1.0);
        hcat->SetMarkerSize(3.0);
        hcat->Draw("COLZ TEXT");

                // Desenhar contornos vermelhos
        int nx = hcat->GetNbinsX();
        int ny = hcat->GetNbinsY();

        for (int ix = 1; ix <= nx; ix++)
        {
            for (int iy = 1; iy <= ny; iy++)
            {
                double x1 = hcat->GetXaxis()->GetBinLowEdge(ix);
                double x2 = hcat->GetXaxis()->GetBinUpEdge(ix);
                double y1 = hcat->GetYaxis()->GetBinLowEdge(iy);
                double y2 = hcat->GetYaxis()->GetBinUpEdge(iy);

                TBox* box = new TBox(x1, y1, x2, y2);
                box->SetFillStyle(0);        // sem preenchimento
                box->SetLineColor(kRed);     // borda vermelha
                box->SetLineWidth(2);        // mais grosso para visualização
                box->Draw("same");
            }
        }

        // ----- MEMBRANE -----
        c->cd(2);

        TH2F* hmem = new TH2F(Form("hmem_%d", rand()), "Membrane;X;Y",
                            1, 0, 1,
                            4, 0, 4);

        // preto
        for(int i=1; i<=4; i++)
            hmem->SetBinContent(1, i, -1);

        // módulos
        hmem->SetBinContent(1, 4, values_q_total[4]);
        hmem->SetBinContent(1, 3, values_q_total[5]);
        hmem->SetBinContent(1, 1, values_q_total[6]);

        hmem->SetMaximum(1.0);
        hmem->SetMinimum(-1.0);

        hmem->SetMarkerSize(3.0);
        hmem->Draw("COLZ TEXT");

        int nx2 = hmem->GetNbinsX();
        int ny2 = hmem->GetNbinsY();

        for (int ix = 1; ix <= nx2; ix++)
        {
            for (int iy = 1; iy <= ny2; iy++)
            {
                double x1 = hmem->GetXaxis()->GetBinLowEdge(ix);
                double x2 = hmem->GetXaxis()->GetBinUpEdge(ix);
                double y1 = hmem->GetYaxis()->GetBinLowEdge(iy);
                double y2 = hmem->GetYaxis()->GetBinUpEdge(iy);

                TBox* box2 = new TBox(x1, y1, x2, y2);
                box2->SetFillStyle(0);
                box2->SetLineColor(kRed);
                box2->SetLineWidth(2);
                box2->Draw("same");
            }
        }

        c->Update();
        gSystem->ProcessEvents();

        std::cin.get();


        /* std::cout << "#---------------------#Cathode#----------------------------#" << std::endl;
        std::cout << "XXXX - " << values_q_total[3] << " - XXXX - XXXX " << std::endl;
        std::cout << "XXXX - XXXX  " << values_q_total[2] << " - XXXX " << std::endl;
        std::cout << values_q_total[1] << " - XXXX - XXXX - XXXX " << std::endl;
        std::cout << "XXXX - XXXX  " << values_q_total[0] << " - XXXX " << std::endl;
        std::cout << "#-------------------------------------------------#" << std::endl;
        std::cout << "#---------------------#Membrane#----------------------------#" << std::endl;
        std::cout << values_q_total[4] << std::endl;
        std::cout << values_q_total[5] << std::endl;
        std::cout << "XXXX" << std::endl;
        std::cout << values_q_total[6] << std::endl;
        std::cout << "#-------------------------------------------------#" << std::endl; */

    }
    
    return n1==3?true:false;
}



int main(int argc, char** argv)
{
    TApplication app("app", &argc, argv);
    std::string file_name;
    //base diretorio
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
            file_name=search_file(base_dir,run_number,std::string("_flash.root"));
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

    //cria um mapa de flashs
    std::map<int,my_Flash> my_map;

    //prenche os flashs por flashID
    {
        TFile* file = TFile::Open(file_name.c_str(), "READ");
        TTree* tree = (TTree*)file->Get("T1");
        my_data* data = nullptr;
        tree->SetBranchAddress("Data", &data);

        for(int i = 0; i < tree->GetEntries(); ++i) 
        {
            tree->GetEntry(i);
            my_map[data->flash].data.push_back(*data);
            my_map[data->flash].n_ch++;
        }

        file->Close();
    }

    // --- Cria arquivo de saida  ---
    std::string file_output;
    if(run_number=="null")
    {
        file_output = fs::path(file_name).replace_extension("").string() + "_filtered.root";
    }
    else
    {   
        file_output=fs::path(file_name).parent_path() / ("data_analysed_" + run_number + "_flash_filtered.root");
    }

    //abre novo arquivo
    TFile* newfile = TFile::Open(file_output.c_str(), "RECREATE");
    std::cout << file_output << std::endl;

    //abre nova TTree
    auto* tree_write = new TTree("T1", "data");
    auto* data_2 = new my_data();
    tree_write->Branch("Data", "my_data", &data_2);

    int n_ch = list_channel.size();

    //carrega arquivo com cargas spe
    std::string file_q_spe="/home/gabriel/Documents/protodune/protodune_vd/light_analysis/charge_spe.txt";
    std::ifstream file_q_spe_ifstream(file_q_spe);
    double x;
    int i_spe=0;
    int j_spe=0;
    if (!file_q_spe_ifstream.is_open()) 
    {
        std::cout << "Erro: não consegui abrir " << file_q_spe << std::endl;
    }
    else
    {
        std::cout << "Abrindo: " << file_q_spe << std::endl;
        while (file_q_spe_ifstream >> x) 
        {
            if(j_spe!=12 && j_spe!=13)
            {
                list_spe_q[i_spe]=x;
                i_spe++;
            }
            j_spe++;
        }
        file_q_spe_ifstream.close();
    }



    //varre os flashs
    for (auto& [flash_id, flash] : my_map)
    {
        bool all_found = true;
        //garente que o flash ao menos tenha o mesmo numero de hits que o esperado
        if (flash.n_ch >= n_ch)
        {
            //varre a lista de canais
            for (int this_ch : list_channel)
            {
                bool found_this = false;
                for (int j = 0; j < flash.data.size(); j++)
                {
                    if (flash.data[j].Channel == this_ch)
                    {
                        //se encontrou o canal, pode ir para o proximo canal da lista
                        found_this = true;
                        break;
                    }
                }
                //se nao encontrou nenhum, pode sair da varrida de canais
                if (!found_this)
                {
                    all_found = false;
                    break; 
                }
            }
        }
        else
        {
            all_found = false;
        }
        // se encontrou todos os canais
        if (all_found)
        {
            if(print_charge(flash.data,print_flash))
            {
                for (int j = 0; j < flash.data.size(); j++)
                {
                    *data_2=flash.data[j];
                    tree_write->Fill();
                }
            }   
            
        }
    }

    tree_write->Write("", TObject::kOverwrite);
    newfile->Close();

    delete data_2;
    app.Run();
    return 0;
}
