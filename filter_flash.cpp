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


typedef struct 
{
    std::vector<my_data> data;
    int n_ch=0;
}my_Flash;

std::vector<int> list_channel = {1050,1051,1060,1061,1070,1071,1080,1081,2030,2031,2040,2041,2080,2081};

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
        file_name = std::string("data_analysed_flash_39510.root");
    }

    std::map<int,my_Flash> my_map;

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

    TFile* newfile = TFile::Open("data_analysed_flash_filtered.root", "RECREATE");
    auto* tree_write = new TTree("T1", "data");
    auto* data_2 = new my_data();
    tree_write->Branch("Data", "my_data", &data_2);

    int n_ch = list_channel.size();
    for (auto& [flash_id, flash] : my_map)
    {
        bool all_found = true;
        if (flash.n_ch >= n_ch)
        {
            for (int this_ch : list_channel)
            {
                bool found_this = false;
                for (int j = 0; j < flash.data.size(); j++)
                {
                    if (flash.data[j].Channel == this_ch)
                    {
                        found_this = true;
                        break;
                    }
                }
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

        if (all_found)
        {
            //std::cout << "um flash\n ";
            for (int j = 0; j < flash.data.size(); j++)
            {
                *data_2=flash.data[j];
                //std::cout<<flash.data[j].Timestamp << std::endl;
                tree_write->Fill();
            }
        }
    }

    tree_write->Write("", TObject::kOverwrite);
    newfile->Close();

    delete data_2;
    return 0;
}
