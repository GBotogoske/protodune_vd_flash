#include "readConfigFit.hh"
#include <TString.h> 

#include <fstream>

readConfigFit::readConfigFit(int ch)
{
    std::string file_name;
    if(ch!=-1)
    {
        file_name = std::string("/home/gabriel/Documents/protodune/protodune_vd/light_analysis/configuration/fit/") + Form("%d.json",ch);
    }
    else
    {
        file_name = std::string("/home/gabriel/Documents/protodune/protodune_vd/light_analysis/configuration/fit/") + Form("default.json");
    }
    std::ifstream file(file_name);
    json j;
    file >> j;

    for (auto &el : j.items())
    {
        my_params[el.key()] = el.value().get<float>();
    }

    file.close();

}

readConfigFit::~readConfigFit()
{
}
