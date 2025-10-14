#ifndef READCONFIGFIT_HH
#define READCONFIGFIT_HH

#include <nlohmann/json.hpp>
using json = nlohmann::json;
#include <iostream>

class readConfigFit
{
    public:
        readConfigFit(int ch);
        ~readConfigFit();

        const std::map<std::string,float> get_params(){return this->my_params;};
        float getParam(const std::string &key) const
        {
            auto it = my_params.find(key);
            if (it != my_params.end())
                return it->second;
            else
                return 0.0;
        }

        void print() const
        {
            for (const auto &p : my_params)
            {
                std::cout << p.first << " = " << p.second << std::endl;
            }
        }
        
    private:
        std::map<std::string,float> my_params;
};



#endif