#include "my_data.hh"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <unordered_map>

ClassImp(my_data)

void my_data::set_parameters(Int_t ch, long time, const std::vector<short>& signal)
{
    this->adcs = signal;
    this->Timestamp = time;
    this->Channel = ch;
}

my_data::my_data()
{
}

my_data::my_data(Int_t ch, long time, std::vector<short> signal)
{
    this->set_parameters(ch,time,signal);
}

void my_data::print_baseline()
{
    std::cout << "Baseline: " << this->baseline <<std::endl;
}

void my_data::print_noise()
{
    std::cout << "Noise: " << this->noise <<std::endl;
}

void my_data::print_t0()
{
    std::cout << "t0: " << this->t0 <<std::endl;
}

void my_data::print_tend()
{
    std::cout << "tend: " << this->tend <<std::endl;
}

void my_data::print_integral()
{
    std::cout << "Integral: " << this->integral <<std::endl;
}

void my_data::print_amplitude()
{
    std::cout << "Amplitude: " << this->amplitude <<std::endl;
}

void my_data::print_all()
{
    this->print_baseline();
    this->print_noise();
    this->print_t0();
    this->print_tend();
    //this->print_zero_crossing();
    this->print_integral();
    this->print_amplitude();
}

void my_data::calc_baseline(const int index)
{
    this->baseline=0;
    for(int i = 0; i<index; i++)
    {
        this->baseline+=this->adcs[i];
    }
    this->baseline/=index;
}

void my_data::calc_baseline(const int index1, const int index2, const int index3, const int index4)
{
    this->baseline=0;
    int cont=0;
    for(int i = index1; i<index2; i++)
    {
        this->baseline+=this->adcs[i];
        cont++;
    }
    for(int i = index3; i<index4; i++)
    {
        this->baseline+=this->adcs[i];
        cont++;
    }
    this->baseline/=cont;
}

void my_data::calc_baseline(const double sigma_cut, const double skip_window)
{
    if (this->adcs.empty())
    {
        this->baseline = 0;
        return;
    }

    std::unordered_map<int, int> freq;
    for (auto val : this->adcs) 
    {
        freq[val]++;
    }
    auto max_it = std::max_element(freq.begin(), freq.end(),
                                   [](auto& a, auto& b){ return a.second < b.second; });
    double mode = static_cast<double>(max_it->first);

    std::vector<double> valid_points;
    double low = mode - sigma_cut;
    double high = mode + sigma_cut;

    for (size_t i = 0; i < this->adcs.size(); ++i)
    {
        double val = this->adcs[i];
        if (val >= low && val <= high) 
        {
            valid_points.push_back(val);
        } 
        else 
        {
            // se encontrou pico, pula N bins
            i += static_cast<size_t>(skip_window);
        }
    }

    // 4️⃣ Média dos pontos válidos
    if (!valid_points.empty())
    {
        this->baseline = std::accumulate(valid_points.begin(), valid_points.end(), 0.0) / valid_points.size();
    } else 
    {
        this->baseline = mode;
    }
}

void my_data::calc_noise(const int index)
{
    this->noise = 0;
    for(int i = 0; i<index; i++)
    {
        double delta = this->adcs[i] - this->baseline;
        this->noise += delta * delta;
    }
    this->noise = std::sqrt(this->noise/index);
}

void my_data::calc_t0(const int index,const int index_lim)
{
    float threshold = this->baseline + 5 * this->noise;

    for (int i = index + 1; i < index_lim ; ++i) 
    {
        if (this->adcs[i] > threshold && this->adcs[i - 1] <= threshold) 
        {
            
            float x1 = (i - 1);
            float y1 = this->adcs[i - 1];
            float x2 = i;
            float y2 = this->adcs[i];
            float m = (float)((y2 - y1) / (x2 - x1));
            float t0 = x1 + (float)((threshold - y1) / m);
            this->t0 = t0;
            return;
        }
    }

    this->t0 = -1000; //error
    return; 
}

void my_data::calc_tend(const int index, const int index_lim)
{
    float threshold = this->baseline + 5 * this->noise;
    if(index>0)
    {
        for (int i = index + 1; i < index_lim ; ++i) 
        {
            if (this->adcs[i] < threshold && this->adcs[i - 1] >= threshold) 
            {
                
                float x1 = (i - 1);
                float y1 = this->adcs[i - 1];
                float x2 = i;
                float y2 = this->adcs[i];
                float m = (float)((y2 - y1) / (x2 - x1));
                float tend = x1 + (float)((threshold - y1) / m);
                this->tend = tend;
                return;
            }
        }
    }
    this->tend = -1000; //error
    return; 
}

void my_data::calc_integral()
{
    this->integral = 0;
    if(this->tend != -1000 && this->t0 != -1000)
    {
        if(this->tend >= this->t0)
        {
            int start = (int)this->t0;
            int end = (int)this->tend;

            for(int i = start; i<=end ; i++)
            {
                this->integral+=(this->adcs[i]-this->baseline);
            }
            return;

        }
    }
    this->integral=-1000;
    return;
}

void my_data::calc_integral(int start, int end)
{
     this->integral = 0;
    if(start != -1000 && end != -1000)
    {
        if(end >= start)
        {
            for(int i = start; i<=end ; i++)
            {
                this->integral+=(this->adcs[i]-this->baseline);
            }
            return;

        }
    }
    this->integral=-1000;
    return;
}

void my_data::calc_amplitude()
{
    this->amplitude = -1000;
    if(this->tend != -1000 && this->t0 != -1000)
    {
        if(this->tend >= this->t0)
        {
            int start = (int)this->t0;
            int end = (int)this->tend;

            for(int i = start; i<=end ; i++)
            {
                double value = this->adcs[i]-this->baseline;
                if(this->amplitude < value)
                {
                    this->amplitude = value;
                }
            }
        }
    }
    return;
}

void my_data::calc_preamplitude(int start, int end)
{
    this->preamplitude = -1000;
    if(start != -1000 && end != -1000)
    {
        if(end >= start)
        {
            for(int i = start; i<=end ; i++)
            {
                double value = this->adcs[i]-this->baseline;
                if(this->preamplitude < value)
                {
                    this->preamplitude = value;
                }
            }
        }
        
    }
    return;
}

void my_data::calc_preminamplitude(int start, int end)
{
    this->preamplitudemin = 1000;
    if(start != -1000 && end != -1000)
    {
        if(end >= start)
        {

            for(int i = start; i<=end ; i++)
            {
                double value = this->adcs[i]-this->baseline;
                if(this->preamplitudemin > value)
                {
                    this->preamplitudemin = value;
                }
            }
        }
    }
    return;
}

void my_data::calc_postamplitude(int start, int end)
{
    this->postamplitude = -1000;
    if(start != -1000 && end != -1000)
    {
        if(end >= start)
        {
            for(int i = start; i<=end ; i++)
            {
                double value = this->adcs[i]-this->baseline;
                if(this->postamplitude < value)
                {
                    this->postamplitude = value;
                }
            }
        }
        
    }
    return;
}

void my_data::calc_postminamplitude(int start, int end)
{
    this->postamplitudemin = 1000;
    if(start != -1000 && end != -1000)
    {
        if(end >= start)
        {

            for(int i = start; i<=end ; i++)
            {
                double value = this->adcs[i]-this->baseline;
                if(this->postamplitudemin > value)
                {
                    this->postamplitudemin = value;
                }
            }
        }
    }
    return;
}

void my_data::calc_amplitude(int start, int end)
{
    this->amplitude = -1000;
    if(start != -1000 && end != -1000)
    {
        if(end >= start)
        {

            for(int i = start; i<=end ; i++)
            {
                double value = this->adcs[i]-this->baseline;
                if(this->amplitude < value)
                {
                    this->amplitude = value;
                }
            }
        }
    }
    return;
}

void my_data::calc_minamplitude(int start, int end)
{
    this->amplitudemin = 1000;
    if(start != -1000 && end != -1000)
    {
        if(end >= start)
        {

            for(int i = start; i<=end ; i++)
            {
                double value = this->adcs[i]-this->baseline;
                if(this->amplitudemin > value)
                {
                    this->amplitudemin = value;
                }
            }
        }
    }
    return;
}

void my_data::calc_prompt(const int index)
{
    this->prompt = -1000;
    if(this->integral!=0 && this->integral!=-1000)
    {
        if(this->tend >= this->t0)
        {
            double diff= this->tend - this->t0;
            if(index <= diff)
            {
                double integral_fast=0.0;
                int start = (int)this->t0;
                int end = start + index;

                for(int i = start; i<=end ; i++)
                {
                    integral_fast+=(this->adcs[i]-this->baseline);
                }
                this->prompt=integral_fast/this->integral;
                return;
            }
        }     
    }
     this->prompt = -1000;
}
