#include "my_data.hh"
#include <iostream>
#include <cmath>

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
