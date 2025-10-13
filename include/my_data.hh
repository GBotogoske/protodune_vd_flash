#ifndef MY_DATA_HH
#define MY_DATA_HH

#include <vector>
#include "TObject.h"

class my_data: public TObject
{
    public:
        Int_t Channel;
        long Timestamp;
        std::vector<short> adcs; 
        int flash;
        float integral = 0;
        float amplitude = 0;
        float baseline = 0;
        float noise = 0;
        int n_peaks = 0;
        float t0 = 0;
        float tend = 0;
        float prompt=0;
        

        void set_parameters(Int_t ch, long time, const std::vector<short>& signal);
        my_data();
        my_data(Int_t ch, long time,  std::vector<short> signal);

        void print_baseline();
        void print_noise();
        void print_t0();
        void print_tend();
        void print_integral();
        void print_amplitude();
        void print_all();

        void calc_baseline(const int index);
        void calc_baseline(const int index1,const int index2,const int index3,const int index4);
        void calc_noise(const int index);
        void calc_t0(const int index,const int index_lim);
        void calc_tend(const int index,const int index_lim);
        
        void calc_integral();
        void calc_integral(int start,int end);
        void calc_amplitude();
        void calc_amplitude(int start,int end);
        void calc_prompt(const int index);

        ClassDef(my_data, 1);
};


#endif
