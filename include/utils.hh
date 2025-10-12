#ifndef UTILS_HH
#define UTILS_HH

#include "my_data.hh"
#include <TCanvas.h>

void print_waveform(TCanvas* c, my_data* data, int i);

typedef struct  
{ double A, mu, s; }Peak ;
bool intersect_two(const Peak& a, const Peak& b, double& xstar);


#endif