#ifndef UTILS_HH
#define UTILS_HH

#include "my_data.hh"
#include <TCanvas.h>


void print_waveform(TCanvas* c, my_data* data, int i);
std::string search_file(std::string path, std::string run,std::string find);

typedef struct  
{ double A, mu, s; }Peak ;
bool intersect_two(const Peak& a, const Peak& b, double& xstar);

std::map<std::string,std::vector<int>> get_map_spe();
std::map<std::string,std::string> get_map_spe2();


#endif