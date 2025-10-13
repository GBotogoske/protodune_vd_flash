#include "utils.hh"

#include <TGraph.h>
#include <TLine.h>
#include "TAxis.h"

#include <TLegend.h>
#include <TApplication.h>
#include <TSystem.h>

#include <iostream>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;

void print_waveform(TCanvas* c, my_data* data, int i)
{
    std::vector<short> wf = data->adcs;
    std::vector<double> x(wf.size()); 
    std::vector<double> y(wf.size()); 
    for (size_t j = 0; j < wf.size(); j++) 
    { 
        x[j] = static_cast<double>(j); 
        y[j] = static_cast<double>(wf[j]); 
    }
    
    auto g = new TGraph(x.size(), x.data(), y.data()); 
    g->SetTitle("Waveform;Time [ticks];ADC"); 
    g->SetLineColor(kBlue); 
    g->SetLineWidth(2);

    c->cd();
    c->Clear();
    g->Draw("AL");

    auto line_start = new TLine(data->t0, g->GetYaxis()->GetXmin(), data->t0, g->GetYaxis()->GetXmax());
    auto line_end = new TLine(data->tend, g->GetYaxis()->GetXmin(), data->tend, g->GetYaxis()->GetXmax());
    line_start->SetLineColor(kRed);
    line_start->SetLineStyle(2);
    line_start->SetLineWidth(2);
    line_start->Draw();
    line_end->SetLineColor(kGreen);
    line_end->SetLineStyle(2);
    line_end->SetLineWidth(2);
    line_end->Draw();

    auto legend = new TLegend(0.7, 0.7, 0.95, 0.95);
    legend->SetTextSize(0.03);   
    legend->SetTextFont(42);
    legend->Clear();
    legend->AddEntry((TObject*)0, Form("Event %d", i), "");
    legend->AddEntry((TObject*)0, Form("Channel %d", data->Channel), "");
    legend->AddEntry((TObject*)0, Form("Baseline %.2f", data->baseline), "");
    legend->AddEntry((TObject*)0, Form("Noise %.2f", data->noise), "");
    legend->AddEntry((TObject*)0, Form("t0 %.2f", data->t0), "");
    legend->AddEntry((TObject*)0, Form("tend %.2f", data->tend), "");
    legend->AddEntry((TObject*)0, Form("Amplitude %.2f", data->amplitude), "");
    legend->AddEntry((TObject*)0, Form("Integral %.2f", data->integral), "");
    legend->AddEntry((TObject*)0, Form("Prompt %.2f", data->prompt), "");
    legend->Draw();

    gSystem->ProcessEvents();
    c->Update();
    c->Draw();
    
    delete legend;
    delete line_start;
    delete line_end;
    delete g;

}

std::string search_file(std::string path, std::string run,std::string find)
{
    std::string file_name;
    if (!run.empty())
    {
       
        bool found = false;
        // Percorre todos os arquivos na pasta procurando o run desejado
        for (const auto& entry : fs::directory_iterator(path))
        {
            std::string fname = entry.path().filename().string();
            if (fname.find(run) != std::string::npos && fname.find(find) != std::string::npos)
            {
                file_name = entry.path().string();
                found = true;
                break;
            }
        }

    }
    return file_name;
}

bool intersect_two(const Peak& a, const Peak& b, double& xstar)
{
    // garante mu_a < mu_b
    Peak left=a, right=b;
    if (left.mu > right.mu) std::swap(left,right);

    const double sa2 = left.s*left.s;
    const double sb2 = right.s*right.s;

    const double R = (right.A/right.s) / (left.A/left.s); // (A_j/σ_j)/(A_i/σ_i)
    const double lnR = std::log(R);

    const double eps = 1e-12;
    const double Acoef = (1.0/sb2) - (1.0/sa2);
    const double Bcoef = (-2.0*right.mu/sb2) + (2.0*left.mu/sa2);
    const double Ccoef = (right.mu*right.mu/sb2) - (left.mu*left.mu/sa2) - 2.0*lnR;

    if (std::fabs(Acoef) < eps) {
        // σ_i ≈ σ_j → solução linear
        // x*(μ_i-μ_j) + (μ_j^2-μ_i^2)/2 - σ^2 lnR = 0, com σ^2 ≈ sa2 ≈ sb2
        double s2 = 0.5*(sa2+sb2);
        double denom = (left.mu - right.mu);
        if (std::fabs(denom) < eps) return false; // meios praticamente iguais
        xstar = ( (left.mu*left.mu - right.mu*right.mu) + 2.0*s2*lnR ) / (2.0*denom);
    } else {
        // solução quadrática
        double disc = Bcoef*Bcoef - 4.0*Acoef*Ccoef;
        if (disc < 0) return false; // sem raiz real
        double sqrtD = std::sqrt(std::max(0.0,disc));
        double x1 = (-Bcoef + sqrtD)/(2.0*Acoef);
        double x2 = (-Bcoef - sqrtD)/(2.0*Acoef);
        // escolha a raiz entre os centros
        double lo = left.mu, hi = right.mu;
        bool in1 = (x1 >= lo && x1 <= hi);
        bool in2 = (x2 >= lo && x2 <= hi);
        if (in1 && !in2) xstar = x1;
        else if (!in1 && in2) xstar = x2;
        else if (in1 && in2) {
            // ambas dentro: pegue a mais próxima do meio
            double mid = 0.5*(lo+hi);
            xstar = (std::fabs(x1-mid) < std::fabs(x2-mid)) ? x1 : x2;
        } else {
            // nenhuma dentro: pegue a mais próxima do intervalo (opcional) ou sinalize falha
            double mid = 0.5*(lo+hi);
            xstar = (std::fabs(x1-mid) < std::fabs(x2-mid)) ? x1 : x2;
        }
    }
    return true;
}

std::map<std::string, std::vector<int>> get_map_spe()
{
    std::map<std::string, std::vector<int>> this_map;
    this_map["39357"] = std::vector<int>{1050,1051};
    this_map["39358"] = std::vector<int>{1060};
    this_map["39359"] = std::vector<int>{1061};
    this_map["39360"] = std::vector<int>{1070,1071,2030,2031,2040,2041};
    this_map["39361"] = std::vector<int>{1080,1081};
    this_map["39365"] = std::vector<int>{2070,2071};
    this_map["39366"] = std::vector<int>{2080,2081};

    return this_map;

}

std::map<std::string,std::string> get_map_spe2()
{
    std::map<std::string,std::string> this_map;
    this_map["1050"]="39357";
    this_map["1051"]="39357";
    this_map["1060"]="39358";
    this_map["1061"]="39359";
    this_map["1070"]="39360";
    this_map["1071"]="39360";
    this_map["1080"]="39361";
    this_map["1081"]="39361";
    this_map["2030"]="39360";
    this_map["2031"]="39360";
    this_map["2040"]="39360";
    this_map["2041"]="39360";
    this_map["2070"]="39365";
    this_map["2071"]="39365";
    this_map["2080"]="39366";
    this_map["2071"]="39366";
    return this_map;
}
