#include "utils.hh"

#include <TGraph.h>
#include <TLine.h>
#include "TAxis.h"

#include <TLegend.h>
#include <TApplication.h>
#include <TSystem.h>

#include <iostream>
#include <string>

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