#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TCanvas.h>
#include <TApplication.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLegend.h>

#include <my_data.hh>

#define T 0

int main(int argc, char** argv)
{
    int Channel = 2030;

    std::ifstream infile("../list_file.txt");
    std::string line;
    std::vector<std::string> list_file;
    while (std::getline(infile, line))
    {
        list_file.push_back(line);
    }
   
    TApplication app("app", &argc, argv);
    std::vector<double> EF = {0,10.0/360,80.0/360,155.7/360};
    // --- histogramas ---
    TCanvas* c = new TCanvas("c", "Histogramas", 800, 600);
    auto legend = new TLegend(0.65, 0.7, 0.9, 0.9);
    std::vector<TH1D*> h_list;

    for(int k=0; k<list_file.size(); k++)
    {
        long min_t,max_t;
        h_list.push_back(new TH1D(Form("Integral %.2f kV/cm", EF[k]), Form("Integral %.2f kV/cm", EF[k]), 100, 0, 2e6));
        auto h  = h_list.back();
        h->SetLineColor(k+1);
        h->SetLineWidth(2);
        {
            TFile* file1 = TFile::Open(list_file[k].c_str(), "READ");
            TTree* tree1 = (TTree*)file1->Get("T1");
            my_data* data = nullptr;
            tree1->SetBranchAddress("Data", &data);

            for(int i = 0; i < tree1->GetEntries(); ++i) 
            {
                tree1->GetEntry(i);
                if(i==0)
                {
                    max_t=data->Timestamp;
                    min_t=data->Timestamp;
                }
                else
                {
                    if(data->Timestamp>max_t)
                    {   
                        max_t=data->Timestamp;
                    }
                    if(data->Timestamp<min_t)
                    {   
                        min_t=data->Timestamp;
                    }
                }
                if (data && data->Channel == Channel && data->integral != -1000 && !data->adcs.empty())
                {
                    if(data->integral >=T)
                    {
                        h->Fill(data->integral);
                    }
                }
                    
            }
            file1->Close();
        }
        long deltaT=(max_t-min_t)*16e-9;
        h->Scale(1.0/deltaT);

        std::cout << "Total Rate for " << EF[k] << " kV/cm " << h->Integral() <<  " events/s" << std::endl;

        if(k==0)
        {
            h->Draw();
        }
        else
        {
            h->Draw("SAME");
        }
        legend->AddEntry(h, h->GetName(), "l");   
    }
    legend->Draw();

    /* auto h1 = h_list[0];
    auto h2 = h_list.back();

    auto func = [h1](double* x, double* par)
    {
        double A = par[0];
        double B = par[1];
        double xprime = B * (x[0]);
        if (xprime < h1->GetXaxis()->GetXmin() || xprime > h1->GetXaxis()->GetXmax() || xprime < T) 
        {
                return 0.0;
        }
        return A * h1->Interpolate(xprime);
    };

    TF1* fitFunc = new TF1("fitFunc", func, h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax(), 2);
    fitFunc->SetParNames("A", "B");
    fitFunc->SetParameters(1.0, 0.6);
    h2->Fit(fitFunc, "R");  

    fitFunc->SetLineColor(kRed);
    fitFunc->Draw("SAME");

    std::cout << "\n--- Fit Results ---" << std::endl;
    std::cout << "A = " << fitFunc->GetParameter(0) << std::endl;
    std::cout << "B = " << fitFunc->GetParameter(1) << std::endl; */

    app.Run();

    for (auto h : h_list) 
    {
        delete h;
    }   
    delete legend;
    delete c;

    return 0;
}
