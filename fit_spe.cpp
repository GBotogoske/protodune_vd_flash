#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TTree.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TApplication.h>
#include <TMath.h>
#include <TLine.h>
#include <TLatex.h>

#include <iostream>
#include <my_data.hh>
#include <utils.hh>

int Npeaks = 6;
bool calc_Intersection = true;

TH1D* hdata = nullptr;
double xmin, xmax;

double model(double* x, const double* par)
{
    double result = 0;
    for (int i = 0; i < Npeaks; ++i)
    {
        double amp,mean,sigma;
        if(i==0)
        {
            amp = par[0];
            mean = par[1];
            sigma = par[2];
        }
        else if(i==1)
        {
            mean = par[3]+par[1];
            amp = par[4];
            sigma = par[5];
        }
        else
        {
            mean = i*par[3]+par[1];
            amp = par[6+2*(i-2)];
            sigma = par[6+2*(i-2)+1];
        }
        result += amp * TMath::Gaus(x[0], mean, sigma, true);
    }
    return result;
}

// --- função χ² ---
void fcn(Int_t& npar, Double_t* grad, Double_t& fval, Double_t* par, Int_t flag)
{
    double chi2 = 0;
    for (int i = 1; i <= hdata->GetNbinsX(); ++i)
    {
        double x = hdata->GetBinCenter(i);
        double y = hdata->GetBinContent(i);
        double e = hdata->GetBinError(i);
        if (e == 0) continue;
        double xx = x;
        double yfit = model(&xx, par);

        chi2 += ((y - yfit) * (y - yfit)) / (e * e);
    }
    fval = chi2;
}

int main(int argc, char** argv)
{
    std::string filename = "/home/gabriel/Documents/protodune/protodune_vd/light_analysis/build/data_analysed.root";
    TFile* file1 = TFile::Open(filename.c_str(),"READ");
    TTree* tree1 = (TTree*)file1->Get("T1");
    my_data* data = nullptr;

    TApplication app("app", &argc, argv);

    tree1->SetBranchAddress("Data", &data);

    TH1D* hist = new TH1D("hist","hist",200,-280,2100);
    for(int i = 0; i < tree1->GetEntries(); ++i)
    {
        tree1->GetEntry(i);
        if(data->baseline>4100 && data->baseline<4700 && data->noise<7.5)
            hist->Fill(data->integral);
    }
    hist->Sumw2();
    hdata = hist;

    xmin = hist->GetXaxis()->GetXmin();
    xmax = hist->GetXaxis()->GetXmax();

    int npars = 3*Npeaks;
    TMinuit minuit(npars);
    minuit.SetFCN(fcn);

    std::vector<double> p(npars); //A0,mean0,std0,meanspe,A1,sigma1,A2,sigma2....
    std::vector<double> A_i = {7700, 150, 100, 100,100,100};
    for(int i=0;i<Npeaks;i++)
    {
        if(i==0)
        {
            p[0]=A_i[i];
            p[1]=0;
            p[2]=200;
        
        }
        else if(i==1)
        {
            p[3]=400;
            p[4]=A_i[i];
            p[5]=200;
        }
        else
        {
            p[6+2*(i-2)]=A_i[i];
            p[6+2*(i-2)+1]=200*sqrt(i);
        }        
    }

    double step=0.1;
    for(int i=0;i<Npeaks;i++)
    {
        if(i==0)
        {
            minuit.DefineParameter(0, Form("Amp_%d",i),  p[0], step, 1000, 40000);
            minuit.DefineParameter(1, Form("Mean_%d",i), p[1], step, -200, 10);
            minuit.DefineParameter(2, Form("Std_%d",i),  p[2], step, 50, 120);
        }
        else if(i==1)
        {
            minuit.DefineParameter(3, Form("Mean_%d",i), p[3], step, 100, 700);
            minuit.DefineParameter(4, Form("Amp_%d",i),  p[4], step, 0, 0);
            minuit.DefineParameter(5, Form("Std_%d",i),  p[5], step, 0, 0);
        }
        else
        {
            minuit.DefineParameter(6+2*(i-2), Form("Amp_%d",i),  p[4+2*(i-1)], step, 0, 0);
            minuit.DefineParameter(6+2*(i-2)+1, Form("Std_%d",i),  p[4+2*(i-1)], step, 0, 0);
        }     
   
    }

    minuit.Migrad();
    minuit.mnmnos();

    std::vector<double> p_fit(npars), e_fit(npars);
    for (int i = 0; i < npars; i++)
    {
        minuit.GetParameter(i, p_fit[i], e_fit[i]);
    }
      

    TF1* ffit = new TF1("ffit",[](double* x, double* p){ return model(x, p); },xmin, xmax, npars);    
    ffit->SetParameters(p_fit.data());

    TCanvas* c1 = new TCanvas("c1","Multi-Gaussian Fit",900,600);
    //hist->SetMarkerStyle(20);
    //hist->SetMarkerSize(0.7);
    hist->SetTitle("Multi-Gaussian Fit;Integral [ADC];Entries");
    hist->Draw("");
    ffit->SetLineColor(kRed);
    ffit->SetLineWidth(2);
    ffit->Draw("SAME");

    const int colors[6] = {kBlue, kGreen+2, kOrange+7, kViolet, kMagenta, kCyan};
    std::vector<TF1*> single_gaus;

    for (int i = 0; i < Npeaks; ++i)
    {
        // cria nome único pra cada gaussiana
        TString name = Form("gaus_%d", i);
        
        TF1* fgaus = nullptr;

        if(i==0)
        {
            fgaus = new TF1(name,
            [=](double* x, double*){
                return p_fit[0] * TMath::Gaus(x[0], p_fit[1], p_fit[2], true);
            },
            xmin, xmax, 0);
        }
        else if(i==1)
        {
            fgaus = new TF1(name,
            [=](double* x, double*){
                return p_fit[4] * TMath::Gaus(x[0], p_fit[3]+p_fit[1], p_fit[5], true);
            },
            xmin, xmax, 0);
        }
        else
        {
            fgaus = new TF1(name,
            [=](double* x, double*){
                return p_fit[6+2*(i-2)] * TMath::Gaus(x[0], i*p_fit[3]+p_fit[1], p_fit[6+2*(i-2)+1], true);
            },
            xmin, xmax, 0);
        }

        fgaus->SetLineColor(colors[i % 6]);
        fgaus->SetLineWidth(2);
        fgaus->SetLineStyle(2); // linha tracejada
        fgaus->Draw("SAME");
        single_gaus.push_back(fgaus);
    }

    // --- legenda bonitinha ---
    auto legend = new TLegend(0.65, 0.7, 0.9, 0.9);
    legend->AddEntry(ffit, "Total Fit", "l");
    for (int i = 0; i < Npeaks; ++i)
    {
        legend->AddEntry(single_gaus[i], Form("Peak %d", i), "l");
    }    
    legend->Draw();

    std::cout << "\n===== Fit Results =====\n";
    for (int i = 0; i < Npeaks; i++)
    {
        if(i==0)
        {
            std::cout << "Amp[" << i << "] = " << p_fit[0]   << " ± " << e_fit[0]   << "\n";
            std::cout << "Mean[" << i << "] = " << p_fit[1] << " ± " << e_fit[1] << "\n";
            std::cout << "Std[" << i << "] = " << p_fit[2]  << " ± " << e_fit[2]  << "\n";
        }
        else if(i==1)
        {
            std::cout << "Mean[" << i << "] = " << p_fit[3] << " ± " << e_fit[3] << "\n";
            std::cout << "Amp[" << i << "] = " << p_fit[4]   << " ± " << e_fit[4]   << "\n";
            std::cout << "Std[" << i << "] = " << p_fit[5]  << " ± " << e_fit[5]  << "\n";
        }
        else
        {
            std::cout << "Amp[" << i << "] = " << p_fit[6+2*(i-2)]   << " ± " << e_fit[6+2*(i-2)]   << "\n";
            std::cout << "Std[" << i << "] = " << p_fit[6+2*(i-2)+1]  << " ± " << e_fit[6+2*(i-2)+1]  << "\n";
        }
        
    }

    if(calc_Intersection)
    {
        std::vector<Peak> peaks;
        peaks.reserve(Npeaks);

        Peak p0 { p_fit[0], p_fit[1], p_fit[2] };
        peaks.push_back(p0);

        double gain = p_fit[3];
        Peak p1 { p_fit[4], p_fit[1] + gain, p_fit[5] };
        peaks.push_back(p1);

        for (int i = 2; i < Npeaks; ++i) {
            double Ai = p_fit[6 + 2*(i-2)];
            double si = p_fit[6 + 2*(i-2) + 1];
            Peak pi { Ai, p_fit[1] + i*gain, si };
            peaks.push_back(pi);
        }

        std::vector<double> xcuts;
        for (int i = 0; i < Npeaks-1; ++i) {
            double xstar;
            if (intersect_two(peaks[i], peaks[i+1], xstar)) {
                xcuts.push_back(xstar);
                // desenhar linha vertical
                auto line = new TLine(xstar, 0, xstar, hist->GetMaximum()*1.05);
                line->SetLineStyle(3);
                line->SetLineColor(kGray+2);
                line->Draw("SAME");
                auto lat = new TLatex(xstar, hist->GetMaximum()*1.07, Form("x_{%d|%d}", i, i+1));
                lat->SetTextAlign(21);
                lat->SetTextSize(0.03);
                lat->Draw("SAME");
            }
        }

        std::cout << "\nInterseções vizinhas (x*):\n";
        for (int i = 0; i < (int)xcuts.size(); ++i)
            std::cout << "Peak " << i << " | " << (i+1) << " : x* = " << xcuts[i] << "\n";
    }


    app.Run();
    return 0;
}
