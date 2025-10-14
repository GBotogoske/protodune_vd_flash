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
#include <fstream>
#include <my_data.hh>
#include <readConfigFit.hh>
#include <utils.hh>


//numero de picos no fit
int Npeaks = 5;
// se deve calcular intersecao
bool calc_Intersection = true;

TH1D* hdata = nullptr;
double xmin, xmax;

//modelo do fit (soma de gaussianas)
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

// --- função χ² --- (funcao de minimizacao)
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
    std::string file_name;
    //base diretorio
    std::string base_dir = "/home/gabriel/Documents/protodune/data/VD/";
    std::string run_number = "null";
    std::string this_ch="39215";
    
    readConfigFit* fitConfig = nullptr;
    if (argc > 1) 
    {
        bool run_found = false;
        for (int i = 1; i < argc; ++i)
        {
            std::string arg = argv[i];
            if (arg == "-ch" && i + 1 < argc) // se o canal foi dado para calcular o spe
            {
                this_ch = argv[i + 1];
                run_number = get_map_spe2()[this_ch]; // busca a run adequada desse canal
                run_found = true;
                break;
            }
        }
        if(run_found)
        {
            file_name=search_file(base_dir,run_number,std::string("_spe.root")); // abre o arquivo
            fitConfig = new readConfigFit(std::stoi(this_ch));
        }
        else
        {
            file_name = std::string(argv[1]);
            fitConfig = new readConfigFit(-1);
        }
    } 
    else 
    {
        file_name = std::string("/home/gabriel/Documents/protodune/data/VD/np02vd_raw_run039357_0000_df-s04-d0_dw_0_20250915T151645_myAnalyser.root");
        fitConfig = new readConfigFit(-1);
    }

    std::cout << "FILE NAME: " << file_name << std::endl;
    std::cout << "RUN: " << run_number << std::endl;
    std::cout << "CH: " << this_ch << std::endl;

    fitConfig->print();

    //abre o arquivo
    TFile* file1 = TFile::Open(file_name.c_str(),"READ");
    //abre a arvore
    TTree* tree1 = (TTree*)file1->Get("T1");
    //criar variavel my_data para ler os dados salvos
    my_data* data = nullptr;

    TApplication app("app", &argc, argv);
    //aponta variavel de leitura ao branch correto do arquivo aberto
    tree1->SetBranchAddress("Data", &data);

    Npeaks = (int) fitConfig->getParam("Npeaks");
    int min_bin = (int) fitConfig->getParam("fit_min");
    int max_bin = (int) fitConfig->getParam("fit_max");
    int bins = (int) fitConfig->getParam("bins");

    //cria histograma
    TH1D* hist = new TH1D("hist","hist",bins,min_bin,max_bin);

    int my_channel = std::stoi(this_ch);
    float baseline_min = fitConfig->getParam("baseline_min");
    float baseline_max = fitConfig->getParam("baseline_max");
    float noise_max = fitConfig->getParam("noise_max");

    //preenche histograma com as variaveis adequadas e canal correto
    for(int i = 0; i < tree1->GetEntries(); ++i)
    {
        tree1->GetEntry(i);
        if (data->Channel==my_channel)
        {
            if(data->baseline>baseline_min && data->baseline<baseline_max && data->noise<noise_max)
                hist->Fill(data->integral);
        }
    }
    //calcula os erros, assumindo cada bin uma distribuicao de possion
    hist->Sumw2();
    hdata = hist;
    //pega o minimo e o maximo
    xmin = hist->GetXaxis()->GetXmin();
    xmax = hist->GetXaxis()->GetXmax();

    //seta numero de paraemetros para o fit
    int npars = 0;
    if(Npeaks<=2)
    {
        npars=3*Npeaks;
    }
    else
    {
        npars=6+(Npeaks-2)*2;
    }
    
    //cria a classe minuiti
    TMinuit minuit(npars);
    //seta a funcao de minimizacao
    minuit.SetFCN(fcn);

    //define os parametros inicias
    std::vector<double> p(npars); //A0,mean0,std0,meanspe,A1,sigma1,A2,sigma2....
    std::vector<double> A_i;
    for(int i=0;i<Npeaks;i++)
    {
        A_i.push_back(fitConfig->getParam(Form("A%d",i)));
    }
    double mean0=fitConfig->getParam("mean0");
    double meanspe=fitConfig->getParam("meanspe");
    double std0=fitConfig->getParam("std0");
    double std1=fitConfig->getParam("std1");
    for(int i=0;i<Npeaks;i++)
    {
        if(i==0)
        {
            p[0]=A_i[i];
            p[1]=mean0;
            p[2]=std0;
        
        }
        else if(i==1)
        {
            p[3]=meanspe;
            p[4]=A_i[i];
            p[5]=std1;
        }
        else
        {
            p[6+2*(i-2)]=A_i[i];
            p[6+2*(i-2)+1]=std1*sqrt(i);
        }        
    }

    //defeni os parametros inicias, steps e valores maximo e minimos    
    double step=fitConfig->getParam("step");
    for(int i=0;i<Npeaks;i++)
    {
        if(i==0)
        {
            minuit.DefineParameter(0, Form("Amp_%d",i),  p[0], step, fitConfig->getParam("A0min"), fitConfig->getParam("A0max"));
            minuit.DefineParameter(1, Form("Mean_%d",i), p[1], step,  fitConfig->getParam("mean0min"), fitConfig->getParam("mean0max"));
            minuit.DefineParameter(2, Form("Std_%d",i),  p[2], step, fitConfig->getParam("std0min"), fitConfig->getParam("std0max"));
        }
        else if(i==1)
        {
            minuit.DefineParameter(3, Form("Mean_%d",i), p[3], step, fitConfig->getParam("mean1min"), fitConfig->getParam("mean1max"));
            minuit.DefineParameter(4, Form("Amp_%d",i),  p[4], step, fitConfig->getParam("A1min"), fitConfig->getParam("A1max"));
            minuit.DefineParameter(5, Form("Std_%d",i),  p[5], step, fitConfig->getParam("std1min"), fitConfig->getParam("std1max"));
        }
        else
        {
            minuit.DefineParameter(6+2*(i-2), Form("Amp_%d",i),  p[4+2*(i-1)], step, fitConfig->getParam("Anmin"), fitConfig->getParam("Anmax"));
            minuit.DefineParameter(6+2*(i-2)+1, Form("Std_%d",i),  p[4+2*(i-1)], step, fitConfig->getParam("stdnmin"), fitConfig->getParam("stdnmax"));
        }     
   
    }
    
    //faz o fit
    minuit.Migrad();
    //calcula incerteza
    minuit.mnmnos();
    //pega os parametros de fit
    std::vector<double> p_fit(npars), e_fit(npars);
    for (int i = 0; i < npars; i++)
    {
        minuit.GetParameter(i, p_fit[i], e_fit[i]);
    }

    //cria histograma utilizando uma initicializao lambda
    TF1* ffit = new TF1("ffit",[](double* x, double* p){ return model(x, p); },xmin, xmax, npars);
    //seta os parametros fitados
    ffit->SetParameters(p_fit.data());
    //cria canvas e opcoes graficas
    TCanvas* c1 = new TCanvas("c1","Multi-Gaussian Fit",900,600);
    hist->SetTitle("Multi-Gaussian Fit;Integral [ADC];Entries");
    hist->Draw("");
    ffit->SetLineColor(kRed);
    ffit->SetLineWidth(2);
    //desenha fit
    ffit->Draw("SAME");

    //agora desenha cada gaussiana
    const int colors[6] = {kBlue, kGreen+2, kOrange+7, kViolet, kMagenta, kCyan};
    std::vector<TF1*> single_gaus;
    for (int i = 0; i < Npeaks; ++i)
    {
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
        //desenha gaussiana
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

    //calcula intersecao das gaussianas visinhas e desenha no plot
    if(calc_Intersection)
    {
        std::vector<Peak> peaks;
        peaks.reserve(Npeaks);

        Peak p0 { p_fit[0], p_fit[1], p_fit[2] };
        peaks.push_back(p0);

        double gain = p_fit[3];
        Peak p1 { p_fit[4], p_fit[1] + gain, p_fit[5] };
        peaks.push_back(p1);

        for (int i = 2; i < Npeaks; ++i) 
        {
            double Ai = p_fit[6 + 2*(i-2)];
            double si = p_fit[6 + 2*(i-2) + 1];
            Peak pi { Ai, p_fit[1] + i*gain, si };
            peaks.push_back(pi);
        }

        std::vector<double> xcuts;
        for (int i = 0; i < Npeaks-1; ++i) 
        {
            double xstar;
            if (intersect_two(peaks[i], peaks[i+1], xstar)) 
            {
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


        std::ofstream output_file(Form("/home/gabriel/Documents/protodune/protodune_vd/light_analysis/data/FIT_SPE/%d.txt", std::stoi(this_ch)));
        std::cout << "\nInterseções vizinhas (x*):\n";
        for (int i = 0; i < (int)xcuts.size(); ++i)
        {
            output_file << xcuts[i] << "\n";
            std::cout << "Peak " << i << " | " << (i+1) << " : x* = " << xcuts[i] << "\n";
        }
        c1->SaveAs(Form("/home/gabriel/Documents/protodune/protodune_vd/light_analysis/data/FIT_SPE/%d.png", std::stoi(this_ch)));
        output_file.close();
    }


    app.Run();
    return 0;
}
