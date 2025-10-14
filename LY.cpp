
#include <TCanvas.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TAxis.h>

#include <iostream>
#include <fstream> 

#define drift_distance 350


std::vector<int> runs = {39510,39511,39512,39514,39515,39516,39517,39518
    ,39519,39521,39522,39523,39525,39526,39526,39528,39529};

std::vector<float> Voltage = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0,155.7};
std::vector<int> list_channel = {1050,1051,1060,1061,1070,1071,1080,1081,2030,2031,2040,2041,2070,2071,2080,2081};

int main(int argc, char** argv)
{
    std::string home_dir = "/home/gabriel/Documents/protodune/protodune_vd/light_analysis/waveforms";    
    //pegar o numero de entradas de campoe eletrico e canais
    int n_runs = runs.size();
    int n_ch = list_channel.size();

    //lista de sinais [run][ch]
    std::vector<std::vector<std::vector<double>>> signals(n_runs,std::vector<std::vector<double>>(n_ch,std::vector<double>(1024, 0.0)));

    //varre todos os arquivos
    int i=0;
    for(int run:runs)
    {
        int j=0;
        for(int ch:list_channel)
        {
            std::string file_wf = home_dir + Form("/%d/waveform_%d.txt",run,ch);             
            double s,t;
            std::ifstream data_wf(file_wf);
            int k=0;
            while(data_wf >> t >> s)
            {
                signals[i][j][k]=s;
                k++; 
            }
            data_wf.close();
            j++;
        }
        i++;
    }
    
    //lista de integrais [ch][run]
    std::vector<std::vector<double>> integral(n_ch,std::vector(n_runs,0.0));
    i=0;
    for(int run:runs)
    {
        int j=0;
        for(int ch:list_channel)
        {
            double sum=0.0;
            int init=50;
            int end=400;
            for(int k=init;k<=end;k++)
            {
                integral[j][i]+=signals[i][j][k];
            }
            j++;
        }
        i++;
    }

    //normalizar a integral
    for(int ix=0;ix<n_ch;ix++)
    {
        double ref = integral[ix][0]; 
        for(int iy=0;iy<n_runs;iy++)
        {
            integral[ix][iy]=integral[ix][iy]/ref;
        }
    }
    
    
    //hora de plotar
    std::string out_folder = "/home/gabriel/Documents/protodune/protodune_vd/light_analysis/data/LY";
    std::vector<double> hv;
    for(int ix=0;ix<n_runs;ix++)
    {
        hv.push_back(Voltage[ix]/drift_distance);
    }

    for(int ix=0;ix<n_ch;ix++)
    {
        
        TCanvas* c = new TCanvas("c", "LY X EF", 800, 600);
        int ch = list_channel[ix];
        TGraph* g = new TGraph(n_runs, hv.data(), integral[ix].data());
        g->SetTitle(Form("LY X EF %d", ch));
        g->GetXaxis()->SetTitle("Electric Field [kV/cm]");
        g->GetYaxis()->SetTitle("Integral");
        g->SetMarkerStyle(20);    // bolinha sólida
        g->SetMarkerSize(1.2);    // tamanho das bolinhas
        g->SetMarkerColor(kBlue); // azul
        g->SetLineColor(0);
        g->Draw("AP");

        std::string fig_name = out_folder+Form("/LY_%d.png", ch);
        c->SaveAs(fig_name.c_str());

        delete g;
        delete c;
 
        
    }

    
}