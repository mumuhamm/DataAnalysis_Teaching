#include <TRandom3.h>
#include <TFile.h>
#include <TNtuple.h>

TRandom3 rnd;

void mksignal()
{
    TFile *fout = new TFile("signal.root","recreate");
    
    TNtuple *tree = new TNtuple("tree","tree","mass:pt:eta");
    
    int count = 0;
    do {
        float var[3];
        
        var[0] = rnd.Gaus(10.,0.1);
        var[1] = rnd.Exp(20.);
        var[2] = rnd.Gaus(0.,1.2);
        
        if (fabs(var[0]-10.)<1. && var[1]>8. && fabs(var[2])<2.8) {
            tree->Fill(var);
            count++;
        }
    }while(count<10000);
    
    fout->Write();
    fout->Close();
}

void mkdata()
{
    TFile *fout = new TFile("data.root","recreate");
    
    TNtuple *tree = new TNtuple("tree","tree","mass:pt:eta");
    
    int count[2] = {0,0};
    do {
        float var[3];
        
        var[0] = rnd.Gaus(10.,0.1);
        var[1] = rnd.Exp(20.);
        var[2] = rnd.Gaus(0.,1.2);
        
        if (fabs(var[0]-10.)<1. && var[1]>8. && fabs(var[2])<2.8) {
            tree->Fill(var);
            count[0]++;
        
            count[1] = 0;
            do {
                var[0] = rnd.Uniform(9.,11.);
                var[1] = rnd.Exp(8.);
                var[2] = rnd.Gaus(0.,5.);
            
                if (fabs(var[0]-10.)<1. && var[1]>8. && fabs(var[2])<2.8) {
                    tree->Fill(var);
                    count[1]++;
                }
            }while (count[1]<20);
        }
    }while(count[0]<2500);
    
    fout->Write();
    fout->Close();
}

void mksample()
{
    mksignal();
    mkdata();
}