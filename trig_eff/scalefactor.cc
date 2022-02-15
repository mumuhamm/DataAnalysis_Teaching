#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
using namespace RooFit ;
using namespace std;

void scalefactor(){

     
Double_t mc_eff[8];
Double_t mc_effErr[8];
Double_t data_eff[8];
Double_t  data_effErr[8];
Double_t sf[8];
Double_t sfErr[8];
Double_t pT[8];
Double_t pTErr[8];
 const Int_t n = 8;
 
    TFile *fIn1 = new TFile("/home/cms/GGToTauTau/GGToTauTau/muonstudy/PBPB_MuonSF15DataMC/Data/results/Mu0_notOpen_pt/Mu0_notOpen_ptMuon_sf_nominal_efficiency.root");
    TFile *fIn2 = new TFile("/home/cms/GGToTauTau/GGToTauTau/muonstudy/PBPB_MuonSF15DataMC/MonteCarlo/results/Mu0_notOpen_pt/Mu0_notOpen_ptMuon_sf_nominal_efficiency.root");   
     	TGraphAsymmErrors * num_data = (TGraphAsymmErrors *) fIn1->Get("Mu0_notOpen_ptMuon_sf_nominal_efficiency");
     	TGraphAsymmErrors * deno_mc = (TGraphAsymmErrors *) fIn2->Get("Mu0_notOpen_ptMuon_sf_nominal_efficiency");

   auto c = new TCanvas("c", "c", 0, 0, 1000, 600);
    c->Divide(2,1);
  c->cd(1);
   num_data->SetLineColor(kBlue);
   num_data->SetLineWidth (2); 
   num_data->SetMarkerStyle(8);
   num_data->SetMarkerColor(kBlue);
   num_data->SetMarkerSize(1.5);
   num_data->Draw();
   deno_mc->SetLineColor(kRed);
   deno_mc->SetLineWidth (2);
   deno_mc->SetMarkerStyle(8);
   deno_mc->SetMarkerColor(kRed);
   deno_mc->SetMarkerSize(1.5);
   deno_mc->Draw( "SAME");
   auto legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry(num_data,"Data");
   legend->AddEntry(deno_mc,"J/#psi MC");
   legend->Draw();
 
   
   Int_t numpoints = num_data->GetN();
   std::cout<<" points in numpoints" << numpoints <<"\n";
   Int_t denopoints = deno_mc->GetN();
   std::cout<<" points in denopoints" << denopoints <<"\n";
   
  
   for (int j=0; j < numpoints-1; j++){

      double x,y, ex, ey;
      num_data->GetPoint(j,x,y);
      ex = num_data->GetErrorX(j);
      ey = num_data->GetErrorY(j);
      std::cout<<"the i "<<j<<"\t"<<x<<"\t"<<y<<"\t"<<ex<<"\t"<<ey<<"\n";
      data_eff[j] = y;
      data_effErr[j]= ey;
      pT[j]= x;
      pTErr[j]=ex;
      

  }
   
   
   
     for (int j=0; j < numpoints-1; j++){ 

       double x,y, ex, ey;
      deno_mc->GetPoint(j,x,y);
      ex = deno_mc->GetErrorX(j);
      ey = deno_mc->GetErrorY(j);
      std::cout<<"the i "<<j<<"\t"<<x<<"\t"<<y<<"\t"<<ex<<"\t"<<ey<<"\n";
      mc_eff[j] = y;
      mc_effErr[j]= ey;

  }
    
   for (int i=0; i<8; ++i )
    {
        
        sf[i] = data_eff[i]/mc_eff[i];
        sfErr[i] =sqrt(pow(data_effErr[i]/data_eff[i], 2)+ pow(mc_effErr[i]/mc_eff[i], 2));
        std::cout<<sf[i]<<"pm"<<sfErr[i]<<"\n";
        
    } 
  
    
   
  c->cd(2);
   auto gr = new TGraphErrors(n,pT,sf,pTErr,sfErr);
   gr->SetTitle("#mu scale factor");
   gr->SetMarkerColor(kBlack);
   gr->SetMarkerStyle(8);
   gr->Draw("AP");
     



}          
