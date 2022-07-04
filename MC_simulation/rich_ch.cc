#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMath.h"
using namespace RooFit ;
using namespace std;


#define mc2_kaon  0.49368  // GeV / c^2
#define mc2_pion  0.13957  // GeV / c^2
#define mc2_muon  0.10565 // GeV / c^2
#define n1        1.03
#define n2        1.04
#define n3        1.05
float beta( float E, float mc2){
   float gamma = E/mc2 + 1;
   return (1-(1/(gamma*gamma)));
}

float theta ( float beta, float ri){
   return TMath::ACos(1/(beta*ri));
}

float RandomFloat(float a, float b) {
   float random = ((float) rand()) / (float) RAND_MAX;
   float diff = b - a;
   float r = random * diff;
   return a + r;
}

void rich_ch(){
   
   auto TP1 = new TProfile("TP1", "TP1; E (GeV); #theta_{ch} (rad)", 100, 0, 20, 0, TMath::Pi());
   auto TP2 = new TProfile("TP2", "TP2; E (GeV); #theta_{ch} (rad)", 100, 0, 20, 0, TMath::Pi());
   auto TP3 = new TProfile("TP3", "TP3; E (GeV); #theta_{ch} (rad)", 100, 0, 20, 0, TMath::Pi());
   
   auto TP11 = new TProfile("TP11", "TP11; E (GeV); #theta_{ch} (rad)", 100, 0, 20, 0, TMath::Pi());
   auto TP21 = new TProfile("TP21", "TP21; E (GeV); #theta_{ch} (rad)", 100, 0, 20, 0, TMath::Pi());
   auto TP31 = new TProfile("TP31", "TP31; E (GeV); #theta_{ch} (rad)", 100, 0, 20, 0, TMath::Pi());
   
   auto TP12 = new TProfile("TP12", "TP12; E (GeV); #theta_{ch} (rad)", 100, 0, 20, 0, TMath::Pi());
   auto TP22 = new TProfile("TP22", "TP22; E (GeV); #theta_{ch} (rad)", 100, 0, 20, 0, TMath::Pi());
   auto TP32 = new TProfile("TP32", "TP32; E (GeV); #theta_{ch} (rad)", 100, 0, 20, 0, TMath::Pi());
   
   
   for(unsigned i =0 ; i<=10e6; ++i){
      
      float E = RandomFloat(0, 20);
      
      float kaon_m = mc2_kaon;
      float muon_m = mc2_muon;
      float pion_m = mc2_pion;
      float kaon_beta = beta(E, kaon_m);
      float pion_beta = beta(E, pion_m);
      float muon_beta = beta(E, muon_m);
      
      TP1->Fill(E, theta(kaon_beta, n1), 1);
      TP2->Fill(E, theta(pion_beta, n1), 1);
      TP3->Fill(E, theta(muon_beta, n1), 1);
      
      TP11->Fill(E, theta(kaon_beta, n2), 1);
      TP21->Fill(E, theta(pion_beta, n2), 1);
      TP31->Fill(E, theta(muon_beta, n2), 1);
      
      TP12->Fill(E, theta(kaon_beta, n3), 1);
      TP22->Fill(E, theta(pion_beta, n3), 1);
      TP32->Fill(E, theta(muon_beta, n3), 1);
     // gr1 = new TGraphErrors(i,E,theta(kaon_beta, n1),E_err,0.1*theta(kaon_beta, n1));
      
   }
   
   
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   
   TCanvas *c = new TCanvas("c", "c", 600,600);
   TP1->SetMarkerStyle(34);TP1->SetMarkerColor(12);TP1->SetLineColor(12);TP1->Draw("E");
   TP2->SetMarkerStyle(34);TP2->SetMarkerColor(28);TP2->SetLineColor(28);TP2->Draw("E SAME");
   TP3->SetMarkerStyle(34);TP3->SetMarkerColor(38);TP3->SetLineColor(38);TP3->Draw("E SAME");
   auto legend1 = new TLegend(0.1,0.4,0.3,0.6);
   legend1->SetHeader("#theta_{ch} (E), n = 1.03","C");
   legend1->AddEntry(TP1,"K^{#pm}","p");
   legend1->AddEntry(TP2,"#pi^{#pm}","p");
   legend1->AddEntry(TP3,"#mu^{#pm}","p");
   legend1->Draw();
   float margin1 = 2.5;
   float margin2 = 15;
   float ymin = TP1->GetMinimum();
   float ymax = TP1->GetMaximum();
   TLine *line1 = new TLine(margin1,ymin,margin1,ymax);
   line1->SetLineWidth(2);
   line1->SetLineColor(kRed);
   line1->SetLineStyle(10);
   line1->Draw();
   TLine *line2 = new TLine(margin2,ymin,margin2,ymax);
   line2->SetLineWidth(2);
   line2->SetLineColor(kBlack);
   line2->SetLineStyle(10);
   line2->Draw();
   c->SetFillColor(0);
   c->SetBorderSize(2);
   c->SetLeftMargin(0.1422222);
   c->SetRightMargin(0.04444445);
   c->Update();
   
   TCanvas *c1 = new TCanvas("c1", "c1", 600,600);
   TP11->SetMarkerStyle(34);TP11->SetMarkerColor(12);TP11->SetLineColor(12);TP11->Draw("E");
   TP21->SetMarkerStyle(34);TP21->SetMarkerColor(28);TP21->SetLineColor(28);TP21->Draw("E SAME");
   TP31->SetMarkerStyle(34);TP31->SetMarkerColor(38);TP31->SetLineColor(38);TP31->Draw("E SAME");
   auto legend11 = new TLegend(0.1,0.4,0.3,0.6);
   legend11->SetHeader("#theta_{ch} (E), n=1.04","C");
   legend11->AddEntry(TP11,"K^{#pm}","p");
   legend11->AddEntry(TP21,"#pi^{#pm}","p");
   legend11->AddEntry(TP31,"#mu^{#pm}","p");
   legend11->Draw();
   float margin11 = 2.5;
   float margin21 = 15;
   float ymin1 = TP11->GetMinimum();
   float ymax1 = TP11->GetMaximum();
   TLine *line11 = new TLine(margin11,ymin1,margin11,ymax1);
   line11->SetLineWidth(2);
   line11->SetLineColor(kRed);
   line11->SetLineStyle(10);
   line11->Draw();
   TLine *line21 = new TLine(margin21,ymin1,margin21,ymax1);
   line21->SetLineWidth(2);
   line21->SetLineColor(kBlack);
   line21->SetLineStyle(10);
   line21->Draw();
   c1->SetFillColor(0);
   c1->SetBorderSize(2);
   c1->SetLeftMargin(0.1422222);
   c1->SetRightMargin(0.04444445);
   c1->Update();
   
   TCanvas *c2 = new TCanvas("c2", "c2", 600,600);
   TP12->SetMarkerStyle(34);TP12->SetMarkerColor(12);TP12->SetLineColor(12);TP12->Draw("E");
   TP22->SetMarkerStyle(34);TP22->SetMarkerColor(28);TP22->SetLineColor(28);TP22->Draw("E SAME");
   TP32->SetMarkerStyle(34);TP32->SetMarkerColor(38);TP32->SetLineColor(38);TP32->Draw("E SAME");
   auto legend12 = new TLegend(0.1,0.4,0.3,0.6);
   legend12->SetHeader("#theta_{ch} (E), n=1.05","C");
   legend12->AddEntry(TP12,"K^{#pm}","p");
   legend12->AddEntry(TP22,"#pi^{#pm}","p");
   legend12->AddEntry(TP32,"#mu^{#pm}","p");
   legend12->Draw();
   float margin12 = 2.5;
   float margin22 = 15;
   float ymin2 = TP12->GetMinimum();
   float ymax2 = TP12->GetMaximum();
   TLine *line12 = new TLine(margin12,ymin2,margin12,ymax2);
   line12->SetLineWidth(2);
   line12->SetLineColor(kRed);
   line12->SetLineStyle(10);
   line12->Draw();
   TLine *line22 = new TLine(margin2,ymin,margin2,ymax2);
   line22->SetLineWidth(2);
   line22->SetLineColor(kBlack);
   line22->SetLineStyle(10);
   line22->Draw();
   c2->SetFillColor(0);
   c2->SetBorderSize(2);
   c2->SetLeftMargin(0.1422222);
   c2->SetRightMargin(0.04444445);
   c2->Update();
   
   
}
