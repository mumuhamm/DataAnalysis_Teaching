//This is .cc file

//#include "TFile.h"
#include <sstream>
#include <iostream>
#include <string>
#include <fstream> 
#include <vector>
#include <map>
#include <cmath>
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TBuffer.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"
#include "TRandom.h"
#include "practice.h"

using namespace std;



int main(){
   TRandom *r = new TRandom();
   int n =10E2;
   float eta_bins[9] = {-2.4,-1.8,-1.2,-0.6,0,0.6, 1.2,1.8,2.4};
   float pt_bins[10] = {0,2,4,6,8,10,12,14,17,20};
   TH1F *h_muon_pt = new TH1F("h_muon_pt","muon_pt",9,pt_bins);
   TH1F *h_muon_eta = new TH1F("h_muon_eta","muon_eta",8,eta_bins);
   
   
   analysis pt_plot;
   float min = 0.0, max = 25.0;
   float pt_val, px, py;
   for(int i = 0; i<n; ++i)
      {
      px = r->Gaus(0.8, 2.78);
      py = r->Gaus(0.9,3.75);
      std::cout<<px <<": \t <= px & py =>  \t : "<<py<<"\n";
      pt_val = pt_plot.PT(px, py);
      std::cout<<" Transverse mom values: \t " << pt_val << "\n";
      h_muon_pt->Fill(pt_val);
      }
   TCanvas *c = new TCanvas();
   h_muon_pt->SetMarkerStyle(8);
   h_muon_pt->SetMarkerColor(kBlue);
   h_muon_pt->GetXaxis()->SetTitle("#mu_{p_{T}} (GeV)");
   h_muon_pt->GetYaxis()->SetTitle("Number of Events (a.u.)");
   h_muon_pt->Draw("E1");
   c->SaveAs("pt_plot.png","png");
   
}

