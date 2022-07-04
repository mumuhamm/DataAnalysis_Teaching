#include "RooWorkspace.h"
#include "RooCategory.h"
#include "TKey.h"
#include "RooExponential.h"
#include <map>
#include "TCut.h"
#include "RooHistPdf.h"
#include "RooHist.h"
#include "TVirtualPad.h"
#include "RooDataHist.h"
#include <string>
#include "TEventList.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "TFile.h"
#include "RooDataSet.h"


void cutflow_test(){

  TH1D * h1 = new TH1D("h1", "h1", 5, 0, 5);//you can chnage the bin by using h1->Rebin(binnumber);
  TRandom * r = new TRandom();
  double px, py, pz, pT;
  int index1;
  for (int i =0; i<10e5;i++){
    px = r->Uniform(0,100);
    py = r->Uniform(0,100);
    pT = sqrt(px*px+py*py); //sqrt(10000+10000)=141.42136, so pT ranges from 0 to 141.42GeV
     
    if (pT > 120)continue;
    index1++;
    h1->SetBinContent(1, index1++);
    if (pT > 100)continue;
    index1++;
    h1->SetBinContent(2, index1++);
    if (pT > 80)continue;
     index1++;
    h1->SetBinContent(3, index1++);
    if (pT > 60)continue;
     index1++;
    h1->SetBinContent(4, index1++);
    if (pT > 40)continue;
     index1++;
    h1->SetBinContent(5, index1++);
    
  }
  TCanvas *c = new TCanvas();
   TCanvas *cut = new TCanvas();
   h1->GetXaxis()->SetBinLabel(1,"pT>120");
   h1->GetXaxis()->SetBinLabel(2,"pT>100" );
   h1->GetXaxis()->SetBinLabel(3,"pT>80");
   h1->GetXaxis()->SetBinLabel(4,"pT>60");
   h1->GetXaxis()->SetBinLabel(5,"pT>40");
  h1->Draw();
}
TCanvas *cut = new TCanvas();
h1->GetXaxis()->SetBinLabel(1,"trhits>12");
h1->GetXaxis()->SetBinLabel(2,"trchi2<10" );
h1->GetXaxis()->SetBinLabel(3,"idhits>6");
h1->GetXaxis()->SetBinLabel(4,"idchi2<5");
h1->GetXaxis()->SetBinLabel(5,"ddg0<8");
h1->GetXaxis()->SetBinLabel(5,"dg0<20");
h1->Draw();
