// Author: Alberto Bragagnolo
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"

void batchVarComparison(TString filename1, TString filename2, TString nameVar)
{
  TFile *file1 = new TFile(filename1);
  TFile *file2 = new TFile(filename2);
  TTree *tree1 = (TTree*) file1->Get("OutTree");
  TTree *tree2 = (TTree*) file2->Get("OutTree");

  cout<<"nentries1 = "<<tree1->GetEntries()<<"  nentries2 = "<<tree2->GetEntries()<<endl;

  TString cut = "isBs == 1 && B_MassFromSV > 5.25 && B_MassFromSV < 5.55";
  cut += " && Mum_Pt > 3.5 && Mup_Pt > 3.5 && B_Pt > 11 && KmK_Pt > 1.2 && KpPi_Pt > 1.2";
  cut += " && abs(Mum_Eta) < 2.5 && abs(Mup_Eta) < 2.5 && abs(KmK_Eta) < 2.5 && abs(KpPi_Eta) < 2.5";
  cut += " && ctau > 0.007 && abs(Jpsi_Mass-3.0969) < 0.150 && abs(PhiKstar_Mass-1.01946) < 0.010";
  cut += " && B_VProb > 0.001 && HLT_JpsiMu == 1";
  cut += " && (((int(KpPi_Hits)/100)%10000)%100) >= 4";
  cut += " && (((int(KmK_Hits)/100)%10000)%100) >= 4";
  cut += " && ctauErr < 0.015";


  tree1->Draw(nameVar + ">>htmp1", cut);
  cout<<"tree1 drawn"<<endl;
  TH1 *h1 = (TH1*)gDirectory->Get("htmp1");
  TH1 *h2 = (TH1*)h1->Clone("h2");
  h2->Reset();

  tree2->Draw(nameVar + ">>h2", cut);
  cout<<"tree2 drawn"<<endl;

  Double_t Int1 = h1->Integral();
  Double_t Int2 = h2->Integral();

  cout<<"Integral histo 1 "<<Int1<<endl;
  cout<<"Integral histo 2 "<<Int2<<endl;

  Double_t scale1 = 1./h1->Integral();
  Double_t scale2 = 1./h2->Integral();
  h1->Scale(scale1);
  h2->Scale(scale2);

  cout<<"Drawing..."<<endl;

  TCanvas *c1 = new TCanvas();
  h1->SetTitle(nameVar);

  h1->Draw("hist");
  h2->SetLineColor(2);
  h2->Draw("same hist");

  gSystem->mkdir("plots");

  c1->Print("plots/" + nameVar + ".png");

  return;
}
