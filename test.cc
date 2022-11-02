#include <iostream>
#include<vector>
#include <list>
#include <cstdlib>
#include <math.h>
#include "TH1D.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TMath.h"
#include "TRandom.h"
#define PI TMath::Pi()
#define N

#define epsilon 0.1
#define epoch 2000

using namespace std;

void test() {
  gROOT->Reset();
  gStyle->SetOptStat(0);
  TRandom * num = new TRandom();
  TCanvas *c1= new TCanvas;

  TH1D* t0 = new TH1D("t0",";x (a.u.); Events",100,-0.5,0.5);
  for (int m=0; m<10000; ++m) {
    double random = num->Gaus(0,0.07);
    t0->Fill(random);
     
   
  }
   t0->Fit("gaus");
   int fitStatus =  t0->Fit("gaus");
   TFitResultPtr r = t0->Fit("gaus","S");
   TMatrixDSym cov = r->GetCovarianceMatrix();
   Double_t chi2   = r->Chi2();
   Double_t par0   = r->Parameter(0);
   Double_t err0   = r->ParError(0);
   r->Print("V");
  
   t0->SetLineColor(kBlack);
   t0->SetFillColor(kGreen );
   t0->SetFillStyle(3353);
  t0->Draw("COLZ");
  c1->Update();
  TLine *l=new TLine(0,c1->GetUymin(),0.0,c1->GetUymax());
  l->SetLineColor(kBlue);
  l->Draw();
}




#include <iostream>
#include <string>
#include "RootFileReader.h"

   

RootFileReader::RootFileReader(const std::string& FileName, const std::string& TreeName) :
m_FileName(FileName),
m_TreeName(TreeName),
m_Fin(m_FileName.c_str()),
m_Iter(0)
{std::cout<<"RootFileReader constructed"<<std::endl;}

bool RootFileReader::OpenFile()
{
   if (!m_Fin.IsOpen()){return false;}
   m_Tin = (TTree *)m_Fin.Get(m_TreeName.c_str());
   if (!m_Tin){return false;}
   return true;
}

void RootFileReader::ListTreeBranches()
{
   for (int ibranch=0;
        ibranch<m_Tin->GetListOfBranches()->GetEntries();
        ibranch++){
      std::cout << m_Tin->GetListOfBranches()->At(ibranch)->GetName() << std::endl;
   }
}

void RootFileReader::FillBranchNames()
{
   for (int ibranch=0;
        ibranch<m_Tin->GetListOfBranches()->GetEntries();
        ibranch++){
      BranchNames.insert(std::string(m_Tin->GetListOfBranches()->At(ibranch)->GetName()));
   }
}

int RootFileReader::GetEntries()
{
   if (!m_Tin){return -1;}
   return m_Tin->GetEntries();
}

  

int RootFileReader::GetEntry(int Entry)
{
   if (!m_Tin){return -1;}
   return m_Tin->GetEntry(Entry);
   m_Tin->Show(0,1);
}

  

void RootFileReader::StartEntry()
{
   m_IsValidEntry=false;
   if (!m_Tin){return;}
   m_Iter = 0;
   if (GetEntry(m_Iter)){m_IsValidEntry=true;}
}

int RootFileReader::EndEntry()
{
   return m_IsValidEntry;
}

void RootFileReader::NextEntry()
{
   m_IsValidEntry=false;
   if (!m_Tin){return;}
   m_Iter++;
   if (GetEntry(m_Iter)){m_IsValidEntry=true;}
}




#ifndef ROOTFILEREADER_H__
#define ROOTFILEREADER_H__

#include <set>
#include <string>
#include <fstream>

#include "TFile.h"
#include "TTree.h"

class RootFileReader
{
   RootFileReader() : m_FileName("") {;}
   
protected:
   const std::string m_FileName;
   const std::string m_TreeName;
   TFile m_Fin;
   TTree * m_Tin;
   
   int m_Iter;
   bool m_IsValidEntry;
   
   std::set<std::string> BranchNames;
   
   void FillBranchNames();
   
public:
   
   RootFileReader(const std::string& FileName, const std::string& TreeName);
   bool OpenFile();
   bool CheckFile()                         {return m_Fin.IsOpen();}
   void ListTreeBranches();
   
      
   int GetEntries();
   int GetEntry(int);
   
     
   void NextEntry();
   void StartEntry();
   int EndEntry();
   
};

#endif



ROOT::EnableImplicitMT();
ROOT::RDataFrame df("T", "data/425287.root");
std::cout << df.GetColumnType("SingleMuons.DDG0") << std::endl;
auto myHisto = df.Histo1D("SingleMuons.DDG0");
TCanvas *c = new Tcanvas();
myHisto->Draw();
