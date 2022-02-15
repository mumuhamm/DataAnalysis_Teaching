#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <TH1.h>
#include <string>
#include <TMath.h>
#include <vector>

#include <TStyle.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TBuffer.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"
using namespace RooFit;
using namespace std;

void plot_efficiency(TString rootFiles, TString triggerFile)

{

  ifstream inRootFile;
  inRootFile.open(rootFiles.Data());
  while(!inRootFile.eof())
    {
    TString rootFile;
    inRootFile>>rootFile;
    if(inRootFile.eof()) break;
    
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("1.2f");
    
    TFile *f1 = new TFile(rootFile.Data());
    ifstream inFile;
    inFile.open(triggerFile.Data());
    TString trigger_name;
    TMultiGraph *mg = new TMultiGraph();
    TCanvas *c = new TCanvas("c","",200,10,600,500);
       //c->SetGridx();
       //c->SetGridy();
    TString SysType = rootFile.ReplaceAll("efficiency","");
    SysType.ReplaceAll(".root","");
    cout<<"systType : "<<SysType.Data()<<"  ,rootFile : "<<rootFile.Data()<<endl;

    while(!inFile.eof())
       {
       TH1F *h_total;
       TH1F *h_pass;
       TH2F *h2_total;
       TH2F *h2_pass;
       TEfficiency* grl=0 ;
       TGraphAsymmErrors * gra=0 ;
       inFile>>trigger_name;
       if(inFile.eof()) break;
       TString total = trigger_name.Data();
       total += "_total";
       if(total.Contains("_SingleMuOpen") ) total.ReplaceAll("_SingleMuOpen","");
       TString pass = trigger_name.Data();
       pass += "_pass";
       
       if(total.Contains("pt_eta"))
          {
          //c->SetLogy();
          TString textFileName = trigger_name.Data();
          textFileName += SysType;
          textFileName+="_efficiency.txt";
          h2_total = (TH2F*)f1->Get(total.Data());
          h2_pass = (TH2F*)f1->Get(pass.Data());
          h2_total->GetYaxis()->SetTitle("pt");
          h2_total->GetXaxis()->SetTitle("eta");
          grl = new TEfficiency(*h2_pass,*h2_total);
          grl->SetTitle(h2_total->GetTitle());
          grl->Draw("colztext");
    
          
          ofstream outfile(textFileName.Data());
          for(int ieta = 1; ieta <= h2_pass->GetNbinsX(); ieta++){
             for(int ipt = 1; ipt <= h2_pass->GetNbinsY(); ipt++){
                int globalBin = grl->GetGlobalBin(ieta,ipt);
                outfile << std::fixed;
                outfile << setprecision(2)<< h2_pass->GetXaxis()->GetBinLowEdge(ieta)<<"\t" << setprecision(2)<< h2_pass->GetXaxis()->GetBinUpEdge(ieta)<< "\t";
                outfile << setprecision(1)<<h2_pass->GetYaxis()->GetBinLowEdge(ipt) <<setprecision(1)<< "\t" << h2_pass->GetYaxis()->GetBinUpEdge(ipt) << "\t" << setprecision(3)<<"\t";
                outfile << std::fixed;
                outfile<< grl->GetEfficiency(globalBin)<<setprecision(3)<<"\t"<<grl->GetEfficiencyErrorLow(globalBin)<<setprecision(3)<<"\t"<<grl->GetEfficiencyErrorUp(globalBin)<<setprecision(3) <<endl;
                outfile << std::fixed;
             }
             
          }
          outfile.close();
          }
       else {
          
          c->SetLogy(0);
          h_total = (TH1F*)f1->Get(total.Data());
          h_pass = (TH1F*)f1->Get(pass.Data());
          if(total.Contains("pt")){h_total->GetXaxis()->SetTitle("pT in (GeV)");h_pass->GetXaxis()->SetTitle("pT in (GeV)");}
          if(total.Contains("eta")){h_total->GetXaxis()->SetTitle("eta");h_total->GetXaxis()->SetTitle("eta");}
          h_total->GetYaxis()->SetRangeUser(0,1);
          h_total->GetYaxis()->SetTitle("efficiency");
          h_pass->GetYaxis()->SetTitle("efficiency");
          h_pass->GetYaxis()->SetRangeUser(0,1);
             //        grl = new TEfficiency(*h_pass,*h_total);
             //   grl->SetLineColor(1);
             //   grl->SetMarkerStyle(20);
             //   grl->SetMarkerColor(1);
             //        grl->SetTitle(h_total->GetTitle());
             //      grl->Draw("APE");
          
          TString graph_name = trigger_name.Data();
          graph_name += SysType;
          graph_name += "_efficiency";
          gra = new TGraphAsymmErrors (h_pass,h_total);
          gra->SetLineColor(1);
          gra->SetName(graph_name.Data());
          gra->SetMarkerStyle(20);
          gra->SetMarkerColor(1);
          gra->SetTitle(h_total->GetTitle());
          gra->Draw("APE");
          gra->GetYaxis()->SetRangeUser(0,1);
          
       }
       
       TString pngFileName = trigger_name.Data() ;
       pngFileName += SysType;
       pngFileName +=  "_efficiency.pdf";
       TString rootFileName = trigger_name.Data() ;
       rootFileName += SysType;
       rootFileName +=  "_efficiency.root";
       c->SaveAs(pngFileName.Data());
       if(!rootFileName.Contains("pt_eta"))
          {   TFile file_out(rootFileName.Data(),"RECREATE");
             file_out.cd();
             gra->Write();
          }
       else
          {c->SaveAs(rootFileName.Data());}
       c->Clear();
       gSystem->Exec("mkdir -p results");
       TString dirName = "results/";
       dirName+=trigger_name.Data();
       gSystem->Exec("mkdir -p  "+ dirName);
       gSystem->Exec("mv " + trigger_name+"*  "+dirName);
       delete grl;
      // delete gra;
       }
    delete c;
    //delete mg;
    }

}
