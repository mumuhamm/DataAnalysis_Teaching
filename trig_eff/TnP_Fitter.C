#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH2.h"
#include "RooFitResult.h"
#include "TStyle.h"
#include "TDirectory.h"

using namespace RooFit;

TH1* massHisto_All(){
TFile *fIn1 = new TFile("Muon_sfefficiency_nominal.root");
TH1F* all_mass = (TH1F*) fIn1->Get("h_Jpsi_notOpen_total");
return all_mass;
}
TH1* massHisto_Pass(){
TFile *fIn2 = new TFile("Muon_sfefficiency_nominal.root");
TH1F* pass_mass = (TH1F*) fIn2->Get("h_Jpsi_notOpen_pass");
return pass_mass;
}


void TnP_Fitter(){
    TH1* all_jpsi = massHisto_All() ;
    TH1* pass_jpsi = massHisto_Pass();
    //double mean_pass = pass_jpsi->GetMean();
   // double meanerr_pass = pass_jpsi->GetMeanError();
   // std::cout<<" mean  "<< mean_pass << " mean error" << meanerr_pass<<"\n";
    RooRealVar *Jpsi_All = new RooRealVar("Jpsi_All"," ; m_{J/#psi} (#mu^{+}#mu^{-}) ) (GeV); Events",2.9, 3.3) ;
    RooDataHist * jpsia = new RooDataHist("jpsia","jpsia", *Jpsi_All , Import(*all_jpsi)) ;
    RooRealVar *Jpsi_Pass = new RooRealVar("Jpsi_Pass"," ; m_{J/#psi} (#mu^{+}#mu^{-}) ) (GeV); Events",2.9, 3.3) ;
    RooDataHist * jpsip = new RooDataHist("jpsip","jpsip", *Jpsi_Pass , Import(*pass_jpsi)) ;
    
    
    
    
    
    

   Jpsi_All->setRange("signalRange",2.9,3.3) ;
   RooRealVar *cbs = new RooRealVar("cbs","cbs",0.05,0.0,3.5) ;
   RooRealVar *cba = new RooRealVar("cba","cba",1,0,15) ;
   RooRealVar *cbn = new RooRealVar("cbn","cbn",3,0,47) ;
   RooRealVar *excon = new RooRealVar("excon","excon",1.4,0.0, 4.5);
   RooExponential *expo = new RooExponential("expo","expo", *Jpsi_All, *excon);
   RooRealVar *mean = new RooRealVar("mean","mean",3.09,2.9,3.1) ;
   //RooRealVar *sigma = new RooRealVar("sigma","sigma",0.12,0.,3.5) ;
   //RooGaussian *gauss = new RooGaussian("gauss","gauss",*Jpsi_All,  *mean, *sigma) ;
   RooRealVar *ns = new RooRealVar("ns","number signal events",0, all_jpsi->GetEntries());
   RooCBShape *cball = new RooCBShape("cball","cball",*Jpsi_All,*mean, *cbs,*cba,*cbn) ;
   RooRealVar *frac = new RooRealVar("frac1","frac1",0.4,0.0,1.0); 
   RooAddPdf *sigcbexpo = new RooAddPdf("sigcbexpo","sigcbexpo",RooArgList(*cball , *expo),RooArgList(*frac));
   RooExtendPdf *esig_all = new RooExtendPdf("esig_all","extended signal p.d.f",*sigcbexpo, *ns,"signalRange") ;
    RooPlot* frame_all = Jpsi_All->frame(Title("All m_{J/#psi} (#mu^{+}#mu^{-}) ) (GeV)")) ;
    jpsia->plotOn(frame_all) ; 
    esig_all->fitTo(*jpsia , Extended(1));
    esig_all->plotOn(frame_all) ;
    esig_all->plotOn(frame_all, LineColor(kBlue), LineWidth(2));
    Double_t chisquare = frame_all->chiSquare();
    cout<<"Chi square of fit is :"<< chisquare<< endl;
    esig_all->plotOn(frame_all, Components(*cball), LineColor(3), LineWidth(2), LineStyle(4));
    esig_all->plotOn(frame_all, Components(*expo), LineColor(2), LineWidth(2), LineStyle(2));
    
    
   Jpsi_Pass->setRange("signalRange_pass",2.9,3.3) ;
   RooRealVar *cbs_pass = new RooRealVar("cbs_pass","cbs",0.05,0.0,3.5) ;
   RooRealVar *cba_pass = new RooRealVar("cba_pass","cba",1,0,6) ;
   RooRealVar *cbn_pass = new RooRealVar("cbn_pass","cbn",2,0,20) ;
   RooRealVar *mean_pass = new RooRealVar("mean_pass","mean",3.09,2.9,3.1) ;
   RooRealVar *excon_pass = new RooRealVar("excon_pass","excon_pass",1.4,0.0, 4.5);
   RooRealVar *np = new RooRealVar("np","number signal events",0, pass_jpsi->GetEntries());
   RooExponential *expo_pass = new RooExponential("expo_pass","expo_pass", *Jpsi_Pass, *excon_pass);
   RooCBShape *cball_pass = new RooCBShape("cball_pass","cball_pass",*Jpsi_Pass,*mean_pass, *cbs_pass,*cba_pass,*cbn_pass) ;
   RooRealVar *frac2 = new RooRealVar("frac2","frac2",0.4,0.0,1.0); 
   RooAddPdf *sigcbexpo_pass = new RooAddPdf("sigcbgau_pass","sigcbgau_pass",RooArgList(*cball_pass , *expo_pass),RooArgList(*frac2));
   RooExtendPdf *esig_pass = new RooExtendPdf("esig_pass","extended signal p.d.f",*sigcbexpo_pass, *np,"signalRange_pass") ;
    RooPlot* frame_pass = Jpsi_Pass->frame(Title("Pass m_{J/#psi} (#mu^{+}#mu^{-}) ) (GeV)")) ;
    jpsip->plotOn(frame_pass) ; 
    esig_pass->fitTo(*jpsip , Extended(1));
    esig_pass->plotOn(frame_pass) ;
    esig_pass->plotOn(frame_pass, LineColor(kBlue), LineWidth(2));
    Double_t chisquare_p = frame_pass->chiSquare();
    cout<<"Chi square of fit is :"<< chisquare_p<< endl;
    esig_pass->plotOn(frame_pass, Components(*cball_pass), LineColor(3), LineWidth(2), LineStyle(4));
    esig_pass->plotOn(frame_pass, Components(*expo_pass), LineColor(2), LineWidth(2), LineStyle(2));
    
    TCanvas *c= new TCanvas("c", "c", 1800, 800);
    c->Divide(3,1);
    c->cd(1);
    frame_all->Draw();
    c->cd(2);
    frame_pass->Draw();
    
}
