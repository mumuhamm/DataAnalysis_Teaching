// Author :  Muhammad Alibordi



#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
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
#include "TTree.h"
#include "TH2D.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooAbsPdf.h"
#include "RooPolynomial.h"
#include "RooGaussModel.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooDecay.h"
#include "RooDataHist.h"
#include "RooAddModel.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooPlot.h"
#include "TH1D.h"
#include "TRandom.h"
#include "RooMinuit.h"
#include "RooExtendPdf.h"
#include "RooChi2Var.h"
#include "Math/Functor.h"
#include "TRandom3.h"
#include "Math/DistFunc.h"
#include "RooClassFactory.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooRealConstant.h"
#include "RooConstVar.h"
#include "Roo1DTable.h"
#include "RooBDecay.h"
#include "RooFormulaVar.h"
#include "RooRealSumPdf.h"
#include "Math/SpecFunc.h"
#include "RooBMixDecay.h"
#include "RooBCPEffDecay.h"
#include "Riostream.h"
#include "RooRandom.h"
#include "TMath.h"
#include "RooMCStudy.h"
#include "RooArgSet.h"
#include "RooLegendre.h"
#include "RooSpHarmonic.h"
#include "RooBifurGauss.h"
#include "complex.h"
#include <iostream>

using namespace RooFit; 
using namespace std;

void likelihood_gl(){
   Int_t nbins = 100;
  //RooRealVar *ReVL = new RooRealVar("ReVL","ReVL",-2.489, 0.504);
  RooRealVar *RegL = new RooRealVar("RegL","RegL",-1.36, 0.83);//constrcutor
  RooRealVar *ImgL = new RooRealVar("ImgL","ImgL",-0.14, 0.49);

// float myparvalue =   RegL->getValV();
  // std::cout<<myparvalue<<"\n";
   
  RooRealVar * mean1 = new RooRealVar("mean1","Re(g_{R}^{V})", -1.3, .8);
  RooRealVar * sigma1 = new RooRealVar("sigma1","sigma1", 0, 1.5);
  RooGaussian *gaus1 = new RooGaussian("gaus1","gaus1", *RegL, *mean1,*sigma1);
    
  RooRealVar * mean2 = new RooRealVar("mean2","Im(g_{R}^{V})", -0.1, 0.4);
  RooRealVar * sigma2 = new RooRealVar("sigma2","sigma2", 0, 0.3);
  RooGaussian *gaus2 = new RooGaussian("gaus2","gaus2", *ImgL, *mean2,*sigma2);
    
  RooRealVar *frac = new RooRealVar("frac", "frac", 0,1);
    
  RooAddPdf * genmodel = new RooAddPdf("genmodel", "genmodel", RooArgList(*gaus1, *gaus2), RooArgList(*frac), true);
  RooDataSet* data = genmodel->generate(RooArgSet(*RegL,*ImgL),100000);
   
   //Extracting histograms
   auto glr_h =(TH1F*)data->createHistogram("glr_h", *RegL);
   auto gli_h =(TH1F*)data->createHistogram("gli_h", *ImgL);
   
  //Plotting histograms
   RooPlot* glr = RegL->frame(Title("#italic{Re}(g^{V}_{L})"),Bins(nbins));
   data->plotOn(glr,DataError(RooAbsData::SumW2));
   RooPlot* gli = ImgL->frame(Title("#italic{Im}(g^{V}_{L})"),Bins(nbins));
   data->plotOn(gli,DataError(RooAbsData::SumW2));
   TCanvas *lz = new TCanvas();
   lz->Divide(2,1);
   lz->cd(1);glr->Draw();lz->cd(2);gli->Draw();
   
   
   
   RooRealVar * mean11 = new RooRealVar("mean11","Re(g_{R}^{V})", -1.3, .8);
   RooRealVar * sigma11 = new RooRealVar("sigma11","sigma1", 0, 1.5);
   RooGaussian *gaus11 = new RooGaussian("gaus11","gaus1", *RegL, *mean11,*sigma11);
   
   RooRealVar * mean22 = new RooRealVar("mean22","Im(g_{R}^{V})", -0.1, 0.4);
   RooRealVar * sigma22 = new RooRealVar("sigma22","sigma2", 0, 0.3);
   RooGaussian *gaus22 = new RooGaussian("gaus22","gaus2", *ImgL, *mean22,*sigma22);
   
   
   
   //Fit results
   
   RooAddPdf * fitmodel = new RooAddPdf("fitmodel", "fitmodel", RooArgList(*gaus11, *gaus22), RooArgList(*frac), true);
   RooFitResult* fitRes = fitmodel->fitTo(*data,Save());//data_SigReg
   fitRes->Print("v");
   
   
  //============================================================real part of the gl
   RooPlot* glr1 = RegL->frame(Title("#italic{Re}(g^{V}_{L})"),Bins(nbins));
   data->plotOn(glr1,DataError(RooAbsData::SumW2));
   fitmodel->plotOn(glr1) ;
   fitmodel->paramOn(glr1);
   fitmodel->plotOn(glr1, LineColor(kBlue), LineWidth(1));
   
   RooPlot* pullglr = RegL->frame(RooFit::Title("glr pull"));
   RooHist* hpull1 = glr1->pullHist();
   
   pullglr->addPlotable(hpull1,"P0") ;
   pullglr->SetMinimum(-3) ;
   pullglr->SetMaximum(+3) ;
   pullglr->SetYTitle("pull");
   pullglr->SetMarkerStyle(20);
   pullglr->SetNdivisions(10);
   Double_t chisquare_glr1 = glr1->chiSquare();
   cout<<"Chi square of glr fit is :"<< chisquare_glr1<< endl;
   
   
   
   TCanvas *c = new TCanvas("c", "c",0,0,600,600);
   TPad *pad1 = new TPad("pad1","pad1",0,0.36,1,1);
   TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.25);
   pad1->SetBottomMargin(0.00001);
   pad1->SetBorderMode(0);
   pad2->SetTopMargin(0.00001);
   pad2->SetBottomMargin(0.1);
   pad2->SetBorderMode(0);
   pad1->Draw();
   pad2->Draw();
   pad1->cd();
   gStyle->SetOptTitle(0);
   c->SetFillColor(0);
   c->SetBorderSize(2);
   c->SetLeftMargin(0.1422222);
   c->SetRightMargin(0.04444445);
   glr1->SetStats(0);
   glr1->Draw();
   pad2->cd();
   pullglr->SetStats(0);
   pullglr->Draw();
   c->cd();
   
 
   
      //============================================================imaginary part of the gl
   RooPlot* gli1 = ImgL->frame(Title("#italic{Im}(g^{V}_{L})"),Bins(nbins));
   data->plotOn(gli1,DataError(RooAbsData::SumW2));
   fitmodel->plotOn(gli1) ;
   fitmodel->paramOn(gli1);
   fitmodel->plotOn(gli1, LineColor(kBlue), LineWidth(1));
   
   RooPlot* pullgli = ImgL->frame(RooFit::Title("gli pull"));
   RooHist* hpull2 = gli1->pullHist();
   
   pullgli->addPlotable(hpull2,"P0") ;
   pullgli->SetMinimum(-3) ;
   pullgli->SetMaximum(+3) ;
   pullgli->SetYTitle("pull");
   pullgli->SetMarkerStyle(20);
   pullgli->SetNdivisions(10);
   Double_t chisquare_gli1 = gli1->chiSquare();
   cout<<"Chi square of gli fit is :"<< chisquare_gli1<< endl;
   
   
   
   TCanvas *cc = new TCanvas("cc", "cc",0,0,600,600);
   TPad *pad11 = new TPad("pad11","pad11",0,0.36,1,1);
   TPad *pad22 = new TPad("pad22","pad22",0,0,1,0.25);
   pad11->SetBottomMargin(0.00001);
   pad11->SetBorderMode(0);
   pad22->SetTopMargin(0.00001);
   pad22->SetBottomMargin(0.1);
   pad22->SetBorderMode(0);
   pad11->Draw();
   pad22->Draw();
   pad11->cd();
   gStyle->SetOptTitle(0);
   cc->SetFillColor(0);
   cc->SetBorderSize(2);
   cc->SetLeftMargin(0.1422222);
   cc->SetRightMargin(0.04444445);
   gli1->SetStats(0);
   gli1->Draw();
   pad22->cd();
   pullgli->SetStats(0);
   pullgli->Draw();
   cc->cd();
   
    //Likelihood
   
   
   Double_t amin,edm,errdef;
   Int_t nvpar,nparx;
   RooRealVar* glr_var = (RooRealVar*)fitRes->floatParsFinal().find("mean11");
   RooRealVar* gli_var = (RooRealVar*)fitRes->floatParsFinal().find("mean22");
   
   
   Float_t glr_mean, gli_mean;
   Float_t glr_err, gli_err;
   
   glr_mean = glr_var->getVal();
   glr_err = glr_var->getAsymErrorHi();
   
   gli_mean = gli_var->getVal();
   gli_err = gli_var->getAsymErrorHi();
   
  
   
   std::cout<<"fitpar : "<<glr_mean<<"\t"<<gli_mean<<"\n";
   std::cout<<"fitpar_err : "<<glr_err<<"\t"<<gli_err<<"\n";
   
   glr_err /= 4.;
   gli_err /= 4.;
  
   
   
   
   RooAbsReal* nll = fitmodel->createNLL(*data,NumCPU(8));
   nll->SetName("nll");
   RooMinuit minuit(*nll);
   minuit.migrad();
   
   
      //result = minuit.save();
   const double nllMin = nll->getVal();
   TString tmp;
   tmp += nllMin;
   tmp += " )";
   RooFormulaVar * likl = new RooFormulaVar("likl", "L", "exp(-nll  + "+tmp, RooArgList(*nll));
   TH2D * manL = new TH2D("manL", "likelihood contour",40,glr_mean-(20*glr_err),glr_mean+(20*glr_err),40,gli_mean-(20*gli_err),gli_mean+(20*gli_err));
   manL->SetXTitle(mean11->GetTitle());
   manL->SetYTitle(mean22->GetTitle());
   Float_t glr_m;
   Float_t gli_m;
   glr_m = glr_mean-(20*glr_err);
   for(int i = 1; i <= manL->GetNbinsX();++i,glr_m+=glr_err)
      {
      gli_m = gli_mean-(20*gli_err);
      for(int j = 1; j <= manL->GetNbinsY();++j,gli_m+=gli_err)
         {
         
         double x = manL->GetXaxis()->GetBinCenter(i);
         mean11->setVal(x);
         double y = manL->GetYaxis()->GetBinCenter(j);
         mean22->setVal(y);
         
         TVirtualFitter *fitter = TVirtualFitter::Fitter(glr_h);
         fitter->GetStats(amin,edm,errdef,nvpar,nparx);
         manL->SetBinContent(i, j, amin-likl->getVal());
         
         }
      }
   TCanvas *con = new TCanvas("con", "con", 0, 0, 800,600);
   manL->Draw("cont4z");
   //con->SaveAs("2Dcontour.root","root");
   
   RooMinuit m(*nll) ;
   m.migrad() ;
   m.hesse() ;
   TCanvas *con1 = new TCanvas("con1", "con1", 0, 0, 800,600);
   RooPlot* p1 = m.contour(*mean11,*mean22,1,2,3) ;
   p1->Draw("CONT4Z") ;
   
   
   
   
   
   
 /*
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
  c->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;
  c->cd(4) ; gPad->SetLeftMargin(0.15) ; frame4->GetYaxis()->SetTitleOffset(1.4) ; frame4->Draw() ;
*/

}
