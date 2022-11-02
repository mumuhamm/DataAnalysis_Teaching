#include "RooFit.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooCmdArg.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooCategory.h"
#include "RooWorkspace.h"
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TClonesArray.h"
#include "TChain.h"

#include "RooFit.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"

#include <TCanvas.h>
#include "TH1F.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include <TVector3.h>
#include <math.h>

using namespace RooFit ;


void UnbinnedFit_BsOnly(string filename){


  RooWorkspace *ws = new RooWorkspace("ws");

  TFile* f=new TFile(filename.c_str());
  TTree* tr=(TTree*) f->Get("OutTree_PlusCuts"); //tree_b

  //RooRealVar deltaM("deltaM","deltaM",0,0.6) ;
  //RooDataSet dataEx("dataEx","dataset with deltaM",tr,deltaM) ; 

  RooRealVar Bs_MassFromSV("Bs_MassFromSV","BsMassFromSV",5.2,5.65) ;
  RooDataSet dataB("dataB","dataset with Bs",tr,Bs_MassFromSV) ; 

  //TCanvas *c1 = new TCanvas("c","MassFit",1200,800);
  //  TCanvas *c1= new TCanvas;
  TCanvas *c1 = new TCanvas("c1", "c1",382,157,776,605);
  c1->Range(5.201544,-2000.826,5.533448,7073.785);
  c1->SetLeftMargin(0.1459948);
  c1->SetBottomMargin(0.2204861);

 ws->import(dataB);
 ws->factory("Exponential::expFunct(Bs_MassFromSV,coefExp[-4.,-10.,10])");

 //ws->factory("Chebychev::expFunct(Bs_MassFromSV,{CoefPol0[-0.05,-1.,1.],CoefPol1[-.4,-1.,1.]})");

 //ws->factory("Exponential::expFunct2(Bs_MassFromSV,coefExp1[-40.,-50.,10])");
 //ws->factory("SUM::expFunct(a[0.5,0,1]*expFunct1,b[0.1,0,1]*expFunct2)");

 // ws->factory("Chebychev::expFunct(Bs_MassFromSV,{CoefPol0[0.,-1.,1.],CoefPol1[0.,-1.,1.]})");

 ws->factory("Gaussian::signalG1(Bs_MassFromSV,meanSig[5.36,5.3,5.5],sigmaSig[0.02,0.01,.03])");
 ws->factory("Gaussian::signalG2(Bs_MassFromSV,meanSig,sigmaSig2[0.005,0.003,.02])");

 ws->factory("SUM::signalG(signalG1,Gfrac[0.5,0.,1.]*signalG2)");

 //RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.,1.) ;
 //RooAddPdf signalG("signalG", "double Gaussian", RooArgList(ws->pdf("signalG1"),ws->pdf("signalG2")), sig1frac);

 ws->factory("SUM::totPDF(NSig[1000000.,50000.,5000000.]*signalG,NBkg[160000.,10000.,5000000.]*expFunct)");

 ws->pdf("totPDF")->fitTo(dataB,Extended(1),Save(1),Minos(0),SumW2Error(kTRUE));

 RooPlot* Mframe = Bs_MassFromSV.frame() ;
 Mframe->GetXaxis()->SetTitle("#mu^{+} #mu^{-} K^{+} K^{-} Inv. Mass. (GeV/c^{2})");
 Mframe->SetTitle("");

 dataB.plotOn(Mframe,Binning(50)) ;
 ws->pdf("totPDF")->plotOn(Mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
 ws->pdf("totPDF")->plotOn(Mframe,Components("expFunct"),LineColor(kBlue),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected));
 ws->pdf("totPDF")->plotOn(Mframe,Components("signalG1"),LineColor(kRed),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected));
 ws->pdf("totPDF")->plotOn(Mframe,Components("signalG2"),LineColor(kGreen),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected));

 Mframe->Draw() ;

 //write results on plot

 char sig[80], ms[80],SovB[80];
 double errsig = ws->var("NSig")->getError();
 double Nsig   = ws->var("NSig")->getVal();

 sprintf(sig,"N_{Sig}= %5.0f #pm %3.0f", Nsig,errsig);

 TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(13);
  t->SetTextFont(63);
  t->SetTextSizePixels(18);
  t->DrawLatex(0.55,0.85,sig);

  double Nbkg = ws->var("NBkg")->getVal();
  double eNbkg = ws->var("NBkg")->getError();  

  double sigma = 1000*(sqrt(ws->var("Gfrac")->getVal()*ws->var("sigmaSig")->getVal()*ws->var("sigmaSig")->getVal()+(1-ws->var("Gfrac")->getVal())*ws->var("sigmaSig2")->getVal()*ws->var("sigmaSig2")->getVal()));
  
  double sigma1=1000*ws->var("sigmaSig")->getVal();
  double sigma2=1000*ws->var("sigmaSig2")->getVal();

  double esigma1=1000*ws->var("sigmaSig")->getError();
  double esigma2=1000*ws->var("sigmaSig2")->getError();

  esigma1=esigma1*sigma1/sigma*ws->var("Gfrac")->getVal();  
  esigma2=esigma2*sigma2/sigma*(1-ws->var("Gfrac")->getVal());  

  double esigma=sqrt(pow(esigma1,2)+pow(esigma2,2));

  //double esigma = 1000*(sqrt(ws->var("Gfrac")->getVal()*ws->var("sigmaSig")->getError()*ws->var("sigmaSig")->getError()+(1-ws->var("Gfrac")->getVal())*ws->var("sigmaSig2")->getError()*ws->var("sigmaSig2")->getError()));
  
  double m0 = ws->var("meanSig")->getVal();
  double em0 = ws->var("meanSig")->getError();
  
  sprintf(ms,"m_{0} = %6.4f #pm %5.4f GeV/c^{2}", m0,em0);
  t->DrawLatex(0.55,0.81,ms);

  sprintf(ms,"#sigma = %4.1f #pm %3.1f MeV/c^{2}",sigma,esigma);
  t->DrawLatex(0.55,0.77,ms);

  //sprintf(ms,"N_{bkg} = %5.0f #pm %3.0f",Nbkg,eNbkg);
  //t->DrawLatex(0.55,0.73,ms);

  //cout << "AAAAAAAAAAAA" << esigma << endl;

  ws->var("Bs_MassFromSV")->setRange("integral",5.2, 5.65);//m0-3*(sigma/1000),m0+3*(sigma/1000));

  RooAbsReal* N= ws->pdf("expFunct")->createIntegral(RooArgSet(*ws->var("Bs_MassFromSV")),NormSet(RooArgSet(*ws->var("Bs_MassFromSV"))),Range("integral"));
  RooAbsReal* S= ws->pdf("signalG")->createIntegral(RooArgSet(*ws->var("Bs_MassFromSV")),NormSet(RooArgSet(*ws->var("Bs_MassFromSV"))),Range("integral"));
  
  cout << "B= " << Nbkg*N->getVal() << " S= " << Nsig*S->getVal() << endl;

  ws->var("Bs_MassFromSV")->setRange("sinistra",5.2, 5.28);
  ws->var("Bs_MassFromSV")->setRange("destra",5.45, 5.65);
  ws->var("Bs_MassFromSV")->setRange("centro",5.28, 5.45);

  RooAbsReal* Nsin= ws->pdf("expFunct")->createIntegral(RooArgSet(*ws->var("Bs_MassFromSV")),NormSet(RooArgSet(*ws->var("Bs_MassFromSV"))),Range("sinistra"));
  RooAbsReal* Ndes= ws->pdf("expFunct")->createIntegral(RooArgSet(*ws->var("Bs_MassFromSV")),NormSet(RooArgSet(*ws->var("Bs_MassFromSV"))),Range("destra"));
  RooAbsReal* Ncent= ws->pdf("expFunct")->createIntegral(RooArgSet(*ws->var("Bs_MassFromSV")),NormSet(RooArgSet(*ws->var("Bs_MassFromSV"))),Range("centro"));


  double CountSin= Nbkg*Nsin->getVal();
  double CountDes= Nbkg*Ndes->getVal();
  double CountCent= Nbkg*Ncent->getVal();

  double fracSB= Nbkg*N->getVal()/(CountSin+CountDes);
  double fracSB1= CountCent/(CountSin+CountDes);

  std::cout << "NBkgTot/NBkgSB " << fracSB << "  Nunder/Nsb= " << fracSB1 << std::endl;

  double SN=(Nsig*S->getVal())/sqrt((Nbkg*N->getVal())+Nsig*S->getVal());

  sprintf(SovB,"S/#sqrt{S+B}= %5.2f",SN);
 //t->DrawLatex(0.55,0.73,SovB);
  
  string on=filename.substr(9,6);

  string oname="BsFit" + on + ".png";

  c1->SaveAs(oname.c_str());
  //  c1->SaveAs("BsFit.png");

 //TCanvas c1;

 /* RooPlot* deltaMframe = deltaM.frame() ;
 dataEx.plotOn(deltaMframe,Binning(20)) ;
 deltaMframe->Draw(); 
 c1->SaveAs("DeltaMPlot.pdf");
 */
}

