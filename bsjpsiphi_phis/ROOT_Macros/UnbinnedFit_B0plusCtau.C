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


void UnbinnedFit_B0plusCtau(string filename){


  RooWorkspace *ws = new RooWorkspace("ws");

  TFile* f=new TFile(filename.c_str());
  TTree* tr=(TTree*) f->Get("OutTreeB0_PlusCuts"); //tree_b

  //RooRealVar deltaM("deltaM","deltaM",0,0.6) ;
  //RooDataSet dataEx("dataEx","dataset with deltaM",tr,deltaM) ; 

  RooRealVar Lxy("Lxy","Lxy",0.03,0.3);
  RooRealVar B0_MassFromSV("B0_MassFromSV","B0_MassFromSV",5.15,5.45) ;
  RooDataSet dataB("dataB","dataset with B0",tr,RooArgSet(B0_MassFromSV,Lxy)) ; 
  //RooDataSet dataB("dataBct","dataset with B0ct",tr,Lxy) ; 


  //TCanvas *c1 = new TCanvas("c","MassFit",1200,800);
  //  TCanvas *c1= new TCanvas;
  TCanvas *c1 = new TCanvas("c1", "c1",382,157,1552/*776*/,605);
  c1->Range(5.201544,-2000.826,5.533448,7073.785);
  c1->SetLeftMargin(0.2459948);
  c1->SetBottomMargin(0.8204861);
  c1->Divide(2,1);

 ws->import(dataB);
 //ws->factory("Exponential::expFunct(B0_MassFromSV,coefExp[-4.,-10.,10])");

 ws->factory("Chebychev::expFunct(B0_MassFromSV,{CoefPol0[-0.435,-0.6,-0.2],CoefPol1[-.34,-0.5,0.1]})"); //bkg inv mass

 //ws->factory("Exponential::expFunct2(B0_MassFromSV,coefExp1[-40.,-50.,10])");
 //ws->factory("SUM::expFunct(a[0.5,0,1]*expFunct1,b[0.1,0,1]*expFunct2)");

 // ws->factory("Chebychev::expFunct(B0_MassFromSV,{CoefPol0[0.,-1.,1.],CoefPol1[0.,-1.,1.]})");

 ws->factory("Gaussian::signalG1(B0_MassFromSV,meanSig[5.279,5.26,5.30],sigmaSig[0.01,0.008,0.02])");// signal Inv Mass
 ws->factory("Gaussian::signalG2(B0_MassFromSV,meanSig,sigmaSig2[0.027,0.02,.05])");

 ws->factory("SUM::signalG(signalG1,Gfrac[0.5,0.,1.]*signalG2)");

 ws->factory("Exponential::ctau(Lxy,coefExp[-25.,-50.,-15])"); //signal ctau

 ws->factory("Exponential::ctaubkg1(Lxy,coefExpctau1[-30.,-50.,-10])"); // bkg ctau as sum of two exponentials
 ws->factory("Exponential::ctaubkg2(Lxy,coefExpctau2[-10.,-50.,-5])");

 ws->factory("SUM::ctaubkg(Frac[.5,0,1]*ctaubkg1,ctaubkg2)");

 //tot signal pdf = product of signal PDFs
 ws->factory("PROD::totPDFSig(signalG, ctau)");

 //tot PDF for bkg = product of Bkg
 ws->factory("PROD::totPDFBkg(expFunct,ctaubkg)");

 ws->factory("SUM::totPDF(NSig[100000.,10000.,2000000.]*totPDFSig,NBkg[1000.,100.,1000000.]*totPDFBkg)");


 //RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.,1.) ;
 //RooAddPdf signalG("signalG", "double Gaussian", RooArgList(ws->pdf("signalG1"),ws->pdf("signalG2")), sig1frac);


 ws->pdf("totPDF")->fitTo(dataB,Extended(1),Save(1),Minos(0),SumW2Error(kTRUE));

 RooPlot* Mframe = B0_MassFromSV.frame() ;
 Mframe->GetXaxis()->SetTitle("#mu^{+} #mu^{-} K^{#pm} #pi^{#pm} Inv. Mass. (GeV/c^{2})");
 Mframe->SetTitle("");

 dataB.plotOn(Mframe,Binning(120)) ;
 ws->pdf("totPDF")->plotOn(Mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
 ws->pdf("totPDF")->plotOn(Mframe,Components("expFunct"),LineColor(kBlue),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected));
 ws->pdf("totPDF")->plotOn(Mframe,Components("signalG1"),LineColor(kRed),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected));
 ws->pdf("totPDF")->plotOn(Mframe,Components("signalG2"),LineColor(kGreen),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected));

 c1->cd(1);
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

  //  double sigma = 1000*ws->var("sigmaSig")->getVal();
  
  double sigma1=1000*ws->var("sigmaSig")->getVal();
  double sigma2=1000*ws->var("sigmaSig2")->getVal();

  double esigma1=1000*ws->var("sigmaSig")->getError();
  double esigma2=1000*ws->var("sigmaSig2")->getError();

  //esigma1=esigma1*sigma1/sigma*ws->var("Gfrac")->getVal();  
  //esigma2=esigma2*sigma2/sigma*(1-ws->var("Gfrac")->getVal());  

  double esigma=sqrt(pow(esigma1,2)+pow(esigma2,2));
  
  double m0 = ws->var("meanSig")->getVal();
  double em0 = ws->var("meanSig")->getError();
  
  sprintf(ms,"m_{0} = %6.4f #pm %5.4f GeV/c^{2}", m0,em0);
  t->DrawLatex(0.55,0.81,ms);

  sprintf(ms,"#sigma = %4.1f #pm %3.1f MeV/c^{2}",sigma,esigma);
  t->DrawLatex(0.55,0.77,ms);

  //sprintf(ms,"N_{bkg} = %5.0f #pm %3.0f",Nbkg,eNbkg);
  //t->DrawLatex(0.55,0.73,ms);

  //cout << "AAAAAAAAAAAA" << esigma << endl;

  ws->var("B0_MassFromSV")->setRange("integral",5.25, 5.5);//m0-3*(sigma/1000),m0+3*(sigma/1000));

  RooAbsReal* N= ws->pdf("expFunct")->createIntegral(RooArgSet(*ws->var("B0_MassFromSV")),NormSet(RooArgSet(*ws->var("B0_MassFromSV"))),Range("integral"));
  RooAbsReal* S= ws->pdf("signalG")->createIntegral(RooArgSet(*ws->var("B0_MassFromSV")),NormSet(RooArgSet(*ws->var("B0_MassFromSV"))),Range("integral"));
  
  cout << "B= " << N->getVal() << " S= " << S->getVal() << endl;

  double SN=(Nsig*S->getVal())/sqrt((Nbkg*N->getVal())+Nsig*S->getVal());

  sprintf(SovB,"S/#sqrt{S+B}= %5.2f",SN);
 //t->DrawLatex(0.55,0.73,SovB);
  
  string on=filename.substr(9,6);

  RooPlot* Ctauframe = Lxy.frame() ;
  Ctauframe->GetXaxis()->SetTitle("ct [cm]");
  Ctauframe->SetTitle("");

  dataB.plotOn(Ctauframe,Binning(100)) ;
  ws->pdf("totPDF")->plotOn(Ctauframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  ws->pdf("totPDF")->plotOn(Ctauframe,Components("ctau"),LineColor(kBlue),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected));
  ws->pdf("totPDF")->plotOn(Ctauframe,Components("ctaubkg"),LineColor(kRed),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected));

  c1->cd(2);
  
  Ctauframe->Draw() ;
  //Ctauframe->GetYaxis()->SetMinimum(0.5);
  //Ctauframe->SetLogy(1);
  gPad->SetLogy(1);

  double ctau=ws->var("coefExp")->getVal();
  double ctauerr=ws->var("coefExp")->getError();

  cout << ctauerr << endl;

  sprintf(ms,"raw c#tau=%4.4f #pm %4.4f" ,-1/ctau, 1/(ctau*ctau)*ctauerr);
  t->DrawLatex(0.55,0.81,ms);

  string oname="B0Fit" + on + ".png";

  c1->SaveAs(oname.c_str());
  
 
 //TCanvas c1;

 /* RooPlot* deltaMframe = deltaM.frame() ;
 dataEx.plotOn(deltaMframe,Binning(20)) ;
 deltaMframe->Draw(); 
 c1->SaveAs("DeltaMPlot.pdf");
 */
}

