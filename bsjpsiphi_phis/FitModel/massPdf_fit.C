// Authors: Md Alibordi, Alberto Bragagnolo, Enrico Lusiani
#pragma once

#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooJohnsonLocal.cxx"

#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TH2.h"

#include "utils.h"


#include "definitions.C"

TDirectory* massPdf_fit(RooDataSet& data, RooRealVar& svmass, bool isData) {
    using namespace RooFit;
    using namespace std;


    RooRealVar mu("mass_mu", "mass_mu", 5.36679, 5.35, 5.37);
    RooRealVar lambda("mass_lambda", "mass_lambda", 0.01, 0, 1);
    RooRealVar gamma("mass_gamma", "mass_gamma", 0., -1, 1);
    RooRealVar delta("mass_delta", "mass_delta", 1., 0, 10);
    
    RooRealVar bkgSlope("mass_bkgSlope", "bkg slope", -100, 100);
    //RooRealVar bkgC1("mass_bkgC1", "", -10, 10);
    //RooRealVar bkgC2("mass_bkgC2", "", -10, 10);

    //===============================================================================================
    RooJohnsonLocal sigPdf("mass_sig", "mass sig", svmass, mu, lambda, gamma, delta);
    RooExponential bkgPdf("mass_bkg", "mass bkg", svmass, bkgSlope);
    //RooChebychev bkgPdf("mass_bkg", "mass_bkg", svmass, RooArgList(bkgC1, bkgC2));

    //==============================Construct decay(t) X Gaussian resolution fn============================================================================ 

    RooRealVar nSig("nSig", "Number of Signal Events in SIGNAL MC", 3./4*data.sumEntries(), 0., 1.5*data.sumEntries());
    RooRealVar nBkg("nBkg", "Number of BG events in SIGNAL MC", 1./4*data.sumEntries(), 0., 1.5*data.sumEntries());

    RooArgList pdfList(sigPdf);
    RooArgList coefList(nSig);
    if (isData) {
        pdfList.add(bkgPdf);
        coefList.add(nBkg);
    }
    RooAddPdf massPdf("mass_pdf", "Total pdf", pdfList, coefList);
    
    RooFitResult *fitRes = massPdf.fitTo(data, Save(), Extended(1));
    fitRes->Print("v");
    
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    TH2 *hcorr = fitRes->correlationHist("massCorrHist");
    TCanvas cCorr("Correlation Matrix", "Correlation Matrix", 1000, 600);
    hcorr->GetYaxis()->SetTitleOffset(1.4);
    hcorr->Draw("colz");
    saveCanvas(cCorr, "massCorrMatrix");
    gStyle->SetPalette(kBird);
    
    RooPlot *massFrame = svmass.frame(Title("M_{B_{s}} (GeV/c^{2})"), Bins(100));
    data.plotOn(massFrame, DataError(RooAbsData::SumW2), Name("data"));
    if (isData) {
        massPdf.plotOn(massFrame, LineColor(kRed), LineStyle(kDashed), Name("bkg"), Components("mass_bkg"));
        massPdf.plotOn(massFrame, LineColor(kGreen), LineStyle(kDashed), Name("sig"), Components("mass_sig"));
    }
    massPdf.plotOn(massFrame, LineColor(kBlue), Name("fit"));
    massPdf.paramOn(massFrame, Layout(0.6));
    Double_t chisquare_mass = massFrame->chiSquare("fit", "data", fitRes->floatParsFinal().getSize());
    cout << "Chi square of mass fit is :" << chisquare_mass << endl;
    
    gStyle->SetOptTitle(0);
    TCanvas c;
    c.SetFillColor(0);
    c.SetBorderSize(2);
    roo_pulls(c, massFrame, "data", "fit", fitRes->floatParsFinal().getSize());
    
    saveCanvas(c, "massPlot");
    
    auto prevDir = gDirectory;
    TString dirPath = "outputs/";
    if (ensureDir(dirPath)) return nullptr;
    auto massDir = TFile::Open(dirPath + "massPdfDir.root", "RECREATE");
    
    RooWorkspace ws("mass");
    ws.import(massPdf);
    
    ws.Write();
    
    massDir->ReOpen("READ");
    
    prevDir->cd();
    return massDir;
}

TDirectory* massPdf_fit(vector<const char*> filenames, bool isData) {
    auto recoTree = getTree(filenames, recoTreeName);
    if (!recoTree) {
        return nullptr;
    }
    
    auto mainWs = getMainWs();
    
    auto BsCt2DMC = mainWs->var(ctName);
    auto BsCt2DMCErr = mainWs->var(ctErrName);
    auto BscosthetaMC = mainWs->var(cosThetaName);
    auto BscospsiMC = mainWs->var(cosPsiName);
    auto BsphiMC = mainWs->var(phiName);
    auto svmass = mainWs->var(massName);
    auto mistag = mainWs->var(mistagName);
    auto tag = mainWs->cat(tagCatName);
    
    RooDataSet recoData("recoData", "raw data1", RooArgSet(*BscosthetaMC, *BscospsiMC, *BsphiMC, *BsCt2DMC, *BsCt2DMCErr, *svmass, *tag, *mistag),
                                      RooFit::Import(*recoTree));
    
    return massPdf_fit(recoData, *svmass, isData);
}
