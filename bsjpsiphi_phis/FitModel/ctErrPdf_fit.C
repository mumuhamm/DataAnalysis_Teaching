// Authors: Md Alibordi, Alberto Bragagnolo, Enrico Lusiani
#pragma once

#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGamma.h"
#include "RooPlot.h"
#include "RooRealVar.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TPad.h"
#include "TStyle.h"
#include "TTree.h"
#include "TLatex.h"

#include "utils.h"

#include "definitions.C"
#include "bkgSigSplitWSB.C"

TDirectory* ctErrPdf_fit(std::map<string, RooDataSet*> splitData, RooRealVar& BsCt2DMCErr, bool isData, TString suffix = "") {
    using namespace RooFit;
    using namespace std;

    BsCt2DMCErr.setRange("bgFitRange", 0.0008, 0.005);
    
    auto& data = *splitData["full"];

    vector<TString> labels;
    labels.push_back("sig");
    if (isData) {
        labels.push_back("bkg");
    }
    
    RooArgList pdfList;
    RooArgList coefList;
    for (TString label: labels) {
        // gamma shape parameter
        auto shape1 = new RooRealVar("ctErr_gamma1_" + label, "ct error pdf, first gamma shape", 8.0440, 0, 16);
        auto shape2 = new RooRealVar("ctErr_gamma2_" + label, "ct error pdf, second gamma shape", 3.798, 0, 8);
        // gamma scale parameter
        auto scale1 = new RooRealVar("ctErr_beta1_" + label, "ct error pdf, first gamma scale", 0.00014473, 0, 0.0003);
        auto scale2 = new RooRealVar("ctErr_beta2_" + label, "ct error pdf, second gamma scale", 0.000541, 0, 0.001);
        // gamma displacement parameter
        auto displ1 = new RooRealVar("ctErr_mu1_" + label, "ct error pdf, first gamma displacement", 0.0000, 0, 0.0002);
        auto displ2 = new RooRealVar("ctErr_mu2_" + label, "ct error pdf, second gamma displacement", 0.0000, 0, 0.0002);

        auto gamma1 = new RooGamma("ctErr_gammaPdf1_" + label, "ct error pdf, first gamma", BsCt2DMCErr, *shape1, *scale1, *displ1);
        auto gamma2 = new RooGamma("ctErr_gammaPdf2_" + label, "ct error pdf, second gamma", BsCt2DMCErr, *shape2, *scale2, *displ2);

        auto frac_uncert = new RooRealVar("cterr_gamma1frac_" + label, "frac", 0.5, 0, 1);

        auto ctErrPdf = new RooAddPdf("ctErr_pdf_" + label, "ct error pdf", RooArgList(*gamma1, *gamma2), *frac_uncert);
    
        const char* fitRange = nullptr;
        if (label == "bkg") {
            fitRange = "bgFitRange";
        }
        cerr << "Now fitting the " << label << " sample" << endl;
        
        RooFitResult *restrictedFitRes = ctErrPdf->fitTo(*splitData[label.Data()], Save(), SumW2Error(false), Range(fitRange));
        restrictedFitRes->Print("V");
        
        auto color0 = (label == "bkg")?kRed:kGreen;
        auto color1 = (label == "bkg")?kOrange+1:kCyan;
        auto color2 = (label == "bkg")?kOrange+7:kTeal+1;
        
        auto frame = BsCt2DMCErr.frame();
        splitData[label.Data()]->plotOn(frame, DataError(RooAbsData::SumW2), Name("data"));
        ctErrPdf->plotOn(frame, Components("ctErr_gammaPdf1_*"), LineColor(color1), LineStyle(2), Name("fitP1"));
        ctErrPdf->plotOn(frame, Components("ctErr_gammaPdf2_*"), LineColor(color2), LineStyle(2), Name("fitP2"));
        ctErrPdf->plotOn(frame, LineColor(color0), Name("fit"), Range("allOfIt"));
        TCanvas c;
        roo_pulls(c, frame, "data", "fit");
        saveCanvas(c, "ctErrPrefit_" + label + "_" + suffix);
    
        displ1->setConstant();
        displ2->setConstant();
        
        shape1->setConstant();
        shape2->setConstant();
       
        pdfList.add(*ctErrPdf);
        coefList.add(*new RooRealVar("ctErr_N" + label, "Number of " + label + "events", data.sumEntries()/labels.size(), 0, data.sumEntries()*2));
    }
    
    RooAddPdf ctErrPdf("ctErr_pdf", "ct error pdf", pdfList, coefList);
    
//    for (auto varBase: *ctErrPdf.getParameters(data)) {
//        auto var = (RooRealVar*)varBase;
//        var->setConstant();
//    }
    
    auto prevDir = gDirectory;
    TString dirPath = "outputs/";
    if (ensureDir(dirPath)) return nullptr;
    auto ctErrPdfDir = TFile::Open(dirPath + "ctErrPdfDir_" + suffix + ".root", "RECREATE");
    
    RooWorkspace ws("ctErr");
    ws.import(ctErrPdf);
    
    ws.Write();
    
    ctErrPdfDir->ReOpen("READ");
    
    prevDir->cd();
    return ctErrPdfDir;
}

TDirectory* ctErrPdf_fit(TString filename, TString massModelFileName, bool isData, TString suffix = "") {
    auto recoTree = getTree(filename, recoTreeName);
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
    
    auto sbStatus = mainWs->catfunc(sbStatusName);
    
    RooDataSet recoData("recoData", "raw data1", RooArgSet(*BscosthetaMC, *BscospsiMC, *BsphiMC, *BsCt2DMC, *BsCt2DMCErr, *svmass, *tag, *mistag),
                                      RooFit::Import(*recoTree));
    
    auto massModelFile = openFileForReading(massModelFileName);
    auto massWs = (RooWorkspace*)massModelFile->Get("mass");
    
    auto splitData = bkgSigSplitWSB(recoData, *massWs, *sbStatus, isData);
    
    return ctErrPdf_fit(splitData, *BsCt2DMCErr, isData, suffix);
    //return ctErrBgPdf_fit(*getSideband(recoData), *BsCt2DMCErr, suffix);
}
