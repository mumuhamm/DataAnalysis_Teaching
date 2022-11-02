#pragma once

#include "TStyle.h"

#include "RooKeysPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "Roo1DTable.h"
#include "RooSimultaneous.h"

#include "RooTagPdf.cxx"
#include "RooDelta.cxx"

#include "utils.h"
#include "definitions.C"
#include "bkgSigSplitWSB.C"

TDirectory* Bs_tagMistagPdf_fit(RooDataSet& data, RooRealVar& mistag, RooCategory& tag, TString suffix) {
    using namespace RooFit;
    data.table(tag)->Print("V");
    auto p = data.sumEntries("tag == tag::Bs");
    auto n = data.sumEntries("tag == tag::Bsbar");
    auto z = data.sumEntries("tag == tag::untag");
    auto t = data.sumEntries();
    clog << "Asym: " << (p-n)/(p+n) << endl;
    clog << "Eff: " << 1-z/t << endl;
  
    auto dsNoUntag = (RooDataSet*)data.reduce(Cut("tag != tag::untag"));

    auto eff = new RooRealVar("eff", "Eff", 0.5, 0, 1);
    auto asym = new RooRealVar("asym", "Asym", 0, -1, 1);
    auto tagPdf = new RooTagPdf("tagPdf", "Tag Pdf", tag, *asym, *eff);
    
    std::clog << "Starting creation of the mistag pdf" << endl;
    auto mistTagPdf = new RooKeysPdf("mistTagPdf", "Tagged Pdf", mistag, *dsNoUntag);
    auto mistUntagPdf = new RooDelta("mistUntagPdf", "Untagged Pdf", mistag, 0.5);
    auto mistagPdf = new RooSimultaneous("mistagPdf", "Mistag Pdf", {{"Bs", mistTagPdf}, {"Bsbar", mistTagPdf}, {"untag", mistUntagPdf}}, tag);

    auto fullPdf = new RooProdPdf("pdf", "", *tagPdf, Conditional(*mistagPdf, mistag));

    tagPdf->fitTo(data, Save(true))->Print("V");
    
    //============================Draw=====================================================
    gStyle->SetOptTitle(false);
    
    auto mistagFrame = mistag.frame();
    dsNoUntag->plotOn(mistagFrame, DataError(RooAbsData::SumW2), Name("data"));
    mistTagPdf->plotOn(mistagFrame, Name("fit"), Precision(1e-5));
    
    TCanvas c;
    c.SetFillColor(0);
    c.SetBorderSize(2);
    roo_pulls(c, mistagFrame, "data", "fit");
    
    saveCanvas(c, "mistagTaggedPdf_" + suffix);
    
    auto mistagFullFrame = mistag.frame();
    data.plotOn(mistagFullFrame, DataError(RooAbsData::SumW2), Name("data"), Binning(5000));
    fullPdf->plotOn(mistagFullFrame, Name("fit"), Precision(1e-5));
    mistagFullFrame->SetMinimum(0.1);
    
    TCanvas cFull;
    cFull.SetLogy();
    cFull.SetFillColor(0);
    cFull.SetBorderSize(2);
    roo_pulls(cFull, mistagFullFrame, "data", "fit");
    
    saveCanvas(cFull, "mistagFullPdf_" + suffix);
    
    //============================Save===================================================
    auto prevDir = gDirectory;
    TString dirPath = "outputs/";
    if (ensureDir(dirPath)) return nullptr;
    auto mistTagDir = TFile::Open(dirPath + "tagMistagPdfDir_" + suffix + ".root", "RECREATE");
    
    RooWorkspace mistWs("mistWs");
    mistWs.import(*mistUntagPdf);
    mistWs.import(*mistTagPdf);
    mistWs.import(*tagPdf);
    
    mistWs.Write();
    
    mistTagDir->ReOpen("READ");
    
    prevDir->cd();
    
    return mistTagDir;
}

TDirectory* Bs_tagMistagPdf_fit(TString filename, TString portion = "full", TString massModelFileName = nullptr, TString suffix = "") {
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
    
    RooDataSet* usedData = &recoData;
    if (portion != "full") {
        if (not massModelFileName) {
            cerr << "Requested sideband subtraction without providing a mass model" << endl;
            return nullptr;
        }
        auto massModelFile = openFileForReading(massModelFileName);
        auto massWs = (RooWorkspace*)massModelFile->Get("mass");
        
        auto splitData = bkgSigSplitWSB(recoData, *massWs, *sbStatus, true);
    
        usedData = splitData[portion.Data()];
    
        suffix = portion + "_" + suffix;
    }
    
    return Bs_tagMistagPdf_fit(*usedData, *mistag, *tag, suffix);
}
