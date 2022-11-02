//Author: Md. Alibordi
#pragma once

#include "TString.h"
#include "TAxis.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"

#include "utils.h"

#include "definitions.C"

void Pwave_Model_plot(RooDataSet& data, RooRealVar& var, RooAbsPdf& model, RooCategory* dsCat, int npars, const char* plotName, 
    const RooArgSet& showParam = RooArgSet(), bool logScale = false) {
    
    using namespace RooFit;
    
    RooPlot *frame = var.frame(Bins(100));
    
    data.plotOn(frame, DataError(RooAbsData::SumW2), Name("data"));
    if (dsCat) {
        model.plotOn(frame, Name("sig"), Components("Angular_Model_*"), LineColor(kGreen), LineStyle(2), ProjWData(*dsCat, data));
        model.plotOn(frame, Name("bkg"), Components("bkg_model_*"), LineColor(kRed), LineStyle(2), ProjWData(*dsCat, data));
        model.plotOn(frame, Name("fit"), ProjWData(*dsCat, data));
    }
    else {
        model.plotOn(frame, Name("sig"), Components("Angular_Model_*"), LineColor(kGreen), LineStyle(2));
        model.plotOn(frame, Name("bkg"), Components("bkg_model_*"), LineColor(kRed), LineStyle(2));
        model.plotOn(frame, Name("fit"));
    }
    
    if (showParam.getSize() > 0) model.paramOn(frame, RooFit::Parameters(showParam));
    
    Double_t chisquare = frame->chiSquare(npars);
    
    TCanvas c;
    if (logScale) {
        frame->SetMinimum(0.01);
        c.SetLogy();
    }
    roo_pulls(c, frame, "data", "fit", npars);
    
    saveCanvas(c, plotName);
    cerr << "Chi2:" << chisquare << "\n";
}

/*
void Pwave_Model_plot(const char* dataFileName, const char* modelFileName, const char* varName, const char* plotName, bool plotParam = false) {
    auto recoTree = getTree(dataFileName, recoTreeName);
    if (!recoTree) {
        return;
    }
    
    auto prevDir = gDirectory;
    auto wsFile = TFile::Open(modelFileName);
    prevDir->cd();
    
    auto ws = (RooWorkspace*)wsFile->Get("main");
    
    auto var = ws->var(varName);
    if (not var){
        cerr << "Variable not found in the dataset" << endl;
        return;
    }
    
    auto BsCt2DMC = ws->var(ctName);
    auto BsCt2DMCErr = ws->var(ctErrName);
    auto BscosthetaMC = ws->var(cosThetaName);
    auto BscospsiMC = ws->var(cosPsiName);
    auto BsphiMC = ws->var(phiName);
    auto svmass = ws->var(massName);
    auto mistag = ws->var(mistagName);
    auto tag = ws->cat(tagCatName);
    
    RooDataSet data("recoData", "raw data1", RooArgSet(*BscosthetaMC, *BscospsiMC, *BsphiMC, *BsCt2DMC, *BsCt2DMCErr, *svmass, *tag, *mistag),
                                      RooFit::Import(*recoTree));
    
    auto model = ws->pdf("Angular_Model");
    if (not model){
        cerr << "Model not found in the dataset" << endl;
        return;
    }
    
    auto fitRes = (RooFitResult*)wsFile->Get("fitRes");
    if (not fitRes){
        cerr << "fit result not found in file" << endl;
        return;
    }
    
    Pwave_Model_plot(data, *var, *model, fitRes->floatParsFinal().getSize(), plotName, plotParam);
}
*/
