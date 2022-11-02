// Author: Enrico Lusiani

#include "definitions.C"
#include "utils.h"

TDirectory* ctBg_fit(RooDataSet& bgData, TString suffix = "") {

    using namespace RooFit;
    
    gStyle->SetOptTitle(false);
    
    auto& ct = *(RooRealVar*)bgData.get()->find(ctName);
    auto& ctErr = *(RooRealVar*)bgData.get()->find(ctErrName);
    
    RooRealVar sigmaScale("ctBg_sigmaScale", "ctBg_sigmaScale", 1);
    //RooGaussModel rm("ctBg_gaussModel", "ctBg_gaussModel", ct, RooConst(0), ctErr, RooConst(0), sigmaScale);
    RooTruthModel rm("ctBg_truthModel", "ctBg_truthModel", ct);
    RooRealVar ctau1("ctBg_ctau1", "c#tau_{1}", 0.004, 0, 0.5);
    RooRealVar ctau2("ctBg_ctau2", "c#tau_{2}", 0.04, 0, 0.5);
    RooDecay d1("ctBg_decay1", "ctBg_decay1", ct, ctau1, rm, RooDecay::SingleSided);
    RooDecay d2("ctBg_decay2", "ctBg_decay2", ct, ctau2, rm, RooDecay::SingleSided);
    
    //RooRealVar frac0("ctBG_frac0", "ctBG_frac0", 0, 1);
    RooRealVar frac1("ctBG_frac1", "model fraction", 0.4, 0, 1);
    RooAddPdf ctBgPdf("ctBg_pdf", "ctBg_pdf", RooArgList(d1, d2), RooArgList(frac1), true);
    
    // conditional while waiting for cterr model
    auto fitRes = ctBgPdf.fitTo(bgData, Save(true), ConditionalObservables(ctErr));
    fitRes->Print("V");
    
    TCanvas c;
    c.SetLogy();
    auto frame = ct.frame();
    bgData.plotOn(frame, Name("data"));
    ctBgPdf.plotOn(frame, Name("fitP1"), Components("ctBg_decay1"), LineColor(kOrange+7), LineStyle(2));
    ctBgPdf.plotOn(frame, Name("fitP2"), Components("ctBg_decay2"), LineColor(kOrange+1), LineStyle(2));
    ctBgPdf.plotOn(frame, Name("fit"), LineColor(kRed));
    ctBgPdf.paramOn(frame, Layout(0.6, 0.99, 0.95));
    roo_pulls(c, frame, "data", "fit", 3);
    saveCanvas(c, "ctBg" + suffix);
    
    auto prevDir = gDirectory;
    TString dirPath = "outputs/";
    if (ensureDir(dirPath)) return nullptr;
    auto file = TFile::Open(dirPath + "ctBgPdfDir" + suffix + ".root", "RECREATE");
    
    RooWorkspace ctBgWs("ctBgWs");
    ctBgWs.import(ctBgPdf, Silence());
    auto constrPdf = fitRes->createHessePdf(*ctBgPdf.getParameters(bgData));
    constrPdf->SetName("ctBg_constr");
    ctBgWs.import(*constrPdf);
    ctBgWs.Write();
    
    file->ReOpen("READ");
    prevDir->cd();

    return file;
}

TDirectory* ctBg_fit(vector<const char*> recoFileNames, TString suffix = "") {
    
    auto recoTree = getTree(recoFileNames, recoTreeName);
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
    
    auto bgData = (RooDataSet*)getRawSideband(recoData);
    cerr << bgData->sumEntries() << endl;
    
    return ctBg_fit(*bgData, suffix);
}
