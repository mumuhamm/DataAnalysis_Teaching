// Authors: Md Alibordi
#pragma once

#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"

#include "Math/MinimizerOptions.h"

#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"

#include "utils.h"
#include "definitions.C"


TDirectory* Bs_ctEff_fit(TTree* recoTree, TTree* genTree, RooRealVar& ct, TString suffix = "") {

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 10000 );
    gStyle->SetOptFit(111);
    gStyle->SetOptStat(0);
    
    std::clog << "Computing the ct efficiency function" << std::endl;
    // register the required histograms/functions
    auto prevDir = gDirectory;
    TString dirPath = "outputs/";
    if (ensureDir(dirPath)) return nullptr;
    auto ctEffDir = TFile::Open(dirPath + "ctEffDir_" + suffix + ".root", "RECREATE");
    
    vector<double> bins = {0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025,
                           0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100,
                           0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.275, 0.300, 0.350, 0.400, 0.500};

    Int_t binnum = bins.size() - 1;

    auto hctaucut = new TH1D("hctaucut", "ct;ct(cm);Events", binnum, bins.data());
    auto hctaunocut = new TH1D("hctaunocut", "ct;ct(cm);Events", binnum, bins.data());
    auto resobulk = new TH1D("resobulk", "ct Resolution;ct_{GEN} - ct_{RECO} [cm];Events ", 100, -0.01, 0.01);
    auto pullbulk = new TH1D("pullbulk", "Pullbulk;Pull;Events", 50, -5.0, 5.0);

    // fill histograms
    Long64_t recoEntries = recoTree->GetEntries();
    std::clog << "Start Processing " << recoEntries << " events from reco tree" << std::endl;
    
    // TODO switch to TTreeReader?
    Float_t ctreco, ctgencut, cterr;
    recoTree->SetBranchAddress(ctName, &ctreco);
    recoTree->SetBranchAddress(genCtName, &ctgencut);
    recoTree->SetBranchAddress(ctErrName, &cterr);
    
    for (Int_t i = 0; i < recoEntries; i++) {
        recoTree->GetEntry(i);
        hctaucut->Fill(ctreco);
        resobulk->Fill(ctgencut - ctreco);
        pullbulk->Fill((ctgencut - ctreco) / cterr);
    }
    
    Int_t genEntries = genTree->GetEntries();
    std::clog << "Start Processing " << genEntries << " events from gen tree" << std::endl;
    
    Float_t ctaugen;
    genTree->SetBranchAddress(ctNameGen, &ctaugen);
    
    auto step20 = genEntries / 20;
    for (Int_t i = 0; i < genEntries; i++) {
        if (i % step20 == 0) std::clog << i << "/" << genEntries << endl;
        genTree->GetEntry(i);
        // TODO what if resolution is a function of ct?
        Double_t actualctgen = ctaugen + resobulk->GetRandom();
        hctaunocut->Fill(actualctgen);
    }
    
    // make the ratio histogram
    // TODO possibly a graph is better?
    auto ctEffHist = (TH1*)hctaucut->Clone("ctEffHist");
    ctEffHist->SetYTitle("#epsilon (ct)");
    ctEffHist->Sumw2();
    ctEffHist->Divide(hctaunocut);
    
    // FIXME We shouldn't need this, but I'm unable to fit with pure RooFit, normalization problems likely'
    TF1 *ctform = new TF1("ctEffFn","expo(0)*ROOT::Math::Chebyshev4(x,[2],[3],[4],[5],[6])",0.007,0.5);
    ctform->SetParameter(1,-1);
    ctform->FixParameter(2,1);

    auto fitRes = ctEffHist->Fit(ctform, "MNS");
    fitRes->SetName("fitRes");
    fitRes->Write();
    
    RooRealVar eff("ctEffVar", "", 0, 1);
    RooArgSet varset(ct, eff);
    RooDataSet ctEffDs("ctEffDs", "", varset, RooFit::StoreError(varset));
    for (int i = 1; i <= ctEffHist->GetNbinsX(); i++) {
        ct.setVal(ctEffHist->GetBinCenter(i));
        eff.setVal(ctEffHist->GetBinContent(i));
        
        ct.setError(ctEffHist->GetBinWidth(i)/2);
        eff.setError(ctEffHist->GetBinError(i));
        
        ctEffDs.add(varset);
    }
    

    RooRealVar ctp0("ctp0", "ctp0", ctform->GetParameter(0));
    ctp0.setError(ctform->GetParError(0));
    RooRealVar ctp1("ctp1", "ctp1", ctform->GetParameter(1), "cm^{-1}");
    ctp1.setError(ctform->GetParError(1));
    RooRealVar ctp2("ctp2", "ctp2", ctform->GetParameter(2));
    ctp2.setError(ctform->GetParError(2));
    RooRealVar ctp3("ctp3", "ctp3", ctform->GetParameter(3), "cm^{-1}");
    ctp3.setError(ctform->GetParError(3));
    RooRealVar ctp4("ctp4", "ctp4", ctform->GetParameter(4), "cm^{-2}");
    ctp4.setError(ctform->GetParError(4));
    RooRealVar ctp5("ctp5", "ctp5", ctform->GetParameter(5), "cm^{-3}");
    ctp5.setError(ctform->GetParError(5));
    RooRealVar ctp6("ctp6", "ctp6", ctform->GetParameter(6), "cm^{-4}");
    ctp6.setError(ctform->GetParError(6));
    
    
    RooGenericPdf ctEffFuncAsPdf("ctEffFuncAsPdf", "exp(@0+@1*BsCt2DMC)*ROOT::Math::Chebyshev4(BsCt2DMC,@2,@3,@4,@5,@6)",
                                                 RooArgList(ctp0, ctp1, ctp2, ctp3, ctp4, ctp5, ctp6, ct));
                                                 
    // just for plotting
    // FIXME We shouldn't need this...
    RooFormulaVar ctEffFuncAsFn("ctEffFunc", "exp(@0+@1*BsCt2DMC)*ROOT::Math::Chebyshev4(BsCt2DMC,@2,@3,@4,@5,@6)",
                                                 RooArgList(ctp0, ctp1, ctp2, ctp3, ctp4, ctp5, ctp6, ct));

    //auto res = ctEffFunc.chi2FitTo(ctEffDs, RooFit::Save(true), RooFit::YVar(eff));
    //res->Print("V");
    
    TCanvas cEff;
    auto frame = ct.frame(RooFit::Name("eff"), RooFit::Title(""));
    ctEffDs.plotOnXY(frame, RooFit::YVar(eff), RooFit::Name("data"), RooFit::MarkerStyle(7));
    ctEffFuncAsFn.plotOn(frame, RooFit::Name("fit"), RooFit::Precision(1e-5));
    ctEffFuncAsPdf.paramOn(frame, 
        RooFit::Parameters(RooArgSet(ctp0, ctp1, ctp3, ctp4, ctp5, ctp6)), 
        RooFit::ShowConstants(true), 
        RooFit::Format("NELU", RooFit::AutoPrecision(1)), 
        RooFit::Layout(0.7, 0.99, 0.95)
    );
    frame->SetYTitle("#epsilon (ct)");
    frame->SetTitle("");
    frame->SetMaximum(ctEffHist->GetBinContent(1)*1.1);
    roo_pulls(cEff, frame, "data", "fit", 6);
    
    saveCanvas(cEff, "ctEffPlot_" + suffix);
    
    TCanvas cRes;
    resobulk->SetMarkerStyle(20);
    resobulk->SetMarkerSize(.5);
    TF1 *resofit = new TF1("resofit", "gaus(0)+gaus(3)", -0.01, 0.01);
    gStyle->SetOptStat(1100);
    resobulk->Draw("PL");
    saveCanvas(cRes, "ctResPlot_" + suffix);
    gStyle->SetOptStat(0);
    
    TCanvas cPull;
    pullbulk->Draw();
    saveCanvas(cPull, "ctPullPlot_" + suffix);
    
    RooWorkspace ctEffWs("ctEff");
    ctEffWs.import(ctEffFuncAsFn);
    ctEffWs.import(ctEffFuncAsPdf);
    ctEffWs.Write();
    
    ctEffDir->Write();
    ctEffDir->ReOpen("READ");
    
    prevDir->cd();
    
    return ctEffDir;
}

TDirectory* Bs_ctEff_fit(TString recoFileName, TString genFileName, TString suffix = "") {
    // reco tree
    auto recoTree = getTree(recoFileName, recoTreeName);
    if (!recoTree) {
        return nullptr;
    }
    
    // gen tree
    auto genTree = getTree(genFileName, genTreeName);
    if (!genTree) {
        return nullptr;
    }
    
    auto mainWs = getMainWs();
    
    return Bs_ctEff_fit(recoTree, genTree, *mainWs->var(ctName), suffix);
}
