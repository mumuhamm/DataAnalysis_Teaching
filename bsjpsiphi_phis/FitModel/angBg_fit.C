// Authors: Enrico Lusiani

#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"

#include "RooDataSet.h"
#include "RooLegendre.h"
#include "RooFormulaVar.h"
#include "RooRealSumPdf.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"

#include "definitions.C"
#include "utils.h"

TDirectory* angBg_fit(RooDataSet& bgData, TString suffix  = "") {
    
    gStyle->SetOptStat(0);
    
    auto& costheta = *(RooRealVar*)bgData.get()->find("BscosthetaMC");
    auto& cospsi = *(RooRealVar*)bgData.get()->find("BscospsiMC");
    auto& phi = *(RooRealVar*)bgData.get()->find("BsphiMC");
    
    // terms and coefficients for cospsi pdf
    RooArgList cpsiTerms;
    cpsiTerms.add(RooFit::RooConst(0.5));
    cpsiTerms.add(*new RooLegendre("angBg_P2cpsi", "angBg_P2cpsi", cospsi, 2, 0));
    cpsiTerms.add(*new RooLegendre("angBg_P4cpsi", "angBg_P4cpsi", cospsi, 4, 0));
    RooArgList cpsiCoefs;
    cpsiCoefs.add(RooFit::RooConst(1));
    cpsiCoefs.add(*new RooRealVar("angBg_k1", "angBg_k1", -0.5, 0.5));
    cpsiCoefs.add(*new RooRealVar("angBg_k2", "angBg_k2", -0.5, 0.5));
    
    // terms and coefficients for costheta and phi pdf
    auto L2cth = new RooLegendre("angBg_cthP2", "angBg_P2cth", costheta, 2, 0);
    auto L4cth = new RooLegendre("angBg_cthP4", "angBg_P4cth", costheta, 4, 0);
    auto L6cth = new RooLegendre("angBg_cthP6", "angBg_P6cth", costheta, 6, 0);
    auto c2ph = new RooFormulaVar("angBg_c2ph", "angBg_c2ph", "cos(2*@0)", phi);
    auto c4ph = new RooFormulaVar("angBg_c4ph", "angBg_c4ph", "cos(4*@0)", phi);
    auto c6ph = new RooFormulaVar("angBg_c6ph", "angBg_c6ph", "cos(6*@0)", phi);
    auto sph2_m_0p5 = new RooFormulaVar("angBg_sph2_m_0p5", "angBg_sph2_m_0p5", "sin(@0)*sin(@0) - 0.5", phi);
    auto sph4_m_3o8 = new RooFormulaVar("angBg_sph4_m_3o8", "angBg_sph4_m_3o8", "sin(@0)*sin(@0)*sin(@0)*sin(@0) - 3./8", phi);
    RooArgList cthphTerms;
    cthphTerms.add(RooFit::RooConst(0.5));
    cthphTerms.add(*L2cth);
    cthphTerms.add(*L6cth);
    cthphTerms.add(*c2ph);
    cthphTerms.add(*c4ph);
    cthphTerms.add(*c6ph);
    cthphTerms.add(*new RooProduct("angBg_P2cth_x_ap_sph2_m_0p5_cp", "angBg_P2cth_x_ap_sph2_m_0p5_cp", 
        RooArgList(
            *L2cth,
            *sph2_m_0p5
        ))
    );
    cthphTerms.add(*new RooProduct("angBg_P4cth_x_ap_sph2_m_0p5_cp", "angBg_P4cth_x_ap_sph2_m_0p5_cp", 
        RooArgList(
            *L4cth,
            *sph2_m_0p5
        ))
    );
    cthphTerms.add(*new RooProduct("angBg_P2cth_x_ap_sph4_m_3o8_cp", "angBg_P2cth_x_ap_sph4_m_3o8_cp", 
        RooArgList(
            *L2cth,
            *sph4_m_3o8
        ))
    );
    RooArgList cthphCoefs;
    cthphCoefs.add(RooFit::RooConst(1));
    cthphCoefs.add(*new RooRealVar("angBg_k3", "angBg_k3", -0.5, 0.5));
    cthphCoefs.add(*new RooRealVar("angBg_k4", "angBg_k4", -0.5, 0.5));
    cthphCoefs.add(*new RooRealVar("angBg_k5", "angBg_k5", -0.5, 0.5));
    cthphCoefs.add(*new RooRealVar("angBg_k6", "angBg_k6", -0.5, 0.5));
    cthphCoefs.add(*new RooRealVar("angBg_k7", "angBg_k7", -0.5, 0.5));
    cthphCoefs.add(*new RooRealVar("angBg_k8", "angBg_k8", -2, 2));
    cthphCoefs.add(*new RooRealVar("angBg_k9", "angBg_k9", -2, 2));
    cthphCoefs.add(*new RooRealVar("angBg_k10", "angBg_k10", -2, 2));
    
    // create the sums
    RooRealSumPdf* cpsiPdf = new RooRealSumPdf("angBg_cpsiPdf", "angBg_cpsiPdf", cpsiTerms, cpsiCoefs);
    RooRealSumPdf* cthphPdf = new RooRealSumPdf("angBg_cthphPdf", "angBg_cthphPdf", cthphTerms, cthphCoefs);
    
    // full pdf
    auto angBgPdf = new RooProdPdf("angBg_pdf", "angBg_pdf", RooArgList(*cpsiPdf, *cthphPdf));
    
    auto fitRes = angBgPdf->fitTo(bgData, RooFit::Save(true));
    fitRes->Print("V");
    
    vector<RooRealVar*> varVect = {&costheta, &cospsi, &phi};
    for (int i = 0; i < 3; i++) {
        auto& var = *varVect[i];
        TString varName = var.GetName();
        auto frame = var.frame();
        bgData.plotOn(frame, RooFit::Name("data"));
        angBgPdf->plotOn(frame, RooFit::Name("fit"), RooFit::LineColor(kRed));
        TCanvas c;
        roo_pulls(c, frame, "data", "fit");
        saveCanvas(c, varName + "bg");
        for (int j = 0; j < i; j++) {
            auto& var2 = *varVect[j];
            TString var2Name = var2.GetName();
            // extend the pdf to match the histogram in height
            auto ndivisions = 10;
            auto nbins = ndivisions*ndivisions;
            RooRealVar ext("ext", "ext", bgData.sumEntries()*(var.getMax() - var.getMin())*(var2.getMax() - var2.getMin())/nbins);
            RooExtendPdf extendBgPdf("extendedPdf", "", *angBgPdf, ext);
            auto pdfHist = extendBgPdf.createHistogram(
                "bgPdfHist" + var2Name + varName + suffix, 
                var2, RooFit::Binning(50, var2.getMin(), var2.getMax()), 
                RooFit::YVar(
                    var, RooFit::Binning(50, var.getMin(), var.getMax())
                ),
                RooFit::Extended(true)
            );
            std::cerr << pdfHist->Integral() << std::endl;
            auto dataHist = bgData.createHistogram(
                "bgDataHist" + var2Name + varName + suffix, 
                var2, RooFit::Binning(ndivisions, var2.getMin(), var2.getMax()), 
                RooFit::YVar(
                    var, RooFit::Binning(ndivisions, var.getMin(), var.getMax())
                )
            );
            std::cerr << dataHist->Integral() << std::endl;
            
            dataHist->GetXaxis()->SetTitleOffset(1.6);
            dataHist->GetYaxis()->SetTitleOffset(1.6);
            dataHist->GetZaxis()->SetTitleOffset(1.4);
            
            pdfHist->GetXaxis()->SetTitleOffset(1.6);
            pdfHist->GetYaxis()->SetTitleOffset(1.6);
            pdfHist->GetZaxis()->SetTitleOffset(1.4);
            if (i == 1) {
                // cospsi needs mode space on y
                dataHist->GetYaxis()->SetTitleOffset(2);
                pdfHist->GetYaxis()->SetTitleOffset(2);
            }
            
            TCanvas c2;
            c2.Divide(2,1);
            c2.GetPad(1)->cd();
            auto dataMax = dataHist->GetMaximum();
            dataHist->GetZaxis()->SetRangeUser(0, dataMax*1.05);
            dataHist->Draw("LEGO2 E");
            c2.GetPad(2)->cd();
            pdfHist->GetZaxis()->SetRangeUser(0, dataMax*1.05);
            pdfHist->Draw("surf2");
            saveCanvas(c2, var2Name+varName+"bg");
        }
    }
    
    
    int npars = angBgPdf->getParameters(RooArgSet(costheta, cospsi, phi))->getSize();
    int nbins = 1000;
    
    auto bgDataAsTH = bgData.createHistogram("bgHist" + suffix, costheta, RooFit::Binning(10, -1, 1), RooFit::YVar(cospsi, RooFit::Binning(10, -1, 1)),
                                                     RooFit::ZVar(phi, RooFit::Binning(10, -TMath::Pi(), TMath::Pi())));
    // extend the pdf to match the histogram in height
    RooRealVar ext("ext", "ext", bgData.sumEntries()*8*TMath::Pi()/nbins);
    RooExtendPdf extendBgPdf("extendedPdf", "", *angBgPdf, ext);
    auto bgPdfAsTH = extendBgPdf.createHistogram("bgPdfHist" + suffix, costheta, RooFit::Binning(10, -1, 1), RooFit::YVar(cospsi, RooFit::Binning(10, -1, 1)),
                                                     RooFit::ZVar(phi, RooFit::Binning(10, -TMath::Pi(), TMath::Pi())), RooFit::Extended(true));
    
    
    double chi2 = 0;
    for (int xBin = 1; xBin <= bgDataAsTH->GetNbinsX(); ++xBin) {
        for (int yBin = 1; yBin <= bgDataAsTH->GetNbinsY(); ++yBin) {
            for (int zBin = 1; zBin <= bgDataAsTH->GetNbinsZ(); ++zBin) {
                auto n = bgDataAsTH->GetBinContent(xBin, yBin, zBin);
                auto funcVal = bgPdfAsTH->GetBinContent(xBin, yBin, zBin);
                if (funcVal < 0) cerr << "wrong" << endl;
                auto diff = n - funcVal;
                double err = 0;
                if (diff < 0) {
                    err = 0.5 + sqrt(n + 0.25);
                }
                else {
                    err = -0.5 + sqrt(n + 0.25);
                }
                diff /= err;
                
                chi2 += diff*diff;
            }
        }
    }
//    std::cerr << "data avg: " << bgAngHist->Integral()/nbins << endl;
//    std::cerr << "func avg: " << angBgFuncHist->Integral()/nbins << endl;
    std::cerr << "npars: " << npars << endl;
    std::cerr << "func vs hist: " << chi2 << "/" << nbins - npars << " = " << chi2/(nbins - npars) << std::endl;
    
    auto prevDir = gDirectory;
    TString dirPath = "outputs/";
    if (ensureDir(dirPath)) return nullptr;
    auto file = TFile::Open(dirPath + "angBgPdfDir" + suffix + ".root", "RECREATE");
    
    RooWorkspace angBgWs("angBgWs");
    angBgWs.import(*angBgPdf, RooFit::Silence());
    auto constrPdf = fitRes->createHessePdf(*angBgPdf->getParameters(bgData));
    constrPdf->SetName("angBg_constr");
    angBgWs.import(*constrPdf);
    angBgWs.Write();
    
    file->ReOpen("READ");
    prevDir->cd();

    return file;
}

TDirectory* angBg_fit(vector<const char*> recoFileNames, TString suffix = "") {
    // reco tree
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
    
    return angBg_fit(*bgData, suffix);
}
