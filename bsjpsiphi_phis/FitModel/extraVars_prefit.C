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

void extraVars_prefit(RooDataSet& data, RooWorkspace& mainWs, RooCategory& dsCat) {
    using namespace RooFit;
    using namespace std;
  
    vector<TString> labels;
    labels.push_back("sig");
    labels.push_back("bkg");
    
    auto massPdf = (RooAddPdf*)mainWs.pdf("mass_pdf");
    
    map<string, RooAbsPdf*> simPdfs;
    for(auto catId: RooCatRange(dsCat)) {
        
        TString dsLabel = catId->GetName();
        
        auto ctErrPdf = (RooAddPdf*)mainWs.pdf("ctErr_pdf_" + dsLabel);
        
        RooArgList pdfList;
        RooArgList coefList;
        for (int i = 0; i <= 1; i++) {
            auto ctErrComponent = ctErrPdf->pdfList().at(i);
            auto massComponent = massPdf->pdfList().at(i);
            
            auto joinedPdf = new RooProdPdf(Form("joined_%s_%s", labels[i].Data(), dsLabel.Data()), "", RooArgSet(*ctErrComponent, *massComponent));
            
            pdfList.add(*joinedPdf);
            coefList.add(*new RooRealVar("N" + labels[i] + "_" + dsLabel, "Number of " + labels[i] + "events", data.sumEntries()/labels.size(), 0, data.sumEntries()*2));
        }
        
        simPdfs[dsLabel.Data()] = new RooAddPdf("joined_pdf_" + dsLabel, "Total pdf", pdfList, coefList);
    }
    
    RooSimultaneous joinedPdf("joined_pdf_sim","", simPdfs, dsCat);
    
    RooFitResult *fitRes = joinedPdf.fitTo(data, Save(), Extended(1));
    fitRes->Print("v");
    
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    TH2 *hcorr = fitRes->correlationHist();
    TCanvas cCorr("Correlation Matrix", "Correlation Matrix", 1000, 600);
    hcorr->GetYaxis()->SetTitleOffset(1.4);
    hcorr->Draw("colz");
    saveCanvas(cCorr, "joinedFitCorrMatrix");
    gStyle->SetPalette(kBird);
    
    auto svmass = (RooRealVar*)data.get()->find(massName);
    RooPlot *massFrame = svmass->frame();
    data.plotOn(massFrame, DataError(RooAbsData::SumW2), Name("data"));
    joinedPdf.plotOn(massFrame, LineColor(kRed), LineStyle(kDashed), Name("bkg"), Components("joined_bkg_*"), ProjWData(dsCat, data));
    joinedPdf.plotOn(massFrame, LineColor(kGreen), LineStyle(kDashed), Name("sig"), Components("joined_sig_*"), ProjWData(dsCat, data));
    joinedPdf.plotOn(massFrame, LineColor(kBlue), Name("fit"), ProjWData(dsCat, data));
    joinedPdf.paramOn(massFrame, Layout(0.6), Parameters(*(RooArgSet*)joinedPdf.getParameters(data)->selectByName("mass_*,Nsig_*,Nbkg_*")));
    
    gStyle->SetOptTitle(0);
    TCanvas cM;
    cM.SetFillColor(0);
    cM.SetBorderSize(2);
    roo_pulls(cM, massFrame, "data", "fit", fitRes->floatParsFinal().getSize());
    
    saveCanvas(cM, "massPlotJoinedFit");
    
    auto ctErr = (RooRealVar*)data.get()->find(ctErrName);
    for(auto catId: RooCatRange(dsCat)) {
        
        TString dsLabel = catId->GetName();
        
        RooPlot *ctErrFrame = ctErr->frame();
        data.plotOn(ctErrFrame, DataError(RooAbsData::SumW2), Name("data"), Cut("dsLabel == dsLabel::" + dsLabel));
        joinedPdf.plotOn(ctErrFrame, LineColor(kRed), LineStyle(kDashed), Name("bkg"), Components("joined_bkg_" + dsLabel), Slice(dsCat, catId->GetName()), ProjWData(dsCat, data));
        joinedPdf.plotOn(ctErrFrame, LineColor(kGreen), LineStyle(kDashed), Name("sig"), Components("joined_sig_" + dsLabel), Slice(dsCat, catId->GetName()), ProjWData(dsCat, data));
        joinedPdf.plotOn(ctErrFrame, LineColor(kBlue), Name("fit"), Slice(dsCat, catId->GetName()), ProjWData(dsCat, data));
        joinedPdf.paramOn(ctErrFrame, Parameters(*(RooArgSet*)joinedPdf.getParameters(data)->selectByName(Form("ctErr_*_%1$s,Nsig_%1$s,Nbkg_%1$s", dsLabel.Data()))));
        
        gStyle->SetOptTitle(0);
        TCanvas cW;
        cW.SetFillColor(0);
        cW.SetBorderSize(2);
        roo_pulls(cW, ctErrFrame, "data", "fit", fitRes->floatParsFinal().getSize());
        
        saveCanvas(cW, "ctErrPlotJoinedFit_" + dsLabel);
    }
    mainWs.import(*(RooArgSet*)joinedPdf.getParameters(data)->selectByName("Nsig_*,Nbkg_*"));
}

//void extraVars_prefit(const char* filename, TString massModelFileName, TString ctErrModelFileName) {
//    auto recoTree = getTree(filename, recoTreeName);
//    if (!recoTree) {
//        return nullptr;
//    }
//    
//    auto mainWs = getMainWs();
//    
//    auto BsCt2DMC = mainWs->var(ctName);
//    auto BsCt2DMCErr = mainWs->var(ctErrName);
//    auto BscosthetaMC = mainWs->var(cosThetaName);
//    auto BscospsiMC = mainWs->var(cosPsiName);
//    auto BsphiMC = mainWs->var(phiName);
//    auto svmass = mainWs->var(massName);
//    auto mistag = mainWs->var(mistagName);
//    auto tag = mainWs->cat(tagCatName);
//    
//    RooDataSet recoData("recoData", "raw data1", RooArgSet(*BscosthetaMC, *BscospsiMC, *BsphiMC, *BsCt2DMC, *BsCt2DMCErr, *svmass, *tag, *mistag),
//                                      RooFit::Import(*recoTree));
//    
//    auto massModelFile = openFileForReading(massModelFileName);
//    auto ctErrModelFile = openFileForReading(ctErrModelFileName);
//    
//    auto massWs = (RooWorkspace*)massModelFile->Get("mass");
//    auto ctErrWs = (RooWorkspace*)ctErrModelFile->Get("ctErr");
//    
//    extraVars_prefit(recoData, *massWs, *ctErrWs);
//}
