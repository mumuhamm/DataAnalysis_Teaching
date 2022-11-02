// Authors: Author: Md. Alibordi
#pragma once

#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"

#include "RooWorkspace.h"
#include "RooDataSet.h"

#include "definitions.C"

#include "utils.h"

void Bs_angEff_closure(TTree* recoTree, TTree* genTree, RooWorkspace* effWs) {
    using namespace RooFit;
    
    auto mainWs = getMainWs();
    
    auto BsCt2DMC = mainWs->var(ctName);
    auto BsCt2DMCErr = mainWs->var(ctErrName);
    auto BscosthetaMC = mainWs->var(cosThetaName);
    auto BscospsiMC = mainWs->var(cosPsiName);
    auto BsphiMC = mainWs->var(phiName);
    auto svmass = mainWs->var(massName);
    auto mistag = mainWs->var(mistagName);
    auto tag = mainWs->cat(tagCatName);
    
    mistag->setVal(0);
    mistag->setConstant();
    
    mainWs->import(*effWs->function("angEffFunc"), RecycleConflictNodes());
    
    auto eff = mainWs->function("angEffFunc");
    
    clog << "Importing reco data" << endl;
    RooDataSet recoData("recoData", "raw data1", RooArgSet(*BscosthetaMC, *BscospsiMC, *BsphiMC, *BsCt2DMC, *BsCt2DMCErr, *svmass, *tag, *mistag),
                                      Import(*recoTree));
    
    clog << "Importing gen data" << endl;
    
    Float_t costhetagen, cospsigen, phigen;
    genTree->SetBranchAddress("angle_costheta_GEN", &costhetagen);
    genTree->SetBranchAddress("angle_cospsi_GEN", &cospsigen);
    genTree->SetBranchAddress("angle_phi_GEN", &phigen);
    
    RooRealVar weight("eff", "", 0, 1);
    
    RooArgSet angVars(*BscosthetaMC, *BscospsiMC, *BsphiMC, weight);
    RooDataSet genDataNoW("genData", "GEN distribution", angVars);
    auto genEntries = genTree->GetEntries();
    auto step20 = genEntries / 20;
    
    for (int iCand = 0; iCand < genEntries; ++iCand) {
        if (iCand % step20 == 0) clog << iCand << "/" << genEntries << endl;
        genTree->GetEntry(iCand);
        
        BscosthetaMC->setVal(costhetagen);
        BscospsiMC->setVal(cospsigen);
        BsphiMC->setVal(phigen);
        weight.setVal(eff->getVal());
        
        genDataNoW.add(angVars);
    }
    
    clog << "Importing with weights" << endl;
    RooDataSet genData("genData", "gen data", RooArgSet(*BscosthetaMC, *BscospsiMC, *BsphiMC, weight), Import(genDataNoW), WeightVar(weight));
    
    clog << "Plotting" << endl;
    TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);
    RooPlot* costhetaframe = BscosthetaMC->frame(Title("cos(#theta) distributions"));
    RooPlot* cospsiframe = BscospsiMC->frame(Title("cos(#psi) distributions"));
    RooPlot* phiframe = BsphiMC->frame(Title("#phi distributions"));
    genData.plotOn(costhetaframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),DataError(RooAbsData::SumW2),Name("plDenDist"));
    genData.plotOn(cospsiframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),DataError(RooAbsData::SumW2));
    genData.plotOn(phiframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),DataError(RooAbsData::SumW2));
    recoData.plotOn(costhetaframe,Binning(30),Name("plNumDist"));
    recoData.plotOn(cospsiframe,Binning(30));
    recoData.plotOn(phiframe,Binning(30));
    costhetaframe->GetYaxis()->SetTitleOffset(1.6);
    cospsiframe->GetYaxis()->SetTitleOffset(1.6);
    phiframe->GetYaxis()->SetTitleOffset(1.6);
    costhetaframe->SetMaximum(costhetaframe->GetMaximum()*1.15);
    cospsiframe->SetMaximum(cospsiframe->GetMaximum()*1.15);
    phiframe->SetMaximum(phiframe->GetMaximum()*1.15);
    leg->SetTextSize(0.03);
    leg->AddEntry(costhetaframe->findObject("plNumDist"),"Post-selection RECO distribution" ,"lep");
    leg->AddEntry(costhetaframe->findObject("plDenDist"),"Efficiency-corrected GEN distribution" ,"lep");
    
    TCanvas c;
    c.Divide(3,1);
    c.cd(1);
    gPad->SetLeftMargin(0.17); 
    costhetaframe->Draw();
    leg->Draw("same");
    c.cd(2);
    gPad->SetLeftMargin(0.17); 
    cospsiframe->Draw();
    leg->Draw("same");
    c.cd(3);
    gPad->SetLeftMargin(0.17); 
    phiframe->Draw();
    leg->Draw("same");
    
    saveCanvas(c, "ang_closure");
}

void Bs_angEff_closure(const char* recoFileName, const char* genFileName, const char* effFileName) {
    auto recoTree = getTree(recoFileName, recoTreeName);
    if (!recoTree) {
        return;
    }
    
    auto genTree = getTree(genFileName, genTreeName);
    if (!genTree) {
        return;
    }
    
    auto prevDir = gDirectory;
    auto effModelFile = TFile::Open(effFileName);
    prevDir->cd();
    
    auto effWs = (RooWorkspace*)effModelFile->Get("ang");
    
    Bs_angEff_closure(recoTree, genTree, effWs);
}
