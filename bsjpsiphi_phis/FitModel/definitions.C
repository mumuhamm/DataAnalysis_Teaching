// Authors: Enrico Lusiani
#pragma once

#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooCategory.h"

#include "TMath.h"

#include <vector>

auto recoTreeName = "treeFit";
auto genTreeName = "OutTree";

const char* ctNameGEN = "ctau_GEN";
const char* cosThetaNameGEN = "angle_costheta_GEN";
const char* cosPsiNameGEN = "angle_cospsi_GEN";
const char* phiNameGEN = "angle_phi_GEN";
const char* tagCatNameGEN = "MC_Flavour";

#ifndef GENONLY
const char* ctName = "BsCt2DMC";
const char* cosThetaName = "BscosthetaMC";
const char* cosPsiName = "BscospsiMC";
const char* phiName = "BsphiMC";
const char* tagCatName = "tag";
#else
const char* ctName = ctNameGEN;
const char* cosThetaName = cosThetaNameGEN;
const char* cosPsiName = cosPsiNameGEN;
const char* phiName = phiNameGEN;
const char* tagCatName = tagCatNameGEN;
#endif

const char* genCtName = "BsCt2DMC_GEN";
const char* ctErrName = "BsCt2DMCErr";
const char* massName = "svmass";
const char* mistagName = "mistag";
const char* dsCatName = "dsLabel";
const char* sbStatusName = "sbStatus";

auto ctNameGen = "ctau_GEN";
auto cosThetaNameGen = "angle_costheta_GEN";
auto cosPsiNameGen = "angle_cospsi_GEN";
auto phiNameGen = "angle_phi_GEN";

// sidebands bounds
auto sbLowerBound = 5.28;
auto sbUpperBound = 5.45;

RooWorkspace* getMainWs() {
    auto mainWs = new RooWorkspace("main");
    
    auto minCt = 0.007;
    auto maxCt = 0.5;
#ifdef GENONLY
    minCt = 0;
    maxCt = 0.7;
#endif
    auto ct = new RooRealVar(ctName, "Bs ct", minCt, maxCt, "cm");
    mainWs->import(*ct);
    auto ctErr = new RooRealVar(ctErrName, "Bs ct Err", 0.0002, 0.005, "cm");
    mainWs->import(*ctErr);
    auto cosTheta = new RooRealVar(cosThetaName, "cos(#theta_{T})", -1, 1);
    mainWs->import(*cosTheta);
    auto cosPsi = new RooRealVar(cosPsiName, "cos(#psi_{T})", -1, 1);
    mainWs->import(*cosPsi);
    auto phi = new RooRealVar(phiName, "#phi_{T}", -TMath::Pi(), TMath::Pi(), "rad");
    mainWs->import(*phi);
    
    auto minMass = 5.24;
    auto maxMass = 5.49;
    auto svmass = new RooRealVar(massName, "M_{B_{s}}", minMass, maxMass, "GeV/c^{2}");
    svmass->setRange("sbLeft", minMass, sbLowerBound);
    svmass->setRange("sbRight", sbUpperBound, maxMass);
    svmass->setRange("center", sbLowerBound, sbUpperBound);
    mainWs->import(*svmass);
    auto mistag = new RooRealVar(mistagName, "Mistag fraction of original B and Bbar", 0, 0, 1.0);
    mainWs->import(*mistag);
        
    auto tag = new RooCategory(tagCatName, "Flavour tag of the B meson");
    tag->defineType("Bsbar", -1);
    tag->defineType("Bs", +1);
    tag->defineType("untag", 0);
//    auto tag = new RooRealVar(tagCatName, "tag", -2, 2);
    mainWs->import(*tag);
    
    RooThresholdCategory sbStatus(sbStatusName, "sbStatus", *svmass, "sb");
    sbStatus.addThreshold(sbLowerBound, "sb");
    sbStatus.addThreshold(sbUpperBound, "cent");
    mainWs->import(sbStatus);
    
    return mainWs;
}

RooCategory createDsCat(std::vector<const char*> dsLabels) {
    RooCategory dsCat(dsCatName, " Dataset label");
    for (auto label: dsLabels) {
        std::clog << "Defining " << label << endl;
        dsCat.defineType(label);
    }
    return dsCat;
}

RooAbsData* getRawSideband(RooAbsData& data) {
    auto cutExpr = TString::Format("%1$s < %2$g || %1$s > %3$g", massName, sbLowerBound, sbUpperBound);
    return data.reduce(RooFit::Cut(cutExpr));
}

void shiftToMCValues(RooArgSet& variables) {
    
    ((RooRealVar*)variables.find("A_0"))->setVal(0.523);
    ((RooRealVar*)variables.find("A_pe"))->setVal(0.249);
    ((RooRealVar*)variables.find("A_S"))->setVal(0.012);
    ((RooRealVar*)variables.find("ctau"))->setVal(0.04527);
    ((RooRealVar*)variables.find("deltaPa"))->setVal(3.19);
    ((RooRealVar*)variables.find("deltaPe"))->setVal(3.2);
    ((RooRealVar*)variables.find("deltaSPe"))->setVal(0.188);
    ((RooRealVar*)variables.find("dGam"))->setVal(3.00);
    ((RooRealVar*)variables.find("dm"))->setVal(592.3);
    ((RooRealVar*)variables.find("phi_s"))->setVal(-0.055);
    ((RooRealVar*)variables.find("lambda"))->setVal(1);
}
