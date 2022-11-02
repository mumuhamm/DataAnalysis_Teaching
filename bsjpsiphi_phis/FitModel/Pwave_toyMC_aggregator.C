// Author Enrico Lusiani

#include "RooFitResult.h"

#include "RooJohnsonLocal.cxx"
#include "RooTagPdf.cxx"

#include "definitions.C"
#include "utils.h"


// TODO remove, this is because of a change during toy run
const char* transformAngBgNames(TString name, bool needChange) {
    if (not needChange) {
        return name;
    }

    if (not name.Contains("angBg")) {
        return name;
    }
    
    static const map<TString, const char*> nameMap = [](){
        map<TString, const char*> m;
        m["angBg_k1"] = "angBg_k5";
        m["angBg_k2"] = "angBg_k6";
        m["angBg_k3"] = "angBg_k1";
        m["angBg_k4"] = "angBg_k2";
        m["angBg_k5"] = "angBg_k10";
        m["angBg_k6"] = "angBg_k11";
        m["angBg_k7"] = "angBg_k12";
        m["angBg_k8"] = "angBg_k3";
        m["angBg_k9"] = "angBg_k13";
        m["angBg_k10"] = "angBg_k4";
        
        return m;
    }();
    
    return nameMap.at(name);
}

void Pwave_toyMC_aggregator(TString dirPrefix, int nGen, int samplesPerGen, RooWorkspace* main) {
    TString resName = "fitresult_Angular_Model_Angular_Model_adjData";
    
    TString dirPath = "outputs/";
    if (ensureDir(dirPath)) return;
    auto outFile = TFile::Open(dirPath + "toyMCPars.root", "RECREATE");
    auto parTree = new TTree("parTree", "parTree");
    map<string, double> bufferpull;
    map<string, double> buffermean;
    map<string, double> buffererror;
    
    int covQual;
    int hesseStatus;
    int migradStatus;
    parTree->Branch("covQual", &covQual);
    parTree->Branch("hesseStat", &hesseStatus);
    parTree->Branch("migrStat", &migradStatus);
    
    
    for (int iGen = 0; iGen < nGen; iGen++) {
        clog << "Processing generation " << iGen << endl;
        for (int iSample = 0; iSample < samplesPerGen; iSample++) {
            clog << "\tsample " << iSample << endl;
            auto prevDir = gDirectory;
            auto file = TFile::Open(Form("%s_%d/fit.%d.root", dirPrefix.Data(), iGen, iSample));
            auto dataFile = TFile::Open(Form("%s_%d/data.%d.root", dirPrefix.Data(), iGen, iSample));
            prevDir->cd();
            auto res = (RooFitResult*)file->Get(resName);
            auto& pars = res->floatParsFinal();
            auto dataWs = (RooWorkspace*)dataFile->Get("genWs");
            auto data = dataWs->data("Angular_Model_adjData");
            
            bool newAngBg = pars.find("angBg_k7");
            bool newAngBgMain = main->var("angBg_k7");
            
            if (buffermean.size() == 0) {
                for (auto par: pars) {
                    auto varname = par->GetName();
                    
                    auto normVarName = transformAngBgNames(varname, newAngBg);
                    
                    parTree->Branch(Form("%sMean",normVarName),&buffermean[normVarName],Form("%sMean/D",normVarName));
                    parTree->Branch(Form("%sError",normVarName),&buffererror[normVarName],Form("%sError/D",normVarName));
                    parTree->Branch(Form("%sPull",normVarName),&bufferpull[normVarName],Form("%sPull/D",normVarName));
                }
            }
            
            for (auto parBase: pars) {
                auto par = (RooRealVar*)parBase;
                auto varname = par->GetName();
                
                auto normVarName = transformAngBgNames(varname, newAngBg and not newAngBgMain);
                
                auto realVal = main->var(normVarName)->getVal();
                TString varnameTStr = varname;
                if (varnameTStr.Contains("Nsig_") or varnameTStr.Contains("Nbkg_")) {
                    bool isSig = varnameTStr.Contains("sig");
                    TString label = varnameTStr(5, 100);
                    
                    TString cutString = TString::Format("sigOrBkg == sigOrBkg::%s && dsLabel == dsLabel::%s", isSig?"sig":"bkg", label.Data());
                    
                    realVal = data->sumEntries(cutString);
                }
                
                
                buffermean[normVarName] = par->getVal();
                buffererror[normVarName] = par->getError();
                bufferpull[normVarName] = (realVal - buffermean[normVarName])/buffererror[normVarName];
            }
            covQual = res->covQual();
            migradStatus = res->statusCodeHistory(0);
            hesseStatus = res->statusCodeHistory(1);
            parTree->Fill();
            
            delete file;
            delete dataFile;
        }
    }
    outFile->cd();
    parTree->Write();
    
    outFile->Close();
}

void Pwave_toyMC_aggregator(TString dirPrefix, int nGen, int samplesPerGen, TString mainFileName) {
    auto mainFile = TFile::Open(mainFileName);
    
    auto mainWs = (RooWorkspace*)mainFile->Get("postFit");
    
    mainWs->cat("blindCat")->setLabel("unblind");
    
    auto allVars = mainWs->allVars();
    shiftToMCValues(allVars);
    
    Pwave_toyMC_aggregator(dirPrefix, nGen, samplesPerGen, mainWs);
}
