//Authors: Md Alibordi, Alberto Bragagnolo, Giacomo Fedi, Enrico Lusiani

#include "TTree.h"
#include "TFile.h"
#include "TString.h"

#include "RooProdPdf.h"
#include "RooRealSumPdf.h"
#include "RooGenericPdf.h"
#include "RooGaussModel.h"
#include "Roo1DTable.h"

#include <map>
#include <vector>
#include <fstream>
#include <sstream>

#include "utils.h"

#include "definitions.C"

#include "Bs_ctEff_fit.C"
#include "Bs_angEff_proj.C"
#include "massPdf_fit.C"
#include "ctErrPdf_fit.C"
#include "ctBg_fit.C"
#include "angBg_fit.C"
#include "extraVars_prefit.C"

#include "Pwave_Model.C"

#include "Pwave_Model_plot.C"

// TODO general issues found by Enrico
// fix some names

void Pwave_Reco_Model_recotag(std::map<string, TTree*> recoTrees, 
                              TTree* genTree, 
                              // mass pdf dir
                              TDirectory* massPdfDir,
                              // efficiencies
                              std::map<string, TDirectory*> ctSigEffDirs, 
                              std::map<string, TDirectory*> angSigEffDirs, 
                              // ctErr (joined)
                              std::map<string, TDirectory*> ctErrPdfDirs,
                              // misc flags
                              bool pdgPars,
                              bool isData,
                              bool plotAngFlag,
                              int numCpu) {
    
    using namespace RooFit;
    using namespace std;
    
    vector<const char*> dsLabels;
    for (auto& entry: recoTrees) {
        cerr << "Label \"" << entry.first << "\"" << endl;
        dsLabels.push_back(entry.first.c_str());
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
    
    auto dsCat = createDsCat(dsLabels);
    mainWs->import(dsCat);
    
    cerr << "Importing reco datasets" << endl;
    
    map<string, RooDataSet*> recoDataByLabel;
    for (auto treeEntry: recoTrees) {
        auto& label = treeEntry.first;
        auto& recoTree = treeEntry.second;
        
        cerr << "Importing for label \"" << label << "\""<< endl;
        
        recoDataByLabel[label] = new RooDataSet(Form("recoData_%s", label.c_str()), "raw data1", 
            RooArgSet(*BscosthetaMC, *BscospsiMC, *BsphiMC, *BsCt2DMC, *BsCt2DMCErr, *svmass, *tag, *mistag),
            Import(*recoTree)
        );
    }
    
    //RooRealVar weight("weight", "", 0, std::numeric_limits<double>::max());
    cerr << "Merging datasets" << endl;
    auto recoData = new RooDataSet("recoData", "raw data", RooArgSet(*BscosthetaMC, *BscospsiMC, *BsphiMC, *BsCt2DMC, *BsCt2DMCErr, *svmass, *tag, *mistag),//, weight), 
        Index(dsCat),
        Import(recoDataByLabel)
    );
    
    auto btable = recoData->table(*tag);
    btable->Print("V");
    
    RooArgSet angVars(*BscosthetaMC, *BscospsiMC, *BsphiMC, *BsCt2DMC);
    auto genData = new RooDataSet("genData", "GEN distribution", angVars);
    
    bool needGen = false;
    for (auto label: dsLabels)
    {
        if (not angSigEffDirs[label]) {
            needGen = true;
            break;
        }
    }
    if (needGen) {
        cerr << "Importing gen dataset" << endl;
    
        Float_t costhetagen, cospsigen, phigen, ctaugen;
        genTree->SetBranchAddress("angle_costheta_GEN", &costhetagen);
        genTree->SetBranchAddress("angle_cospsi_GEN", &cospsigen);
        genTree->SetBranchAddress("angle_phi_GEN", &phigen);
        genTree->SetBranchAddress(ctNameGen, &ctaugen);
    
        auto genEntries = genTree->GetEntries();
        auto step20 = genEntries / 20;
        
        for (int iCand = 0; iCand < genEntries; ++iCand) {
            if (iCand % step20 == 0) cerr << iCand << "/" << genEntries << endl;
            genTree->GetEntry(iCand);
            
            BscosthetaMC->setVal(costhetagen);
            BscospsiMC->setVal(cospsigen);
            BsphiMC->setVal(phigen);
            BsCt2DMC->setVal(ctaugen);
            
            genData->add(angVars);
        }
    }
    
    // =====================================CT EFF======================================================================
    
    for (auto label: dsLabels) {
    
        if (ctSigEffDirs.count(label) < 1) {
            cerr << "No ct efficiency entry for label " << label << endl;
            cerr << "This is a BUG" << endl;
            cerr << "ABORTING" << endl;
        }
        
        auto& ctEffDir = ctSigEffDirs[label];
        
        if (isData and not ctEffDir) {
            cerr << "No ctefficiency specified for a data sample" << endl;
            cerr << "This is a BUG" << endl;
            cerr << "ABORTING" << endl;
            
            return;
        }
        
        if (not ctEffDir) {
            ctEffDir = Bs_ctEff_fit(recoTrees[label], genTree, *BsCt2DMC, label);
        }
        
        if (!ctEffDir) {
            cerr << "Unable to compute the ct efficiency for " << label << endl;
            return;
        }
        auto ctPullHist = (TH1*)ctEffDir->Get("pullbulk");
        
        auto kappa_val = ctPullHist->GetRMS();
        
        auto kappa = new RooRealVar(Form("kappa%s", label), "kappa", kappa_val);
        mainWs->import(*kappa);
    
        auto ctEffWs = (RooWorkspace*)ctEffDir->Get("ctEff");
        
        auto ctEff = ctEffWs->function("ctEffFunc");
        for (auto varBase: *ctEff->getParameters(*recoData)) {
            auto var = (RooRealVar*)varBase;
            var->setConstant();
        }
        
        //RooFormulaVar weight("weight", "", "1/ctEffFunc", *ctEff);
        
        //recoDataByLabel[label]->addColumn(weight);
        
        // mixed style for formula because the ct variable has always the same name, the parameters do not
        mainWs->import(*ctEff, 
            RenameAllVariablesExcept(label, ctName), 
            RenameAllNodes(label), 
            Silence()
        );
    }
    
    // ====================================MASS=============================================================================
    
    if (!massPdfDir) {
        massPdfDir = massPdf_fit(*recoData, *svmass, isData);
    }
    
    if (!massPdfDir) {
        cerr << "Unable to fit the mass pdf" << endl;
        return;
    }
    
    auto massWs = (RooWorkspace*)massPdfDir->Get("mass");
    
    auto massPdf = massWs->pdf("mass_pdf");
    
//    for (auto varBase: *massPdf->getParameters(*recoData)) {
//        auto var = (RooRealVar*)varBase;
//        var->setConstant();
//    }
    
    mainWs->import(*massPdf, Silence());
    
    map<string, map<string, RooDataSet*>> sbSplitDataByLabel;
    for (auto label: dsLabels) {
        sbSplitDataByLabel[label] = bkgSigSplitWSB(*recoDataByLabel[label], *mainWs, *sbStatus, isData);
    }
    
    // ====================================ANG EFF======================================================================
    for (auto label: dsLabels) {
    
        if (angSigEffDirs.count(label) < 1) {
            cerr << "No angular efficiency entry for label " << label << endl;
            cerr << "This is a BUG" << endl;
            cerr << "ABORTING" << endl;
        }
        
        auto& angEffDir = angSigEffDirs[label];
        
        if (isData and not angEffDir) {
            cerr << "No angular efficiency specified for a data sample" << endl;
            cerr << "This is a BUG" << endl;
            cerr << "ABORTING" << endl;
            return;
        }
        
        if (!angEffDir) {
            angEffDir = Bs_angEff_proj(*recoDataByLabel[label], *genData, angVars, label);
        }
        
        if (!angEffDir) {
            cerr << "Unable to compute the angular function" << endl;
            return;
        }
        
        auto angWs = (RooWorkspace*)angEffDir->Get("ang");
        
        auto angEffFn = angWs->function("angEffFunc");
        
        mainWs->import(*angEffFn, 
            RenameAllVariablesExcept(label, Form("%s,%s,%s", cosThetaName, cosPsiName, phiName)),
            RenameAllNodes(label), 
            Silence()
        );
    }
    
    if (isData){
    
        auto rawSbData = (RooDataSet*)getRawSideband(*recoData);
        
        RooArgSet bgConstr;
        
        // =================================CT BKG PDF =================================================================
        
        auto ctBkgPdfDir = ctBg_fit(*rawSbData);
        
        auto ctBkgPdfWs = (RooWorkspace*)ctBkgPdfDir->Get("ctBgWs");
        
        auto ctBgPdf = ctBkgPdfWs->pdf("ctBg_pdf");
    
//        for (auto varBase: *ctBgPdf->getParameters(*recoData)) {
//            auto var = (RooRealVar*)varBase;
//            
//            var->setConstant();
//        }
        
        mainWs->import(*ctBgPdf, Silence());
        bgConstr.add(*ctBkgPdfWs->pdf("ctBg_constr"));
        
        // =================================ANG BKG PDF ================================================================
        
        auto angBkgPdfDir = angBg_fit(*rawSbData);
        
        auto angBkgPdfWs = (RooWorkspace*)angBkgPdfDir->Get("angBgWs");
        
        auto angBgPdf = angBkgPdfWs->pdf("angBg_pdf");
    
//        for (auto varBase: *angBgPdf->getParameters(*recoData)) {
//            auto var = (RooRealVar*)varBase;
//            
//            var->setConstant();
//        }
        
        mainWs->import(*angBgPdf, Silence());
        bgConstr.add(*angBkgPdfWs->pdf("angBg_constr"));
        
        mainWs->defineSet("bgConstr", bgConstr, true);
    }
    
    // =====================================CT ERR======================================================================
    
    for (auto label: dsLabels) {
    
        if (ctErrPdfDirs.count(label) < 1) {
            cerr << "No ctErr pdf entry for label " << label << endl;
            cerr << "This is a BUG" << endl;
            cerr << "ABORTING" << endl;
        }
        
        auto& ctErrPdfDir = ctErrPdfDirs[label];
        
        if (!ctErrPdfDir) {
            ctErrPdfDir = ctErrPdf_fit(sbSplitDataByLabel[label], *BsCt2DMCErr, isData, label);
        }
        
        if (!ctErrPdfDir) {
            cerr << "Unable to fit the ct error pdf" << endl;
            return;
        }
        
        auto ctErrWs = (RooWorkspace*)ctErrPdfDir->Get("ctErr");
        
        mainWs->import(*ctErrWs->function("ctErr_pdf"), 
            RenameAllVariablesExcept(label, ctErrName), 
            RenameAllNodes(label), 
            Silence()
        );
    }
    
    // ======================================REFIT MASS+CTERR==============================================================
    
    if (isData) {
        extraVars_prefit(*recoData, *mainWs, dsCat);
    }
    else {
        for (auto label: dsLabels) {
            auto nsigByLabel = new RooRealVar(Form("Nsig_%s", label), "", recoDataByLabel[label]->sumEntries(), 0, 1.5*recoDataByLabel[label]->sumEntries());
            nsigByLabel->setConstant();
            mainWs->import(*nsigByLabel);
        }
    }
    
    for (auto label: dsLabels) {
        auto ctErrPdf = mainWs->pdf(Form("ctErr_pdf_%s", label));
        
        for(auto varBase: *ctErrPdf->getParameters(*recoData)) {
            auto var = (RooRealVar*)varBase;
            
            var->setConstant();
        }
    }
    
    // ====================================================================================================================
    
    Pwave_Model(*recoData, mainWs, false, pdgPars, isData, plotAngFlag, numCpu);
}

//void Pwave_Reco_Model_recotag(TString recoFileName, TString genFileName,
//    bool pdgPars, 
//    const char* label,
//    const char* ctEffDirName = nullptr, 
//    const char* angEffDirName = nullptr, 
//    const char* massPdfDirName = nullptr,
//    const char* ctErrPdfDirName = nullptr,
//    bool plotAngFlag = true) {
//    
//    // reco tree
//    auto recoTree = getTree(recoFileName, recoTreeName);
//    if (!recoTree) {
//        return;
//    }
//    
//    // gen tree
//    auto genTree = getTree(genFileName, genTreeName);
//    if (!genTree) {
//        return;
//    }
//    
//    TDirectory* ctEffDir = nullptr;
//    if (ctEffDirName) {
//        auto prevDir = gDirectory;
//        ctEffDir = TFile::Open(ctEffDirName);
//        prevDir->cd();
//    }
//    
//    TDirectory* angEffDir = nullptr;
//    if (angEffDirName) {
//        auto prevDir = gDirectory;
//        angEffDir = TFile::Open(angEffDirName);
//        prevDir->cd();
//    }
//    
//    TDirectory* massPdfDir = nullptr;
//    if (massPdfDirName) {
//        auto prevDir = gDirectory;
//        massPdfDir = TFile::Open(massPdfDirName);
//        prevDir->cd();
//    }
//    
//    TDirectory* ctErrPdfDir = nullptr;
//    if (ctErrPdfDirName) {
//        auto prevDir = gDirectory;
//        ctErrPdfDir = TFile::Open(ctErrPdfDirName);
//        prevDir->cd();
//    }
//    
//    string theOnlyLabel = label;
//    
//    Pwave_Reco_Model_recotag({{theOnlyLabel, recoTree}}, genTree, {{theOnlyLabel, ctEffDir}}, {{theOnlyLabel, angEffDir}}, massPdfDir, {{theOnlyLabel, ctErrPdfDir}}, pdgPars, plotAngFlag, maxOrder, xbins_angEff, ybins_angEff, zbins_angEff);
//}

void Pwave_Reco_Model_recotag(TString recoDatasets, TString extraConfig) {
    
    using namespace std;
    
    ifstream extraConfigIn(extraConfig);
    string line;
    
    map<string, string> extraConfigMap;
    while(getline(extraConfigIn, line)) {
        stringstream lineReader(line);
        
        string key;
        lineReader >> key;
        
        string val;
        lineReader >> val;
        
        extraConfigMap[key] = val;
    }
    
    auto isData = false;
    if (extraConfigMap.count("isData") > 0){
        isData = extraConfigMap["isData"] == "true";
    }
    else {
        cerr << "Missing configuration: isData" << endl;
        return;
    }
    
    const char* massPdfDirName = nullptr;
    if (extraConfigMap.count("massPdfDir") > 0){
        massPdfDirName = extraConfigMap["massPdfDir"].c_str();
    }
    TDirectory* massPdfDir = nullptr;
    if (massPdfDirName) {
        auto prevDir = gDirectory;
        massPdfDir = TFile::Open(massPdfDirName);
        prevDir->cd();
    }
    
    auto usePdgPars = false;
    if (extraConfigMap.count("usePdgPars") > 0){
        usePdgPars = extraConfigMap["usePdgPars"] == "true";
    }
    else {
        cerr << "Missing configuration: usePdgPars" << endl;
        return;
    }
    
    auto doPlot = false;
    if (extraConfigMap.count("doPlot") > 0){
        doPlot = extraConfigMap["doPlot"] == "true";
    }
    else {
        // it's fine, we just skip it
    }
    
    int numCpu = 2;
    if (extraConfigMap.count("numCpu") > 0){
        numCpu = stoi(extraConfigMap["numCpu"]);
    }
    else {
        // it's fine, we just skip it
    }
    
    ifstream recoConfigFile(recoDatasets);
    
    bool needgen = false;
    map<string, TTree*> recoTrees;
    map<string, TDirectory*> ctEffDirs;
    map<string, TDirectory*> angEffDirs;
    map<string, TDirectory*> ctErrDirs;
    while(getline(recoConfigFile, line)) {
        stringstream lineReader(line);
        
        string dsLabel;
        lineReader >> dsLabel;
        
        // reco tree
        TString recoFileName;
        lineReader >> recoFileName;
        recoTrees[dsLabel] = getTree(recoFileName, recoTreeName);
        if (!recoTrees[dsLabel]) {
            return;
        }
        
        TString ctEffDirName;
        lineReader >> ctEffDirName;
        if (ctEffDirName != "null") {
            auto prevDir = gDirectory;
            ctEffDirs[dsLabel] = TFile::Open(ctEffDirName);
            prevDir->cd();
            if (not ctEffDirs[dsLabel]) {
                cerr << "Unable to open the ct eff file " << ctEffDirName << endl;
                return;
            }
        }
        else {
            if (isData)
            {
                cerr << "No ct efficiency was provided but we are running in data" << endl
                    << "ABORTING" << endl;
                return;
            }
            cerr << "Scheduling computation for ct eff for label " << dsLabel << endl;
            ctEffDirs[dsLabel] = nullptr;
            needgen = true;
        }
        
        TString angEffDirName;
        lineReader >> angEffDirName;
        if (angEffDirName != "null") {
            auto prevDir = gDirectory;
            angEffDirs[dsLabel] = TFile::Open(angEffDirName);
            prevDir->cd();
            if (not angEffDirs[dsLabel]) {
                cerr << "Unable to open the ang eff file " << angEffDirName << endl;
                return;
            }
        }
        else {
            if (isData)
            {
                cerr << "No angular efficiency was provided but we are running in data" << endl
                    << "ABORTING" << endl;
                return;
            }
            cerr << "Scheduling computation for ang eff for label " << dsLabel << endl;
            angEffDirs[dsLabel] = nullptr;
            needgen = true;
        }
        
        TString ctErrDirName;
        lineReader >> ctErrDirName;
        if (ctErrDirName != "null") {
            auto prevDir = gDirectory;
            ctErrDirs[dsLabel] = TFile::Open(ctErrDirName);
            prevDir->cd();
            if (not ctErrDirs[dsLabel]) {
                cerr << "Unable to open the ct err pdf file " << ctErrDirName << endl;
                return;
            }
        }
        else {
            cerr << "Scheduling computation for ct err for label " << dsLabel << endl;
            ctErrDirs[dsLabel] = nullptr;
        }
    }
    
    // gen tree
    TTree* genTree = nullptr;
    if (needgen) {
        if (extraConfigMap.count("genFile") < 1) {
            cerr << "Missing at least one efficiency but no gen file provided" << endl;
            return;
        }
        auto genFileName = extraConfigMap["genFile"].c_str();
        
        genTree = getTree(genFileName, genTreeName);
        
        if (!genTree) {
            return;
        }
    }
    
    Pwave_Reco_Model_recotag(recoTrees, genTree, massPdfDir, ctEffDirs, angEffDirs, ctErrDirs, usePdgPars, isData, doPlot, numCpu);
}
