// Author Enrico Lusiani

#include "TFile.h"

#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooWorkspace.h"

#include "RooJohnsonLocal.cxx"
#include "RooTagPdf.cxx"

#include "definitions.C"

#include "utils.h"

#include "massPdf_fit.C"
#include "ctErrPdf_fit.C"
#include "ctBg_fit.C"
#include "angBg_fit.C"
#include "extraVars_prefit.C"

void runPrefitAndFit(RooWorkspace* mainWs, RooWorkspace* modelWs, RooDataSet& data, RooCategory& dsCat, int numCpu, bool isData, TDirectory* outDir) {
    
    using namespace RooFit;
    using namespace std;
    
    auto BsCt2DMC = mainWs->var(ctName);
    auto BsCt2DMCErr = mainWs->var(ctErrName);
    auto BscosthetaMC = mainWs->var(cosThetaName);
    auto BscospsiMC = mainWs->var(cosPsiName);
    auto BsphiMC = mainWs->var(phiName);
    auto svmass = mainWs->var(massName);
    auto mistag = mainWs->var(mistagName);
    auto tag = mainWs->cat(tagCatName);
    
    auto sbStatus = mainWs->catfunc(sbStatusName);
    
    auto splitData = data.split(dsCat);
    
    map<string, RooDataSet*> recoDataByLabel;
    for (auto ds: TRangeDynCast<RooDataSet>(splitData)) {
        recoDataByLabel[ds->GetName()] = ds;
    }
    
    auto massPdfDir = massPdf_fit(data, *svmass, isData);
    
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
    
    mainWs->import(*massPdf);
    
    map<string, map<string, RooDataSet*>> sbSplitDataByLabel;
    
    for (auto dsCatEntry: RooCatRange(dsCat)) {
        auto label = dsCatEntry->GetName();
        sbSplitDataByLabel[label] = bkgSigSplitWSB(*recoDataByLabel[label], *mainWs, *sbStatus, isData);
    }
    
    if (isData){
    
        auto rawSbData = (RooDataSet*)getRawSideband(data);
        
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
        
        mainWs->import(*ctBgPdf);
        
        // =================================ANG BKG PDF ================================================================
        
        auto angBkgPdfDir = angBg_fit(*rawSbData);
        
        auto angBkgPdfWs = (RooWorkspace*)angBkgPdfDir->Get("angBgWs");
        
        auto angBgPdf = angBkgPdfWs->pdf("angBg_pdf");
    
//        for (auto varBase: *angBgPdf->getParameters(*recoData)) {
//            auto var = (RooRealVar*)varBase;
//            
//            var->setConstant();
//        }
        
        mainWs->import(*angBgPdf);
        
        angBgPdf->fitTo(data);
    }
    
    // =====================================CT ERR======================================================================
    
    for (auto dsCatEntry: RooCatRange(dsCat)) {
        auto label = dsCatEntry->GetName();
        
        auto ctErrPdfDir = ctErrPdf_fit(sbSplitDataByLabel[label], *BsCt2DMCErr, isData, label);
        
        if (!ctErrPdfDir) {
            cerr << "Unable to fit the ct error pdf" << endl;
            return;
        }
        
        auto ctErrWs = (RooWorkspace*)ctErrPdfDir->Get("ctErr");
        
        mainWs->import(*ctErrWs->function("ctErr_pdf"), 
            RenameAllVariablesExcept(label, ctErrName), 
            RenameAllNodes(label)
        );
    }
    
    // ======================================REFIT MASS+CTERR==============================================================
    
    if (isData) {
        extraVars_prefit(data, *mainWs, dsCat);
    }
    
    for (auto dsCatEntry: RooCatRange(dsCat)) {
        auto label = dsCatEntry->GetName();
        
        auto ctErrPdf = mainWs->pdf(Form("ctErr_pdf_%s", label));
        
        for(auto varBase: *ctErrPdf->getParameters(data)) {
            auto var = (RooRealVar*)varBase;
            
            var->setConstant();
        }
    }
    
    // ====================================================================================================================
    
    modelWs->cat("blindCat")->setLabel("unblind");
    auto model = modelWs->pdf("Angular_Model");
    
    mainWs->import(*model, RecycleConflictNodes());
    model = mainWs->pdf("Angular_Model");
    mainWs->cat("blindCat")->setLabel("unblind");
    
    auto res = model->fitTo(data, NumCPU(numCpu), PrintLevel(3), Verbose(true), Save(true), Offset(true));
    res->Print("V");
    
    RooWorkspace fitWs("fitWs");
    fitWs.import(*model);
    
    fitWs.Write();
    res->Write();
}

void Pwave_toyMC_fitterV2(TString modelFilename, TString dataFileDir, int sampleIndex, int numCpu, bool isData) {

    auto modelFile = TFile::Open(modelFilename);
    auto prefitWs = (RooWorkspace*)modelFile->Get("preFit");
    
    auto mainWs = getMainWs();
    
    auto filename = TString::Format("/data.%d.root", sampleIndex);
    auto dataFile = TFile::Open(dataFileDir + filename);
    auto dataWs = (RooWorkspace*)dataFile->Get("genWs");
    auto data = (RooDataSet*)dataWs->data("Angular_Model_adjData");
    
    auto dsCat = dataWs->cat(dsCatName);
    
    filename = TString::Format("/fit.%d.root", sampleIndex);
    auto outDir = TFile::Open(dataFileDir + filename, "CREATE");
    if (!outDir) {
        clog << "Unable to open file " << dataFileDir << " fit.root for writing" << endl;
        return;
    }
    
    runPrefitAndFit(mainWs, prefitWs, *data, *dsCat, numCpu, isData, outDir);
    
    delete outDir;
}
