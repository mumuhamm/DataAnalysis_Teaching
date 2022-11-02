// Author Enrico Lusiani

#include "TFile.h"
#include "TF1.h"

#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "RooSimultaneous.h"
#include "RooEffProd.h"
#include "RooProduct.h"
#include "RooRealSumPdf.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"

#include "RooJohnsonLocal.cxx"
#include "RooTagPdf.cxx"
#include "RooDelta.cxx"

#include "definitions.C"

std::pair<RooSimultaneous*, RooSimultaneous*> spliceModel(RooWorkspace& ws);

void Pwave_toyMC_generator(RooWorkspace* main, std::map<string, TDirectory*> tagDirs, TString dataDir, int sampleSize, int nGen, int seed) {
    RooRandom::randomGenerator()->SetSeed(seed);
    using namespace RooFit;
    
    auto BsCt2DMC = main->var(ctName);
    auto BsCt2DMCErr = main->var(ctErrName);
    auto BscosthetaMC = main->var(cosThetaName);
    auto BscospsiMC = main->var(cosPsiName);
    auto BsphiMC = main->var(phiName);
    auto svmass = main->var(massName);
    auto mistag = main->var(mistagName);
    auto tag = main->cat(tagCatName);
    
    auto fittedModel = main->pdf("Angular_Model");
    
    // fix model for generation
    auto adjustedModelAndSigBkg = spliceModel(*main);
    auto adjustedModel = adjustedModelAndSigBkg.first;
    auto sigBkgPdf = adjustedModelAndSigBkg.second;
    auto& indexCat = (RooSuperCategory&)adjustedModel->indexCat();
    auto& dsLabel = *(RooAbsCategoryLValue*)indexCat.inputCatList().find("dsLabel")->clone();
    auto& sigOrBkg = *(RooAbsCategoryLValue*)indexCat.inputCatList().find("sigOrBkg")->clone();
    
    auto variables = adjustedModel->getVariables();
    shiftToMCValues(*variables);
    
    // =================================== initialize mistag ================================================
    std::map<string, RooAbsPdf*> tagPdfMap;
    std::map<string, RooAbsPdf*> mistPdfMap;
    for (auto tagDirEntry: tagDirs) {
        auto label = tagDirEntry.first.c_str();
        auto mistDir = tagDirEntry.second;
        
        auto mistWs = (RooWorkspace*)mistDir->Get("mistWs");
        
        tagPdfMap[label] = mistWs->pdf("tagPdf");
        tagPdfMap[label]->SetName(Form("%s_%s", tagPdfMap[label]->GetName(), label));
        auto mistTagPdf = mistWs->pdf("mistTagPdf");
        mistTagPdf->SetName(Form("%s_%s", mistTagPdf->GetName(), label));
        auto mistUntagPdf = new RooDelta(Form("mistUntagPdf_%s", label), "Untagged Pdf", *mistag, 0.5);
        
        // RooCatConditionalPdf is needed for the nested generation which RooSimultaneous does not support 
        // (I suspect it's a bug, I will report it after the end of the analysis)
        mistPdfMap[label] = new RooSimultaneous(Form("mistagPdf_%s", label), "Mistag Pdf", {{"Bs", mistTagPdf}, {"Bsbar", mistTagPdf}, {"untag", mistUntagPdf}}, *tag);
    }
    auto tagPdf = new RooSimultaneous("tagPdf", "tag pdf", tagPdfMap, dsLabel);
    auto mistPdf = new RooSimultaneous("mistPdf", "mistag pdf", mistPdfMap, dsLabel);
    
    // =================================== generate prototypes =========================================================
    
    //RooMsgService::instance().addStream(DEBUG, Topic(Generation));
    
    std::clog << "Generating tag/mistag..." << std::endl;
    auto protoDs = new RooDataSet("protoDs", "proto ds", dsLabel);
    // sqrt(2) to avoid repetition at all
    for (int i = 0; i < sampleSize*2*sqrt(2.); i++) {
        dsLabel.randomize();
        protoDs->add(dsLabel);
    }
    protoDs->table(dsLabel)->Print("V");
    auto protoSigBkg = sigBkgPdf->generate(RooArgList(sigOrBkg), Name("protoSigBkg"), ProtoData(*protoDs));
    protoSigBkg->table(dsLabel)->Print("V");
    protoSigBkg->table(sigOrBkg)->Print("V");
    auto protoTag = tagPdf->generate(RooArgList(*tag), Name("protoTag"), ProtoData(*protoDs));
    protoTag->table(*tag)->Print("V");
    protoTag->table(dsLabel)->Print("V");
    auto protoTagMist = mistPdf->generate(RooArgList(*mistag), Name("protoTagMistag"), ProtoData(*protoTag));
    protoTagMist->table(*tag)->Print("V");
    protoTagMist->table(dsLabel)->Print("V");
    
    protoTagMist->merge(protoSigBkg);
    protoTagMist->table(*tag)->Print("V");
    protoTagMist->table(sigOrBkg)->Print("V");
    protoTagMist->table(dsLabel)->Print("V");
    
    // =================================== initialize generator =======================================================
    
    RooArgSet varToGen(*BsCt2DMC, *BscosthetaMC, *BscospsiMC, *BsphiMC, *BsCt2DMCErr, *svmass);
    auto genSpec = adjustedModel->prepareMultiGen(varToGen, NumEvents(sampleSize), Name("data"), Verbose(true), ProtoData(*protoTagMist, true));
    
    // =================================== generate ===================================================================
    
    for (int iGen = 0; iGen < nGen; iGen++ ) {
        
        std::clog << "Generating main model..." << std::endl;
        
        auto data = adjustedModel->generate(*genSpec);
    
        auto filename = TString::Format("/data.%d.root", iGen);
        auto file = TFile::Open(dataDir + filename, "RECREATE");
        
        RooWorkspace ws("genWs");
        ws.import(*data);
        ws.Write();
    
        delete file;
        delete data;
    }
}


#include "utils.h"

RooAbsArg* selectFirst(const RooArgList& list, const char* name) {
    return ((RooArgList*)list.selectByName(name))->at(0);
}


// the model used for fitting is not apt for generation: efficiencies are applied too soon, which is great for integration but terrible for generation. 
// Furthermore the extra correlation that we have inside the final pdf forces RooFit to use accept/reject for pretty much every variable
// Result is that the original model takes 8h to be INITIALIZED at 1%, which is not viable
// This function applies some transformation to the model to optimize for generation, all while remaining mathematically consistent with the theorical model
std::pair<RooSimultaneous*, RooSimultaneous*> spliceModel(RooWorkspace& ws) {
    using namespace RooFit;
    auto origModel = (RooSimultaneous*)ws.pdf("Angular_Model")->clone("model_clone");
    auto& indexCat = (RooAbsCategoryLValue&)origModel->indexCat();
    
    ws.cat("blindCat")->setLabel("unblind");
    
    auto ct = ws.var("BsCt2DMC");
    auto costheta = ws.var("BscosthetaMC");
    auto cospsi = ws.var("BscospsiMC");
    auto phi = ws.var("BsphiMC");
    auto mistag = ws.var("mistag");
    
    auto sigOrBkg = new RooCategory("sigOrBkg", "");
    sigOrBkg->defineType("sig", 1);
    sigOrBkg->defineType("bkg", 0);
    
    map<string, RooAbsPdf*> pdfStates;
    map<string, RooAbsPdf*> sigBkgPdfStates;
    for (auto indexState: RooCatRange(indexCat)) {
        TString indexStateLabel = indexState->GetName();
        
        // the pdf associated to this sim state
        auto currentPdf = (RooAddPdf*)origModel->getPdf(indexStateLabel);
        
        // signal model (without extended term)
        auto sigModel = (RooProdPdf*)currentPdf->pdfList().at(0);
        auto bkgModel = (RooProdPdf*)currentPdf->pdfList().at(1);
        
        // and now we hack it up
        auto& prodList = sigModel->pdfList();
        
        auto massPdf = (RooAbsPdf*)prodList.find("mass_sig");
        auto ctErrPdf = (RooAbsPdf*)selectFirst(prodList, "ctErr_pdf_*");
        auto pdfDef = (RooRealSumPdf*)selectFirst(prodList, "PDFdefP_*");
        
        auto& funcList = pdfDef->funcList();
        auto& coefList = pdfDef->coefList();
        
        RooAbsReal* angEffFunc = nullptr;
        RooAbsReal* ctEffFunc = nullptr;
        
        RooArgList newFuncList;
        
        for (auto eqObj: funcList) {
            auto eq = dynamic_cast<RooProduct*>(eqObj);
            
            auto eqComponents = eq->components();
            
            RooArgList newComponents;
            
            for (auto compObj: eqComponents) {
                auto comp = (RooAbsReal*)compObj;
                
                TString compName = comp->GetName();
                
                if (compName.Contains("angEff")) {
                    if (angEffFunc and angEffFunc != comp) {
                        throw 1;
                    }
                    angEffFunc = comp;
                    continue;
                }
                if (compName.Contains("ctEff")) {
                    if (ctEffFunc and ctEffFunc != comp) {
                        throw 1;
                    }
                    ctEffFunc = comp;
                    continue;
                }
                newComponents.add(*comp);
            }
            newFuncList.add(*new RooProduct(eq->GetName(), eq->GetTitle(), newComponents));
        }
        auto newPdfDef = new RooRealSumPdf(pdfDef->GetName(), pdfDef->GetTitle(), newFuncList, coefList);
        
        RooArgList newProdList(*massPdf, *ctErrPdf);
        auto newNormModel = new RooProdPdf(sigModel->GetName(), sigModel->GetTitle(), newProdList, Conditional(*newPdfDef, RooArgSet(*ct, *costheta, *cospsi, *phi)));
        
        auto angEffTF = angEffFunc->asTF(RooArgList(*costheta, *cospsi, *phi));
        auto angEffMax = angEffTF->GetMaximum();
        std::clog << "ang eff maximum is " << angEffMax << endl;
        
        auto ctEffTF = ctEffFunc->asTF(RooArgList(*ct));
        auto ctEffMax = ctEffTF->GetMaximum();
        std::clog << "ct eff maximum is " << ctEffMax << endl;
        
        auto fullEff = new RooProduct(Form("eff_%s", indexStateLabel.Data()), "", RooArgList(*angEffFunc, *ctEffFunc, RooConst(0.9/(ctEffMax*angEffMax))));
        
        map<string, RooAbsPdf*> models;
        
        models["sig"] = new RooEffProd(Form("%sSig_wCtEff", currentPdf->GetName()), currentPdf->GetTitle(), 
            *newNormModel, *fullEff
        );
        models["bkg"] = bkgModel;
        
        pdfStates[indexStateLabel.Data()] = new RooSimultaneous(currentPdf->GetName(), currentPdf->GetTitle(), models, *sigOrBkg);
        
        
        auto sigNum = (RooAbsReal*)currentPdf->coefList().at(0);
        auto bkgNum = (RooAbsReal*)currentPdf->coefList().at(1);
        
        auto bkgFrac = bkgNum->getVal()/(sigNum->getVal() + bkgNum->getVal());
        auto bkgFracVar = new RooRealVar(Form("bkgFrac_%s", indexStateLabel.Data()), "", bkgFrac);
        
        sigBkgPdfStates[indexStateLabel.Data()] = new RooGenericPdf(Form("sigVsBkg_%s", indexStateLabel.Data()), "", "sigOrBkg == sigOrBkg::bkg?@1:(1-@1)", RooArgList(*sigOrBkg, *bkgFracVar));
    }
    
    auto newFullModel = new RooSimultaneous("Angular_Model_adj", origModel->GetTitle(), pdfStates, indexCat);
    auto sigBkgPdf = new RooSimultaneous("sisVsBkg", "", sigBkgPdfStates, indexCat);
    
    return {newFullModel, sigBkgPdf};
}

void Pwave_toyMC_generator(TString filename, TString mistagConfigFileName, TString dataDir, int sampleSize, int nGen, int seed) {
    auto file = TFile::Open(filename);
    
    auto mainWs = (RooWorkspace*)file->Get("postFit");
    
    ifstream mistConfigFile(mistagConfigFileName);
    string line;
    
    std::map<string, TDirectory*> mistFiles;
    
    while(getline(mistConfigFile, line)) {
        stringstream lineReader(line);
        
        string dsLabel;
        lineReader >> dsLabel;
        
        TString mistagDirName;
        lineReader >> mistagDirName;
        
        auto prevDir = gDirectory;
        mistFiles[dsLabel] = TFile::Open(mistagDirName);
        prevDir->cd();
    }
    
    Pwave_toyMC_generator(mainWs, mistFiles, dataDir, sampleSize, nGen, seed);
}





