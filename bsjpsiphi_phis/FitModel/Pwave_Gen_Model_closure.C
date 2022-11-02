// Authors: Enrico Lusiani
#include "TTree.h"
#include "TFile.h"
#include "TString.h"

#include "RooProdPdf.h"
#include "RooRealSumPdf.h"
#include "RooGenericPdf.h"
#include "RooGaussModel.h"
#include "Roo1DTable.h"

#include "utils.h"

#define GENONLY
#include "definitions.C"

#include "Pwave_Model.C"

#include "Pwave_Model_plot.C"

void Pwave_Gen_Model_closure(TTree* genTree, bool plotAngFlag = true) {
    
    using namespace RooFit;
    using namespace std;
    // ====================================VARIABLE AND DATASET DEFINITIONS=============================================
    
    auto mainWs = getMainWs();
    
    auto BsCt2DMC = mainWs->var(ctName);
    auto BscosthetaMC = mainWs->var(cosThetaName);
    auto BscospsiMC = mainWs->var(cosPsiName);
    auto BsphiMC = mainWs->var(phiName);
    auto tag = mainWs->cat(tagCatName);

    mainWs->var(mistagName)->setVal(0);
    mainWs->var(mistagName)->setConstant();    

    auto genData = new RooDataSet("genData", "genData", RooArgSet(*BscosthetaMC, *BscospsiMC, *BsphiMC, *BsCt2DMC, *tag),
                                      Import(*genTree));
                                      

    //====================================Fit===============================================================================
    
    Pwave_Model(*genData, mainWs, true, false, false, plotAngFlag, 8);
}

void Pwave_Gen_Model_closure(TString genFileName, bool plotAngFlag = true) {
    // gen tree
    auto genTree = getTree(genFileName, genTreeName);
    if (!genTree) {
        return;
    }
    
    Pwave_Gen_Model_closure(genTree, plotAngFlag);
}
