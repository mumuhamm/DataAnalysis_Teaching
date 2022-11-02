//Author: Md. Alibordi


#pragma once

#include "RooDataSet.h"
#include "RooFormulaVar.h"
#include "RooAbsPdf.h"

#include "definitions.C"

#include "RooJohnsonLocal.cxx"

std::map<string, RooDataSet*> bkgSigSplitWSB(RooDataSet& data, RooWorkspace& massWs, RooAbsCategory& sbStatus, bool isData) {
    if (!isData) {
        return {{"sig", &data}, {"full", &data}};
    }
    
    auto svmass = (RooRealVar*)data.get()->find(massName);
    
    auto massPdf = massWs.pdf("mass_pdf");
    massPdf->attachDataSet(data);
    //massPdf->fitTo(data);
    
    auto sbSum = data.sumEntries(nullptr, "sbLeft,sbRight");
    auto centSum = data.sumEntries(nullptr, "center");

    auto& mass_sigPdf = *massWs.pdf("mass_sig");
    auto& mass_bkgPdf = *massWs.pdf("mass_bkg");

    auto bkgInSbInt = mass_bkgPdf.createIntegral(*svmass, *svmass, "sbLeft,sbRight")->getVal();
    auto bkgInCentInt = mass_bkgPdf.createIntegral(*svmass, *svmass, "center")->getVal();

    auto sigInSbInt = mass_sigPdf.createIntegral(*svmass, *svmass, "sbLeft,sbRight")->getVal();
    auto sigInCentInt = mass_sigPdf.createIntegral(*svmass, *svmass, "center")->getVal();
    cerr << bkgInSbInt << " " << bkgInCentInt << endl;
    cerr << sigInSbInt << " " << sigInCentInt << endl;
    
//    if (bkgInSbInt*sbSum > bkgInCentInt*centSum) {
//        bkgInSbInt = -bkgInSbInt;
//        bkgInCentInt = -bkgInCentInt;
//    }
//    if (sigInCentInt*centSum > sigInSbInt*sbSum) {
//        sigInSbInt = -sigInSbInt;
//        sigInCentInt = -sigInCentInt;
//    }

    RooFormulaVar sigWeight("sig_w", "", Form("sbStatus == sbStatus::sb?-%g:%g", bkgInCentInt, bkgInSbInt), sbStatus);
    RooFormulaVar bkgWeight("bkg_w", "", Form("sbStatus == sbStatus::sb?%g:-%g", sigInCentInt, sigInSbInt), sbStatus);
    
    data.addColumn(sigWeight);
    data.addColumn(bkgWeight);
    
    auto sigData = new RooDataSet(Form("%s_sig", data.GetName()), data.GetTitle(), &data, *data.get(), 0, "sig_w");
    auto bkgData = new RooDataSet(Form("%s_bkg", data.GetName()), data.GetTitle(), &data, *data.get(), 0, "bkg_w");
    
    return {{"sig", sigData}, {"bkg", bkgData}, {"full", &data}};
}
