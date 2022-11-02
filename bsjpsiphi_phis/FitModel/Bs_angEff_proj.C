// Authors: Alessio Boletti, Enrico Lusiani
#pragma once

#include "TTree.h"
#include "TString.h"
#include "TH3.h"
#include "TMath.h"
#include "TStyle.h"
#include "TNtuple.h"

#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "RooAddition.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooLegendre.h"
#include "RooFormulaVar.h"
#include "RooProduct.h"
#include "RooDataSet.h"

#include "utils.h"

#include "definitions.C"

#include "Bs_angEff_slicePlot.C"

void plotProjections(TH3*, TString);
RooAbsReal* buildAngularFunction(TH3* recoAngHist, TH3* genAngHist, const RooArgSet& vars, int maxOrder);
TH3* plotAngFunc(RooAbsReal* angFunc, RooRealVar&, RooRealVar&, RooRealVar&, int maxOrder, int xbins, int ybins, int zbins, TString);

TDirectory* Bs_angEff_proj(RooDataSet& recoData, RooDataSet& genData, const RooArgSet& vars, TString suffix = "", int maxOrder = 6, int xbins = 70, int ybins = 70, int zbins = 30) {
    using namespace std;
    using namespace RooFit;
    
    std::clog << "Projecting the angular efficiency" << endl;
    
    if (maxOrder < 0)
        return nullptr;
    if (ybins < 1)
        ybins = xbins;
    if (zbins < 1)
        zbins = xbins;
    if (xbins < 1)
        return nullptr;
        
    auto BscosthetaMC = (RooRealVar*)vars.find("BscosthetaMC");
    auto BscospsiMC = (RooRealVar*)vars.find("BscospsiMC");
    auto BsphiMC = (RooRealVar*)vars.find("BsphiMC");

    auto genAngHist = (TH3D *)genData.createHistogram("genHist" + suffix, *BscosthetaMC, Binning(xbins, -1, 1), YVar(*BscospsiMC, Binning(ybins, -1, 1)),
                                                      ZVar(*BsphiMC, Binning(zbins, -TMath::Pi(), TMath::Pi())));
    auto recoAngHist = (TH3D *)recoData.createHistogram("recoHist" + suffix, *BscosthetaMC, Binning(xbins, -1, 1), YVar(*BscospsiMC, Binning(ybins, -1, 1)),
                                                     ZVar(*BsphiMC, Binning(zbins, -TMath::Pi(), TMath::Pi())));
    
    
    auto prevDir = gDirectory;
    TString dirPath = "outputs/";
    if (ensureDir(dirPath)) return nullptr;
    auto angEffDir = TFile::Open(dirPath + "angEffDir_" + suffix + ".root", "RECREATE");
    
    TH3D *angEff = (TH3D *)recoAngHist->Clone("angEff" + suffix);
    angEff->Sumw2();
    angEff->Divide(genAngHist);
    
    plotProjections(angEff, suffix);
    
    auto angFunc = (RooAddition*)buildAngularFunction(recoAngHist, genAngHist, vars, maxOrder);
    
    auto angEffFuncHist = plotAngFunc(angFunc, *BscosthetaMC, *BscospsiMC, *BsphiMC, maxOrder, xbins, ybins, zbins, suffix);
    plotSlices(recoAngHist, genAngHist, *angFunc, RooArgSet(*BscosthetaMC, *BscospsiMC, *BsphiMC), maxOrder, xbins, ybins, zbins, suffix);
    
    auto ws = new RooWorkspace("ang");
    ws->import(*angFunc, Silence());
    
    ws->Write();
    
    angEffDir->Write();
    angEffDir->ReOpen("READ");
    
    prevDir->cd();
    auto diff2 = (TH3*)angEff->Clone("diff" + suffix);
    diff2->Add(angEffFuncHist, -1);
    diff2->Multiply(diff2);
    auto chi2 = diff2->Integral();
    std::clog << "avg eff: " << angEff->Integral()/(xbins*ybins*zbins) << endl;
    std::clog << "avg proj: " << angEffFuncHist->Integral()/(xbins*ybins*zbins) << endl;
    std::clog << "func vs hist: " << chi2 << endl;
    std::clog << "average chi2 : " << chi2/(xbins*ybins*zbins) << endl; 
    std::clog << "components: " << angFunc->list().getSize() << endl;
    
    return angEffDir;
}

void plotProjections(TH3* angEff, TString suffix) {

    TCanvas c1("c1", "c1", 0, 0, 1200, 800);
    
    auto proj2d1 = angEff->Project3D("xy");
    auto proj2d2 = angEff->Project3D("xz");
    auto proj2d3 = angEff->Project3D("yz");
    proj2d1->SetTitle("Efficiency cos(#theta_{T})/cos(#psi_{T}) projection ");
    proj2d2->SetTitle("Efficiency cos(#theta_{T})/#phi projection");
    proj2d3->SetTitle("Efficiency cos(#psi_{T})/#phi projection");
    c1.Divide(3, 1);
    c1.cd(1);
    proj2d1->GetXaxis()->SetTitleOffset(1.4);
    proj2d1->GetYaxis()->SetTitleOffset(2);
    proj2d1->SetMinimum(0.0);
    proj2d1->Draw("SURF3");
    c1.cd(2);
    proj2d2->GetXaxis()->SetTitleOffset(1.4);
    proj2d2->GetYaxis()->SetTitleOffset(2);https://old.reddit.com/r/frackinuniverse/comments/dp7dyt/any_tutorials_on_storage_automation/
    proj2d2->SetMinimum(0.0);
    proj2d2->Draw("SURF3");
    c1.cd(3);
    proj2d3->GetXaxis()->SetTitleOffset(1.4);
    proj2d3->GetYaxis()->SetTitleOffset(2);
    proj2d3->SetMinimum(0.0);
    proj2d3->Draw("SURF3");
    
    saveCanvas(c1, "angEffPlot_" + suffix);
}

RooAbsReal* buildAngularFunction(TH3* recoAngHist, TH3* genAngHist, const RooArgSet& vars, int maxOrder) {
    using namespace RooFit;
    using namespace std;

    vector<RooRealVar *> factors;
    vector<double> proj;
    vector<RooProduct *> vectFunc;
    
    auto BsphiMC = (RooRealVar*)vars.find(phiName);
    auto BscospsiMC = (RooRealVar*)vars.find(cosPsiName);
    auto BscosthetaMC = (RooRealVar*)vars.find(cosThetaName);
    double avg = recoAngHist->Integral() / genAngHist->Integral();

    RooArgList facList;
    RooArgList funList;

    // FIXME the vectors are useless
    for (int xOrder = 0; xOrder <= maxOrder; ++xOrder) {
        for (int yOrder = 0; yOrder <= maxOrder; ++yOrder) {
            for (int zOrder = -1 * TMath::Min(xOrder, yOrder); zOrder <= TMath::Min(xOrder, yOrder); ++zOrder) {

                // vector of coefficients for the function basis
                factors.push_back(new RooRealVar(Form("l%i_k%i_m%i", xOrder, yOrder, zOrder), Form("l%i_k%i_m%i", xOrder, yOrder, zOrder), 0));

                RooArgList prodList;

                // phi terms by trigonometric polynomials (degree zOrder)
                if (zOrder > 0) {
                    prodList.add(*new RooFormulaVar(Form("funcPoly%i_%i_%i", xOrder, yOrder, zOrder), Form("funcPoly%i_%i_%i", xOrder, yOrder, zOrder),
                                                             Form("cos(%i*BsphiMC)", zOrder), *BsphiMC));
                }
                if (zOrder < 0) {
                    prodList.add(*new RooFormulaVar(Form("funcPoly%i_%i_%i", xOrder, yOrder, zOrder), Form("funcPoly%i_%i_%i", xOrder, yOrder, zOrder),
                                                             Form("sin(%i*BsphiMC)", -1 * zOrder), *BsphiMC));
                }

                // costhetaT terms by associated Legendre polynomials (degree l=xOrder m=zOrder)
                prodList.add(*new RooLegendre(Form("funcLegctK%i_%i_%i", xOrder, yOrder, zOrder),
                                                               Form("funcLegctK%i_%i_%i", xOrder, yOrder, zOrder), *BscosthetaMC, xOrder, abs(zOrder)));

                // cospsiT terms by associated Legendre polynomials (degree l=yOrder m=zOrder)
                prodList.add(*new RooLegendre(Form("funcLegctL%i_%i_%i", xOrder, yOrder, zOrder),
                                                             Form("funcLegctL%i_%i_%i", xOrder, yOrder, zOrder), *BscospsiMC, yOrder, abs(zOrder)));

                // build member of the basis of 3D functions
                vectFunc.push_back(new RooProduct(Form("func%i_%i_%i", xOrder, yOrder, zOrder), Form("func%i_%i_%i", xOrder, yOrder, zOrder), prodList));

                // coefficients values to be filled later
                proj.push_back(0);

                // preparation of RooArgList objects
                funList.add(*vectFunc.back());
                facList.add(*factors.back());
            }
        }
    }

    std::clog << "Number of parameters used in angular efficiency: " << factors.size() << endl;
    
    int iOrder = -1;
//    auto prevDir = gDirectory;
//    auto file = TFile::Open("coefs.root", "RECREATE");
//    TNtuple coefHist("coefHist", "coef hist", "xOrd:yOrd:zOrd:coef");

    // loop over the coefficients
    for (int xOrder = 0; xOrder <= maxOrder; ++xOrder) {
        for (int yOrder = 0; yOrder <= maxOrder; ++yOrder) {
            for (int zOrder = -1 * TMath::Min(xOrder, yOrder); zOrder <= TMath::Min(xOrder, yOrder); ++zOrder) {

                ++iOrder;

                // project the binned efficiency on the [iOrder] function
                for (int xBin = 1; xBin <= genAngHist->GetNbinsX(); ++xBin) {
                    for (int yBin = 1; yBin <= genAngHist->GetNbinsY(); ++yBin) {
                        for (int zBin = 1; zBin <= genAngHist->GetNbinsZ(); ++zBin) {

                            BscosthetaMC->setVal(genAngHist->GetXaxis()->GetBinCenter(xBin));
                            BscospsiMC->setVal(genAngHist->GetYaxis()->GetBinCenter(yBin));
                            BsphiMC->setVal(genAngHist->GetZaxis()->GetBinCenter(zBin));

                            // contribution of one bin
                            if (genAngHist->GetBinContent(xBin, yBin, zBin) > 0)
                                proj[iOrder] += (recoAngHist->GetBinContent(xBin, yBin, zBin) / genAngHist->GetBinContent(xBin, yBin, zBin)*
                                                 genAngHist->GetXaxis()->GetBinWidth(xBin) * genAngHist->GetYaxis()->GetBinWidth(yBin) *
                                                 genAngHist->GetZaxis()->GetBinWidth(zBin) * vectFunc[iOrder]->getVal(vars));
                            else
                                proj[iOrder] += (avg * genAngHist->GetXaxis()->GetBinWidth(xBin) * genAngHist->GetYaxis()->GetBinWidth(yBin) *
                                                 genAngHist->GetZaxis()->GetBinWidth(zBin) * vectFunc[iOrder]->getVal(vars));
                        }
                    }
                }

                // normalization of 0-degree trigonometric polynomial differs by a factor 2
                if (zOrder == 0)
                    proj[iOrder] = proj[iOrder] / 2.0;

                // set coefficient value, normalised
                factors[iOrder]->setVal(proj[iOrder] * (2 * xOrder + 1) * TMath::Factorial(xOrder - abs(zOrder)) / 2 /
                                        TMath::Factorial(xOrder + abs(zOrder)) // associated legendre poly
                                        * (2 * yOrder + 1) * TMath::Factorial(yOrder - abs(zOrder)) / 2 / TMath::Factorial(yOrder + abs(zOrder)) /
                                        TMath::Pi() // trigonometric polynomial
                );
//                coefHist.Fill(xOrder, yOrder, zOrder, 
//                    proj[iOrder] * sqrt((2 * xOrder + 1) * TMath::Factorial(xOrder - abs(zOrder)) / 2 /
//                        TMath::Factorial(xOrder + abs(zOrder)) // associated legendre poly
//                        * (2 * yOrder + 1) * TMath::Factorial(yOrder - abs(zOrder)) / 2 / TMath::Factorial(yOrder + abs(zOrder)) /
//                        TMath::Pi())
//                );

                std::clog << xOrder << " " << yOrder << " " << zOrder << "\t" << iOrder << "\t" << factors[iOrder]->getValV() << endl;
            }
        }
    }
//    coefHist.Write();
//    file->Close();
//    delete file;
//    prevDir->cd();
    
    RooArgList funListSkimmed;
    RooArgList facListSkimmed;
    iOrder = -1;
    for (int xOrder = 0; xOrder <= maxOrder; ++xOrder) {
        for (int yOrder = 0; yOrder <= maxOrder; ++yOrder) {
            for (int zOrder = -1 * TMath::Min(xOrder, yOrder); zOrder <= TMath::Min(xOrder, yOrder); ++zOrder) {
                ++iOrder;
                auto norm = proj[iOrder] * sqrt((2 * xOrder + 1) * TMath::Factorial(xOrder - abs(zOrder)) / 2 /
                        TMath::Factorial(xOrder + abs(zOrder)) // associated legendre poly
                        * (2 * yOrder + 1) * TMath::Factorial(yOrder - abs(zOrder)) / 2 / TMath::Factorial(yOrder + abs(zOrder)) /
                        TMath::Pi());
                // arbitrary
                if (log10(abs(norm)) > -3.375) {
                    facListSkimmed.add(*facList.at(iOrder));
                    funListSkimmed.add(*funList.at(iOrder));
                }
            }
        }
    }

    // Sum function
    RooAddition *projectedFuncSkimmed = new RooAddition("angEffFunc", "angular efficiency function", funListSkimmed, facListSkimmed);
    RooAddition *projectedFunc = new RooAddition("angEffFunc", "angular efficiency function", funList, facList);
    
    return projectedFuncSkimmed;    // Sum function
}

TH3* plotAngFunc(RooAbsReal* angFunc, RooRealVar& x, RooRealVar& y, RooRealVar& z, int maxOrder, int xbins, int ybins, int zbins, TString suffix) {

    using namespace RooFit;
    auto h3_xyz = (TH3 *)angFunc->createHistogram("angEffFuncHist" + suffix, x, Binning(xbins), YVar(y, Binning(ybins)), ZVar(z, Binning(zbins)), Scaling(false));
    auto h3_x = h3_xyz->ProjectionX();
    auto h3_y = h3_xyz->ProjectionY();
    auto h3_z = h3_xyz->ProjectionZ();
    auto h3_xy = h3_xyz->Project3D("xy");
    auto h3_xz = h3_xyz->Project3D("xz");
    auto h3_yz = h3_xyz->Project3D("yz");
    
    gStyle->SetOptStat(0);

    // 2D projections
    TCanvas c3new("c3new", "c3new", 1200, 800);
    h3_xy->SetTitle("Efficiency cos(#theta_{T})/cos(#psi_{T}) projected function");
    h3_xz->SetTitle("Efficiency cos(#theta_{T})/#phi projected function");
    h3_yz->SetTitle("Efficiency cos(#psi_{T})/#phi projected function");
    h3_xy->GetXaxis()->SetTitleOffset(1.4);
    h3_xy->GetYaxis()->SetTitleOffset(2);
    h3_xz->GetXaxis()->SetTitleOffset(1.4);
    h3_xz->GetYaxis()->SetTitleOffset(2);
    h3_yz->GetXaxis()->SetTitleOffset(1.4);
    h3_yz->GetYaxis()->SetTitleOffset(2);
    h3_xy->SetMinimum(0.0);
    h3_xz->SetMinimum(0.0);
    h3_yz->SetMinimum(0.0);
    c3new.Divide(3, 1);
    c3new.cd(1);
    h3_xy->Draw("SURF3");
    c3new.cd(2);
    h3_xz->Draw("SURF3");
    c3new.cd(3);
    h3_yz->Draw("SURF3");
    saveCanvas(c3new, Form("EffProjectionFunction_%i_%i_%i_2DProj_SpH%iOrder" + suffix + ".pdf", xbins, ybins, zbins, maxOrder));
    
    return h3_xyz;
}

TDirectory* Bs_angEff_proj(TString recoFileName, TString genFileName, TString suffix = "", int maxOrder = 6, int xbins = 70, int ybins = 70, int zbins = 30) {
    // reco tree
    auto recoTree = getTree(recoFileName, recoTreeName);
    if (!recoTree) {
        return nullptr;
    }
    
    // gen tree
    auto genTree = getTree(genFileName, genTreeName);
    if (!genTree) {
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
    
    
    Float_t costhetagen, cospsigen, phigen;
    genTree->SetBranchAddress(cosThetaNameGen, &costhetagen);
    genTree->SetBranchAddress(cosPsiNameGen, &cospsigen);
    genTree->SetBranchAddress(phiNameGen, &phigen);
    
    RooArgSet angVars(*BscosthetaMC, *BscospsiMC, *BsphiMC);
    RooDataSet genData("genData", "GEN distribution", angVars);
    auto genEntries = genTree->GetEntries();
    auto step20 = genEntries / 20;
    
    std::clog << "Importing gen dataset" << endl;
    for (int iCand = 0; iCand < genEntries; ++iCand) {
        if (iCand % step20 == 0) std::clog << iCand << "/" << genEntries << endl;
        genTree->GetEntry(iCand);
        
        BscosthetaMC->setVal(costhetagen);
        BscospsiMC->setVal(cospsigen);
        BsphiMC->setVal(phigen);
        
        genData.add(angVars);
    }
    
    return Bs_angEff_proj(recoData, genData, angVars, suffix, maxOrder, xbins, ybins, zbins);
}
