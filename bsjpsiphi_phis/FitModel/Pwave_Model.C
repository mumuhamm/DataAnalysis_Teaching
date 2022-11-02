// Author Md Alibordi, Giacomo Fedi
#pragma once

#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH2.h"
#include "TStyle.h"

#include "RooProdPdf.h"
#include "RooRealSumPdf.h"
#include "RooGenericPdf.h"
#include "RooGaussModel.h"
#include "RooTruthModel.h"
#include "RooProduct.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooUnblindPrecision.h"

#include "utils.h"

#include "definitions.C"

#include "Pwave_Model_plot.C"
#include "RooTagPdf.cxx"

RooAbsReal& thisCoefIs0 = RooFit::RooConst(0);
RooAbsReal& thisCoefIs1 = RooFit::RooConst(1);

map<int, RooAbsReal*> getFunNT(RooRealVar& costheta, RooRealVar& cospsi, RooRealVar& phi);
map<int, map<int, RooAbsReal*>> getCoefs(RooAbsReal& phi_s, RooAbsReal& lambda, RooCategory& tag, RooRealVar& mistag, RooRealVar& deltaPe, RooRealVar& deltaPa, RooRealVar& deltaSPe);

void Pwave_Model(RooDataSet& data, RooWorkspace* mainWs, bool isClosure, bool pdgPars, bool isData, bool plotAngFlag = true, int numCpu = 2) {
    using namespace std;
    using namespace RooFit;
    
    std::cerr << "Fitting the full model" << endl;
    
    auto BsCt2DMC = mainWs->var(ctName);
    auto BsCt2DMCErr = mainWs->var(ctErrName);
    auto BscosthetaMC = mainWs->var(cosThetaName);
    auto BscospsiMC = mainWs->var(cosPsiName);
    auto BsphiMC = mainWs->var(phiName);
    auto svmass = mainWs->var(massName);
    auto mistag = mainWs->var(mistagName);
    auto tag = mainWs->cat(tagCatName);
    
    auto dsCat = mainWs->cat(dsCatName);

    float A_0_default = pdgPars ? 0.52 : 0.6;
    float A_pe_default = pdgPars ? 0.25 : 0.15;
    float A_S_default = 0.01;
    float deltaPa_default = pdgPars ? 3 : 2.5;
    float deltaPe_default = pdgPars ? 3 : -0.15;
    float phi_s_default = -0.5;
    float dm_default = isData ? 582 : 592;
    float dGam_default = pdgPars ? 3 : 0;
    float dGam_min = pdgPars ? 0 : -1;
    float ctau_default = 0.04413;
    float lambda_default = 1.0;

    // Amplitudes
    RooRealVar *A_0 = new RooRealVar("A_0", "|A_{0}|^{2}", A_0_default, 0.35, 0.7); //DG0 0.600625 -- DG!=0 0.522729
    RooRealVar *A_pe = new RooRealVar("A_pe", "|A_{#perp}|^{2}", A_pe_default, 0., 0.35); //DG0 0.16 -- DG!=0 0.246016
    RooRealVar *A_S = new RooRealVar("A_S", "|A_{S}|^{2}", A_S_default, 0., 0.1); // only s-wave

    // Strong Phases
    RooRealVar *deltaPa = new RooRealVar("deltaPa", "#delta_{#parallel}", deltaPa_default, -2*TMath::Pi(), 2*TMath::Pi()); //DG0 2.5 -- DG!=0 3.14
    RooRealVar *deltaPe = new RooRealVar("deltaPe", "#delta_{#perp}", deltaPe_default, -2*TMath::Pi(), 2*TMath::Pi()); //DG0 -0.17 -- DG!=0 3.14
    RooRealVar *deltaSPe = new RooRealVar("deltaSPe", "#delta_{S}", 0., -2*TMath::Pi(), 2*TMath::Pi()); // only s-wave

    //Physics Parameters
    RooRealVar *phi_s = new RooRealVar("phi_s", "#phi_{s}", phi_s_default, -1, 1); // -0.04
    RooRealVar *dm = new RooRealVar("dm", "#Deltam_{s}", dm_default, 500, 700, "cm^{-1}"); // 593.74409
    RooRealVar *dGam = new RooRealVar("dGam", "#Delta#Gammma_{s}", dGam_default, dGam_min, 10, "cm^{-1}"); //DG0 0 -- DG!=0 unknowm (maybe 2)
    RooRealVar *ctau = new RooRealVar("ctau", "c#tau", ctau_default, 0.040, 0.050, "cm"); //DG0 0.04413  -- DG!=0 unknown
    RooRealVar *lambda = new RooRealVar("lambda", "|#lambda|", lambda_default, 0.0, 2.0);

    auto blindFlag = isData;
    // Blind Parameters
    Bool_t isBlind(blindFlag); // blind flag
    
    TString blindString("BlindStringData");
    TString blind("blind"), unblind("unblind");
    
    RooCategory blindCat("blindCat","blind state Category");
    blindCat.defineType(unblind, 0);
    blindCat.defineType(blind, 1);
    
    if(isBlind) blindCat.setLabel(blind);
    else        blindCat.setLabel(unblind);

    RooUnblindPrecision *A_0_ub = new RooUnblindPrecision("A_0_ub", "|A_{0}|^{2} (unblind)", blindString, A_0_default, 0.02, *A_0, blindCat, 0);
    RooUnblindPrecision *A_pe_ub = new RooUnblindPrecision("A_pe_ub", "|A_{#perp}|^{2} (unblind)", blindString, A_pe_default, 0.02, *A_pe, blindCat, 0);
    RooUnblindPrecision *A_S_ub = new RooUnblindPrecision("A_S_ub", "|A_{S}|^{2} (unblind)", blindString, A_S_default, 0.002, *A_S, blindCat, 0);

    RooUnblindPrecision *phi_s_ub = new RooUnblindPrecision("phi_s_ub", "#phi_{s} (unblind)", blindString, phi_s_default, 0.1, *phi_s, blindCat, 0);
    RooUnblindPrecision *dGam_ub = new RooUnblindPrecision("dGam_ub", "#Delta#Gammma_{s} (unblind)", blindString, dGam_default, 0.5, *dGam, blindCat, 0);
    RooUnblindPrecision *ctau_ub = new RooUnblindPrecision("ctau_ub", "c#tau (unblind)", blindString, ctau_default, 0.001, *ctau, blindCat, 0);
    RooUnblindPrecision *lambda_ub = new RooUnblindPrecision("lambda_ub", "|#lambda| (unblind)", blindString, lambda_default, 0.2, *lambda, blindCat, 0);

    
    RooFormulaVar *A_pa = new RooFormulaVar("A_pa", "|A_{#parallel}|^{2}", "1-@0-@1", RooArgList(*A_0_ub, *A_pe_ub));

    map<int, RooAbsReal*> A_N;
    A_N[1] = A_0_ub;
    A_N[2] = A_pa;
    A_N[3] = A_pe_ub;
    A_N[4] = new RooFormulaVar("A_4", "A_{4}", "sqrt(@0)*sqrt(@1)", RooArgList(*A_pa, *A_pe_ub));
    A_N[5] = new RooFormulaVar("A_5", "|A_{0}||A_{pa}|", "sqrt(@0)*sqrt(@1)", RooArgList(*A_0_ub, *A_pa));
    A_N[6] = new RooFormulaVar("A_6", "A_6", "sqrt(@0)*sqrt(@1)", RooArgList(*A_0_ub, *A_pe_ub));
    A_N[7] = A_S_ub;
    A_N[8] = new RooFormulaVar("A_8", "A_8", "sqrt(@0)*sqrt(@1)", RooArgList(*A_S_ub, *A_pa));
    A_N[9] = new RooFormulaVar("A_9", "|A_9", "sqrt(@0)*sqrt(@1)", RooArgList(*A_S_ub, *A_pe_ub));
    A_N[10] = new RooFormulaVar("A_10", "A_10", "sqrt(@0)*sqrt(@1)", RooArgList(*A_S_ub, *A_0_ub));


    map<int, RooFormulaVar*> basis;
    basis[1] = new RooFormulaVar("coshGBasis", "exp(-@0/@1)*cosh(@0*@2/2)", RooArgList(*BsCt2DMC, *ctau_ub, *dGam_ub));
    basis[2] = new RooFormulaVar("sinhGBasis", "exp(-@0/@1)*sinh(@0*@2/2)", RooArgList(*BsCt2DMC, *ctau_ub, *dGam_ub));
    basis[3] = new RooFormulaVar("cosGBasis", "exp(-@0/@1)*cos(@0*@2)", RooArgList(*BsCt2DMC, *ctau_ub, *dm));
    basis[4] = new RooFormulaVar("sinGBasis", "exp(-@0/@1)*sin(@0*@2)", RooArgList(*BsCt2DMC, *ctau_ub, *dm));
    
    map<string, RooResolutionModel*> ctRes;
    if (!isClosure) {
        for (auto dsCatLabel: RooCatRange(*dsCat)) {
            RooRealVar *kappa = mainWs->var(Form("kappa%s", dsCatLabel->GetName()));
            
            ctRes[dsCatLabel->GetName()] = new RooGaussModel(Form("gm_%s", dsCatLabel->GetName()), "gauss model scaled bt per-event error", *BsCt2DMC, RooConst(0), *BsCt2DMCErr, *kappa);
        }
    }
    else {
        ctRes["gen"] = new RooTruthModel("tm", "truth", *BsCt2DMC);
    }
    map<int, map<string, RooAbsReal*>> convBasis;
    for(int iBasis = 1; iBasis <= 4; iBasis++) {
        for (auto& resModel: ctRes) {
            auto dsLabel = resModel.first;
            auto model = resModel.second;
        
            convBasis[iBasis][dsLabel] = model->convolution(basis[iBasis], BsCt2DMC);
        }
    }

    auto funNT = getFunNT(*BscosthetaMC, *BscospsiMC, *BsphiMC);

    auto coefs = getCoefs(*phi_s_ub, *lambda_ub, *tag, *mistag, *deltaPe, *deltaPa, *deltaSPe);
    
    map<int, map<int, map<string, RooAbsReal*>>> eqs;
    for (int iEq = 1; iEq <= 10; iEq++) {
        auto funT = funNT[iEq];
        
        for (int iBasis = 1; iBasis <= 4; iBasis++) {
            auto coef = coefs[iEq][iBasis];
            
            if (coef == &thisCoefIs0) {
                continue;
            }
            RooArgSet prodFn(*funT);
            if (coef != &thisCoefIs1) {
                prodFn.add(*coef);
            }
            
            if (not isClosure) {
                for (auto dsCatLabel: RooCatRange(*dsCat)) {
                
                    auto basis = convBasis[iBasis][dsCatLabel->GetName()];
                    
                    RooArgSet allFn(prodFn);
                    allFn.add(*basis);
                    allFn.add(*mainWs->function(Form("angEffFunc_%s", dsCatLabel->GetName())));
                    allFn.add(*mainWs->function(Form("ctEffFunc_%s", dsCatLabel->GetName())));
                    
                    eqs[iEq][iBasis][dsCatLabel->GetName()] = new RooProduct(Form("eq%d%d_%s", iEq, iBasis, dsCatLabel->GetName()), "amp0*g1", allFn);
                }
            }
            else {
                auto basis = convBasis[iBasis]["gen"];
                
                RooArgSet allFn(prodFn);
                allFn.add(*basis);
                
                eqs[iEq][iBasis]["gen"] = new RooProduct(Form("eq%d%d", iEq, iBasis), "amp0*g1", allFn);
            }
        }
    }

    map<string, RooArgList> funtot;
    RooArgList coeftot;
    map<string, RooArgList> funtotP;
    RooArgList coeftotP;
    
        
    for (int iEq = 1; iEq <= 10; iEq++) {
        auto A = A_N[iEq];
        
        for (int iBasis = 1; iBasis <= 4; iBasis++) {
            
            if (eqs.count(iEq) > 0 and eqs[iEq].count(iBasis) > 0) {
                coeftotP.add(*A);
                if (iEq <= 6) {
                    coeftot.add(*A);
                }
                for (auto& byYear: eqs[iEq][iBasis]) {
                    auto& dsLabel = byYear.first;
                    auto& eq = byYear.second;
                    
                    funtotP[dsLabel].add(*eq);
                    if (iEq <= 6) {
                        funtot[dsLabel].add(*eq);
                    }
                }
            }
        }
    }
    
    RooAbsPdf *Angular_Model;
    if (!isClosure) {
        map<string, RooAbsPdf*> angModelByYear;
        for (auto dsCatLabel: RooCatRange(*dsCat)) {
            auto dsLabel = dsCatLabel->GetName();
            
            RooRealSumPdf *PDFdefP = new RooRealSumPdf(Form("PDFdefP_%s", dsLabel), "Signal PDF Amp_i*funiT", funtotP[dsLabel], coeftotP);
            RooRealSumPdf *PDFdef = new RooRealSumPdf(Form("PDFdef_%s", dsLabel), "Signal PDF Amp_i*funiT", funtot[dsLabel], coeftot);
            
            
            RooAbsPdf* selectedPDFdef = PDFdef;
            if (isData) {
                selectedPDFdef = PDFdefP;
            }
            
            auto massSigPdf = mainWs->pdf("mass_sig");
            auto ctErrSigPdf = mainWs->pdf(Form("ctErr_pdf_sig_%s", dsLabel));
            
            auto asym = new RooRealVar(Form("asym_%s", dsLabel), "", -1, 1);
            auto eff = new RooRealVar(Form("eff_%s", dsLabel), "", 0,1);
            auto tagPdf = new RooTagPdf(Form("tagPdf_%s", dsLabel), "", *tag, *asym, *eff);
            
            auto sigModel = new RooProdPdf(Form("Angular_Model_%s", dsLabel), "PDF_Signaltag", 
                RooArgList(*massSigPdf, *ctErrSigPdf, *selectedPDFdef, *tagPdf),
                Conditional(*selectedPDFdef, RooArgList(*BsCt2DMC, *BscosthetaMC, *BscospsiMC, *BsphiMC, *mistag))
            );
            auto nSig = mainWs->var(Form("Nsig_%s", dsLabel));
            //nSig->setConstant();
            //new RooRealVar(Form("nSig_%s", dsCatLabel->GetName()), "nSig", data.sumEntries()/2, 0, data.sumEntries());

            RooArgList pdfList(*sigModel);
            RooArgList coefList(*nSig);
            
            RooAbsPdf* massBkgPdf = nullptr;
            RooAbsPdf* ctErrBkgPdf = nullptr;
            RooAbsPdf* ctBkgPdf = nullptr;
            RooAbsPdf* angBkgPdf = nullptr;
            
            RooAbsPdf* bkgModel = nullptr;
            RooRealVar* nBkg = nullptr;
            
            if (isData) {
                massBkgPdf = mainWs->pdf("mass_bkg");
                ctErrBkgPdf = mainWs->pdf(Form("ctErr_pdf_bkg_%s", dsLabel));
                ctBkgPdf = mainWs->pdf("ctBg_pdf");
                angBkgPdf = mainWs->pdf("angBg_pdf");
                
                bkgModel = new RooProdPdf(Form("bkg_model_%s", dsLabel), "bkg model",RooArgList(*massBkgPdf, *ctBkgPdf, *ctErrBkgPdf, *angBkgPdf, *tagPdf));
                nBkg = mainWs->var(Form("Nbkg_%s", dsLabel));
                //nBkg->setConstant();
                //new RooRealVar(Form("nBkg_%s", dsCatLabel->GetName()), "nBkg", data.sumEntries()/2, 0, data.sumEntries());
                
                pdfList.add(*bkgModel);
                coefList.add(*nBkg);
            }
            
            angModelByYear[dsLabel] = 
                new RooAddPdf(Form("full_Model_%s", dsCatLabel->GetName()), "Extended model", pdfList, coefList);
        }
        Angular_Model = new RooSimultaneous("Angular_Model", "PDF_Signaltag", angModelByYear, *dsCat);
    }
    else {
        RooRealSumPdf *PDFdef = new RooRealSumPdf("PDFdef", "Signal PDF Amp_i*funiT", funtot["gen"], coeftot);
        
        PDFdef->SetName("Angular_Model");
        Angular_Model = PDFdef;
    }
    
    // Too many warnings
    auto& msg = RooMsgService::instance();
    for(int iStream = 0; iStream < msg.numStreams(); iStream++) {
        auto& stream = msg.getStream(iStream);
        if (stream.minLevel <= RooFit::WARNING) {
            stream.removeTopic(RooFit::Eval);
        }
    }
    
    RooWorkspace preFitWs("preFit");
    preFitWs.import(*Angular_Model, Silence());
    
    mainWs->import(*Angular_Model, Silence());
    
    std::cerr << endl << "  ----- BEGIN FITTING  ----- " << endl << endl;
    
    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
    RooFitResult *PwaveMCResults = Angular_Model->fitTo(data, NumCPU(numCpu), Save(kTRUE), Offset(true), PrintLevel(3), Verbose(true));//, SumW2Error(true));
    //, ConditionalObservables(RooArgSet(*mistag, *tag)), ExternalConstraints(*mainWs->set("bgConstr")));
    PwaveMCResults->Print("v");
    
    RooWorkspace postFitWs("postFit");
    postFitWs.import(*Angular_Model, Silence());

    auto prevDir = gDirectory;
    if (ensureDir("outputs/")) return;
    TFile modelFile("outputs/model.root", "RECREATE");
    mainWs->Write();
    postFitWs.Write();
    preFitWs.Write();
    PwaveMCResults->SetName("fitRes");
    PwaveMCResults->Write();
    modelFile.Close();
    prevDir->cd();

    std::cerr << endl << "  ----- BEGIN PLOTTING  ----- " << endl << endl;

    // Construct 2D color plot of correlation matrix
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    TH2 *hcorr = PwaveMCResults->correlationHist("mainCorrMatrix");
    TCanvas *c = new TCanvas("Correlation Matrix", "Correlation Matrix", 1000, 600);
    hcorr->GetYaxis()->SetTitleOffset(1.4);
    hcorr->Draw("colz");
    if (ensureDir("plots/")) return;
    c->Print("plots/corrMatrix.png", "png");
    
    auto& floatPars = PwaveMCResults->floatParsFinal();
    auto phPars = floatPars.selectByName("A_0,A_pe,A_S,dGam,dm,deltaPe,deltaPa,deltaSPe,phi_s,ctau,lambda");
    TH2 *reduHcorr = reducedCorrelationHist(PwaveMCResults, *phPars, "reduCorrMatrix");
    TCanvas *cRedu = new TCanvas("Correlation Matrix", "Correlation Matrix", 1000, 600);
    reduHcorr->GetYaxis()->SetTitleOffset(1.4);
    reduHcorr->Draw("colz");
    if (ensureDir("plots/")) return;
    cRedu->Print("plots/corrMatrixReduced.png", "png");
    
    gStyle->SetPalette(kBird);
    
    Angular_Model->graphVizTree("outputs/model.dot");
    
    //auto floatPars = PwaveMCResults->floatParsFinal();
    //auto phPars = (RooArgSet*)floatPars.selectByName("A_0,A_pe,dGam,deltaPe,deltaPa,dm,mass_mean,phi_s,ctau");
    
    auto dsCatForPlotting = dsCat;
    if (not dynamic_cast<RooSimultaneous*>(Angular_Model)) {
        dsCat = nullptr;
    }
    
    Pwave_Model_plot(data, *BsCt2DMC, *Angular_Model, dsCat, 0, "finalct", RooArgSet(), true);
    //Pwave_Model_plot(data, *BsCt2DMC, *Angular_Model, dsCat, PwaveMCResults->floatParsFinal().getSize(), "ct", *phPars, true);
    //Pwave_Model_plot(data, *BsCt2DMC, *Angular_Model, dsCat, PwaveMCResults->floatParsFinal().getSize(), "ctOnlyTau", RooArgSet(*floatPars.find("ctau")), true);
    Pwave_Model_plot(data, *svmass, *Angular_Model, dsCat, 0, "finalmass", RooArgSet());
    Pwave_Model_plot(data, *BsCt2DMCErr, *Angular_Model, dsCat, 0, "finalcterr", RooArgSet());
    //RooArgSet(*(RooArgSet*)floatPars.selectByName("mass_*,nSig,nBkg")));
    //Pwave_Model_plot(data, *svmass, *Angular_Model, dsCat, PwaveMCResults->floatParsFinal().getSize(), "finalmass", RooArgSet(*(RooArgSet*)floatPars.selectByName("nSig,nBkg")));
    
//    if (isClosure) {    // not enough points to make this plot on non-GEN
        BsCt2DMC->setRange("Asymme",0.02,0.08) ;
        RooPlot* asymplot = BsCt2DMC->frame(Title("BsCt Asymetry"), Range("Asymme"));
        //data.plotOn(asymplot, LineStyle(kDashed), Binning(300), Asymmetry(*tag), Name("dataAsym"));
        //BsCt2DMCErr->setVal(0.0005);
        //data.plotOn(asymplot, LineStyle(kDashed), Binning(300), Cut("tag == tag::Bs"), CutRange("asymme"), Name("dataP"), MarkerColor(kViolet+1));
        tag->setIndex(1);
        Angular_Model->plotOn(asymplot, Range("Asymme"), ProjWData(*dsCat, data), Name("modelP"), Slice(RooArgSet(*tag)), LineColor(kViolet+1));
        
        //data.plotOn(asymplot, LineStyle(kDashed), Binning(300), Cut("tag == tag::Bsbar"), CutRange("asymme"), Name("dataM"), MarkerColor(kMagenta));
        tag->setIndex(-1);
        Angular_Model->plotOn(asymplot, Range("Asymme"), ProjWData(*dsCat, data), Name("modelM"), Slice(RooArgSet(*tag)), LineColor(kMagenta));
        
        //Angular_Model->plotOn(asymplot, Range("Asymme"), ProjWData(*dsCat, data), Name("model"));
        TCanvas cAsym;
        asymplot->Draw();
        saveCanvas(cAsym, "asym");
//    }
    
    if(!plotAngFlag) return;
    
    Pwave_Model_plot(data, *BscosthetaMC, *Angular_Model, dsCat, 0, "finalcostheta");
    Pwave_Model_plot(data, *BscospsiMC, *Angular_Model, dsCat, 0, "finalcospsi");
    Pwave_Model_plot(data, *BsphiMC, *Angular_Model, dsCat, 0, "finalphi");
}

map<int, RooAbsReal*> getFunNT(RooRealVar& costheta, RooRealVar& cospsi, RooRealVar& phi) {
    map<int, RooAbsReal*> funNT;
    
    funNT[1] = new RooFormulaVar("fun1T", "2*@1*@1*(1-(1-@0*@0)*cos(@2)*cos(@2))",
                                            RooArgSet(costheta, cospsi, phi));
    funNT[2] = new RooFormulaVar("fun2T", "(1-@1*@1)*(1-(1-@0*@0)*sin(@2)*sin(@2))",
                                            RooArgSet(costheta, cospsi, phi));
    funNT[3] = new RooFormulaVar("fun3T", "(1-@1*@1)*(1-@0*@0)", 
                                            RooArgSet(costheta, cospsi));
    funNT[4] = new RooFormulaVar("fun4T", "-(1-@1*@1)*2*@0*sqrt(1-@0*@0)*sin(@2)",
                                            RooArgSet(costheta, cospsi, phi));
    funNT[5] = new RooFormulaVar("fun5T", "2/sqrt(2.)*@1*sqrt(1-@1*@1)*(1-@0*@0)*sin(2*@2)",
                                            RooArgSet(costheta, cospsi, phi));
    funNT[6] = new RooFormulaVar("fun6T", "2/sqrt(2.)*@1*sqrt(1-@1*@1)*2*@0*sqrt(1-@0*@0)*cos(@2)",
                                            RooArgSet(costheta, cospsi, phi));
    funNT[7] = new RooFormulaVar("fun7T", "2/3*(1-(1-@0*@0)*cos(@1)*cos(@1))", 
                                            RooArgSet(costheta, phi));
    funNT[8] = new RooFormulaVar("fun8T", "sqrt(6.)/3*sqrt(1-@1*@1)*(1-@0*@0)*sin(2*@2)",
                                            RooArgSet(costheta, cospsi, phi));
    funNT[9] = new RooFormulaVar("fun9T", "sqrt(6.)/3*sqrt(1-@1*@1)*2*@0*sqrt(1-@0*@0)*cos(@2)",
                                            RooArgSet(costheta, cospsi, phi));
    funNT[10] = new RooFormulaVar("fun10T", "sqrt(3.)*4/3*@1*(1-(1-@0*@0)*cos(@2)*cos(@2))",
                                            RooArgSet(costheta, cospsi, phi));
    
    return funNT;
}

map<int, map<int, RooAbsReal*>> getCoefs(RooAbsReal& phi_s, RooAbsReal& lambda, RooCategory& tag, RooRealVar& mistag, RooRealVar& deltaPe, RooRealVar& deltaPa, RooRealVar& deltaSPe) {
    map<int, map<int, RooAbsReal*>> coefs;
    
    //RooAbsReal* S = new RooFormulaVar("S", "S", "-sin(@0)", phi_s);
    //RooAbsReal* D = new RooFormulaVar("D", "D", "-cos(@0)", phi_s);
    RooAbsReal& dilution = *new RooFormulaVar("dilution", "dilution", "1-2*@0", mistag);
    RooAbsReal& weightedTag = *new RooProduct("weightedTag", "weighted tag", RooArgList(tag, dilution));
    
    auto& C = *new RooFormulaVar("C", "C", "  (1 - @0*@0)/(1 + @0*@0)", RooArgList(lambda));
    auto& S = *new RooFormulaVar("S", "S", "-2*@0*sin(@1)/(1 + @0*@0)", RooArgList(lambda, phi_s));
    auto& D = *new RooFormulaVar("D", "D", "-2*@0*cos(@1)/(1 + @0*@0)", RooArgList(lambda, phi_s));
    
    coefs[1][1] = &thisCoefIs1;
    coefs[1][2] = &D;
    coefs[1][3] = new RooFormulaVar("coef13", "coef13", " C*@0", RooArgList(weightedTag, C));
    coefs[1][4] = new RooFormulaVar("coef14", "coef14", "-S*@0", RooArgList(weightedTag, S));
    
    coefs[2][1] = &thisCoefIs1;
    coefs[2][2] = &D;
    coefs[2][3] = new RooFormulaVar("coef23", "coef23", " C*@0", RooArgList(weightedTag, C));
    coefs[2][4] = new RooFormulaVar("coef24", "coef24", "-S*@0", RooArgList(weightedTag, S));
    
    coefs[3][1] = &thisCoefIs1;
    coefs[3][2] = new RooFormulaVar("coef32", "coef32", "-D   ", RooArgList(D));
    coefs[3][3] = new RooFormulaVar("coef33", "coef33", " C*@0", RooArgList(weightedTag, C));
    coefs[3][4] = new RooFormulaVar("coef34", "coef34", " S*@0", RooArgList(weightedTag, S));
    
    coefs[4][1] = new RooFormulaVar("coef41", "coef41", " C*sin(@0-@1)   ", RooArgList(deltaPe, deltaPa, C));
    coefs[4][2] = new RooFormulaVar("coef42", "coef42", " S*cos(@0-@1)   ", RooArgList(deltaPe, deltaPa, S));
    coefs[4][3] = new RooFormulaVar("coef43", "coef43", "   sin(@0-@1)*@2", RooArgList(deltaPe, deltaPa, weightedTag));
    coefs[4][4] = new RooFormulaVar("coef44", "coef44", " D*cos(@0-@1)*@2", RooArgList(deltaPe, deltaPa, weightedTag, D));
    
    coefs[5][1] = new RooFormulaVar("coef51", "coef51", "   cos(@0)   ", RooArgList(deltaPa));
    coefs[5][2] = new RooFormulaVar("coef52", "coef52", " D*cos(@0)   ", RooArgList(deltaPa, D));
    coefs[5][3] = new RooFormulaVar("coef53", "coef53", " C*cos(@0)*@1", RooArgList(deltaPa, weightedTag, C));
    coefs[5][4] = new RooFormulaVar("coef54", "coef54", "-S*cos(@0)*@1", RooArgList(deltaPa, weightedTag, S));
    
    coefs[6][1] = new RooFormulaVar("coef61", "coef61", " C*sin(@0)   ", RooArgList(deltaPe, C));
    coefs[6][2] = new RooFormulaVar("coef62", "coef62", " S*cos(@0)   ", RooArgList(deltaPe, S));
    coefs[6][3] = new RooFormulaVar("coef63", "coef63", "   sin(@0)*@1", RooArgList(deltaPe, weightedTag));
    coefs[6][4] = new RooFormulaVar("coef64", "coef64", " D*cos(@0)*@1", RooArgList(deltaPe, weightedTag, D));
    
    coefs[7][1] = &thisCoefIs1;
    coefs[7][2] = new RooFormulaVar("coef72", "coef72", "-D   ", RooArgList(D));
    coefs[7][3] = new RooFormulaVar("coef73", "coef73", " C*@0", RooArgList(weightedTag, C));
    coefs[7][4] = new RooFormulaVar("coef74", "coef74", " S*@0", RooArgList(weightedTag, S));
    
    //deltaSPe = deltaS - deltaPe
    coefs[8][1] = new RooFormulaVar("coef81", "coef81", " C*sin(@0-@1-@2)   ", RooArgList(deltaPa, deltaSPe, deltaPe, C));
    coefs[8][2] = new RooFormulaVar("coef82", "coef82", " S*sin(@0-@1-@2)   ", RooArgList(deltaPa, deltaSPe, deltaPe, S));
    coefs[8][3] = new RooFormulaVar("coef83", "coef83", "   cos(@0-@1-@2)*@3", RooArgList(deltaPa, deltaSPe, deltaPe, weightedTag));
    coefs[8][4] = new RooFormulaVar("coef84", "coef84", " D*sin(@0-@1-@2)*@3", RooArgList(deltaPa, deltaSPe, deltaPe, weightedTag, D));
    
    coefs[9][1] = new RooFormulaVar("coef91", "coef91", "   sin(-@0)   ", RooArgList(deltaSPe));
    coefs[9][2] = new RooFormulaVar("coef92", "coef92", "-D*sin(-@0)   ", RooArgList(deltaSPe, D));
    coefs[9][3] = new RooFormulaVar("coef93", "coef93", " C*sin(-@0)*@1", RooArgList(deltaSPe, weightedTag, C));
    coefs[9][4] = new RooFormulaVar("coef94", "coef94", " S*sin(-@0)*@1", RooArgList(deltaSPe, weightedTag, S));
    
    coefs[10][1] = new RooFormulaVar("coef101", "coef101", " C*sin(-@0-@1)   ", RooArgList(deltaSPe, deltaPe, C));
    coefs[10][2] = new RooFormulaVar("coef102", "coef102", " S*sin(-@0-@1)   ", RooArgList(deltaSPe, deltaPe, S));
    coefs[10][3] = new RooFormulaVar("coef103", "coef103", "   cos(-@0-@1)*@2", RooArgList(deltaSPe, deltaPe, weightedTag));
    coefs[10][4] = new RooFormulaVar("coef104", "coef104", " D*sin(-@0-@1)*@2", RooArgList(deltaSPe, deltaPe, weightedTag, D));
    
    return coefs;
}
