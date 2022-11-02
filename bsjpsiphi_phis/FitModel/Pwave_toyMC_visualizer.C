// Author Enrico Lusiani

#include "TH1.h"
#include "TStyle.h"
#include "TTree.h"

#include "RooFitResult.h"

#include "utils.h"

map<string, const char*> titleMap = []() {
    map<string, const char*> m;
    m["A_0"] = "A_{0} pull";
    m["A_pe"] = "A_{#perp}  pull";
    m["A_S"] = "A_{S}  pull";
    m["dGam"] = "#Delta#Gamma_{s} pull";
    m["deltaPa"] = "#delta_{#parallel} pull";
    m["deltaPe"] = "#delta_{#perp}  pull";
    m["deltaSPe"] = "#delta_{S#perp}  pull";
    m["dm"] = "#Deltam_{s} pull";
    m["phi_s"] = "#phi_{s} pull";
    m["ctau"] = "c#tau pull";
    m["lambda"] = "|#lambda| pull";
    m["Nsig_dt17"] = "N_{sig}^{17}";
    m["Nsig_dt18"] = "N_{sig}^{18}";
    m["Nbkg_dt17"] = "N_{bkg}^{17}";
    m["Nbkg_dt18"] = "N_{bkg}^{18}";
    
    return m;
}();


void Pwave_toyMC_visualizer(TTree* parTree, TTree* parTreeRef, bool drawCov) {
    gStyle->SetOptFit();
    gErrorIgnoreLevel = kWarning;
    
    stringstream output;
    
    // we want to reuse RooFit formatting;
    RooRealVar dummyVar("dummyVar", "", -5, 5);
    
    vector<TString> vars = {"A_0", "A_pe", "A_S", "deltaPa", "deltaPe", "deltaSPe", "ctau", "dGam", "dm", "lambda", "phi_s"};//, "Nsig_dt17", "Nsig_dt18", "Nbkg_dt17", "Nbkg_dt18"};
    vector<double> sigmas(vars.size());
    vector<double> means(vars.size());
    
    TMatrixDSym cov(vars.size());
    
    for (int iVar = 0; iVar < vars.size(); iVar++) {
        auto var = vars[iVar];
        
        TCanvas c;
        c.Divide(2,2);
        c.GetPad(1)->cd();
        parTree->Draw(var + "Mean >> mean(50)", "covQual == 3 && migrStat == 0 && hesseStat == 0");
        
        c.GetPad(2)->cd();
        parTree->Draw(var + "Error >> error(50)", "covQual == 3 && migrStat == 0 && hesseStat == 0");
        
        c.GetPad(3)->cd();
        parTree->Draw(var + "Pull >> pull(50)", "covQual == 3 && migrStat == 0 && hesseStat == 0");
        
        saveCanvas(c, "toyMC_" + var);
        
        TCanvas cPull;
        parTree->Draw(Form("%1$sPull >> %1$s_pull(40, -5, 5)", var.Data()), "covQual == 3 && migrStat == 0 && hesseStat == 0", "goff");
        auto h = (TH1*)gDirectory->Get(Form("%s_pull", var.Data()));
        gStyle->SetOptStat(false);
        h->SetMarkerStyle(20);
        h->SetLineColor(kBlack);
        h->Sumw2();
        h->GetXaxis()->SetRangeUser(-5, 5);
        h->SetXTitle(titleMap[var.Data()]);
        h->GetYaxis()->SetRangeUser(0, TMath::Max(h->GetMaximum()*1.1, 120.));
        h->SetTitle("");
        auto f = new TF1("f", "gaus");
        h->Fit(f, "Q", "P E1 X0");
        dummyVar.setVal(f->GetParameter(1));
        dummyVar.setError(f->GetParError(1));
        output << var << "\tsist: " << f->GetParameter(1) << " +/- " << f->GetParError(1);//*dummyVar.format(0, "EUP");
        sigmas[iVar] = f->GetParameter(2);
        means[iVar] = f->GetParameter(1);
        //h->Draw("P E1 X0");
        if (parTreeRef) {
            parTreeRef->Draw(Form("%1$sPull >> %1$s_pull2(40, -5, 5)", var.Data()), "covQual == 3 && migrStat == 0 && hesseStat == 0", "goff");
            auto hRef = (TH1*)gDirectory->Get(Form("%s_pull2", var.Data()));
            auto fRef = new TF1("fRef", "gaus");
            hRef->Scale(h->Integral()/hRef->Integral());
            fRef->SetLineColor(kBlue);
            fRef->SetLineStyle(kDashed);
            hRef->SetMarkerColor(kBlue);
            hRef->SetLineColor(kBlue);
            hRef->SetXTitle(titleMap[var.Data()]);
            hRef->GetYaxis()->SetRangeUser(0, TMath::Max(hRef->GetMaximum()*1.1, 120.));
            hRef->SetTitle("");
            hRef->Fit(fRef, "Q", "FUNC SAME");
            h->Draw("P E1 X0 SAME");
            dummyVar.setVal(fRef->GetParameter(1));
            dummyVar.setError(fRef->GetParError(1));
            output << "\tref: " << fRef->GetParameter(1) << " +/- " << fRef->GetParError(1);//*dummyVar.format(0, "EUP");
            output << "\tsd: " << (f->GetParameter(1) - fRef->GetParameter(1))/hypot(fRef->GetParError(1), sqrt(hRef->GetEntries()*1./h->GetEntries())*f->GetParError(1));
            output << "\tdiff: " << f->GetParameter(1) - fRef->GetParameter(1);
            clog<< hypot(fRef->GetParError(1), sqrt(hRef->GetEntries()*1./h->GetEntries())*f->GetParError(1)) << " " << hypot(fRef->GetParError(1), f->GetParError(1)) << endl;
            clog << sqrt(hRef->GetEntries()*1./h->GetEntries()) <<  " " << fRef->GetParError(1) << " " << f->GetParError(1) << endl;
        }
        output << endl;
        saveCanvas(cPull, "toyMCPresentation_" + var);
        gStyle->SetOptStat();
        
        
    }
    
    if (drawCov) {
        gStyle->SetOptStat(0);
        for (int iVar = 0; iVar < vars.size(); iVar++) {
            auto var = vars[iVar];
            
            for (int jVar = iVar + 1; jVar < vars.size(); jVar++) {
                auto var2 = vars[jVar];
                
                clog << var << " " << var2 << endl;
                
                parTree->Draw(Form("%2$sPull:%1$sPull >> %1$s_%2$s_pull(40, -5, 5, 40, -5, 5)", var.Data(), var2.Data()), "covQual == 3 && migrStat == 0 && hesseStat == 0", "goff");
                auto hXY = (TH2*)gDirectory->Get(Form("%1$s_%2$s_pull", var.Data(), var2.Data()));
                hXY->SetTitle("");
                hXY->SetXTitle(titleMap[var.Data()]);
                hXY->SetYTitle(titleMap[var2.Data()]);
                auto fXY = new TF2("fXY", "[Constant]*exp((-0.5*pow(((x-[MeanX])/[SigmaX]),2 ) - 0.5*pow(((y-[MeanY])/[SigmaY]),2) + (x-[MeanX])*(y-[MeanY])*[CorrXY]/([SigmaX]*[SigmaY]))/((1 + [CorrXY])*(1 - [CorrXY])))");
                
                fXY->SetParameters(100,0, 0, 0, 1, 1);
                
                fXY->SetParLimits(1, -0.9, 0.9);
                fXY->SetParLimits(2, -1, 1);
                fXY->SetParLimits(3, -1, 1);
                fXY->SetParLimits(4, 0, 3);
                fXY->SetParLimits(5, 0, 3);
                
                fXY->FixParameter(2, means[iVar]);
                fXY->FixParameter(3, means[jVar]);
                fXY->FixParameter(4, sigmas[iVar]);
                fXY->FixParameter(5, sigmas[jVar]);
                
                fXY->SetLineWidth(1);
                fXY->SetNpx(1000);
                fXY->SetNpy(1000);
                
                hXY->GetZaxis()->SetRangeUser(0, 20);
                
                TCanvas cXY;
                hXY->Fit(fXY, "L", "COL");
                if (ensureDir("plots/corrPlots")) return;
                saveCanvas(cXY, Form("corrPlots/%1$s_%2$s", var.Data(), var2.Data()));
                cov(iVar, jVar) = fXY->GetParameter(1);
                
                clog << "\n\n=========================================================\n\n" << endl;
            }
            cov(iVar, iVar) = 1;
        }
        
        cov.Print("V");
    }
    
    clog << "Number of valid events: " << parTree->Draw("1", "covQual == 3 && migrStat == 0 && hesseStat == 0", "goff") << endl;
    if (parTreeRef) {
        clog << "(were " << parTreeRef->Draw("1", "covQual == 3 && migrStat == 0 && hesseStat == 0", "goff") << " in reference)" << endl;
    }
    clog << output.str();
}

void Pwave_toyMC_visualizer(TString filename, TString filename2 = nullptr, bool drawCov = true) {
    auto parTree = getTree(filename, "parTree");
    TTree* parTreeRef = nullptr;
    if (not filename2.IsNull()) {
        parTreeRef = getTree(filename2, "parTree");
    }
    
    Pwave_toyMC_visualizer(parTree, parTreeRef, drawCov);
}
