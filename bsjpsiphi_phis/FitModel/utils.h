//Authors: Enrico Lusiani
#ifndef RECOFIT_UTILS
#define RECOFIT_UTILS

#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TChain.h"

#include "RooPlot.h"
#include "RooConstVar.h"
#include "RooHist.h"
#include "RooAbsRealLValue.h"
#include "RooAbsCategory.h"

#include <iostream>

TTree* getTree(const char* filename, const char* treename) {
    auto prevDir = gDirectory;
    auto file = TFile::Open(filename);
    if (!file) {
        std::clog << "Error: Unable to open file, aborting" << std::endl;
        std::clog << "File path: " << filename << std::endl;
        return nullptr;
    }
    
    auto tree = (TTree *)file->Get(treename);
    if (!tree) {
        std::clog << "Error: Unable to find the tree, aborting" << std::endl;
        std::clog << "No tree named " << treename << " in directory " << file->GetPath() << std::endl;
        return nullptr;
    }
    
    prevDir->cd();
    return tree;
}

TTree* getTree(std::vector<const char*> filenames, const char* treename) {
    auto prevDir = gDirectory;
    
    auto tree = new TChain(treename);
    for(auto filename: filenames) {
        tree->Add(filename);
    }
    
    prevDir->cd();
    return tree;
}

int ensureDir(const char* dirPath) {
    if (gSystem->mkdir(dirPath, true))
    {
        if (gSystem->AccessPathName(dirPath))
        {
            std::clog << "Unable to create the directory " << dirPath << std::endl;
            return 1;
        }
    }
    return 0;
}

void saveCanvas(TCanvas& c, TString name) {
    auto dirPath = "plots/";
    if (ensureDir(dirPath)) return;
    for (auto ext: {"png", "pdf", "C", "root"}) {
        auto filename = dirPath + name + "." + ext;
        c.SaveAs(filename);
    }
}

void roo_pulls(TCanvas& c, RooPlot* frame, const char* data, const char* fit, int chi2ndf = -1) {
    c.Divide(1, 2);
    c.GetPad(1)->SetLogy(c.GetLogy());
    c.GetPad(1)->SetPad(0.0, 0.25, 1, 1);
    c.GetPad(1)->SetBottomMargin(0.015);
    c.GetPad(1)->SetRightMargin(0.05);
//    c.GetPad(1)->SetTicks(1, 1);
    c.GetPad(2)->SetPad(0.0, 0.0, 1, 0.25);
    c.GetPad(2)->SetBottomMargin(0.32);
    c.GetPad(2)->SetTopMargin(0.0);
    c.GetPad(2)->SetRightMargin(0.05);
//    c.GetPad(2)->SetTicks(1, 1);
    
    c.GetPad(1)->SetBorderMode(0);
    c.GetPad(2)->SetBorderMode(0);
    
    auto pulls = frame->pullHist(data, fit, true); // true for useAverage
    pulls->SetMarkerStyle(frame->getAttMarker(data)->GetMarkerStyle());
    
    auto pullFrame = frame->emptyClone("pulls");
    auto line0 = RooFit::RooConst(0);
    line0.plotOn(pullFrame, RooFit::Name("line0"), RooFit::LineWidth(1), RooFit::LineColor(kBlack));
    pullFrame->SetTitle("");
    pullFrame->addPlotable(pulls, "PE1");
    
    auto label_scale = 1.;
    
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->CenterTitle();
    pullFrame->GetXaxis()->SetTitleSize(0.18);
    pullFrame->GetYaxis()->SetTitleSize(0.18);
    pullFrame->GetYaxis()->SetTitleOffset(0.28);
    pullFrame->GetXaxis()->SetTitleOffset(.82);
    pullFrame->GetXaxis()->SetLabelSize(0.12 * label_scale);
    pullFrame->GetYaxis()->SetLabelSize(0.12 * label_scale);
    pullFrame->GetYaxis()->SetLabelOffset(0.006);
    
    frame->GetXaxis()->SetLabelOffset(0.006);
    frame->GetXaxis()->SetLabelSize(0.1);
    frame->GetYaxis()->SetTitleOffset(0.8);
    frame->GetYaxis()->SetTitleSize(0.06);
    
    if (pulls->GetMaximum() > 3.5 or pulls->GetMinimum() < -3.5)
    {
        pullFrame->SetMinimum(-5.5);
        pullFrame->SetMaximum(5.5);
    }
    else
    {
        pullFrame->SetMinimum(-3.5);
        pullFrame->SetMaximum(3.5);
    }
    
    c.GetPad(1)->cd();
    frame->Draw();
    
    if (chi2ndf >= 0) {
        Double_t chisquare = frame->chiSquare(fit, data, chi2ndf);
          
        auto& var = *frame->getPlotVar();

        auto range = var.getMax() - var.getMin();
        auto chi2 = new TLatex(var.getMin() + range/1000, frame->GetMaximum()*1.04, Form("#chi^{2}/NDF = %f", chisquare));
        chi2->Draw();
    }
    
    
    c.GetPad(2)->cd();
    pullFrame->Draw();
}

TFile* openFileForReading(const char* name) {
    auto prevDir = gDirectory;
    auto file = TFile::Open(name);
    prevDir->cd();
    return file;
}

// WARNING this is JUST for presentation, it's likely statustically WRONG
// copied and modified from RooFitResult::correlationHist
TH2* reducedCorrelationHist(RooFitResult* res, const RooAbsCollection& selectedPars, const char* name = "correlation_matrix") {
    auto& corr = res->correlationMatrix();
    
    Int_t n = selectedPars.getSize();
    Int_t N = corr.GetNcols();
    
    auto& allPars = res->floatParsFinal();

    TH2D* hh = new TH2D(name,name,n,0,n,n,0,n) ;
    int realI = 0, realJ = 0;
    for (Int_t i = 0; i < N; i++) {
        if (not selectedPars.find(allPars.at(i)->GetName())) {
            continue;
        }
        for (Int_t j = 0; j < N; j++) {
            if (not selectedPars.find(allPars.at(j)->GetName())) {
                continue;
            }
            hh->Fill(realI + 0.5, n - realJ - 0.5, corr(i,j));
            realJ++;
        }
        hh->GetXaxis()->SetBinLabel(realI + 1,allPars.at(i)->GetName());
        hh->GetYaxis()->SetBinLabel(n - realI,allPars.at(i)->GetName());
        realI++;
    }
    hh->SetMinimum(-1) ;
    hh->SetMaximum(+1) ;


    return hh ;
}


// adapter for roocollections to emulate features of newer ROOT versions
#if ROOT_VERSION_CODE < ROOT_VERSION(6, 18, 0)
TIter begin(const RooAbsCollection& coll)
{
	return TIter(new RooLinkedListIter(std::move(coll.iterator())));
}

TIter end(const RooAbsCollection& coll)
{
	return TIter::End();
}
#endif

class RooCatRange
{
    public: 
        RooCatRange(const RooAbsCategory& sim)
            :sim(sim)
        {
        };
        
#if ROOT_VERSION_CODE > ROOT_VERSION(6, 20, 0)
        std::vector< RooCatType * >::const_iterator begin() const
        {
            return sim.begin();
#else
        ROOT::Internal::TRangeDynCastIterator<RooCatType> begin() const
        {
            auto iter = sim.typeIterator();
            iter->Next();
            return ROOT::Internal::TRangeDynCastIterator<RooCatType>(iter);
#endif
        };
        
#if ROOT_VERSION_CODE > ROOT_VERSION(6, 20, 0)
        std::vector< RooCatType * >::const_iterator end() const
        {
            return sim.end();
#else
        ROOT::Internal::TRangeDynCastIterator<RooCatType> end() const
        {
            return ROOT::Internal::TRangeDynCastIterator<RooCatType>(TIter::End());
#endif
        };
    private:
        const RooAbsCategory& sim;
};

#endif
