#include <iostream>
#include <time.h>

#include <TROOT.h>
#include <TMath.h>
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF2.h>
#include <TAxis.h>

#define NBINS               30

TH2F*          template_;

// Converts a TF2 to TH2F, matching the dimentions of hRef
TH2F* funcToHistogram2D(TF2* f, TH2F* hRef) {
  int         nBinsX = hRef->GetNbinsX();
  int         nBinsY = hRef->GetNbinsY();
  double      xMin   = hRef->GetXaxis()->GetXmin();
  double      xMax   = hRef->GetXaxis()->GetXmax();
  double      yMin   = hRef->GetYaxis()->GetXmin();
  double      yMax   = hRef->GetYaxis()->GetXmax();
  const char* name   = Form("%s_fit", hRef->GetName());
  const char* title  = Form("Fit %s", hRef->GetTitle());

  TH2F* h = new TH2F(name, title, nBinsX, xMin, xMax, nBinsY, yMin, yMax);

  for (int i=1; i<=nBinsX; ++i) {
    for (int j=1; j<=nBinsY; ++j) {
      double x = h->GetXaxis()->GetBinCenter(i);
      double y = h->GetYaxis()->GetBinCenter(j);
      double v = f->Eval(x, y);
      double e = 0.0; // TODO: Calculate bin error
      h->SetBinContent(i, j, v);
      h->SetBinError(i, j, e);
    }
  }
  return h;
}

// Projects a 2D hist and function to one dimension
void plotProjection(TH2F* h, TF2* f, const char* axis="X", TVirtualPad* pad = NULL) {
  TVirtualPad* pad1 = NULL;
  TVirtualPad* pad2 = NULL;
  if (pad != NULL) {
    pad->Divide(1,2, 0.01, 0.01);
    pad1 = pad->cd(1);
    pad2 = pad->cd(2);
    pad->cd(1);
  }

  /** Project 2D hist and function ***********/
  TH2D* h_func  = (TH2D*) funcToHistogram2D(f, h);
  TH1D* h_proj  = NULL;
  TH1D* h_fproj = NULL;
  if        (strcmp(axis, "X") == 0 || strcmp(axis, "x") == 0) {
    h_proj  = h->ProjectionX();
    h_fproj = h_func->ProjectionX();
  } else if (strcmp(axis, "Y") == 0 || strcmp(axis, "y") == 0) {
    h_proj = h->ProjectionY();
    h_fproj = h_func->ProjectionY();
  }
  if (h_proj == NULL || h_proj == NULL) return;
  h_proj->SetMinimum(0.0);
  h_fproj->SetMinimum(0.0);

  h_proj->SetTitle(Form("%s %s projection", h->GetTitle(), axis));
  h_proj->Draw();
  h_fproj->SetLineWidth(2);
  h_fproj->SetLineColor(kRed);
  h_fproj->SetFillStyle(4000);
  h_fproj->Draw("SAMEHIST");

  /** Calculate fit rediduals ***********/
  if (pad != NULL) {
    pad->cd(2);

    TH1F* hdiff = (TH1F*) h_proj->Clone(Form("%s_dif", h_proj->GetName()));
    hdiff->SetTitle("Residuals;x;y");
    hdiff->Add(h_fproj, -1.0);

    double yRange = TMath::Max(fabs(hdiff->GetMinimum()), fabs(hdiff->GetMaximum()));
    hdiff->SetMinimum(-yRange);
    hdiff->SetMaximum(yRange);
    hdiff->SetFillColor(kBlue+2);
    hdiff->SetStats(0);
    hdiff->Draw("HIST");
  }
  pad->cd(1);
}

// Plots the fit results
void plotTempFitResults(const char* filename,
			TH2F* data,
			TF2*  fit,
			TH2F* temp) {
  TCanvas* c = new TCanvas("fit", "fit");
  c->Divide(2,2);
  TVirtualPad* pad = NULL;

  pad = c->cd(1);
  data->SetMinimum(0.0);
  data->Draw("COLZ");
  pad = c->cd(2);
  temp->SetMinimum(0.0);
  temp->Draw("COLZ");
  pad = c->cd(3);
  plotProjection(data, fit, "X", pad);
  pad = c->cd(4);
  plotProjection(data, fit, "Y", pad);

  c->Update();
  c->Print(filename);
}

// Fitting function given to TF2
double templateHistFunc(Double_t *x, Double_t* par) {
  Double_t xx   = x[0];
  Double_t yy   = x[1];
  Int_t    binx = template_->GetXaxis()->FindBin(xx);
  Int_t    biny = template_->GetYaxis()->FindBin(yy);
  Double_t t1   = par[0]*template_->GetBinContent(binx, biny);
  return (t1);
}

// Fits a 2D function to given data in a 2D hist
TF2* templateFit2D(const char* name, TH2F* data) {
  TF2* f = new TF2(name, templateHistFunc, -1.0, 1.0, -TMath::Pi(), TMath::Pi(), 1);
  f->SetParameter(0, 2.0);
  f->SetParLimits(0, -2.0, 2.0);
  data->Fit(f, "EM0");
  return f;
}

/*****************************************
 * TempFitStandalone
 *****************************************/
void TempFitStandalone(int nEvents=10000) {
  const char* plotFilename = "fitplot.pdf";

  /**** INITIALIZING TEMPLATE ****/
  TFile* templateFile   = new TFile("template.root", "READ");
  TH2F* tmp = (TH2F*) templateFile->Get("phi_cosTh_detEvts_pt0_rap0_f0");
  if (tmp != NULL) {
    template_ = (TH2F*) tmp->Clone("template");
    template_->SetTitle("Template;x;y");
  }

  /**** GENERATING EVENTS ****/
  // Dummy function used to fill data
  TF2* fTemplate = new TF2("fTemplate", templateHistFunc, -1.0, 1.0, -TMath::Pi(), TMath::Pi(),1);
  fTemplate->SetParameter(0,1.0);
  TH2F* dataHist = new TH2F("dataHist", "Data", NBINS, -1.0, 1.0, NBINS, -TMath::Pi(), TMath::Pi());
  dataHist->FillRandom("fTemplate", nEvents);
  dataHist->SetTitle("Data;x;y");

  /**** FITTING ****/
  TF2* tempFit = templateFit2D("tempFit2D", dataHist);
  tempFit->SetTitle("Fit results;x;y");

  /**** PLOTTING ****/
  plotTempFitResults(plotFilename, dataHist, tempFit, template_);
}
