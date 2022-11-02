//Author: Md. Alibordi
void plotSlices(TH3* numHist, TH3* denHist, RooAbsReal& projectedFunc, const RooArgSet& vars, int maxOrder, int xbins, int ybins, int zbins, TString suffix) {
  using namespace RooFit;
  
  auto& costhetaMC = *(RooRealVar*)vars.find("BscosthetaMC");
  auto& cospsiMC = *(RooRealVar*)vars.find("BscospsiMC");
  auto& phiMC = *(RooRealVar*)vars.find("BsphiMC");

  vector<TEfficiency *> effHistsX;
  vector<TEfficiency *> effHistsY;
  vector<TEfficiency *> effHistsZ;
  vector<RooPlot *> xframes;
  vector<RooPlot *> yframes;
  vector<RooPlot *> zframes;
  auto cx1new = new TCanvas(
      "cx1new", "Projected efficiency Cos(#theta_{T}) slices", 1400, 1400);
  auto cy1new = new TCanvas(
      "cy1new", "Projected efficiency Cos(#psi_{T}) slices", 1400, 1400);
  auto cz1new =
      new TCanvas("cz1new", "Projected efficiency #phi slices", 1400, 1400);
  cx1new->Divide(5, 5);
  cy1new->Divide(5, 5);
  cz1new->Divide(5, 5);

  TLegend *leg = new TLegend(0.35, 0.8, 0.9, 0.9);

  // width of the slices in the hidden variables ("border" is half of it)
  double border = 0.035;
  // variables to be filled with global efficiency maximum
  double maxEffX = 0;
  double maxEffY = 0;
  double maxEffZ = 0;

  // loop over slice grid
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {

      // central values and borders of the slices in the hidden variables
      double centA = -0.8 + 1.6 * i / 4;
      double centB = -0.8 + 1.6 * j / 4;
      double lowA = TMath::Max(centA - border, 1e-4 - 1);
      double lowB = TMath::Max(centB - border, 1e-4 - 1);
      double highA = TMath::Min(centA + border, -1e-4 + 1);
      double highB = TMath::Min(centB + border, -1e-4 + 1);

      // slicing num and den distributions
      auto numProjX = numHist->ProjectionX(
          "numProjX", numHist->GetYaxis()->FindBin(lowA),
          numHist->GetYaxis()->FindBin(highA),
          numHist->GetZaxis()->FindBin(lowB * TMath::Pi()),
          numHist->GetZaxis()->FindBin(highB * TMath::Pi()), "e");
      auto numProjY = numHist->ProjectionY(
          "numProjY", numHist->GetXaxis()->FindBin(lowA),
          numHist->GetXaxis()->FindBin(highA),
          numHist->GetZaxis()->FindBin(lowB * TMath::Pi()),
          numHist->GetZaxis()->FindBin(highB * TMath::Pi()), "e");
      auto numProjZ =
          numHist->ProjectionZ("numProjZ", numHist->GetXaxis()->FindBin(lowA),
                               numHist->GetXaxis()->FindBin(highA),
                               numHist->GetYaxis()->FindBin(lowB),
                               numHist->GetYaxis()->FindBin(highB), "e");
      auto denProjX = denHist->ProjectionX(
          "denProjX", denHist->GetYaxis()->FindBin(lowA),
          denHist->GetYaxis()->FindBin(highA),
          denHist->GetZaxis()->FindBin(lowB * TMath::Pi()),
          denHist->GetZaxis()->FindBin(highB * TMath::Pi()), "e");
      auto denProjY = denHist->ProjectionY(
          "denProjY", denHist->GetXaxis()->FindBin(lowA),
          denHist->GetXaxis()->FindBin(highA),
          denHist->GetZaxis()->FindBin(lowB * TMath::Pi()),
          denHist->GetZaxis()->FindBin(highB * TMath::Pi()), "e");
      auto denProjZ =
          denHist->ProjectionZ("denProjZ", denHist->GetXaxis()->FindBin(lowA),
                               denHist->GetXaxis()->FindBin(highA),
                               denHist->GetYaxis()->FindBin(lowB),
                               denHist->GetYaxis()->FindBin(highB), "e");

      // producing 1D efficiencies from the slices
      effHistsX.push_back(new TEfficiency(*numProjX, *denProjX));
      effHistsX.back()->SetName(Form("effHistX_%i_%i", i, j));
      effHistsX.back()->SetTitle(
          Form("Efficiency - slice cospsiMC=%1.2f "
               "phiMC=%1.2f;cos(#theta_{T});Efficiency",
               centA, centB * TMath::Pi()));

      effHistsY.push_back(new TEfficiency(*numProjY, *denProjY));
      effHistsY.back()->SetName(Form("effHistY_%i_%i", i, j));
      effHistsY.back()->SetTitle(Form("Efficiency - slice costhetaMC=%1.2f "
                                      "phiMC=%1.2f;cos(#psi_{T});Efficiency",
                                      centA, centB * TMath::Pi()));

      effHistsZ.push_back(new TEfficiency(*numProjZ, *denProjZ));
      effHistsZ.back()->SetName(Form("effHistZ_%i_%i", i, j));
      effHistsZ.back()->SetTitle(Form("Efficiency - slice costhetaMC=%1.2f "
                                      "cospsiMC=%1.2f;#phi;Efficiency",
                                      centA, centB));

      // producing 1D slices of efficiency description
      cospsiMC.setVal(centA);
      phiMC.setVal(centB * TMath::Pi());
      xframes.push_back(
          costhetaMC.frame(Name(Form("frameslice_costheta_%i_%i_", i, j))));
      projectedFunc.plotOn(xframes.back(), LineColor(kRed),
                           Name(Form("projectedFuncx_%i_%i", i, j)));

      costhetaMC.setVal(centA);
      phiMC.setVal(centB * TMath::Pi());
      yframes.push_back(
          cospsiMC.frame(Name(Form("frameslice_cospsi_%i_%i_", i, j))));
      projectedFunc.plotOn(yframes.back(), LineColor(kRed),
                           Name(Form("projectedFuncy_%i_%i", i, j)));

      costhetaMC.setVal(centA);
      cospsiMC.setVal(centB);
      zframes.push_back(
          phiMC.frame(Name(Form("frameslice_phi_%i_%i_", i, j))));
      projectedFunc.plotOn(zframes.back(), LineColor(kRed),
                           Name(Form("projectedFuncz_%i_%i", i, j)));

      // plot in canvas
      cx1new->cd(5 * j + i + 1);
      gPad->SetLeftMargin(0.18);
      effHistsX.back()->Draw();
      cx1new->cd(5 * j + i + 1)->Update();
      auto graphx = effHistsX.back()->GetPaintedGraph();
      graphx->SetMinimum(0);
      auto effValsX = graphx->GetY();
      for (int iBin = 0; iBin < graphx->GetN(); ++iBin)
        if (maxEffX < effValsX[iBin])
          maxEffX = effValsX[iBin];
      // if (maxEffX<graphx->GetYaxis()->GetXmax()) maxEffX =
      // graphx->GetYaxis()->GetXmax();
      graphx->GetYaxis()->SetTitleOffset(1.7);
      cx1new->cd(5 * j + i + 1)->Update();
      xframes.back()->Draw("same");

      if (i + j == 0) {
        leg->AddEntry(effHistsX.back(), "Binned efficiency", "lep");
        leg->AddEntry(
            xframes.back()->findObject(Form("projectedFuncx_%i_%i", i, j)),
            "Projected spherical harmonics", "l");
      }

      leg->Draw("same");

      cy1new->cd(5 * j + i + 1);
      gPad->SetLeftMargin(0.18);
      effHistsY.back()->Draw();
      cy1new->cd(5 * j + i + 1)->Update();
      auto graphy = effHistsY.back()->GetPaintedGraph();
      graphy->SetMinimum(0);
      auto effValsY = graphy->GetY();
      for (int iBin = 0; iBin < graphy->GetN(); ++iBin)
        if (maxEffY < effValsY[iBin])
          maxEffY = effValsY[iBin];
      // if (maxEffY<graphy->GetYaxis()->GetXmax()) maxEffY =
      // graphy->GetYaxis()->GetXmax();
      graphy->GetYaxis()->SetTitleOffset(1.7);
      cy1new->cd(5 * j + i + 1)->Update();
      yframes.back()->Draw("same");
      leg->Draw("same");

      cz1new->cd(5 * j + i + 1);
      gPad->SetLeftMargin(0.18);
      effHistsZ.back()->Draw();
      cz1new->cd(5 * j + i + 1)->Update();
      auto graphz = effHistsZ.back()->GetPaintedGraph();
      graphz->SetMinimum(0);
      auto effValsZ = graphz->GetY();
      for (int iBin = 0; iBin < graphz->GetN(); ++iBin)
        if (maxEffZ < effValsZ[iBin])
          maxEffZ = effValsZ[iBin];
      // if (maxEffZ<graphz->GetYaxis()->GetXmax()) maxEffZ =
      // graphz->GetYaxis()->GetXmax();
      graphz->GetYaxis()->SetTitleOffset(1.7);
      cz1new->cd(5 * j + i + 1)->Update();
      zframes.back()->Draw("same");
      leg->Draw("same");
    }
  }

  // set uniform y-axis ranges
  for (int i = 0; i < effHistsX.size(); ++i) {
    (effHistsX[i]->GetPaintedGraph())->GetXaxis()->SetRangeUser(-1, 1);
    (effHistsX[i]->GetPaintedGraph())->SetMaximum(maxEffX * 1.25);
  }
  for (int i = 0; i < effHistsY.size(); ++i) {
    (effHistsY[i]->GetPaintedGraph())->GetXaxis()->SetRangeUser(-1, 1);
    (effHistsY[i]->GetPaintedGraph())->SetMaximum(maxEffY * 1.25);
  }
  for (int i = 0; i < effHistsZ.size(); ++i) {
    (effHistsZ[i]->GetPaintedGraph())->GetXaxis()->SetRangeUser(-TMath::Pi(), TMath::Pi());
    (effHistsZ[i]->GetPaintedGraph())->SetMaximum(maxEffZ * 1.25);
  }
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j) {
      cx1new->cd(5 * j + i + 1)->Update();
      cy1new->cd(5 * j + i + 1)->Update();
      cz1new->cd(5 * j + i + 1)->Update();
    }

  saveCanvas(*cx1new, Form("EffProj_%i_%i_%i_CosThetaslices_SpH%io_dp%i_%s",
                      xbins, ybins, zbins, maxOrder, (int)(border * 200), suffix.Data()));
  saveCanvas(*cy1new, Form("EffProj_%i_%i_%i_CosPsislices_SpH%io_dp%i_%s",
                      xbins, ybins, zbins, maxOrder, (int)(border * 200), suffix.Data()));
  saveCanvas(*cz1new, Form("EffProj_%i_%i_%i_Phislices_SpH%io_dp%i_%s", xbins,
                      ybins, zbins, maxOrder, (int)(border * 200), suffix.Data()));
}
