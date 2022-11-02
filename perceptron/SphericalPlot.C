

Bool_t SphericalPlotWrong ()
{
	Double_t borders[19];
	// Generate bin borders for latitude
	for (UInt_t i=0; i<19; i++)
		borders[i] = TMath::RadToDeg() * TMath::ACos(1. - i*2./18.);
	TH2I* Hist = new TH2I ("SphericalWrong", "Latitude X, Longitude Y", 18, borders, 36, 0, 360);
	for (UInt_t i=1; i<18; i++)
		for (UInt_t j=1; j<36; j++)
			Hist->SetBinContent(i, j, 1.5);
	Hist->Draw("LEGO SPH");
	return 0;
}

Bool_t SphericalPlotCorrect ()
{
	Double_t borders[19];
	// Generate bin borders for latitude
	for (UInt_t i=0; i<19; i++)
		borders[i] = TMath::RadToDeg() * TMath::ACos(1. - i*2./18.);
	TH2I* Hist = new TH2I ("SphericalCorrect", "Latitude Y, Longitude X", 36, 0, 360, 18, borders);
	for (UInt_t i=1; i<18; i++)
		for (UInt_t j=1; j<36; j++)
			Hist->SetBinContent(j, i, 1.5);
	Hist->Draw("LEGO SPH");
	return 0;
}