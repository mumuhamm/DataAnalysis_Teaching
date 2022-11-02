RooRealVar x("x","MM_{(#pi #pi)} GeV/c^{2}",3.04,3.17) ;

RooDataHist dh("dh","dataset",x,h1);
RooPlot* frame = x.frame(Title("@4.26 GeV :: e^{+}e^{-} -> J/#psi #pi^{+} #pi^{-}")) ;
x.setRange("Range1",3.04,3.089);
x.setRange("Range2",3.11,3.17);
x.setRange("fullrange",3.04,3.17);


dh.plotOn(frame,Name("dh")) ;


RooRealVar a0("a0","a0",-2.,2.) ;
RooRealVar a1("a1","a1",-2.,2.) ;
RooRealVar a2("a2","a2",-2.,2.) ;

RooChebychev bkg("bkg","background p.d.f.",x,RooArgList(a0));
RooChebychev bkg1("bkg1","background p.d.f.",x,RooArgList(a0,a1));
RooChebychev bkg2("bkg2","background p.d.f.",x,RooArgList(a0,a1,a2));


RooRealVar nsig("N_{SIG}","signal events",0,1000000);
RooRealVar nbkg("N_{BKG}","signal background even0ts",0,100000000);

RooAddPdf all("all","all",RooArgList(bkg),RooArgList(nbkg));
RooFitResult* r = all.fitTo(dh,Extended(kTRUE),Save()) ;  all.plotOn(frame,Range("Range1,Range2"),NormRange("Range1,Range2"),LineColor(kRed),Components(bkg));
all.paramOn(frame,Layout(0.5,0.90,0.55));

TH1 *hbkg = bkg.createHistogram("hbkg",x,Binning(150,3.04,3.17));
frame->Draw();

c2->Divide(1,2);
c2->cd(1);
hbkg->Scale(nbkg.getVal());
hbkg->Sumw2();
hbkg->Draw();
hbkg->SetTitle("Background histogram");

c2->cd(2);
TH1F *h4230 = (TH1F*)hist4230->Clone("h4230");
h4230->Sumw2();
h4230->Add(hbkg,-1);
h4230->Draw();
h4230->GetXaxis()->SetTitle("MM_{(#pi#pi)} GeV/c^{2}");
TFile *file = new TFile("4230_2RG.root", "RECREATE");
h4230->Write();
file->Close();
