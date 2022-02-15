#define TnPAnalysis_cxx
#include "TnPAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
using namespace RooFit;
using namespace std;

void TnPAnalysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L TnPAnalysis.C
//      root> TnPAnalysis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

   if (fChain == 0) return;
   
              // TCanvas *c1= new TCanvas();
	      // TCanvas *c2= new TCanvas();
	      // TCanvas *c3= new TCanvas();
       bool RunSystematic = true;
       vector<TString> systematicVar = {"nominal"};
       if(RunSystematic)
		     {
			     systematicVar.push_back("tagpt_up");
			     systematicVar.push_back("tagpt_down");
			     systematicVar.push_back("jpsi_up");
			     systematicVar.push_back("jpsi_down");
		     }
		
       for(int i=0; i<systematicVar.size();i++)
       {
	       double ptTag = 3;
	       double JpsiMassL = 2.4;
	       double JpsiMassR = 3.7;
	      
	       if(systematicVar.at(i).Contains("tagPt_up")) ptTag = 3;
	       if(systematicVar.at(i).Contains("tagPt_down")) ptTag =2.5 ;
	       if(systematicVar.at(i).Contains("Jpsi_up")) {JpsiMassL = 2.8; JpsiMassR = 3.7;}
	       if(systematicVar.at(i).Contains("Jpsi_down")) {JpsiMassL = 2.4; JpsiMassR = 3.3;}
	       cout<<systematicVar.at(i).Data()<<endl;
	       cout<<" ptTag : "<<ptTag<<" , JpsiMassL : "<<JpsiMassL<<" , JpsiMassR : "<<JpsiMassR<<endl;
	       TString output = "Muon_sf";
	       output +="efficiency_";
	       output += systematicVar.at(i);
	       output +=".root";
	        Long64_t nentries = fChain->GetEntriesFast();
	        std::cout<<nentries<<"\n";
		double eta_bins[14] = {-2.4,-2.1,-1.6,-1.2,-0.8,-0.3,-0.2,0.2,0.3,0.8,1.2,1.6,2.1,2.4};
	        double pt_bins_Mu0_notOpen[10] = {0,2,4,6,8,10,12,14,17,20};
	        TFile *file = new TFile(output.Data(),"RECREATE");
                
                
              TH1F *h_Jpsi_notOpen_total = new TH1F("h_Jpsi_notOpen_total", "h_Jpsi_notOpen_total", 50, 2.4, 3.7);
               TH1F *h_Jpsi_notOpen_pass = new TH1F("h_Jpsi_notOpen_pass", "h_Jpsi_notOpen_pass", 50, 2.4, 3.7);
	       TH1F *h_Mu0_notOpen_pt_total = new TH1F("Mu0_notOpen_pt_total","Mu0_notOpen_pt",9,pt_bins_Mu0_notOpen);
               TH1F *h_Mu0_notOpen_eta_total = new TH1F("Mu0_notOpen_eta_total","Mu0_notOpen_eta",13,eta_bins);
	       TH2F *h_Mu0_notOpen_pt_eta_total = new TH2F("Mu0_notOpen_pt_eta_total","Mu0_notOpen_pt_eta",13,eta_bins,9,pt_bins_Mu0_notOpen);
	       TH1F *h_Mu0_notOpen_pt_pass = new TH1F("Mu0_notOpen_pt_pass","Mu0_notOpen_pt",9,pt_bins_Mu0_notOpen);
	       TH1F *h_Mu0_notOpen_eta_pass = new TH1F("Mu0_notOpen_eta_pass","Mu0_notOpen_eta",13,eta_bins);
	       TH2F *h_Mu0_notOpen_pt_eta_pass = new TH2F("Mu0_notOpen_pt_eta_pass","Mu0_notOpen_pt_eta",13,eta_bins,9,pt_bins_Mu0_notOpen);
               
	       h_Mu0_notOpen_pt_total->Sumw2();
	       h_Mu0_notOpen_eta_total->Sumw2();
	       h_Mu0_notOpen_pt_eta_total->Sumw2();
	       h_Mu0_notOpen_pt_pass->Sumw2();
	       h_Mu0_notOpen_eta_pass->Sumw2();
	       h_Mu0_notOpen_pt_eta_pass->Sumw2();
	       h_Jpsi_notOpen_total->Sumw2();
	       h_Jpsi_notOpen_pass->Sumw2();
	       




      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry< nentries;jentry++) { 
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry % 50000 ==0)cout<<"Number of events processed : "<<jentry<<endl;
      if(nMu != 2) continue;
         int first  = rand()%2;
         int second = (first+1)%2;
         if(mu_charge->at(first) * mu_charge->at(second)>0)continue;
         bool soft_ID = mu_d0->at(first) < 0.3 && mu_dz->at(first) < 20 && mu_TrkLayers->at(first) > 5 && mu_PixelLayers->at(first) > 0 && mu_TrkQuality->at(first) == 1 ;
         
         bool tag_MuKin = mu_pt->at(first)>ptTag && fabs(mu_eta->at(first))<2.4;
         bool tag_ImpactID = mu_d0->at(first) > 1 && mu_dz->at(first) > 20 ;
         bool tag_samID =  mu_Stations->at(first) > 1; 
         bool tag_TriggerMatch = passFilter_trigNotOpen->at(first);
         if(!( tag_MuKin && tag_TriggerMatch && soft_ID ))continue;
         bool probe_MuId = GGMM_Muon_Def(second, mu_pt->at(second));
         bool probe_MuKin = fabs(mu_eta->at(second))<2.4;    
         if(!(probe_MuKin)) continue;
	 TLorentzVector tag_muLV, probe_muLV, Jpsi_candLV;
	 tag_muLV.SetPtEtaPhiE(mu_pt->at(first), mu_eta->at(first), mu_phi->at(first), mu_energy->at(first));
	 probe_muLV.SetPtEtaPhiE(mu_pt->at(second), mu_eta->at(second), mu_phi->at(second), mu_energy->at(second));
	 Jpsi_candLV = tag_muLV + probe_muLV;
	 if (Jpsi_candLV.M()<JpsiMassL || Jpsi_candLV.M() >JpsiMassR) continue;
         //std::cout<<" Jpsi mass " <<Jpsi_candLV.M()<<"\n";
	      h_Jpsi_notOpen_total->Fill(Jpsi_candLV.M());
	      h_Mu0_notOpen_pt_total->Fill(mu_pt->at(second));
         h_Mu0_notOpen_eta_total->Fill(mu_eta->at(second));
         h_Mu0_notOpen_pt_eta_total->Fill(mu_eta->at(second),mu_pt->at(second));	  
	      if(passFilter_trigNotOpen->at(second)){
        
         h_Jpsi_notOpen_pass->Fill(Jpsi_candLV.M());
         h_Mu0_notOpen_pt_pass->Fill(mu_pt->at(second));
         h_Mu0_notOpen_eta_pass->Fill(mu_eta->at(second));
         h_Mu0_notOpen_pt_eta_pass->Fill(mu_eta->at(second),mu_pt->at(second));
         
        
         //Jpsi_fit_total(h_Jpsi_notOpen_pass, "Pass mass" , c1);
         //Jpsi_fit_total(h_Jpsi_notOpen_total, "Total mass" , c1);
         // if((Jpsi_candLV.M() > SB1_L && Jpsi_candLV.M() < SB1_H) || (Jpsi_candLV.M() > SB2_L && Jpsi_candLV.M() < SB2_H))continue;
         
	 }
	 
	 
   }
   file->Write();
}
}


bool TnPAnalysis::GGMM_Muon_Def(int i, double pt)
{
          bool pass_dz = mu_dz->at(i) < 30;
          if (!pass_dz) return false;
          if (pt <= 1){
      	  if (mu_d0->at(i) > 1) return false;
	  }
	  else {
		  if (mu_d0->at(i) > 1.5) return false;
	  }
	return true;
}
void TnPAnalysis::Jpsi_fit_total(TH1F * Hist, string category, TCanvas *c)
{

 RooRealVar *Jpsi_All = new RooRealVar("Jpsi_All","  m_{J/#psi} (#mu^{+}#mu^{-}) ) (GeV)",2.4, 3.6) ;
 RooDataHist * jpsia = new RooDataHist("jpsia","jpsia", *Jpsi_All , Import(*Hist)) ;
 RooRealVar *mean = new RooRealVar("mean","mean",3.09691,2.9,3.1) ;
 RooRealVar *sigma = new RooRealVar("sigma","sigma",0.12,0.,3.5) ;
 RooGaussian *gauss = new RooGaussian("gauss","gauss",*Jpsi_All,  *mean, *sigma) ;
 RooRealVar *a0 = new RooRealVar("a0","a0",1.0,-1.0,1.0) ;
 RooRealVar *a1 = new RooRealVar("a1","a1",0.2,0,100.) ;
 RooChebychev *cbcv = new RooChebychev("cbcv","chebychev",*Jpsi_All,RooArgSet(*a0,*a1)) ; 
 RooRealVar *ns = new RooRealVar("ns","number signal events",0, Hist->GetEntries());  
 RooRealVar *frac = new RooRealVar("frac","frac",0.4,0.0,1.0); 
 RooAddPdf *siggaucbcv = new RooAddPdf("siggaucbcv","siggaucbcv",RooArgList(*gauss , *cbcv),RooArgList(*frac));
 RooExtendPdf *esig_all = new RooExtendPdf("esig_all","extended signal p.d.f",*siggaucbcv, *ns,"signalRange") ;
 
 RooFitResult* fitRes = esig_all->fitTo(*jpsia,Save(),NumCPU(8));
 fitRes->Print("v");
 gStyle->SetOptStat(0) ;
 gStyle->SetPalette(1) ;
 TH2* hcorr = fitRes->correlationHist() ;
 TCanvas* ccor = new TCanvas("Corr Matrix","M/ct/cterr correaltion matrix",800,400) ;
 gPad->SetLeftMargin(0.15) ; 
 hcorr->GetYaxis()->SetTitleOffset(1.4) ; 
 hcorr->Draw("colz") ;
 ccor->Print(("correlation_Matrix_"+category+"_data.png").c_str(),"png");
 
 
 RooPlot* frame = Jpsi_All->frame(Title("All m_{J/#psi} (#mu^{+}#mu^{-}) ) (GeV)")) ;
 jpsia->plotOn(frame) ; 
 esig_all->fitTo(*jpsia , Extended(1));
 esig_all->plotOn(frame);
 esig_all->paramOn(frame);
 RooPlot* pullframe = Jpsi_All->frame(RooFit::Title("Mass pull"));
 RooHist* hpull1 = frame->pullHist();
 pullframe->addPlotable(hpull1,"P0") ;
 pullframe->SetMinimum(-3) ;
 pullframe->SetMaximum(+3) ;
 pullframe->SetYTitle("pull");
 pullframe->SetMarkerStyle(20);
 pullframe->SetNdivisions(10); 
 Double_t chisquare = frame->chiSquare();
 cout<<"Chi square of total mass distribution fit is :"<< chisquare<< endl;
 esig_all->plotOn(frame, Components(*gauss),RooFit::Name("signal "),  LineColor(3), LineWidth(2), LineStyle(4));
 esig_all->plotOn(frame, Components(*cbcv), RooFit::Name("background "),LineColor(2), LineWidth(2), LineStyle(2));
 TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
 leg->AddEntry(frame->findObject("signal"),"J/#psi (#to #mu#mu)","l");
 leg->AddEntry(frame->findObject("background"),"Background","l");
 c = new TCanvas("c", "c",0,0,600,600);
    TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
    pad1->SetBottomMargin(0.00001);
    pad1->SetBorderMode(0);
    pad2->SetTopMargin(0.00001);
    pad2->SetBottomMargin(0.1);
    pad2->SetBorderMode(0);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    gStyle->SetOptTitle(0);
    c->SetFillColor(0);
    c->SetBorderSize(2);
    c->SetLeftMargin(0.1422222);
    c->SetRightMargin(0.04444445);
    frame->SetStats(0);
    frame->Draw();
    leg->Draw("same");
    auto cms1 = new TLatex(5.15, 1400, "#bf{CMS} #it{Preliminary} 2015, #sqrt{s} = 5.02 TeV");
    cms1->SetNDC(false);
    cms1->SetTextColor(12);
    cms1->SetTextFont(42);
    cms1->SetTextSize(0.055);
    cms1-> Draw();
    pad2->cd();
    pullframe->SetStats(0);
    pullframe->Draw();
    c->cd();
    c->Print(("jpsimass"+category+"_data.png").c_str(), "png");
       
}


