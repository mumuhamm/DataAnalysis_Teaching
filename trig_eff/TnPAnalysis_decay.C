#define TnPAnalysis_cxx
#include "TnPAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void TnPAnalysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L TnPAnalysis_decay.C
//      root> TnPAnalysis_decay t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

   if (fChain == 0) return;
   
   
      
	        Long64_t nentries = fChain->GetEntriesFast();
	        std::cout<<nentries<<"\n";
                 double ptTag = 1;
               const double SB1_L=2.5;
               const double SB1_H=2.8; 
               const double SR_L=3;
               const double SR_H=3.2;
               const double SB2_L=3.3;
               const double SB2_H=3.6;
           
              Double_t Jpsi_mass, Jpsi_pt, Jpsi_eta, Jpsi_phi, muon1_pt, muon1_eta, muon1_phi, muon2_pt, muon2_eta, muon2_phi;   
       
    TFile *f = new TFile("VariablesJpsi_15PbPb_data.root","recreate");
    TTree *fit = new TTree("fit","selected ntcltestle");
    fit->Branch("Jpsi_mass",&Jpsi_mass,"Jpsi_mass/D");
    fit->Branch("Jpsi_eta",&Jpsi_eta,"Jpsi_eta/D");
    fit->Branch("Jpsi_pt",&Jpsi_pt,"Jpsi_pt/D");
    fit->Branch("Jpsi_phi",&Jpsi_phi,"Jpsi_phi/D");
    fit->Branch("muon1_pt",&muon1_pt,"muon1_pt/D");
    fit->Branch("muon1_eta",&muon1_eta,"muon1_eta/D");
    fit->Branch("muon1_phi",&muon1_phi,"muon1_phi/D");
     fit->Branch("muon2_pt",&muon2_pt,"muon2_pt/D");
    fit->Branch("muon2_eta",&muon2_eta,"muon2_eta/D");
    fit->Branch("muon2_phi",&muon2_phi,"muon2_phi/D");
    
	       




      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) { 
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry % 50000 ==0)cout<<"Number of events processed : "<<jentry<<endl;
      if(nMu != 2) continue;
         int first  = rand()%2;
         int second = (first+1)%2;
         if(mu_charge->at(first) * mu_charge->at(second)>0)continue;
         bool tag_MuKin = mu_pt->at(first)>ptTag && fabs(mu_eta->at(first))<2.4;
          if(!( tag_MuKin ))continue;   
         muon1_pt = mu_pt->at(first);
         muon1_eta = mu_eta->at(first);
         muon1_phi = mu_phi->at(first);
          bool probe_MuId = GGMM_Muon_Def(second, mu_pt->at(second));
         bool probe_MuKin = fabs(mu_eta->at(second))<2.4;
         if(!(probe_MuId && probe_MuKin)) continue;
         muon2_pt = mu_pt->at(second);
         muon2_eta = mu_eta->at(second);
         muon2_phi = mu_phi->at(second);
         if(!( probe_MuKin ))continue; 
	 TLorentzVector tag_muLV, probe_muLV, Jpsi_candLV;
	 tag_muLV.SetPtEtaPhiE(mu_pt->at(first), mu_eta->at(first), mu_phi->at(first), mu_energy->at(first));
	 probe_muLV.SetPtEtaPhiE(mu_pt->at(second), mu_eta->at(second), mu_phi->at(second), mu_energy->at(second));
	 Jpsi_candLV = tag_muLV + probe_muLV;
	 if (Jpsi_candLV.M()<2.4 || Jpsi_candLV.M() >3.6) continue;
	 Jpsi_mass = Jpsi_candLV.M();
	 //std::cout<<" Jpsi mass " << Jpsi_candLV.M()<<"\n";
	 Jpsi_pt =  Jpsi_candLV.Pt();
	 Jpsi_eta = Jpsi_candLV.Eta();
         Jpsi_phi = Jpsi_candLV.Phi();
	// if((Jpsi_candLV.M() > SB1_L && Jpsi_candLV.M() < SB1_H) || (Jpsi_candLV.M() > SB2_L && Jpsi_candLV.M() < SB2_H))continue;
	fit->Fill();
         }
         fit->Write();
         fit->Print();
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

