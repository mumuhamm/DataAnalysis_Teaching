//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 14 06:50:05 2021 by ROOT version 6.14/09
// from TTree eff_muon/eff_muon
// found on file: data_hlt_15PbPb.root
//////////////////////////////////////////////////////////

#ifndef TnPAnalysis_h
#define TnPAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class TnPAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TTree          *tree;
   TChain         *Chain; 
   TFile          *file_bin_00_03;
   TFile          *file_bin_03_06_1;
   TFile          *file_bin_03_06_2;
   TFile          *file_bin_03_06_3;
   TFile          *file_bin_06_09;
   TFile          *file_bin_09_12;
   TFile          *file_bin_12_15;
   TFile          *file_bin_15_30;
   TFile          *file_bin_30_Inf;
   RooRealVar * weight_bin_00_03 = new RooRealVar("weight_00_03", "weight_00_03", 0.2469, 0.2469-0.0214, 0.2469+0.0214);
   RooRealVar * weight_bin_03_06 = new RooRealVar("weight_03_06", "weight_03_06", 0.1721, 0.1721-0.0157, 0.1721+0.0157);
   RooRealVar * weight_bin_06_09 = new RooRealVar("weight_06_09", "weight_06_09", 0.0204, 0.0204-0.00202, 0.0204+0.00202);
   RooRealVar * weight_bin_09_12 = new RooRealVar("weight_09_12", "weight_09_12", 0.00403, 0.00403-0.000402, 0.00403+0.000402);
   RooRealVar * weight_bin_12_15 = new RooRealVar("weight_12_15", "weight_12_15", 0.001211, 0.001211-0.000402, 0.001211+0.000402);
   RooRealVar * weight_bin_15_30 = new RooRealVar("weight_15_30", "weight_15_30", 0.0007182, 0.0007182-0.0000718, 0.0007182+0.0000718);
   RooRealVar * weight_bin_30_Inf = new RooRealVar("weight_30_Inf", "weight_30_Inf", 0.0000328, 0.0000328-0.00000328, 0.0000328+0.00000328);

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           pvNTracks;
   Int_t           good_vertices;
   Int_t           nPV;
   Double_t        genWeight;
   vector<string>  *triggerPath;
   vector<bool>    *triggerDecision;
   vector<bool>    *passFilter_trigNotOpen;
   Double_t        muon_soft_trk_pt;
   Double_t        muon_soft_trk_eta;
   Double_t        muon_soft_trk_phi;
   Bool_t          triggerbit_HLT_UPC_SingleMu0;
   Double_t        BSx;
   Double_t        BSy;
   Double_t        BSz;
   Double_t        BSdx;
   Double_t        BSdy;
   Double_t        BSdz;
   Double_t        BSdxdz;
   Double_t        BSsigmaZ;
   Double_t        BSdsigmaZ;
   Double_t        MuonsDCA;
   Double_t        vtxProb_dimuon;
   Double_t        CosAlpha;
   Int_t           nMu;
   vector<double>  *mu_pt;
   vector<double>  *mu_eta;
   vector<double>  *mu_phi;
   vector<double>  *mu_energy;
   vector<int>     *mu_charge;
   vector<int>     *mu_type;
   vector<double>  *mu_d0;
   vector<double>  *mu_dz;
   vector<double>  *mu_SIP;
   vector<double>  *mu_Chi2NDF;
   vector<double>  *mu_InnerD0;
   vector<double>  *mu_InnerDz;
   vector<int>     *mu_TrkLayers;
   vector<int>     *mu_PixelLayers;
   vector<int>     *mu_PixelHits;
   vector<int>     *mu_MuonHits;
   vector<int>     *mu_Stations;
   vector<int>     *mu_Matches;
   vector<int>     *mu_TrkQuality;
   vector<double>  *mu_IsoTrk;
   vector<double>  *mu_PFChIso;
   vector<double>  *mu_PFPhoIso;
   vector<double>  *mu_PFNeuIso;
   vector<double>  *mu_PFPUIso;
   vector<double>  *mu_PFCHIso03;
   vector<double>  *mu_PFPhoIso03;
   vector<double>  *mu_PFNeuIso03;
   vector<double>  *mu_PFPUIso03;
   vector<double>  *mu_InnervalidFraction;
   vector<double>  *mu_segmentCompatibility;
   vector<double>  *mu_chi2LocalPosition;
   vector<double>  *mu_trkKink;
   vector<double>  *mu_BestTrkPtError;
   vector<double>  *mu_BestTrkPt;
   vector<int>     *mu_BestTrkType;
   vector<double>  *genMuon_pt;
   vector<double>  *genMuon_eta;
   vector<double>  *genMuon_phi;
   vector<double>  *genMuon_energy;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_pvNTracks;   //!
   TBranch        *b_good_vertices;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_triggerPath;   //!
   TBranch        *b_triggerDecision;   //!
   TBranch        *b_passFilter_trigNotOpen;   //!
   TBranch        *b_muon_soft_trk_pt;   //!
   TBranch        *b_muon_soft_trk_eta;   //!
   TBranch        *b_muon_soft_trk_phi;   //!
   TBranch        *b_triggerbit_HLT_UPC_SingleMu0;   //!
   TBranch        *b_BSx;   //!
   TBranch        *b_BSy;   //!
   TBranch        *b_BSz;   //!
   TBranch        *b_BSdx;   //!
   TBranch        *b_BSdy;   //!
   TBranch        *b_BSdz;   //!
   TBranch        *b_BSdxdz;   //!
   TBranch        *b_BSsigmaZ;   //!
   TBranch        *b_BSdsigmaZ;   //!
   TBranch        *b_MuonsDCA;   //!
   TBranch        *b_vtxProb_dimuon;   //!
   TBranch        *b_CosAlpha;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_energy;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_type;   //!
   TBranch        *b_mu_d0;   //!
   TBranch        *b_mu_dz;   //!
   TBranch        *b_mu_SIP;   //!
   TBranch        *b_mu_Chi2NDF;   //!
   TBranch        *b_mu_InnerD0;   //!
   TBranch        *b_mu_InnerDz;   //!
   TBranch        *b_mu_TrkLayers;   //!
   TBranch        *b_mu_PixelLayers;   //!
   TBranch        *b_mu_PixelHits;   //!
   TBranch        *b_mu_MuonHits;   //!
   TBranch        *b_mu_Stations;   //!
   TBranch        *b_mu_Matches;   //!
   TBranch        *b_mu_TrkQuality;   //!
   TBranch        *b_mu_IsoTrk;   //!
   TBranch        *b_mu_PFChIso;   //!
   TBranch        *b_mu_PFPhoIso;   //!
   TBranch        *b_mu_PFNeuIso;   //!
   TBranch        *b_mu_PFPUIso;   //!
   TBranch        *b_mu_PFCHIso03;   //!
   TBranch        *b_mu_PFPhoIso03;   //!
   TBranch        *b_mu_PFNeuIso03;   //!
   TBranch        *b_mu_PFPUIso03;   //!
   TBranch        *b_mu_InnervalidFraction;   //!
   TBranch        *b_mu_segmentCompatibility;   //!
   TBranch        *b_mu_chi2LocalPosition;   //!
   TBranch        *b_mu_trkKink;   //!
   TBranch        *b_mu_BestTrkPtError;   //!
   TBranch        *b_mu_BestTrkPt;   //!
   TBranch        *b_mu_BestTrkType;   //!
   TBranch        *b_genMuon_pt;   //!
   TBranch        *b_genMuon_eta;   //!
   TBranch        *b_genMuon_phi;   //!
   TBranch        *b_genMuon_energy;   //!

   TnPAnalysis(TTree *tree=0);
   virtual ~TnPAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   bool    GGMM_Muon_Def(int i, double pt);
   void Jpsi_fit_total(TH1F * Hist, string category, TCanvas *c);
};

#endif

#ifdef TnPAnalysis_cxx
TnPAnalysis::TnPAnalysis(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
/*   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Old_samples/data_hlt_15PbPb.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Old_samples/data_hlt_15PbPb.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("Old_samples/data_hlt_15PbPb.root:/eff");
      dir->GetObject("eff_muon",tree);

   }
   Init(tree);
 */
   
   // Integrate filter efficiency
   file_bin_00_03 = TFile::Open("/eos/user/m/mumuhamm/HLT_data_mc/mc/promptjpsi_00_03_Hydjet.root");
   TTree *t1 = (TTree*)file_bin_00_03->Get("eff/eff_muon");
   Int_t n1 = t1->GetEntries();
   tree = (TTree*)file_bin_00_03->Get("eff/eff_muon");
   tree->SetWeight(weight_bin_00_03->getValV()/n1,"global");
   tree->AutoSave();
   
   file_bin_03_06_1 = TFile::Open("/eos/user/m/mumuhamm/HLT_data_mc/mc/ForwardMC/Pythia8_JpsiMM_ptJpsi_03_06_Hydjet_MB/crab_HLT_Efficiency/promptjpsi_03_06_Hydjet_1st.root");
   TTree *t2 = (TTree*)file_bin_03_06_1->Get("eff/eff_muon");
   Int_t n2 = t2->GetEntries();
   tree = (TTree*)file_bin_03_06_1->Get("eff/eff_muon");
   tree->SetWeight(weight_bin_03_06->getValV()/n2,"global");
   tree->AutoSave();
   
   file_bin_03_06_2 = TFile::Open("/eos/user/m/mumuhamm/HLT_data_mc/mc/ForwardMC/Pythia8_JpsiMM_ptJpsi_03_06_Hydjet_MB/crab_HLT_Efficiency/promptjpsi_03_06_Hydjet_2nd.root");
   TTree *t3 = (TTree*)file_bin_03_06_2->Get("eff/eff_muon");
   Int_t n3 = t3->GetEntries();
   tree = (TTree*)file_bin_03_06_2->Get("eff/eff_muon");
   tree->SetWeight(weight_bin_03_06->getValV()/n3,"global");
   tree->AutoSave();
   
   file_bin_03_06_3 = TFile::Open("/eos/user/m/mumuhamm/HLT_data_mc/mc/ForwardMC/Pythia8_JpsiMM_ptJpsi_03_06_Hydjet_MB/crab_HLT_Efficiency/promptjpsi_03_06_Hydjet_3rd.root");
   TTree *t4 = (TTree*)file_bin_03_06_3->Get("eff/eff_muon");
   Int_t n4 = t4->GetEntries();
   tree = (TTree*)file_bin_03_06_3->Get("eff/eff_muon");
   tree->SetWeight(weight_bin_03_06->getValV()/n4,"global");
   tree->AutoSave();
   
   file_bin_06_09 = TFile::Open("/eos/user/m/mumuhamm/HLT_data_mc/mc/promptjpsi_06_09_Hydjet.root");
   TTree *t5 = (TTree*)file_bin_06_09->Get("eff/eff_muon");
   Int_t n5 = t5->GetEntries();
   tree = (TTree*)file_bin_06_09->Get("eff/eff_muon");
   tree->SetWeight(weight_bin_06_09->getValV()/n5,"global");
   tree->AutoSave();
   
   file_bin_09_12 = TFile::Open("/eos/user/m/mumuhamm/HLT_data_mc/mc/promptjpsi_09_12_Hydjet.root");
   TTree *t6 = (TTree*)file_bin_09_12->Get("eff/eff_muon");
   Int_t n6 = t6->GetEntries();
   tree = (TTree*)file_bin_09_12->Get("eff/eff_muon");
   tree->SetWeight(weight_bin_09_12->getValV()/n6,"global");
   tree->AutoSave();
   
   file_bin_12_15 = TFile::Open("/eos/user/m/mumuhamm/HLT_data_mc/mc/promptjpsi_12_15_Hydjet.root");
   TTree *t7 = (TTree*)file_bin_12_15->Get("eff/eff_muon");
   Int_t n7 = t7->GetEntries();
   tree = (TTree*)file_bin_12_15->Get("eff/eff_muon");
   tree->SetWeight(weight_bin_12_15->getValV()/n7,"global");
   tree->AutoSave();
   
   file_bin_15_30 = TFile::Open("/eos/user/m/mumuhamm/HLT_data_mc/mc/promptjpsi_15_30_Hydjet.root");
   TTree *t8 = (TTree*)file_bin_15_30->Get("eff/eff_muon");
   Int_t n8 = t8->GetEntries();
   tree = (TTree*)file_bin_15_30->Get("eff/eff_muon");
   tree->SetWeight(weight_bin_15_30->getValV()/n8,"global");
   tree->AutoSave();
   
   file_bin_30_Inf = TFile::Open("/eos/user/m/mumuhamm/HLT_data_mc/mc/promptjpsi_30_Inf_Hydjet.root");
   TTree *t9 = (TTree*)file_bin_30_Inf ->Get("eff/eff_muon");
   Int_t n9 = t9->GetEntries();
   tree = (TTree*)file_bin_30_Inf ->Get("eff/eff_muon");
   tree->SetWeight(weight_bin_30_Inf->getValV()/n9,"global");
   tree->AutoSave();
   
   // Add files
   Chain = new TChain("eff/eff_muon");
   Chain->Add("/eos/user/m/mumuhamm/HLT_data_mc/mc/promptjpsi_00_03_Hydjet.root");
   Chain->Add("/eos/user/m/mumuhamm/HLT_data_mc/mc/ForwardMC/Pythia8_JpsiMM_ptJpsi_03_06_Hydjet_MB/crab_HLT_Efficiency/promptjpsi_03_06_Hydjet_1st.root");
   Chain->Add("/eos/user/m/mumuhamm/HLT_data_mc/mc/ForwardMC/Pythia8_JpsiMM_ptJpsi_03_06_Hydjet_MB/crab_HLT_Efficiency/promptjpsi_03_06_Hydjet_2nd.root");
   Chain->Add("/eos/user/m/mumuhamm/HLT_data_mc/mc/ForwardMC/Pythia8_JpsiMM_ptJpsi_03_06_Hydjet_MB/crab_HLT_Efficiency/promptjpsi_03_06_Hydjet_3rd.root");
   Chain->Add("/eos/user/m/mumuhamm/HLT_data_mc/mc/promptjpsi_06_09_Hydjet.root");
   Chain->Add("/eos/user/m/mumuhamm/HLT_data_mc/mc/promptjpsi_09_12_Hydjet.root");
   Chain->Add("/eos/user/m/mumuhamm/HLT_data_mc/mc/promptjpsi_12_15_Hydjet.root");
   Chain->Add("/eos/user/m/mumuhamm/HLT_data_mc/mc/promptjpsi_15_30_Hydjet.root");
   Chain->Add("/eos/user/m/mumuhamm/HLT_data_mc/mc/promptjpsi_30_Inf_Hydjet.root");
    Init(Chain);
}

TnPAnalysis::~TnPAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TnPAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TnPAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TnPAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   triggerPath = 0;
   triggerDecision = 0;
   passFilter_trigNotOpen = 0;
   mu_pt = 0;
   mu_eta = 0;
   mu_phi = 0;
   mu_energy = 0;
   mu_charge = 0;
   mu_type = 0;
   mu_d0 = 0;
   mu_dz = 0;
   mu_SIP = 0;
   mu_Chi2NDF = 0;
   mu_InnerD0 = 0;
   mu_InnerDz = 0;
   mu_TrkLayers = 0;
   mu_PixelLayers = 0;
   mu_PixelHits = 0;
   mu_MuonHits = 0;
   mu_Stations = 0;
   mu_Matches = 0;
   mu_TrkQuality = 0;
   mu_IsoTrk = 0;
   mu_PFChIso = 0;
   mu_PFPhoIso = 0;
   mu_PFNeuIso = 0;
   mu_PFPUIso = 0;
   mu_PFCHIso03 = 0;
   mu_PFPhoIso03 = 0;
   mu_PFNeuIso03 = 0;
   mu_PFPUIso03 = 0;
   mu_InnervalidFraction = 0;
   mu_segmentCompatibility = 0;
   mu_chi2LocalPosition = 0;
   mu_trkKink = 0;
   mu_BestTrkPtError = 0;
   mu_BestTrkPt = 0;
   mu_BestTrkType = 0;
   genMuon_pt = 0;
   genMuon_eta = 0;
   genMuon_phi = 0;
   genMuon_energy = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("pvNTracks", &pvNTracks, &b_pvNTracks);
   fChain->SetBranchAddress("good_vertices", &good_vertices, &b_good_vertices);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("triggerPath", &triggerPath, &b_triggerPath);
   fChain->SetBranchAddress("triggerDecision", &triggerDecision, &b_triggerDecision);
   fChain->SetBranchAddress("passFilter_trigNotOpen", &passFilter_trigNotOpen, &b_passFilter_trigNotOpen);
   fChain->SetBranchAddress("muon_soft_trk_pt", &muon_soft_trk_pt, &b_muon_soft_trk_pt);
   fChain->SetBranchAddress("muon_soft_trk_eta", &muon_soft_trk_eta, &b_muon_soft_trk_eta);
   fChain->SetBranchAddress("muon_soft_trk_phi", &muon_soft_trk_phi, &b_muon_soft_trk_phi);
   fChain->SetBranchAddress("triggerbit_HLT_UPC_SingleMu0", &triggerbit_HLT_UPC_SingleMu0, &b_triggerbit_HLT_UPC_SingleMu0);
   fChain->SetBranchAddress("BSx", &BSx, &b_BSx);
   fChain->SetBranchAddress("BSy", &BSy, &b_BSy);
   fChain->SetBranchAddress("BSz", &BSz, &b_BSz);
   fChain->SetBranchAddress("BSdx", &BSdx, &b_BSdx);
   fChain->SetBranchAddress("BSdy", &BSdy, &b_BSdy);
   fChain->SetBranchAddress("BSdz", &BSdz, &b_BSdz);
   fChain->SetBranchAddress("BSdxdz", &BSdxdz, &b_BSdxdz);
   fChain->SetBranchAddress("BSsigmaZ", &BSsigmaZ, &b_BSsigmaZ);
   fChain->SetBranchAddress("BSdsigmaZ", &BSdsigmaZ, &b_BSdsigmaZ);
   fChain->SetBranchAddress("MuonsDCA", &MuonsDCA, &b_MuonsDCA);
   fChain->SetBranchAddress("vtxProb_dimuon", &vtxProb_dimuon, &b_vtxProb_dimuon);
   fChain->SetBranchAddress("CosAlpha", &CosAlpha, &b_CosAlpha);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
   fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_energy", &mu_energy, &b_mu_energy);
   fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
   fChain->SetBranchAddress("mu_type", &mu_type, &b_mu_type);
   fChain->SetBranchAddress("mu_d0", &mu_d0, &b_mu_d0);
   fChain->SetBranchAddress("mu_dz", &mu_dz, &b_mu_dz);
   fChain->SetBranchAddress("mu_SIP", &mu_SIP, &b_mu_SIP);
   fChain->SetBranchAddress("mu_Chi2NDF", &mu_Chi2NDF, &b_mu_Chi2NDF);
   fChain->SetBranchAddress("mu_InnerD0", &mu_InnerD0, &b_mu_InnerD0);
   fChain->SetBranchAddress("mu_InnerDz", &mu_InnerDz, &b_mu_InnerDz);
   fChain->SetBranchAddress("mu_TrkLayers", &mu_TrkLayers, &b_mu_TrkLayers);
   fChain->SetBranchAddress("mu_PixelLayers", &mu_PixelLayers, &b_mu_PixelLayers);
   fChain->SetBranchAddress("mu_PixelHits", &mu_PixelHits, &b_mu_PixelHits);
   fChain->SetBranchAddress("mu_MuonHits", &mu_MuonHits, &b_mu_MuonHits);
   fChain->SetBranchAddress("mu_Stations", &mu_Stations, &b_mu_Stations);
   fChain->SetBranchAddress("mu_Matches", &mu_Matches, &b_mu_Matches);
   fChain->SetBranchAddress("mu_TrkQuality", &mu_TrkQuality, &b_mu_TrkQuality);
   fChain->SetBranchAddress("mu_IsoTrk", &mu_IsoTrk, &b_mu_IsoTrk);
   fChain->SetBranchAddress("mu_PFChIso", &mu_PFChIso, &b_mu_PFChIso);
   fChain->SetBranchAddress("mu_PFPhoIso", &mu_PFPhoIso, &b_mu_PFPhoIso);
   fChain->SetBranchAddress("mu_PFNeuIso", &mu_PFNeuIso, &b_mu_PFNeuIso);
   fChain->SetBranchAddress("mu_PFPUIso", &mu_PFPUIso, &b_mu_PFPUIso);
   fChain->SetBranchAddress("mu_PFCHIso03", &mu_PFCHIso03, &b_mu_PFCHIso03);
   fChain->SetBranchAddress("mu_PFPhoIso03", &mu_PFPhoIso03, &b_mu_PFPhoIso03);
   fChain->SetBranchAddress("mu_PFNeuIso03", &mu_PFNeuIso03, &b_mu_PFNeuIso03);
   fChain->SetBranchAddress("mu_PFPUIso03", &mu_PFPUIso03, &b_mu_PFPUIso03);
   fChain->SetBranchAddress("mu_InnervalidFraction", &mu_InnervalidFraction, &b_mu_InnervalidFraction);
   fChain->SetBranchAddress("mu_segmentCompatibility", &mu_segmentCompatibility, &b_mu_segmentCompatibility);
   fChain->SetBranchAddress("mu_chi2LocalPosition", &mu_chi2LocalPosition, &b_mu_chi2LocalPosition);
   fChain->SetBranchAddress("mu_trkKink", &mu_trkKink, &b_mu_trkKink);
   fChain->SetBranchAddress("mu_BestTrkPtError", &mu_BestTrkPtError, &b_mu_BestTrkPtError);
   fChain->SetBranchAddress("mu_BestTrkPt", &mu_BestTrkPt, &b_mu_BestTrkPt);
   fChain->SetBranchAddress("mu_BestTrkType", &mu_BestTrkType, &b_mu_BestTrkType);
   fChain->SetBranchAddress("genMuon_pt", &genMuon_pt, &b_genMuon_pt);
   fChain->SetBranchAddress("genMuon_eta", &genMuon_eta, &b_genMuon_eta);
   fChain->SetBranchAddress("genMuon_phi", &genMuon_phi, &b_genMuon_phi);
   fChain->SetBranchAddress("genMuon_energy", &genMuon_energy, &b_genMuon_energy);
   Notify();
}

Bool_t TnPAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TnPAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TnPAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TnPAnalysis_cxx
