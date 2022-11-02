#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include "TProfile.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"
#include "TProfile.h"

void PrintHisto(TH1F &histo){
  TCanvas c ("canv","canv", 900, 600);
  histo.Draw();
  float RMS=histo.GetRMS();

  string oname=string(histo.GetName())+".png";
  if (oname.find("angle")!=string::npos)
    histo.Fit("gaus","","",(-1*RMS)/3  , (1*RMS)/3);
  else histo.Fit("gaus","","",(-1*RMS)/2, (1*RMS)/2);
  c.SaveAs(oname.c_str());
}

void PrintHisto2(TH2F &histo){
  TCanvas c("canv","canv", 900, 600);
  histo.Draw("COLZ");
  string oname=string(histo.GetName())+".png";
  c.SaveAs(oname.c_str());
}


void CutAndPlot_Bs(string filename, bool mc=0, int trig=1, bool sb=0, bool draw=0){
  //trig=0 -> no trigger selection, trig=1 -> JpsiMu trig, trig=2 -> JpsiTkTk trig
  //for sb=1 sidebands are selected
  //if draw=1 some resolution plots are produced and printed (works only if mc==1)
  if (!mc) draw=false;

  TTree          *fChain, *copy;   //!pointer to the analyzed TTree or TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.
  TH1F *hBsPt, *hBsEta, *hBsMass, *hBsPhi, *hBsMassNoRefit,*hBsMassFromSV;
  TH1F *hJpsiPt, *hJpsiEta, *hJpsiMass, *hJpsiPhi,*hJpsiMassNoRefit,*hJpsiMassFromSV;
  TH1F *hPhiPt, *hPhiEta, *hPhiMass, *hPhiPhi,*hPhiMassNoRefit,*hPhiMassFromSV;

  TH1F *hMupPt, *hMupEta, *hMupPhi;
  TH1F *hMumPt, *hMumEta, *hMumPhi;

  TH1F *hKpPt, *hKpEta, *hKpPhi;
  TH1F *hKmPt, *hKmEta, *hKmPhi;

  TH1F *hPix, *hStrip;

  TH1F *hangle_cospsi, *hangle_costheta, *hangle_phi, *hlxy, *hlxye ;

  //MC histograms
  TH1F *hBsPt_GEN, *hBsEta_GEN, *hBsMass_GEN, *hBsPhi_GEN;
  TH1F *hJpsiPt_GEN, *hJpsiEta_GEN, *hJpsiMass_GEN, *hJpsiPhi_GEN;
  TH1F *hPhiPt_GEN, *hPhiEta_GEN, *hPhiMass_GEN, *hPhiPhi_GEN;
  
  TH1F *hMupPt_GEN, *hMupEta_GEN, *hMupPhi_GEN;
  TH1F *hMumPt_GEN, *hMumEta_GEN, *hMumPhi_GEN;
  
  TH1F *hKpPt_GEN, *hKpEta_GEN, *hKpPhi_GEN;
  TH1F *hKmPt_GEN, *hKmEta_GEN, *hKmPhi_GEN;
  
  TH1F *hangle_cospsi_GEN, *hangle_costheta_GEN, *hangle_phi_GEN, *hlxy_GEN;

  //Residuals
  
  TH1F *hBsPt_RES, *hBsEta_RES, *hBsMass_RES, *hBsPhi_RES;
  TH1F *hJpsiPt_RES, *hJpsiEta_RES, *hJpsiMass_RES, *hJpsiPhi_RES;
  TH1F *hPhiPt_RES, *hPhiEta_RES, *hPhiMass_RES, *hPhiPhi_RES;
  
  TH1F *hMupPt_RES, *hMupEta_RES, *hMupPhi_RES;
  TH1F *hMumPt_RES, *hMumEta_RES, *hMumPhi_RES;

  TH1F *hKpPt_RES, *hKpEta_RES, *hKpPhi_RES;
  TH1F *hKmPt_RES, *hKmEta_RES, *hKmPhi_RES;

  TH1F *hangle_cospsi_RES, *hangle_costheta_RES, *hangle_phi_RES, *hlxy_RES;
  TH1F *hBsMassFromSV_RES, *hJpsiMassFromSV_RES, *hPhiMassFromSV_RES; 

  TProfile *hBsMassFromSV_RESvsPT, *hJpsiMassFromSV_RESvsPT, *hPhiMassFromSV_RESvsPT; 
  TProfile *hBsPt_RESvsPT, *hBsEta_RESvsPT, *hBsMass_RESvsPT, *hBsPhi_RESvsPT;
  TProfile *hJpsiPt_RESvsPT, *hJpsiEta_RESvsPT, *hJpsiMass_RESvsPT, *hJpsiPhi_RESvsPT;
  TProfile *hPhiPt_RESvsPT, *hPhiEta_RESvsPT, *hPhiMass_RESvsPT, *hPhiPhi_RESvsPT;
  
  TProfile *hMupPt_RESvsPT, *hMupEta_RESvsPT, *hMupPhi_RESvsPT;
  TProfile *hMumPt_RESvsPT, *hMumEta_RESvsPT, *hMumPhi_RESvsPT;

  TProfile *hKpPt_RESvsPT, *hKpEta_RESvsPT, *hKpPhi_RESvsPT;
  TProfile *hKmPt_RESvsPT, *hKmEta_RESvsPT, *hKmPhi_RESvsPT;

  TProfile *hangle_cospsi_RESvsPT, *hangle_costheta_RESvsPT, *hangle_phi_RESvsPT, *hlxy_RESvsPT;

  TH2F *hBsMassFromSV_RESvsPTTH2, *hJpsiMassFromSV_RESvsPTTH2, *hPhiMassFromSV_RESvsPTTH2; 
  TH2F *hBsPt_RESvsPTTH2, *hBsEta_RESvsPTTH2, *hBsMass_RESvsPTTH2, *hBsPhi_RESvsPTTH2;
  TH2F *hJpsiPt_RESvsPTTH2, *hJpsiEta_RESvsPTTH2, *hJpsiMass_RESvsPTTH2, *hJpsiPhi_RESvsPTTH2;
  TH2F *hPhiPt_RESvsPTTH2, *hPhiEta_RESvsPTTH2, *hPhiMass_RESvsPTTH2, *hPhiPhi_RESvsPTTH2;
  
  TH2F *hMupPt_RESvsPTTH2, *hMupEta_RESvsPTTH2, *hMupPhi_RESvsPTTH2;
  TH2F *hMumPt_RESvsPTTH2, *hMumEta_RESvsPTTH2, *hMumPhi_RESvsPTTH2;

  TH2F *hKpPt_RESvsPTTH2, *hKpEta_RESvsPTTH2, *hKpPhi_RESvsPTTH2;
  TH2F *hKmPt_RESvsPTTH2, *hKmEta_RESvsPTTH2, *hKmPhi_RESvsPTTH2;
   
  TH2F *hangle_cospsi_RESvsPTTH2, *hangle_costheta_RESvsPTTH2, *hangle_phi_RESvsPTTH2, *hlxy_RESvsPTTH2;

  TH1F *hSVrecoMassRES;
  TH2F *hSVrecoMassRESvsPTTH2;
  TProfile *hSVrecoMassRESvsPT;
  
  string nome;

  // Declaration of leaf types
  Float_t run;
  Float_t evt;
  Float_t lumi;

  Float_t HLT_JpsiTkTk;
  Float_t HLT_JpsiTk;
  Float_t HLT_JpsiMu;

  Float_t  Bs_MassNoRefit;
  Float_t  Bs_MassFromSV;

  Float_t  Jpsi_MassNoRefit;
  Float_t  Jpsi_MassFromSV;

  Float_t  Phi_MassNoRefit;
  Float_t  Phi_MassFromSV;
  
  Float_t  HLT_MatchedJpsi;
  Float_t  HLT_MatchedPhi;

  Float_t         angle_cospsi;
  Float_t         angle_costheta;
  Float_t         angle_phi;
  Float_t         PV_cosPoint;
  Float_t         Lxy;
  Float_t         Lxyz;
  Float_t         LxyErr;
  Float_t         LxyzErr;
  Float_t         Bs_VProb;
  Float_t         Jpsi_VProb;
  Float_t         Phi_VProb;
  Float_t         Bs_Mass;
  Float_t         Bs_Pt;
  Float_t         Bs_Eta;
  Float_t         Bs_Phi;
  Float_t         Jpsi_Mass;
  Float_t         Jpsi_Pt;
  Float_t         Jpsi_Eta;
  Float_t         Jpsi_Phi;
  Float_t         Phi_Mass;
  Float_t         Phi_Pt;
  Float_t         Phi_Eta;
  Float_t         Phi_Phi;
  Float_t         Mup_Pt;
  Float_t         Mup_Eta;
  Float_t         Mup_Phi;
  Float_t         Mup_Hits;
  Float_t         Mum_Pt;
  Float_t         Mum_Eta;
  Float_t         Mum_Phi;
  Float_t         Mum_Hits;
  Float_t         Km_Pt;
  Float_t         Km_Eta;
  Float_t         Km_Phi;
  Float_t         Km_Hits;
  Float_t         Kp_Pt;
  Float_t         Kp_Eta;
  Float_t         Kp_Phi;
  Float_t         Kp_Hits;
  Float_t         MC_matched;
  Float_t         MC_Flavour;
  Float_t         angle_cospsi_GEN;
  Float_t         angle_costheta_GEN;
  Float_t         angle_phi_GEN;
  Float_t         Lxy_GEN;
  Float_t         Lxyz_GEN;
  Float_t         Bs_Mass_GEN;
  Float_t         Bs_Pt_GEN;
  Float_t         Bs_Eta_GEN;
  Float_t         Bs_Phi_GEN;
  Float_t         Jpsi_Mass_GEN;
  Float_t         Jpsi_Pt_GEN;
  Float_t         Jpsi_Eta_GEN;
  Float_t         Jpsi_Phi_GEN;
  Float_t         Phi_Mass_GEN;
  Float_t         Phi_Pt_GEN;
  Float_t         Phi_Eta_GEN;
  Float_t         Phi_Phi_GEN;
  Float_t         Mup_Pt_GEN;
  Float_t         Mup_Eta_GEN;
  Float_t         Mup_Phi_GEN;
  Float_t         Mum_Pt_GEN;
  Float_t         Mum_Eta_GEN;
  Float_t         Mum_Phi_GEN;
  Float_t         Km_Pt_GEN;
  Float_t         Km_Eta_GEN;
  Float_t         Km_Phi_GEN;
  Float_t         Kp_Pt_GEN;
  Float_t         Kp_Eta_GEN;
  Float_t         Kp_Phi_GEN;

  // end variable declarations 
  
  TFile* f = new TFile(filename.c_str());

  fChain = new TTree;

  fChain=(TTree*)f->Get("OutTree");
    
  if (!fChain){
    cout << "No TTree found in imput file, returning" << endl;
    return;
  }
  
  fChain->SetBranchAddress("run", &run);
  fChain->SetBranchAddress("evt", &lumi);
  fChain->SetBranchAddress("lumi", &lumi);

  fChain->SetBranchAddress("HLT_JpsiTkTk", &HLT_JpsiTkTk);
  fChain->SetBranchAddress("HLT_JpsiTk", &HLT_JpsiTk);
  fChain->SetBranchAddress("HLT_JpsiMu", &HLT_JpsiMu);

  fChain->SetBranchAddress("Bs_MassNoRefit", &Bs_MassNoRefit);
  fChain->SetBranchAddress("Bs_MassFromSV", &Bs_MassFromSV);

  fChain->SetBranchAddress("Jpsi_MassNoRefit", &Jpsi_MassNoRefit);
  fChain->SetBranchAddress("Jpsi_MassFromSV", &Jpsi_MassFromSV);

  fChain->SetBranchAddress("Phi_MassNoRefit", &Phi_MassNoRefit);
  fChain->SetBranchAddress("Phi_MassFromSV", &Phi_MassFromSV);

  fChain->SetBranchAddress("HLT_MatchedJpsi", &HLT_MatchedJpsi);
  fChain->SetBranchAddress("HLT_MatchedPhi", &HLT_MatchedPhi);

  fChain->SetBranchAddress("angle_cospsi", &angle_cospsi);
  fChain->SetBranchAddress("angle_costheta", &angle_costheta);
  fChain->SetBranchAddress("angle_phi", &angle_phi);
  fChain->SetBranchAddress("PV_cosPoint", &PV_cosPoint);
  fChain->SetBranchAddress("Lxy", &Lxy);
  fChain->SetBranchAddress("Lxyz", &Lxyz);
  fChain->SetBranchAddress("LxyErr", &LxyErr);
  fChain->SetBranchAddress("LxyzErr", &LxyzErr);
  fChain->SetBranchAddress("Bs_VProb", &Bs_VProb);
  fChain->SetBranchAddress("Jpsi_VProb", &Jpsi_VProb);
  fChain->SetBranchAddress("Phi_VProb", &Phi_VProb);
  fChain->SetBranchAddress("Bs_Mass", &Bs_Mass);
  fChain->SetBranchAddress("Bs_Pt", &Bs_Pt);
  fChain->SetBranchAddress("Bs_Eta", &Bs_Eta);
  fChain->SetBranchAddress("Bs_Phi", &Bs_Phi);
  fChain->SetBranchAddress("Jpsi_Mass", &Jpsi_Mass);
  fChain->SetBranchAddress("Jpsi_Pt", &Jpsi_Pt);
  fChain->SetBranchAddress("Jpsi_Eta", &Jpsi_Eta);
  fChain->SetBranchAddress("Jpsi_Phi", &Jpsi_Phi);
  fChain->SetBranchAddress("Phi_Mass", &Phi_Mass);
  fChain->SetBranchAddress("Phi_Pt", &Phi_Pt);
  fChain->SetBranchAddress("Phi_Eta", &Phi_Eta);
  fChain->SetBranchAddress("Phi_Phi", &Phi_Phi);
  fChain->SetBranchAddress("Mup_Pt", &Mup_Pt);
  fChain->SetBranchAddress("Mup_Eta", &Mup_Eta);
  fChain->SetBranchAddress("Mup_Phi", &Mup_Phi);
  fChain->SetBranchAddress("Mup_Hits", &Mup_Hits);
  fChain->SetBranchAddress("Mum_Pt", &Mum_Pt);
  fChain->SetBranchAddress("Mum_Eta", &Mum_Eta);
  fChain->SetBranchAddress("Mum_Phi", &Mum_Phi);
  fChain->SetBranchAddress("Mum_Hits", &Mum_Hits);
  fChain->SetBranchAddress("Km_Pt", &Km_Pt);
  fChain->SetBranchAddress("Km_Eta", &Km_Eta);
  fChain->SetBranchAddress("Km_Phi", &Km_Phi);
  fChain->SetBranchAddress("Km_Hits", &Km_Hits);
  fChain->SetBranchAddress("Kp_Pt", &Kp_Pt);
  fChain->SetBranchAddress("Kp_Eta", &Kp_Eta);
  fChain->SetBranchAddress("Kp_Phi", &Kp_Phi);
  fChain->SetBranchAddress("Kp_Hits", &Kp_Hits);
  fChain->SetBranchAddress("MC_matched", &MC_matched);
  fChain->SetBranchAddress("MC_Flavour", &MC_Flavour);
  fChain->SetBranchAddress("angle_cospsi_GEN", &angle_cospsi_GEN);
  fChain->SetBranchAddress("angle_costheta_GEN", &angle_costheta_GEN);
  fChain->SetBranchAddress("angle_phi_GEN", &angle_phi_GEN);
  fChain->SetBranchAddress("Lxy_GEN", &Lxy_GEN);
  fChain->SetBranchAddress("Lxyz_GEN", &Lxyz_GEN);
  fChain->SetBranchAddress("Bs_Mass_GEN", &Bs_Mass_GEN);
  fChain->SetBranchAddress("Bs_Pt_GEN", &Bs_Pt_GEN);
  fChain->SetBranchAddress("Bs_Eta_GEN", &Bs_Eta_GEN);
  fChain->SetBranchAddress("Bs_Phi_GEN", &Bs_Phi_GEN);
  fChain->SetBranchAddress("Jpsi_Mass_GEN", &Jpsi_Mass_GEN);
  fChain->SetBranchAddress("Jpsi_Pt_GEN", &Jpsi_Pt_GEN);
  fChain->SetBranchAddress("Jpsi_Eta_GEN", &Jpsi_Eta_GEN);
  fChain->SetBranchAddress("Jpsi_Phi_GEN", &Jpsi_Phi_GEN);
  fChain->SetBranchAddress("Phi_Mass_GEN", &Phi_Mass_GEN);
  fChain->SetBranchAddress("Phi_Pt_GEN", &Phi_Pt_GEN);
  fChain->SetBranchAddress("Phi_Eta_GEN", &Phi_Eta_GEN);
  fChain->SetBranchAddress("Phi_Phi_GEN", &Phi_Phi_GEN);
  fChain->SetBranchAddress("Mup_Pt_GEN", &Mup_Pt_GEN);
  fChain->SetBranchAddress("Mup_Eta_GEN", &Mup_Eta_GEN);
  fChain->SetBranchAddress("Mup_Phi_GEN", &Mup_Phi_GEN);
  fChain->SetBranchAddress("Mum_Pt_GEN", &Mum_Pt_GEN);
  fChain->SetBranchAddress("Mum_Eta_GEN", &Mum_Eta_GEN);
  fChain->SetBranchAddress("Mum_Phi_GEN", &Mum_Phi_GEN);
  fChain->SetBranchAddress("Km_Pt_GEN", &Km_Pt_GEN);
  fChain->SetBranchAddress("Km_Eta_GEN", &Km_Eta_GEN);
  fChain->SetBranchAddress("Km_Phi_GEN", &Km_Phi_GEN);
  fChain->SetBranchAddress("Kp_Pt_GEN", &Kp_Pt_GEN);
  fChain->SetBranchAddress("Kp_Eta_GEN", &Kp_Eta_GEN);
  fChain->SetBranchAddress("Kp_Phi_GEN", &Kp_Phi_GEN);

  gStyle->SetPalette(1);

  string name;
  if (trig==1)name = "Plots_JpsiMu_";
  else if (trig==2)name = "Plots_JpsiTkTk_";
  else	name="Plots_";
  if (sb) name=name+"_SB_";
  if (mc) name=name+"_MC_";
  name = name + filename;
	
  TFile *out=new TFile(name.c_str(),"RECREATE");

  copy= new TTree;
  copy=fChain->CloneTree(0);
  copy->SetName("OutTree_PlusCuts");

  hlxy=new TH1F("c #tau","c #tau", 100, -.05, 0.5);
  hlxye=new TH1F("c #tau Err","c #tau Err", 100, 0, 0.01);

  hPix=new TH1F("Pix Hits","Pix hits", 10, 0.5, 10.5);
  hStrip=new TH1F("Strip Hits","Strip hits", 30, 0.5, 30.5);

  hBsPt=new TH1F("BsPt","BsPt", 100, 0, 50);
  hBsEta=new TH1F("BsEta","BsEta", 100, -3, 3);
  hBsMass=new TH1F("BsMass","BsMass", 50, 5.2, 5.6);
  hBsPhi=new TH1F("BsPhi","BsPhi", 50, -3.2, 3.2);

  hBsMassFromSV=new TH1F("BsMassFromSV","BsMassFromSV", 50, 5.2, 5.6);
  hBsMassNoRefit=new TH1F("BsMassNoRefit","BsMassNoRefit", 50, 5.2, 5.6);

  hJpsiPt=new TH1F("JpsiPt","JpsiPt", 100, 0, 50);
  hJpsiEta=new TH1F("JpsiEta","JpsiEta", 100, -3, 3);
  hJpsiMass=new TH1F("JpsiMass","JpsiMass", 50, 2.9, 3.2);
  hJpsiPhi=new TH1F("JpsiPhi","JpsiPhi", 50, -3.2, 3.2);

  hJpsiMassFromSV=new TH1F("JpsiMassFromSV","JpsiMassFromSV", 50, 2.9, 3.2);
  hJpsiMassNoRefit=new TH1F("JpsiMassNoRefit","JpsiMassNoRefit", 50, 2.9, 3.2);

  hPhiPt=new TH1F("PhiPt","PhiPt", 100, 0, 50);
  hPhiEta=new TH1F("PhiEta","PhiEta", 100, -3, 3);
  hPhiMass=new TH1F("PhiMass","PhiMass", 50, 0.98, 1.05);
  hPhiPhi=new TH1F("PhiPhi","PhiPhi", 50, -3.2, 3.2);

  hPhiMassFromSV=new TH1F("PhiMassFromSV","PhiMassFromSV", 50, 0.98, 1.05);
  hPhiMassNoRefit=new TH1F("PhiMassNoRefit","PhiMassNoRefit", 50, 0.98, 1.05);

  hMupPt=new TH1F("MupPt","MupPt", 100, 0, 50);
  hMupEta=new TH1F("MupEta","MupEta", 100, -3, 3);
  hMupPhi=new TH1F("MupPhi","MupPhi", 50, -3.2, 3.2);

  hMumPt=new TH1F("MumPt","MumPt", 100, 0, 50);
  hMumEta=new TH1F("MumEta","MumEta", 100, -3, 3);
  hMumPhi=new TH1F("MumPhi","MumPhi", 50, -3.2, 3.2);

  hKpPt=new TH1F("KpPt","KpPt", 100, 0, 50);
  hKpEta=new TH1F("KpEta","KpEta", 100, -3, 3);
  hKpPhi=new TH1F("KpPhi","KpPhi", 50, -3.2, 3.2);

  hKmPt=new TH1F("KmPt","KmPt", 100, 0, 50);
  hKmEta=new TH1F("KmEta","KmEta", 100, -3, 3);
  hKmPhi=new TH1F("KmPhi","KmPhi", 50, -3.2, 3.2);

  hangle_cospsi=new TH1F("angle_cospsi","angle_cospsi", 50, -1, 1);
  hangle_costheta=new TH1F("angle_costheta","angle_costheta", 50, -1, 1);
  hangle_phi=new TH1F("angle_phi","angle_phi", 50, -3.2, 3.2);

  hSVrecoMassRES= new TH1F("SVvsRECO_Mass_RES","SVvsRECO_Mass_RES",50, -0.15, 0.15);
  hSVrecoMassRESvsPTTH2= new TH2F("SVvsRECO_Mass_RESvsPTTH2","SVvsRECO_Mass_RESvsPTTH2",100,0,50,50, -0.15, 0.15);
  hSVrecoMassRESvsPT= new TProfile("SVvsRECO_Mass_RESvsPT","SVvsRECO_Mass_RESvsPT",100,0,50,-0.15, 0.15,"");

  if (mc){

    hBsPt_GEN=new TH1F("BsPt_GEN","BsPt_GEN", 100, 0, 50);
    hBsEta_GEN=new TH1F("BsEta_GEN","BsEta_GEN", 100, -3, 3);
    hBsMass_GEN=new TH1F("BsMass_GEN","BsMass_GEN", 50, 5.2, 5.6);
    hBsPhi_GEN=new TH1F("BsPhi_GEN","BsPhi_GEN", 50, -3.2, 3.2);

    hJpsiPt_GEN=new TH1F("JpsiPt_GEN","JpsiPt_GEN", 100, 0, 50);
    hJpsiEta_GEN=new TH1F("JpsiEta_GEN","JpsiEta_GEN", 100, -3, 3);
    hJpsiMass_GEN=new TH1F("JpsiMass_GEN","JpsiMass_GEN", 50, 2.9, 3.2);
    hJpsiPhi_GEN=new TH1F("JpsiPhi_GEN","JpsiPhi_GEN", 50, -3.2, 3.2);

    hPhiPt_GEN=new TH1F("PhiPt_GEN","PhiPt_GEN", 100, 0, 50);
    hPhiEta_GEN=new TH1F("PhiEta_GEN","PhiEta_GEN", 100, -3, 3);
    hPhiMass_GEN=new TH1F("PhiMass_GEN","PhiMass_GEN", 50, 0.98, 1.05);
    hPhiPhi_GEN=new TH1F("PhiPhi_GEN","PhiPhi_GEN", 50, -3.2, 3.2);

    hMupPt_GEN=new TH1F("MupPt_GEN","MupPt_GEN", 100, 0, 50);
    hMupEta_GEN=new TH1F("MupEta_GEN","MupEta_GEN", 100, -3, 3);
    hMupPhi_GEN=new TH1F("MupPhi_GEN","MupPhi_GEN", 50, -3.2, 3.2);

    hMumPt_GEN=new TH1F("MumPt_GEN","MumPt_GEN", 100, 0, 50);
    hMumEta_GEN=new TH1F("MumEta_GEN","MumEta_GEN", 100, -3, 3);
    hMumPhi_GEN=new TH1F("MumPhi_GEN","MumPhi_GEN", 50, -3.2, 3.2);

    hKpPt_GEN=new TH1F("KpPt_GEN","KpPt_GEN", 100, 0, 50);
    hKpEta_GEN=new TH1F("KpEta_GEN","KpEta_GEN", 100, -3, 3);
    hKpPhi_GEN=new TH1F("KpPhi_GEN","KpPhi_GEN", 50, -3.2, 3.2);

    hKmPt_GEN=new TH1F("KmPt_GEN","KmPt_GEN", 100, 0, 50);
    hKmEta_GEN=new TH1F("KmEta_GEN","KmEta_GEN", 100, -3, 3);
    hKmPhi_GEN=new TH1F("KmPhi_GEN","KmPhi_GEN", 50, -3.2, 3.2);

    hangle_cospsi_GEN=new TH1F("angle_cospsi_GEN","angle_cospsi_GEN", 50, -1, 1);
    hangle_costheta_GEN=new TH1F("angle_costheta_GEN","angle_costheta_GEN", 50, -1, 1);
    hangle_phi_GEN=new TH1F("angle_phi_GEN","angle_phi_GEN", 50, -3.2, 3.2);

    //residuals 
    hBsPt_RES=new TH1F("BsPt_RES","BsPt_RES", 100, -0.5, 0.5);
    hBsEta_RES=new TH1F("BsEta_RES","BsEta_RES", 100, -0.01, 0.01);
    hBsMass_RES=new TH1F("BsMass_RES","BsMass_RES", 50, -0.15, 0.15);
    hBsPhi_RES=new TH1F("BsPhi_RES","BsPhi_RES", 50, -0.01, 0.01);

    hBsMassFromSV_RES=new TH1F("BsMassFromSV_RES","BsMassFromSV_RES", 100, -0.15, 0.15);
    hJpsiMassFromSV_RES=new TH1F("JpsiMassFromSV_RES","JpsiMassFromSV_RES", 100, -0.15, 0.15);
    hPhiMassFromSV_RES=new TH1F("PhiMassFromSV_RES","PhiMassFromSV_RES", 100, -0.025, 0.025);

    hJpsiPt_RES=new TH1F("JpsiPt_RES","JpsiPt_RES", 100, -0.5, 0.5);
    hJpsiEta_RES=new TH1F("JpsiEta_RES","JpsiEta_RES", 100, -0.15, 0.15);
    hJpsiMass_RES=new TH1F("JpsiMass_RES","JpsiMass_RES", 50, -0.15, 0.15);
    hJpsiPhi_RES=new TH1F("JpsiPhi_RES","JpsiPhi_RES", 50, -0.15, 0.15);

    hPhiPt_RES=new TH1F("PhiPt_RES","PhiPt_RES", 100, -0.5, 0.5);
    hPhiEta_RES=new TH1F("PhiEta_RES","PhiEta_RES", 100, -0.15, 0.15);
    hPhiMass_RES=new TH1F("PhiMass_RES","PhiMass_RES", 50, -0.15, 0.15);
    hPhiPhi_RES=new TH1F("PhiPhi_RES","PhiPhi_RES", 50, -0.15, 0.15);

    hMupPt_RES=new TH1F("MupPt_RES","MupPt_RES", 100, -0.5, 0.5);
    hMupEta_RES=new TH1F("MupEta_RES","MupEta_RES", 100, -0.15, 0.15);
    hMupPhi_RES=new TH1F("MupPhi_RES","MupPhi_RES", 50, -0.15, 0.15);

    hMumPt_RES=new TH1F("MumPt_RES","MumPt_RES", 100, -0.5, 0.5);
    hMumEta_RES=new TH1F("MumEta_RES","MumEta_RES", 100, -0.15, 0.15);
    hMumPhi_RES=new TH1F("MumPhi_RES","MumPhi_RES", 50, -0.15, 0.15);

    hKpPt_RES=new TH1F("KpPt_RES","KpPt_RES", 100, -0.5, 0.5);
    hKpEta_RES=new TH1F("KpEta_RES","KpEta_RES", 100, -0.15, 0.15);
    hKpPhi_RES=new TH1F("KpPhi_RES","KpPhi_RES", 50, -0.15, 0.15);

    hKmPt_RES=new TH1F("KmPt_RES","KmPt_RES", 100, -0.5, 0.5);
    hKmEta_RES=new TH1F("KmEta_RES","KmEta_RES", 100, -0.15, 0.15);
    hKmPhi_RES=new TH1F("KmPhi_RES","KmPhi_RES", 50, -0.15, 0.15);

    hangle_cospsi_RES=new TH1F("angle_cospsi_RES","angle_cospsi_RES", 300, -0.15, 0.15);
    hangle_costheta_RES=new TH1F("angle_costheta_RES","angle_costheta_RES", 300, -0.15, 0.15);
    hangle_phi_RES=new TH1F("angle_phi_RES","angle_phi_RES", 300, -0.15, 0.15);

    hlxy_GEN=new TH1F("Lxy_GEN","Lxy_GEN", 100, -.05, 0.5);

    hlxy_RES=new TH1F("ctau_RES","ctau_RES", 200, -.02, 0.02);

    hBsMassFromSV_RESvsPT=new TProfile("BsMassFromSV_RESvsPT","BsMassFromSV_RESvsPT", 100,0,50, -0.5, 0.5,"");
    hJpsiMassFromSV_RESvsPT=new TProfile("JpsiMassFromSV_RESvsPT","JpsiMassFromSV_RESvsPT", 100,0,50, -0.5, 0.5,"");
    hPhiMassFromSV_RESvsPT=new TProfile("PhiMassFromSV_RESvsPT","PhiMassFromSV_RESvsPT", 100,0,50, -0.025, 0.025,"");

    hBsPt_RESvsPT=new TProfile("BsPt_RESvsPT","BsPt_RESvsPT",100,0,50, -0.5, 0.5,"");
    hBsEta_RESvsPT=new TProfile("BsEta_RESvsPT","BsEta_RESvsPT", 100,0,50, -0.15, 0.15,"");
    hBsMass_RESvsPT=new TProfile("BsMass_RESvsPT","BsMass_RESvsPT",100,0,50, -0.15, 0.15,"");
    hBsPhi_RESvsPT=new TProfile("BsPhi_RESvsPT","BsPhi_RESvsPT",100,0,50, -0.15, 0.15,"");

    hJpsiPt_RESvsPT=new TProfile("JpsiPt_RESvsPT","JpsiPt_RESvsPT", 100,0,50, -0.5, 0.5,"");
    hJpsiEta_RESvsPT=new TProfile("JpsiEta_RESvsPT","JpsiEta_RESvsPT", 100,0,50, -0.15, 0.15,"");
    hJpsiMass_RESvsPT=new TProfile("JpsiMass_RESvsPT","JpsiMass_RESvsPT",100,0,50, -0.15, 0.15,"");
    hJpsiPhi_RESvsPT=new TProfile("JpsiPhi_RESvsPT","JpsiPhi_RESvsPT", 100,0,50, -0.15, 0.15,"");

    hPhiPt_RESvsPT=new TProfile("PhiPt_RESvsPT","PhiPt_RESvsPT", 100,0,50, -0.5, 0.5,"");
    hPhiEta_RESvsPT=new TProfile("PhiEta_RESvsPT","PhiEta_RESvsPT", 100,0,50, -0.15, 0.15,"");
    hPhiMass_RESvsPT=new TProfile("PhiMass_RESvsPT","PhiMass_RESvsPT", 100,0,50, -0.15, 0.15,"");
    hPhiPhi_RESvsPT=new TProfile("PhiPhi_RESvsPT","PhiPhi_RESvsPT", 100,0,50, -0.15, 0.15,"");

    hMupPt_RESvsPT=new TProfile("MupPt_RESvsPT","MupPt_RESvsPT", 100,0,50, -0.5, 0.5,"");
    hMupEta_RESvsPT=new TProfile("MupEta_RESvsPT","MupEta_RESvsPT", 100,0,50, -0.15, 0.15,"");
    hMupPhi_RESvsPT=new TProfile("MupPhi_RESvsPT","MupPhi_RESvsPT", 100,0,50, -0.15, 0.15,"");

    hMumPt_RESvsPT=new TProfile("MumPt_RESvsPT","MumPt_RESvsPT", 100,0,50, -0.5, 0.5,"");
    hMumEta_RESvsPT=new TProfile("MumEta_RESvsPT","MumEta_RESvsPT", 100,0,50, -0.15, 0.15,"");
    hMumPhi_RESvsPT=new TProfile("MumPhi_RESvsPT","MumPhi_RESvsPT", 100,0,50, -0.15, 0.15,"");

    hKpPt_RESvsPT=new TProfile("KpPt_RESvsPT","KpPt_RESvsPT", 100,0,50, -0.5, 0.5,"");
    hKpEta_RESvsPT=new TProfile("KpEta_RESvsPT","KpEta_RESvsPT", 100,0,50, -0.15, 0.15,"");
    hKpPhi_RESvsPT=new TProfile("KpPhi_RESvsPT","KpPhi_RESvsPT", 100,0,50, -0.15, 0.15,"");

    hKmPt_RESvsPT=new TProfile("KmPt_RESvsPT","KmPt_RESvsPT", 100,0,50, -0.5, 0.5,"");
    hKmEta_RESvsPT=new TProfile("KmEta_RESvsPT","KmEta_RESvsPT", 100,0,50, -0.15, 0.15,"");
    hKmPhi_RESvsPT=new TProfile("KmPhi_RESvsPT","KmPhi_RESvsPT", 100,0,50, -0.15, 0.15,"");

    hangle_cospsi_RESvsPT=new TProfile("angle_cospsi_RESvsPT","angle_cospsi_RESvsPT", 100,0,50, -0.15, 0.15,"");
    hangle_costheta_RESvsPT=new TProfile("angle_costheta_RESvsPT","angle_costheta_RESvsPT", 100,0,50, -0.15, 0.15,"");
    hangle_phi_RESvsPT=new TProfile("angle_phi_RESvsPT","angle_phi_RESvsPT", 100,0,50, -0.15, 0.15,"");

    hlxy_RESvsPT=new TProfile("Lxy_RESvsPT","Lxy_RESvsPT", 100,0,50, -.05, 0.05,"");

    hBsMassFromSV_RESvsPTTH2=new TH2F("BsMassFromSV_RESvsPTTH2","BsMassFromSV_RESvsPTTH2", 100,0,50, 100, -0.5, 0.5);
    hJpsiMassFromSV_RESvsPTTH2=new TH2F("JpsiMassFromSV_RESvsPTTH2","JpsiMassFromSV_RESvsPTTH2", 100,0,50, 100, -0.5, 0.5);
    hPhiMassFromSV_RESvsPTTH2=new TH2F("PhiMassFromSV_RESvsPTTH2","PhiMassFromSV_RESvsPTTH2", 100,0,50, 100, -0.025, 0.025);

    hBsPt_RESvsPTTH2=new TH2F("BsPt_RESvsPTTH2","BsPt_RESvsPTTH2",100,0,50, 100, -0.5, 0.5);
    hBsEta_RESvsPTTH2=new TH2F("BsEta_RESvsPTTH2","BsEta_RESvsPTTH2", 100,0,50, 100, -0.01, 0.01);
    hBsMass_RESvsPTTH2=new TH2F("BsMass_RESvsPTTH2","BsMass_RESvsPTTH2",100,0,50, 100, -0.15, 0.15);
    hBsPhi_RESvsPTTH2=new TH2F("BsPhi_RESvsPTTH2","BsPhi_RESvsPTTH2",100,0,50, 100, -0.01, 0.01);

    hJpsiPt_RESvsPTTH2=new TH2F("JpsiPt_RESvsPTTH2","JpsiPt_RESvsPTTH2", 100,0,50, 100, -0.5, 0.5);
    hJpsiEta_RESvsPTTH2=new TH2F("JpsiEta_RESvsPTTH2","JpsiEta_RESvsPTTH2", 100,0,50, 100, -0.15, 0.15);
    hJpsiMass_RESvsPTTH2=new TH2F("JpsiMass_RESvsPTTH2","JpsiMass_RESvsPTTH2",100,0,50, 100, -0.15, 0.15);
    hJpsiPhi_RESvsPTTH2=new TH2F("JpsiPhi_RESvsPTTH2","JpsiPhi_RESvsPTTH2", 100,0,50, 100, -0.15, 0.15);

    hPhiPt_RESvsPTTH2=new TH2F("PhiPt_RESvsPTTH2","PhiPt_RESvsPTTH2", 100,0,50, 100, -0.5, 0.5);
    hPhiEta_RESvsPTTH2=new TH2F("PhiEta_RESvsPTTH2","PhiEta_RESvsPTTH2", 100,0,50, 100, -0.15, 0.15);
    hPhiMass_RESvsPTTH2=new TH2F("PhiMass_RESvsPTTH2","PhiMass_RESvsPTTH2", 100,0,50, 100, -0.15, 0.15);
    hPhiPhi_RESvsPTTH2=new TH2F("PhiPhi_RESvsPTTH2","PhiPhi_RESvsPTTH2", 100,0,50, 100, -0.15, 0.15);

    hMupPt_RESvsPTTH2=new TH2F("MupPt_RESvsPTTH2","MupPt_RESvsPTTH2", 100,0,50, 100, -0.5, 0.5);
    hMupEta_RESvsPTTH2=new TH2F("MupEta_RESvsPTTH2","MupEta_RESvsPTTH2", 100,0,50, 100, -0.15, 0.15);
    hMupPhi_RESvsPTTH2=new TH2F("MupPhi_RESvsPTTH2","MupPhi_RESvsPTTH2", 100,0,50, 100, -0.15, 0.15);

    hMumPt_RESvsPTTH2=new TH2F("MumPt_RESvsPTTH2","MumPt_RESvsPTTH2", 100,0,50, 100, -0.5, 0.5);
    hMumEta_RESvsPTTH2=new TH2F("MumEta_RESvsPTTH2","MumEta_RESvsPTTH2", 100,0,50, 100, -0.15, 0.15);
    hMumPhi_RESvsPTTH2=new TH2F("MumPhi_RESvsPTTH2","MumPhi_RESvsPTTH2", 100,0,50, 100, -0.15, 0.15);

    hKpPt_RESvsPTTH2=new TH2F("KpPt_RESvsPTTH2","KpPt_RESvsPTTH2", 100,0,50, 100, -0.5, 0.5);
    hKpEta_RESvsPTTH2=new TH2F("KpEta_RESvsPTTH2","KpEta_RESvsPTTH2", 100,0,50, 100, -0.15, 0.15);
    hKpPhi_RESvsPTTH2=new TH2F("KpPhi_RESvsPTTH2","KpPhi_RESvsPTTH2", 100,0,50, 100, -0.15, 0.15);

    hKmPt_RESvsPTTH2=new TH2F("KmPt_RESvsPTTH2","KmPt_RESvsPTTH2", 100,0,50, 100, -0.5, 0.5);
    hKmEta_RESvsPTTH2=new TH2F("KmEta_RESvsPTTH2","KmEta_RESvsPTTH2", 100,0,50, 100, -0.15, 0.15);
    hKmPhi_RESvsPTTH2=new TH2F("KmPhi_RESvsPTTH2","KmPhi_RESvsPTTH2", 100,0,50, 100, -0.15, 0.15);

    hangle_cospsi_RESvsPTTH2=new TH2F("angle_cospsi_RESvsPTTH2","angle_cospsi_RESvsPTTH2", 100,0,50, 200, -0.15, 0.15);
    hangle_costheta_RESvsPTTH2=new TH2F("angle_costheta_RESvsPTTH2","angle_costheta_RESvsPTTH2", 100,0,50, 200, -0.15, 0.15);
    hangle_phi_RESvsPTTH2=new TH2F("angle_phi_RESvsPTTH2","angle_phi_RESvsPTTH2", 100,0,50, 200, -0.15, 0.15);

    hlxy_RESvsPTTH2=new TH2F("Lxy_RESvsPTTH2","Lxy_RESvsPTTH2", 100,0,50, 100, -.05, 0.05);

  }

   Long64_t nentries = fChain->GetEntries();
   cout << "Start Processing " << nentries << " events" << endl;
 
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      
      fChain->GetEntry(jentry);
      if (jentry%100000==0) cout << "processing event " << jentry << "/" << nentries << endl;

      if (mc && MC_matched<0) continue;

      //actual code here
      TLorentzVector mup, mum, kp, km, Bs, Jpsi, Phi;

      const float MUMASS=0.105658;
      const float KMASS=0.493677;

      mum.SetPtEtaPhiM(Mum_Pt,Mum_Eta,Mum_Phi,MUMASS);
      mup.SetPtEtaPhiM(Mup_Pt,Mup_Eta,Mup_Phi,MUMASS);

      km.SetPtEtaPhiM(Km_Pt,Km_Eta,Km_Phi,KMASS);
      kp.SetPtEtaPhiM(Kp_Pt,Kp_Eta,Kp_Phi,KMASS);

      Jpsi.SetPtEtaPhiM(Jpsi_Pt,Jpsi_Eta,Jpsi_Phi,Jpsi_MassFromSV);      
      Phi.SetPtEtaPhiM(Phi_Pt,Phi_Eta,Phi_Phi,Phi_MassFromSV);      
      Bs.SetPtEtaPhiM(Bs_Pt,Bs_Eta,Bs_Phi,Bs_MassFromSV);
      
      int mpnpst = ( int(Mup_Hits) / 100 ) % 10000; 
      int mpnpix = mpnpst / 100;
      int mpntrk = mpnpst % 100;

      int mmnpst = ( int(Mum_Hits) / 100 ) % 10000; 
      int mmnpix = mmnpst / 100;
      int mmntrk = mmnpst % 100;

      int kpnpst = ( int(Kp_Hits) / 100 ) % 10000; 
      int kpnpix = kpnpst / 100;
      int kpntrk = kpnpst % 100;

      int kmnpst = ( int(Km_Hits) / 100 ) % 10000; 
      int kmnpix = kmnpst / 100;
      int kmntrk = kmnpst % 100;

      hPix->Fill(mpnpix);
      hStrip->Fill(mpntrk-mpnpix);

      //      if (!HLT_MatchedJpsi || !HLT_MatchedPhi) continue;
      
      if (trig==1 && HLT_JpsiMu<0) continue;
      if (trig==2 && HLT_JpsiTkTk<0) continue;      
      //if (HLT_JpsiTk<0) continue; // this is a third trigger available in the ntuples

      //analysis cuts
      if (abs(mum.Eta()) > 2.1) continue;
      if (mum.Pt() < 4.0) continue;

      if (abs(mup.Eta()) > 2.1) continue;
      if (mup.Pt() < 4.0) continue;

      if (kp.Pt() < 0.7 || km.Pt() < 0.7) continue;

      if (kmntrk<4 || kpntrk<4 ) continue;
      if (abs(km.Eta()) > 2.5 || abs(kp.Eta())>2.5) continue;

      if (trig==2 && Lxy< 0.02) continue;//for JpsiTkTk trigger
      if (trig==2 && Lxy/LxyErr < 3) continue;

      if (trig==1 && Lxy< 0.005) continue; //for JpsiMu trigger
      
      if (Bs.Pt() < 8) continue;
      if (Bs_MassFromSV > 5.65 || Bs_MassFromSV < 5.2) continue;

      if (sb && Bs_MassFromSV > 5.28 && Bs_MassFromSV < 5.45) continue; // Sideband selection
      //if (Bs_MassFromSV < 5.28 || Bs_MassFromSV > 5.45) continue; //Under Signal selection

      if (abs(Bs.Eta()) >2.4) continue;

      if (Jpsi.Pt()<7.) continue;
      if (Bs_VProb < 0.02) continue;

      if (Jpsi_MassFromSV < 2.9469 || Jpsi_MassFromSV > 3.2469) continue;
      if (Phi_MassFromSV < 1.009 || Phi_MassFromSV > 1.029) continue;

      copy->Fill();

      hBsPt->Fill(Bs_Pt);
      hBsEta->Fill(Bs_Eta);
      hBsMass->Fill(Bs_Mass);
      hBsPhi->Fill(Bs_Phi);

      hSVrecoMassRES->Fill(Bs_MassFromSV-Bs_Mass);
      hSVrecoMassRESvsPT->Fill(Bs_Pt,Bs_MassFromSV-Bs_Mass);
      hSVrecoMassRESvsPTTH2->Fill(Bs_Pt,Bs_MassFromSV-Bs_Mass);

      hBsMassFromSV->Fill(Bs_MassFromSV);
      hBsMassNoRefit->Fill(Bs_MassNoRefit);

      hJpsiMassFromSV->Fill(Jpsi_MassFromSV);
      hJpsiMassNoRefit->Fill(Jpsi_MassNoRefit);

      hPhiMassFromSV->Fill(Phi_MassFromSV);
      hPhiMassNoRefit->Fill(Phi_MassNoRefit);

      hJpsiPt->Fill(Jpsi_Pt);
      hJpsiEta->Fill(Jpsi_Eta);
      hJpsiMass->Fill(Jpsi_Mass);
      hJpsiPhi->Fill(Jpsi_Phi);

      hPhiPt->Fill(Phi_Pt);
      hPhiEta->Fill(Phi_Eta);
      hPhiMass->Fill(Phi_Mass);
      hPhiPhi->Fill(Phi_Phi);

      hMupPt->Fill(Mup_Pt);
      hMupEta->Fill(Mup_Eta);
      hMupPhi->Fill(Mup_Phi);

      hMumPt->Fill(Mum_Pt);
      hMumEta->Fill(Mum_Eta);
      hMumPhi->Fill(Mum_Phi);

      hKpPt->Fill(Kp_Pt);
      hKpEta->Fill(Kp_Eta);
      hKpPhi->Fill(Kp_Phi);

      hKmPt->Fill(Km_Pt);
      hKmEta->Fill(Km_Eta);
      hKmPhi->Fill(Km_Phi);

      hangle_cospsi->Fill(angle_cospsi);
      hangle_costheta->Fill(angle_costheta); 
      hangle_phi->Fill(angle_phi);

      hlxy->Fill(Lxy);
      hlxye->Fill(LxyErr);


      if (mc){

	hBsPt_GEN->Fill(Bs_Pt_GEN);
	hBsEta_GEN->Fill(Bs_Eta_GEN);
	hBsMass_GEN->Fill(Bs_Mass_GEN);
	hBsPhi_GEN->Fill(Bs_Phi_GEN);

	hJpsiPt_GEN->Fill(Jpsi_Pt_GEN);
	hJpsiEta_GEN->Fill(Jpsi_Eta_GEN);
	hJpsiMass_GEN->Fill(Jpsi_Mass_GEN);
	hJpsiPhi_GEN->Fill(Jpsi_Phi_GEN);

	hPhiPt_GEN->Fill(Phi_Pt_GEN);
	hPhiEta_GEN->Fill(Phi_Eta_GEN);
	hPhiMass_GEN->Fill(Phi_Mass_GEN);
	hPhiPhi_GEN->Fill(Phi_Phi_GEN);

	hMupPt_GEN->Fill(Mup_Pt_GEN);
	hMupEta_GEN->Fill(Mup_Eta_GEN);
	hMupPhi_GEN->Fill(Mup_Phi_GEN);

	hMumPt_GEN->Fill(Mum_Pt_GEN);
	hMumEta_GEN->Fill(Mum_Eta_GEN);
	hMumPhi_GEN->Fill(Mum_Phi_GEN);
	
	hKpPt_GEN->Fill(Kp_Pt_GEN);
	hKpEta_GEN->Fill(Kp_Eta_GEN);
	hKpPhi_GEN->Fill(Kp_Phi_GEN);
	
	hKmPt_GEN->Fill(Km_Pt_GEN);
	hKmEta_GEN->Fill(Km_Eta_GEN);
	hKmPhi_GEN->Fill(Km_Phi_GEN);
	
	hangle_cospsi_GEN->Fill(angle_cospsi_GEN);
	hangle_costheta_GEN->Fill(angle_costheta_GEN); 
	hangle_phi_GEN->Fill(angle_phi_GEN);

	hlxy_GEN->Fill(Lxy_GEN);

	hBsPt_RES->Fill(Bs_Pt-Bs_Pt_GEN);
	hBsEta_RES->Fill(Bs_Eta-Bs_Eta_GEN);
	hBsMass_RES->Fill(Bs_Mass-Bs_Mass_GEN);
	hBsPhi_RES->Fill(Bs_Phi-Bs_Phi_GEN);

	hJpsiPt_RES->Fill(Jpsi_Pt-Jpsi_Pt_GEN);
	hJpsiEta_RES->Fill(Jpsi_Eta-Jpsi_Eta_GEN);
	hJpsiMass_RES->Fill(Jpsi_Mass-Jpsi_Mass_GEN);
	hJpsiPhi_RES->Fill(Jpsi_Phi-Jpsi_Phi_GEN);

	hPhiPt_RES->Fill(Phi_Pt-Phi_Pt_GEN);
	hPhiEta_RES->Fill(Phi_Eta-Phi_Eta_GEN);
	hPhiMass_RES->Fill(Phi_Mass-Phi_Mass_GEN);
	hPhiPhi_RES->Fill(Phi_Phi-Phi_Phi_GEN);

	hMupPt_RES->Fill(Mup_Pt-Mup_Pt_GEN);
	hMupEta_RES->Fill(Mup_Eta-Mup_Eta_GEN);
	hMupPhi_RES->Fill(Mup_Phi-Mup_Phi_GEN);

	hMumPt_RES->Fill(Mum_Pt-Mum_Pt_GEN);
	hMumEta_RES->Fill(Mum_Eta-Mum_Eta_GEN);
	hMumPhi_RES->Fill(Mum_Phi-Mum_Phi_GEN);
	
	hKpPt_RES->Fill(Kp_Pt-Kp_Pt_GEN);
	hKpEta_RES->Fill(Kp_Eta-Kp_Eta_GEN);
	hKpPhi_RES->Fill(Kp_Phi-Kp_Phi_GEN);
	
	hKmPt_RES->Fill(Km_Pt-Km_Pt_GEN);
	hKmEta_RES->Fill(Km_Eta-Km_Eta_GEN);
	hKmPhi_RES->Fill(Km_Phi-Km_Phi_GEN);
	
	hBsMassFromSV_RES->Fill(Bs_MassFromSV-Bs_Mass_GEN);
	hJpsiMassFromSV_RES->Fill(Jpsi_MassFromSV-Jpsi_Mass_GEN);
	hPhiMassFromSV_RES->Fill(Phi_MassFromSV-Phi_Mass_GEN);

	hangle_cospsi_RES->Fill(angle_cospsi-angle_cospsi_GEN);
	hangle_costheta_RES->Fill(angle_costheta-angle_costheta_GEN); 
	hangle_phi_RES->Fill(angle_phi-angle_phi_GEN);

	hlxy_RES->Fill(Lxy-Lxy_GEN);

	hBsPt_RESvsPT->Fill(Bs_Pt,Bs_Pt-Bs_Pt_GEN);
	hBsEta_RESvsPT->Fill(Bs_Pt,Bs_Eta-Bs_Eta_GEN);
	hBsMass_RESvsPT->Fill(Bs_Pt,Bs_Mass-Bs_Mass_GEN);
	hBsPhi_RESvsPT->Fill(Bs_Pt,Bs_Phi-Bs_Phi_GEN);

	hJpsiPt_RESvsPT->Fill(Jpsi_Pt,Jpsi_Pt-Jpsi_Pt_GEN);
	hJpsiEta_RESvsPT->Fill(Jpsi_Pt,Jpsi_Eta-Jpsi_Eta_GEN);
	hJpsiMass_RESvsPT->Fill(Jpsi_Pt,Jpsi_Mass-Jpsi_Mass_GEN);
	hJpsiPhi_RESvsPT->Fill(Jpsi_Pt,Jpsi_Phi-Jpsi_Phi_GEN);

	hPhiPt_RESvsPT->Fill(Phi_Pt,Phi_Pt-Phi_Pt_GEN);
	hPhiEta_RESvsPT->Fill(Phi_Pt,Phi_Eta-Phi_Eta_GEN);
	hPhiMass_RESvsPT->Fill(Phi_Pt,Phi_Mass-Phi_Mass_GEN);
	hPhiPhi_RESvsPT->Fill(Phi_Pt,Phi_Phi-Phi_Phi_GEN);

	hMupPt_RESvsPT->Fill(Mup_Pt,Mup_Pt-Mup_Pt_GEN);
	hMupEta_RESvsPT->Fill(Mup_Pt,Mup_Eta-Mup_Eta_GEN);
	hMupPhi_RESvsPT->Fill(Mup_Pt,Mup_Phi-Mup_Phi_GEN);

	hMumPt_RESvsPT->Fill(Mum_Pt,Mum_Pt-Mum_Pt_GEN);
	hMumEta_RESvsPT->Fill(Mum_Pt,Mum_Eta-Mum_Eta_GEN);
	hMumPhi_RESvsPT->Fill(Mum_Pt,Mum_Phi-Mum_Phi_GEN);
	
	hKpPt_RESvsPT->Fill(Kp_Pt,Kp_Pt-Kp_Pt_GEN);
	hKpEta_RESvsPT->Fill(Kp_Pt,Kp_Eta-Kp_Eta_GEN);
	hKpPhi_RESvsPT->Fill(Kp_Pt,Kp_Phi-Kp_Phi_GEN);
	
	hKmPt_RESvsPT->Fill(Km_Pt,Km_Pt-Km_Pt_GEN);
	hKmEta_RESvsPT->Fill(Km_Pt,Km_Eta-Km_Eta_GEN);
	hKmPhi_RESvsPT->Fill(Km_Pt,Km_Phi-Km_Phi_GEN);
	
	hangle_cospsi_RESvsPT->Fill(Bs_Pt,angle_cospsi-angle_cospsi_GEN);
	hangle_costheta_RESvsPT->Fill(Bs_Pt,angle_costheta-angle_costheta_GEN); 
	hangle_phi_RESvsPT->Fill(Bs_Pt,angle_phi-angle_phi_GEN);

	hlxy_RESvsPT->Fill(Bs_Pt,Lxy-Lxy_GEN);

	hBsPt_RESvsPTTH2->Fill(Bs_Pt,Bs_Pt-Bs_Pt_GEN);
	hBsEta_RESvsPTTH2->Fill(Bs_Pt,Bs_Eta-Bs_Eta_GEN);
	hBsMass_RESvsPTTH2->Fill(Bs_Pt,Bs_Mass-Bs_Mass_GEN);
	hBsPhi_RESvsPTTH2->Fill(Bs_Pt,Bs_Phi-Bs_Phi_GEN);

	hBsMassFromSV_RESvsPT->Fill(Bs_Pt,Bs_MassFromSV-Bs_Mass_GEN);
	hJpsiMassFromSV_RESvsPT->Fill(Jpsi_Pt,Jpsi_MassFromSV-Jpsi_Mass_GEN);
	hPhiMassFromSV_RESvsPT->Fill(Phi_Pt,Phi_MassFromSV-Phi_Mass_GEN);

	hBsMassFromSV_RESvsPTTH2->Fill(Bs_Pt,Bs_MassFromSV-Bs_Mass_GEN);
	hJpsiMassFromSV_RESvsPTTH2->Fill(Jpsi_Pt,Jpsi_MassFromSV-Jpsi_Mass_GEN);
	hPhiMassFromSV_RESvsPTTH2->Fill(Phi_Pt,Phi_MassFromSV-Phi_Mass_GEN);

	hJpsiPt_RESvsPTTH2->Fill(Jpsi_Pt,Jpsi_Pt-Jpsi_Pt_GEN);
	hJpsiEta_RESvsPTTH2->Fill(Jpsi_Pt,Jpsi_Eta-Jpsi_Eta_GEN);
	hJpsiMass_RESvsPTTH2->Fill(Jpsi_Pt,Jpsi_Mass-Jpsi_Mass_GEN);
	hJpsiPhi_RESvsPTTH2->Fill(Jpsi_Pt,Jpsi_Phi-Jpsi_Phi_GEN);

	hPhiPt_RESvsPTTH2->Fill(Phi_Pt,Phi_Pt-Phi_Pt_GEN);
	hPhiEta_RESvsPTTH2->Fill(Phi_Pt,Phi_Eta-Phi_Eta_GEN);
	hPhiMass_RESvsPTTH2->Fill(Phi_Pt,Phi_Mass-Phi_Mass_GEN);
	hPhiPhi_RESvsPTTH2->Fill(Phi_Pt,Phi_Phi-Phi_Phi_GEN);

	hMupPt_RESvsPTTH2->Fill(Mup_Pt,Mup_Pt-Mup_Pt_GEN);
	hMupEta_RESvsPTTH2->Fill(Mup_Pt,Mup_Eta-Mup_Eta_GEN);
	hMupPhi_RESvsPTTH2->Fill(Mup_Pt,Mup_Phi-Mup_Phi_GEN);

	hMumPt_RESvsPTTH2->Fill(Mum_Pt,Mum_Pt-Mum_Pt_GEN);
	hMumEta_RESvsPTTH2->Fill(Mum_Pt,Mum_Eta-Mum_Eta_GEN);
	hMumPhi_RESvsPTTH2->Fill(Mum_Pt,Mum_Phi-Mum_Phi_GEN);
	
	hKpPt_RESvsPTTH2->Fill(Kp_Pt,Kp_Pt-Kp_Pt_GEN);
	hKpEta_RESvsPTTH2->Fill(Kp_Pt,Kp_Eta-Kp_Eta_GEN);
	hKpPhi_RESvsPTTH2->Fill(Kp_Pt,Kp_Phi-Kp_Phi_GEN);
	
	hKmPt_RESvsPTTH2->Fill(Km_Pt,Km_Pt-Km_Pt_GEN);
	hKmEta_RESvsPTTH2->Fill(Km_Pt,Km_Eta-Km_Eta_GEN);
	hKmPhi_RESvsPTTH2->Fill(Km_Pt,Km_Phi-Km_Phi_GEN);
	
	hangle_cospsi_RESvsPTTH2->Fill(Bs_Pt,angle_cospsi-angle_cospsi_GEN);
	hangle_costheta_RESvsPTTH2->Fill(Bs_Pt,angle_costheta-angle_costheta_GEN); 
	hangle_phi_RESvsPTTH2->Fill(Bs_Pt,angle_phi-angle_phi_GEN);

	hlxy_RESvsPTTH2->Fill(Bs_Pt,Lxy-Lxy_GEN);

      }

      // if (Cut(ientry) < 0) continue;
   }// end loop on entries

   out->cd();

   copy->Write();

   hlxy->Write();
   hlxye->Write();


   hBsMassNoRefit->Write();

   hBsPt->Write();
   hBsEta->Write();
   hBsMass->Write();
   hBsPhi->Write();

   hJpsiPt->Write();
   hJpsiEta->Write();
   hJpsiMass->Write();
   hJpsiPhi->Write();

   hPhiPt->Write();
   hPhiEta->Write();
   hPhiMass->Write();
   hPhiPhi->Write();

   hMupPt->Write();
   hMupEta->Write();
   hMupPhi->Write();

   hMumPt->Write();
   hMumEta->Write();
   hMumPhi->Write();

   hKpPt->Write();
   hKpEta->Write();
   hKpPhi->Write();

   hKmPt->Write();
   hKmEta->Write();
   hKmPhi->Write();

   hangle_cospsi->Write();
   hangle_costheta->Write(); 
   hangle_phi->Write();

   hBsMassFromSV->Write();
   hJpsiMassFromSV->Write();
   hPhiMassFromSV->Write();

   hSVrecoMassRES->Write();
   hSVrecoMassRESvsPT->Write();
   hSVrecoMassRESvsPTTH2->Write();

   hPix->Write();
   hStrip->Write();

   if (mc){

     hBsMassFromSV_RES->Write();
     hBsMassFromSV_RESvsPT->Write();
     hBsMassFromSV_RESvsPTTH2->Write();

     hJpsiMassFromSV_RES->Write();
     hJpsiMassFromSV_RESvsPT->Write();
     hJpsiMassFromSV_RESvsPTTH2->Write();

     hPhiMassFromSV_RES->Write();
     hPhiMassFromSV_RESvsPT->Write();
     hPhiMassFromSV_RESvsPTTH2->Write();

     hlxy_GEN->Write();

     hBsPt_GEN->Write();
     hBsEta_GEN->Write();
     hBsMass_GEN->Write();
     hBsPhi_GEN->Write();
     
     hJpsiPt_GEN->Write();
     hJpsiEta_GEN->Write();
     hJpsiMass_GEN->Write();
     hJpsiPhi_GEN->Write();

     hPhiPt_GEN->Write();
     hPhiEta_GEN->Write();
     hPhiMass_GEN->Write();
     hPhiPhi_GEN->Write();

     hMupPt_GEN->Write();
     hMupEta_GEN->Write();
     hMupPhi_GEN->Write();
     
     hMumPt_GEN->Write();
     hMumEta_GEN->Write();
     hMumPhi_GEN->Write();

     hKpPt_GEN->Write();
     hKpEta_GEN->Write();
     hKpPhi_GEN->Write();

     hKmPt_GEN->Write();
     hKmEta_GEN->Write();
     hKmPhi_GEN->Write();

     hangle_cospsi_GEN->Write();
     hangle_costheta_GEN->Write(); 
     hangle_phi_GEN->Write();

     hlxy_RES->Write();

     hBsPt_RES->Write();
     hBsEta_RES->Write();
     hBsMass_RES->Write();
     hBsPhi_RES->Write();
     
     hJpsiPt_RES->Write();
     hJpsiEta_RES->Write();
     hJpsiMass_RES->Write();
     hJpsiPhi_RES->Write();

     hPhiPt_RES->Write();
     hPhiEta_RES->Write();
     hPhiMass_RES->Write();
     hPhiPhi_RES->Write();

     hMupPt_RES->Write();
     hMupEta_RES->Write();
     hMupPhi_RES->Write();
     
     hMumPt_RES->Write();
     hMumEta_RES->Write();
     hMumPhi_RES->Write();

     hKpPt_RES->Write();
     hKpEta_RES->Write();
     hKpPhi_RES->Write();

     hKmPt_RES->Write();
     hKmEta_RES->Write();
     hKmPhi_RES->Write();

     hangle_cospsi_RES->Write();
     hangle_costheta_RES->Write(); 
     hangle_phi_RES->Write();

     hlxy_RESvsPT->Write();

     hBsPt_RESvsPT->Write();
     hBsEta_RESvsPT->Write();
     hBsMass_RESvsPT->Write();
     hBsPhi_RESvsPT->Write();
     
     hJpsiPt_RESvsPT->Write();
     hJpsiEta_RESvsPT->Write();
     hJpsiMass_RESvsPT->Write();
     hJpsiPhi_RESvsPT->Write();

     hPhiPt_RESvsPT->Write();
     hPhiEta_RESvsPT->Write();
     hPhiMass_RESvsPT->Write();
     hPhiPhi_RESvsPT->Write();

     hMupPt_RESvsPT->Write();
     hMupEta_RESvsPT->Write();
     hMupPhi_RESvsPT->Write();
     
     hMumPt_RESvsPT->Write();
     hMumEta_RESvsPT->Write();
     hMumPhi_RESvsPT->Write();

     hKpPt_RESvsPT->Write();
     hKpEta_RESvsPT->Write();
     hKpPhi_RESvsPT->Write();

     hKmPt_RESvsPT->Write();
     hKmEta_RESvsPT->Write();
     hKmPhi_RESvsPT->Write();

     hangle_cospsi_RESvsPT->Write();
     hangle_costheta_RESvsPT->Write(); 
     hangle_phi_RESvsPT->Write();

     hlxy_RESvsPTTH2->Write();

     hBsPt_RESvsPTTH2->Write();
     hBsEta_RESvsPTTH2->Write();
     hBsMass_RESvsPTTH2->Write();
     hBsPhi_RESvsPTTH2->Write();
     
     hJpsiPt_RESvsPTTH2->Write();
     hJpsiEta_RESvsPTTH2->Write();
     hJpsiMass_RESvsPTTH2->Write();
     hJpsiPhi_RESvsPTTH2->Write();

     hPhiPt_RESvsPTTH2->Write();
     hPhiEta_RESvsPTTH2->Write();
     hPhiMass_RESvsPTTH2->Write();
     hPhiPhi_RESvsPTTH2->Write();

     hMupPt_RESvsPTTH2->Write();
     hMupEta_RESvsPTTH2->Write();
     hMupPhi_RESvsPTTH2->Write();
     
     hMumPt_RESvsPTTH2->Write();
     hMumEta_RESvsPTTH2->Write();
     hMumPhi_RESvsPTTH2->Write();

     hKpPt_RESvsPTTH2->Write();
     hKpEta_RESvsPTTH2->Write();
     hKpPhi_RESvsPTTH2->Write();

     hKmPt_RESvsPTTH2->Write();
     hKmEta_RESvsPTTH2->Write();
     hKmPhi_RESvsPTTH2->Write();

     hangle_cospsi_RESvsPTTH2->Write();
     hangle_costheta_RESvsPTTH2->Write(); 
     hangle_phi_RESvsPTTH2->Write();

     if (draw){

       PrintHisto(*hBsMassFromSV_RES);
       PrintHisto2(*hBsMassFromSV_RESvsPTTH2);

       PrintHisto(*hJpsiMassFromSV_RES);
       PrintHisto2(*hJpsiMassFromSV_RESvsPTTH2);

       PrintHisto(*hPhiMassFromSV_RES);
       PrintHisto2(*hPhiMassFromSV_RESvsPTTH2);

       PrintHisto(*hlxy_RES);

       PrintHisto(*hBsPt_RES);
       PrintHisto(*hBsEta_RES);
       PrintHisto(*hBsPhi_RES);

       PrintHisto(*hangle_cospsi_RES);
       PrintHisto(*hangle_costheta_RES);
       PrintHisto(*hangle_phi_RES);

       PrintHisto2(*hlxy_RESvsPTTH2);

       PrintHisto2(*hBsPt_RESvsPTTH2);
       PrintHisto2(*hBsEta_RESvsPTTH2);
       PrintHisto2(*hBsPhi_RESvsPTTH2);

       PrintHisto2(*hangle_cospsi_RESvsPTTH2);
       PrintHisto2(*hangle_costheta_RESvsPTTH2);
       PrintHisto2(*hangle_phi_RESvsPTTH2);
     }
   }
   out->Close();
   //delete copy;
   //delete out;
}

