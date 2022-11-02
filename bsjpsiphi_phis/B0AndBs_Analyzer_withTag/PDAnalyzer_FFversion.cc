#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>
//
#include "PDAnalyzer.h"
//
#include "TMatrixF.h"
#include "TVectorF.h"
#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "Math/LorentzVector.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TVector3.h"
#include "TMatrixTBase.h"
#include "Math/SMatrix.h"
#include "PDMuonVar.cc"
#include "PDSoftMuonMvaEstimator.cc"
#include "AlbertoUtil.cc"
#include "OSMuonMvaTag.cc"

// #include "PDSecondNtupleWriter.h"

using namespace std;

PDAnalyzer::PDAnalyzer() {

  std::cout << "new PDAnalyzer" << std::endl;

  // user parameters are set as names associated to a string, 
  // default values can be set in the analyzer class contructor

  setUserParameter( "verbose", "f" );
  setUserParameter( "ptCut", "4.0" );
  setUserParameter( "maxAngle", "1.0" );
  setUserParameter( "genOnly", "f" ); //to be set to true to run on GEN Only samples
  setUserParameter( "debug", "f" ); // set to true to activate debugging printout
}


PDAnalyzer::~PDAnalyzer() {
}



void PDAnalyzer::beginJob() {

  PDAnalyzerUtil::beginJob();

  // user parameters are retrieved as strings by using their names;
  // numeric parameters ( int, float or whatever ) can be directly seture
  // by passing the corresponding variable,
  // e.g. getUserParameter( "name", x )

  getUserParameter( "verbose", verbose );
  getUserParameter( "ptCut", ptCut );
  getUserParameter( "maxAngle", maxAngle );
  getUserParameter( "genOnly", genOnly );
  getUserParameter( "debug", debug );

  inizializeMuonMvaReader();
  inizializeOSMuonMvaReader();
  bool osInit = inizializeOSMuonCalibration("BuJPsiKData2018", "BuJPsiKMC2018", "BsJPsiPhiMC2018");
  if(!osInit) cout<<endl<<"!!! FAILED TO INIZIALIZED TAG CALIBRATION"<<endl<<endl;

  return;
}


void PDAnalyzer::book() {

  if (!genOnly) outTree = new TTree("OutTree","OutTree");
  else outTree = new TTree("OutTreeGEN","OutTreeGEN");

  if (!genOnly){
  // putting "autoSavedObject" in front of the histo creation 
  // it's automatically marked for saving on file; the option 
  // is uneffective when not using the full utility

    autoSavedObject =
      hPhi          = new TH1F( "Phi mass", "Phi mass", 240, 0.95, 1.1 );
    
    autoSavedObject =
      hJpsi         = new TH1F( "Jpsi mass", "Jpsi mass", 240, 2.8, 3.3 );

    autoSavedObject =
      hPhiHLT       = new TH1F( "Phi mass HLT", "Phi mass HLT", 240, 0.95, 1.1 );
    
    autoSavedObject =
      hJpsiHLT      = new TH1F( "Jpsi mass HLT", "Jpsi mass HLT", 240, 2.8, 3.3 );
    
    autoSavedObject =
      hBs           = new TH1F( "Bs mass", "Bs mass", 50, 5.2, 5.5 );
    
    autoSavedObject =
      hAnglePhi     = new TH1F( "Phi angle", "Phi Angle", 50, -3.2, 3.2 );
    
    autoSavedObject =
      hCosPsi       = new TH1F( "CosPsi", "CosPsi", 51, -1, 1 );
    
    autoSavedObject =
      hCosTheta     = new TH1F( "CosTheta", "CosTheta", 51, -1, 1 );

    outTree->Branch("run",&Trun,"F");
    outTree->Branch("evt",&Tevt,"evt/l");
    outTree->Branch("lumi",&Tlumi,"F");

    outTree->Branch("Tag",&Ttag,"Tag/I");
    outTree->Branch("MisTag",&Tmtag,"F");
    outTree->Branch("MisTagCal",&TmtagCal,"F");
    outTree->Branch("MisTagCalBs",&TmtagCalBs,"F");

    outTree->Branch("isBs",&TisBs,"F");
    
    outTree->Branch("HLT_JpsiTkTk",&THLT_Jtktk,"F");
    outTree->Branch("HLT_JpsiTk",&THLT_Jtk,"F");
    outTree->Branch("HLT_JpsiMu",&THLT_Jmu,"F");

    outTree->Branch("HLT_MatchedJpsi",&TmatchHLTmu,"F");
    outTree->Branch("HLT_MatchedTracks",&TmatchHLTk,"F");

    //polarization angles
    outTree->Branch("angle_cospsi",&Tangle_cospsi,"F");
    outTree->Branch("angle_costheta",&Tangle_costheta,"F");
    outTree->Branch("angle_phi",&Tangle_phi,"F");

    outTree->Branch("PV_cosPoint",&TcosPoint,"F");
    
    outTree->Branch("ctau",&TLxy,"F");
    outTree->Branch("ctauPD",&TLxyPD,"F");
    outTree->Branch("ctauANBS",&TLxyANBS,"F");
    outTree->Branch("ctauAT",&TLxyAT,"F");
    outTree->Branch("ctauCNBS",&TLxyCNBS,"F");
    outTree->Branch("ctauCO",&TLxyCO,"F");
    outTree->Branch("ctau3D",&TLxyz,"F");
  
    outTree->Branch("ctauErr",&TLxyErr,"F");
    outTree->Branch("ctau3DErr",&TLxyzErr,"F");
    
    outTree->Branch("B_VProb",&TVprob,"F");
    outTree->Branch("Jpsi_VProb",&TVprobJpsi,"F");
    outTree->Branch("PhiKstar_VProb",&TVprobPhi,"F");
    
    outTree->Branch("B_Mass",&TBsMass,"F");
    outTree->Branch("B_Pt",&TBsPt,"F");
    outTree->Branch("B_Eta",&TBsEta,"F");
    outTree->Branch("B_Phi",&TBsPhi,"F");

    outTree->Branch("B_MassNoRefit",&TBsMassNoRefit,"F");
    outTree->Branch("Jpsi_MassNoRefit",&TJpsiMassNoRefit,"F");
    outTree->Branch("PhiKstar_MassNoRefit",&TPhiMassNoRefit,"F");

    outTree->Branch("B_MassFromSV",&TBsMassFromSV,"F");
    outTree->Branch("Jpsi_MassFromSV",&TJpsiMassFromSV,"F");
    outTree->Branch("PhiKstar_MassFromSV",&TPhiMassFromSV,"F");

    outTree->Branch("Jpsi_Mass",&TJpsiMass,"F");
    outTree->Branch("Jpsi_Pt",&TJpsiPt,"F");
    outTree->Branch("Jpsi_Eta",&TJpsiEta,"F");
    outTree->Branch("Jpsi_Phi",&TJpsiPhi,"F");

    outTree->Branch("PhiKstar_Mass",&TPhiMass,"F");
    outTree->Branch("PhiKstar_Pt",&TPhiPt,"F");
    outTree->Branch("PhiKstar_Eta",&TPhiEta,"F");
    outTree->Branch("PhiKstar_Phi",&TPhiPhi,"F");
    
    outTree->Branch("Mup_Pt",&TMupPt,"F");
    outTree->Branch("Mup_Eta",&TMupEta,"F");
    outTree->Branch("Mup_Phi",&TMupPhi,"F");
    outTree->Branch("Mup_Hits",&TMupHits,"F");
    
    outTree->Branch("Mum_Pt",&TMumPt,"F");
    outTree->Branch("Mum_Eta",&TMumEta,"F");
    outTree->Branch("Mum_Phi",&TMumPhi,"F");
    outTree->Branch("Mum_Hits",&TMumHits,"F");
    
    outTree->Branch("KmK_Pt",&TKmPt,"F");
    outTree->Branch("KmK_Eta",&TKmEta,"F");
    outTree->Branch("KmK_Phi",&TKmPhi,"F");
    outTree->Branch("KmK_Hits",&TKmHits,"F");
    
    outTree->Branch("KpPi_Pt",&TKpPt,"F");
    outTree->Branch("KpPi_Eta",&TKpEta,"F");
    outTree->Branch("KpPi_Phi",&TKpPhi,"F");
    outTree->Branch("KpPi_Hits",&TKpHits,"F");

    outTree->Branch("Pi_Pt",&TPiPt,"F");
    outTree->Branch("Pi_Eta",&TPiEta,"F");
    outTree->Branch("Pi_Phi",&TPiPhi,"F");
    outTree->Branch("Pi_Hits",&TPiHits,"F");
    

    outTree->Branch("Kaon_Deltadz",&TKaon_Deltadz,"F");
    outTree->Branch("Kaon_DeltaR",&TKaon_DeltaR,"F");

    outTree->Branch("N_PV",&TNpv,"F");

    outTree->Branch("Bc_Mass",&TBcMass,"F");

  }

  if (use_gen || genOnly){

    outTree->Branch("MC_Flavour",&TgenFlavour,"I");
    outTree->Branch("MC_matched",&Tmatched,"F");

    outTree->Branch("angle_cospsi_GEN",&Tangle_cospsiGen,"F");
    outTree->Branch("angle_costheta_GEN",&Tangle_costhetaGen,"F");
    outTree->Branch("angle_phi_GEN",&Tangle_phiGen,"F");

    outTree->Branch("ctau_GEN",&TLxyGen,"F");
    outTree->Branch("ctau3D_GEN",&TLxyzGen,"F");

    outTree->Branch("B_Mass_GEN",&TBsMassGen,"F");
    outTree->Branch("B_Pt_GEN",&TBsPtGen,"F");
    outTree->Branch("B_Eta_GEN",&TBsEtaGen,"F");
    outTree->Branch("B_Phi_GEN",&TBsPhiGen,"F");

    outTree->Branch("Jpsi_Mass_GEN",&TJpsiMassGen,"F");
    outTree->Branch("Jpsi_Pt_GEN",&TJpsiPtGen,"F");
    outTree->Branch("Jpsi_Eta_GEN",&TJpsiEtaGen,"F");
    outTree->Branch("Jpsi_Phi_GEN",&TJpsiPhiGen,"F");

    outTree->Branch("PhiKstar_Mass_GEN",&TPhiMassGen,"F");
    outTree->Branch("PhiKstar_Pt_GEN",&TPhiPtGen,"F");
    outTree->Branch("PhiKstar_Eta_GEN",&TPhiEtaGen,"F");
    outTree->Branch("PhiKstar_Phi_GEN",&TPhiPhiGen,"F");

    outTree->Branch("Mup_Pt_GEN",&TMupPtGen,"F");
    outTree->Branch("Mup_Eta_GEN",&TMupEtaGen,"F");
    outTree->Branch("Mup_Phi_GEN",&TMupPhiGen,"F");

    outTree->Branch("Mum_Pt_GEN",&TMumPtGen,"F");
    outTree->Branch("Mum_Eta_GEN",&TMumEtaGen,"F");
    outTree->Branch("Mum_Phi_GEN",&TMumPhiGen,"F");

    outTree->Branch("KmK_Pt_GEN",&TKmPtGen,"F");
    outTree->Branch("KmK_Eta_GEN",&TKmEtaGen,"F");
    outTree->Branch("KmK_Phi_GEN",&TKmPhiGen,"F");

    outTree->Branch("KpPi_Pt_GEN",&TKpPtGen,"F");
    outTree->Branch("KpPi_Eta_GEN",&TKpEtaGen,"F");
    outTree->Branch("KpPi_Phi_GEN",&TKpPhiGen,"F");

  }

  return;
}


void PDAnalyzer::reset() {
  autoReset();
  return;
}


bool PDAnalyzer::analyze( int entry, int event_file, int event_tot ) {

  // flag to be set "true" or "false" for events accepted or rejected
  bool flag;

  if (!genOnly){

    if ( (!(event_tot%10) && event_tot<100 ) || 
    (!(event_tot %100) && event_tot<1000 ) || 
    (!(event_tot %1000)&& event_tot<10000 ) || 
    (!(event_tot %10000) && event_tot<100000 ) || 
    (!(event_tot %100000) && event_tot<1000000 ) || 
    (!(event_tot %1000000) && event_tot<10000000 ) )
        cout << " == at event " << event_file << " " << event_tot << endl;

    computeMuonVar();   //compute muon variable for soft id
    convSpheCart(jetPt, jetEta, jetPhi, jetPx, jetPy, jetPz);
    convSpheCart(muoPt, muoEta, muoPhi, muoPx, muoPy, muoPz);
    convSpheCart(trkPt, trkEta, trkPhi, trkPx, trkPy, trkPz);
    convSpheCart(pfcPt, pfcEta, pfcPhi, pfcPx, pfcPy, pfcPz);

    //FF
    //  int mu1_idx, mu2_idx, k1_idx, k2_idx;
    
    const float MUMASS=0.105658;
    const float KMASS=0.493677;
    const float BSMASS=5.36689;
    const float PIMASS=0.139570;
    const float B0MASS=5.27963;

    if (debug) cout << "PDG values for BS Mass " << BSMASS << " B0 Mass " << B0MASS << " Mu Mass " << MUMASS << " Pi Mass " << PIMASS << " K Mass " << KMASS << endl;

    int iSV, nsv=0, nCandidates=0, MuMatchHLT=0, KMatchHLT=0, SelectedSV=-999;
    float selVprob=-999;
    
    TmatchHLTmu=-999;
    TmatchHLTk=-999;
    TisBs=-999;//these need to be initialized here and not in "resetTree()"

    for ( iSV = 0; iSV < nSVertices; ++iSV ){
      if ( svtType->at( iSV ) != PDEnumString::svtBsJPsiPhi) continue;
      nsv++;

      MuMatchHLT=0; KMatchHLT=0;
      float tmpVprob=TMath::Prob(svtChi2->at(iSV),svtNDOF->at(iSV)); 

      const vector<int>& sub = subVtxFromSV( iSV );
      
      if (sub.size() !=2) continue;   
      
      for (uint k=0; k<sub.size(); k++){
        if (svtType->at(sub[k]) == PDEnumString::svtJPsi){
          const vector<int>& mutks = tracksFromSV( sub[k] );
        
          if (mutks.size()!=2) continue;

          int m1=-999, m2=-999;
          if (use_hlto) {
            
            m1=HLTMatch(mutks[0],1);
            m2=HLTMatch(mutks[1],1);
            
            if (m1 != m2 && m1>=0 && m2>=0) {
              MuMatchHLT=1;
              TLorentzVector hltm1, hltm2;
              hltm1.SetPtEtaPhiM(hltPt->at(m1),hltEta->at(m1),hltPhi->at(m1),MUMASS);
              hltm2.SetPtEtaPhiM(hltPt->at(m2),hltEta->at(m2),hltPhi->at(m2),MUMASS);
              hJpsiHLT->Fill((hltm1+hltm2).M());
            }
          }
        }
        else if (svtType->at(sub[k]) == PDEnumString::svtPhi){
          const vector<int>& ktks = tracksFromSV( sub[k] );

          if (ktks.size()!=2) continue;

          int k1=-999, k2=-999;
          if (use_hlto) {
            //cout << "matching" << endl; 
            k1=HLTMatch(ktks[0],0);
            k2=HLTMatch(ktks[1],0);
            if (k1 != k2 && k1>=0 && k2>=0) {
              KMatchHLT=1;
              //draw Psi Mass at HLT
              TLorentzVector hltk1, hltk2;
              hltk1.SetPtEtaPhiM(hltPt->at(k1),hltEta->at(k1),hltPhi->at(k1),KMASS);
              hltk2.SetPtEtaPhiM(hltPt->at(k2),hltEta->at(k2),hltPhi->at(k2),KMASS);
              hPhiHLT->Fill((hltk1+hltk2).M());
            }
          }
        }
      }
      if (tmpVprob>selVprob && MuMatchHLT==1){
        SelectedSV=iSV; 
        selVprob=tmpVprob; 
        TisBs=1;
        TmatchHLTmu=MuMatchHLT;
        TmatchHLTk=KMatchHLT;
        if (debug) cout << "Selected Bs SV "  << SelectedSV << endl;
      } 
    }
      
    //Put the B0 code here, to be executed only if the Bs has not been found
    if (TisBs!=1){
      
      selVprob=-999;

      for ( iSV = 0; iSV < nSVertices; ++iSV ) {
        if ( svtType->at( iSV ) != PDEnumString::svtBdJPsiKx) continue;
    
        float tmpVprob=TMath::Prob(svtChi2->at(iSV),svtNDOF->at(iSV)); 
        MuMatchHLT=0; KMatchHLT=0;

        const vector<int>& sub = subVtxFromSV( iSV );
          
        if (sub.size() !=2) continue;   
      
        for (uint k=0; k<sub.size(); k++){
        
          if (svtType->at(sub[k]) == PDEnumString::svtJPsi){
            const vector<int>& mutks = tracksFromSV( sub[k] );

            if (mutks.size()!=2) continue;

            int m1=-999,m2=-999;
            if (use_hlto) {
               m1=HLTMatch(mutks[0],1);
               m2=HLTMatch(mutks[1],1);
              if (m1 != m2 && m1>=0 && m2>=0) {
            MuMatchHLT=1;
              }
            }
          }
          else if (svtType->at(sub[k]) == PDEnumString::svtKx0){
            const vector<int>& ktks = tracksFromSV( sub[k] );

            if (ktks.size()!=2) continue;
            
            int k1=-999,k2=-999;
            if (use_hlto) {
              k1=HLTMatch(ktks[0],0);
              k2=HLTMatch(ktks[1],0);
              if (k1 != k2 && k1>=0 && k2>=0) {
            KMatchHLT=1;
              }
            }
          }
        }

        if (tmpVprob>selVprob && MuMatchHLT==1) {
          SelectedSV=iSV; 
          selVprob=tmpVprob; 
          TisBs=0;
          TmatchHLTmu=MuMatchHLT;
          TmatchHLTk=KMatchHLT;
          if (debug) cout << "Selected B0 SV "  << SelectedSV << " Vprob " << selVprob << " evt " << eventNumber << endl;
        }   
      }
    }
    //end code on the B0
      
    if (TisBs==-999){
      if (debug) cout << "no Bs or B0 vertex found in the event, returning" << endl;
      return false;
    }//return if no Bs or B0 candidate is found
    //if (nsv>1 && use_gen) {
    //  if (debug) cout << nsv << " Secondary vertices, returning " << endl; 
    //   return false;
    //} //in MC avoid events with a double candidate

    /*int nBs=0;
    if (use_gen){
      for (int p=0; p< nGenP; ++p){
    if (abs(genId->at(p))==531 && (abs(genId->at(genMother->at(p)))==533 || abs(genId->at(genMother->at(p)))==21 || abs(genId->at(genMother->at(p)))<10)) nBs++;
      }
    
      if (nBs>1){if (debug) cout << nBs << " GenBs, returning " << endl; return false;} // consider only MC events with a single Bs
    } */ //why do remove the events with two Bs?
     
      //Beginning of the reconstruction code
      
    if (debug) cout << " NEW EVENT---------------------------------------" << endl;
      
    for ( iSV = 0; iSV < nSVertices; ++iSV ){
      
      if (iSV!=SelectedSV) continue;
      
      //Make sure the vertex association is right
      if ( TisBs==1 && svtType->at( iSV ) != PDEnumString::svtBsJPsiPhi) {cout<<"Something bad happened"<<endl; continue;}
      if ( TisBs==0 && svtType->at( iSV ) != PDEnumString::svtBdJPsiKx) {cout<<"Something bad happened"<<endl; continue;}

      if (debug) cout << "Selected Vertex " << iSV << " IsBs " << TisBs << endl;
      
      resetTree();//reset the Tree branches variables
      
      Trun=runNumber;
      Tevt=ULong64_t(eventNumber);
      Tlumi=lumiSection;
      TNpv=nPVertices;

      if (debug) cout << "+++++++++++ Identifying Objects ++++++++++++" << endl;
      
      if (use_hlts){//filling trigger info

        if (debug) cout << "Checking HLT decision" << endl;

        for (int path=0; path<nHLTStatus; path++){
          
          if (!hltRun->at(path)) continue;
          
          if (!(hltPath->at(path)==PDEnumString::HLT_DoubleMu4_JpsiTrkTrk_Displaced_v || hltPath->at(path)==PDEnumString::HLT_Dimuon0_Jpsi_Muon_v ||  hltPath->at(path)==PDEnumString::HLT_DoubleMu4_JpsiTrk_Displaced_v ||  hltPath->at(path)==PDEnumString::HLT_Dimuon0_Jpsi3p5_Muon2_v)) continue; 
          
          if (hltPath->at(path)==PDEnumString::HLT_DoubleMu4_JpsiTrkTrk_Displaced_v) THLT_Jtktk=hltAccept->at(path);
          if (hltPath->at(path)==PDEnumString::HLT_Dimuon0_Jpsi_Muon_v || hltPath->at(path)==PDEnumString::HLT_Dimuon0_Jpsi3p5_Muon2_v) THLT_Jmu=hltAccept->at(path);
          if (hltPath->at(path)==PDEnumString::HLT_DoubleMu4_JpsiTrk_Displaced_v) THLT_Jtk=hltAccept->at(path);

          if (debug) cout << "TkTk trigger " << THLT_Jtktk << endl;

        }
      }

      if (TisBs==1){
    
        std::vector<int> mu_idx, k_idx, mu_idx_refit, k_idx_refit;
          
        TLorentzVector mu_p, mu_m, k_p, k_m, mu_p_norefit, mu_m_norefit, k_p_norefit, k_m_norefit, Bs, Jpsi, Phi, BsNoRefit, JpsiNoRefit, PhiNoRefit, pi, Bc;
        TVector3 SVpos, PVpos;
          
        float VProb=-999, VProbJpsi=-999, VProbPhi=-999;
        float mupH=-999,mumH=-999,kpH=-999, kmH=-999;
          
        double angle_costheta;
        double angle_phi;
        double angle_cospsi;
          
        const vector<int>& tks = tracksFromSV( iSV );
        int n = tks.size();

        TMatrixF covP_Track0(3,3);      
        TMatrixF covP_Track1(3,3);      

        if (!n) continue;
      
        SVpos.SetXYZ(svtX->at( iSV ),svtY->at( iSV ),svtZ->at( iSV ));

        //filling some histos
        hBs->Fill( svtMass->at( iSV ) );
        VProb=TMath::Prob(svtChi2->at(iSV),svtNDOF->at(iSV));    

        TBsMassFromSV=svtMass->at( iSV );

        //Select the sub-secondary vertices (Jpsi and Phi)
        const vector<int>& sub = subVtxFromSV( iSV );
          
        if (sub.size() !=2) continue;   
      
        for (uint k=0; k<sub.size(); k++){
        
          if (svtType->at(sub[k]) == PDEnumString::svtJPsi){
            const vector<int>& mutks = tracksFromSV( sub[k] );
         
            for (uint t=0; t< mutks.size(); ++t) {
              //cout << " m " << mutks[t];
              mu_idx.push_back(mutks[t]); 
              mu_idx_refit.push_back(vpIndex(mutks[t],iSV));//sub[k]));
            }
          
            hJpsi->Fill(svtMass->at(sub[k]));
            TJpsiMassFromSV=svtMass->at(sub[k]);
          
            VProbJpsi=TMath::Prob(svtChi2->at(sub[k]),svtNDOF->at(sub[k]));
          
            if (debug) {
              if (debug) cout << "Size of kaon collection " << mutks.size() << endl;
              cout << "J/psi mass (from vertex fit) " << svtMass->at(sub[k]) << " n tracks in SV " << mutks.size()  << " n tracks saved " << mu_idx.size() << " n tracks refitted " << mu_idx_refit.size() << endl;
              cout << "Muon indices norefit ";
              for (uint t=0; t< mutks.size(); ++t) cout << mutks[t] << " ";
              cout << endl;
              cout << "Muon indices refit ";
              for (uint t=0; t< mutks.size(); ++t) cout << mu_idx_refit[t] << " ";
              cout << endl;
              cout << endl;
            }
          }else if (svtType->at(sub[k]) == PDEnumString::svtPhi){
            const vector<int>& ktks = tracksFromSV( sub[k] );
      
            uint track0id=99999;
            uint track1id=99999;         
            for (uint t=0; t<tppTrk->size();t++) {
               if (tppTrk->at(t)==ktks[0]) track0id=t;
               if (tppTrk->at(t)==ktks[1]) track1id=t;
            }
 

 
            for (uint t=0; t< ktks.size(); ++t) {
              //cout << " k " << ktks[t];
              k_idx.push_back(ktks[t]); 
              k_idx_refit.push_back(vpIndex(ktks[t],iSV));
            }

            float trk0Qop = tppQop->at(track0id);
            float trk0Q = trk0Qop/TMath::Abs(trk0Qop);
            float trk0P = TMath::Abs(1/trk0Qop);
            float trk0Lam = tppLam->at(track0id);
            float trk0Phi = tppPhi->at(track0id);
            float trk0Dxy = tppDxy->at(track0id);
            float trk0Dsx = tppDsx->at(track0id);
            float trk1Qop = tppQop->at(track1id);
            float trk1Q = trk1Qop/TMath::Abs(trk1Qop);
            float trk1P = TMath::Abs(1/trk1Qop);
            float trk1Lam = tppLam->at(track1id);
            float trk1Phi = tppPhi->at(track1id);
            float trk1Dxy = tppDxy->at(track1id);
            float trk1Dsx = tppDsx->at(track1id);
            TKaon_Deltadz=tppDsx->at(track0id)/cos(tppLam->at(track0id))-tppDsx->at(track1id)/cos(tppLam->at(track1id));
            TLorentzVector track0par;
            TLorentzVector track1par;
            track0par.SetPtEtaPhiM(trkPt->at(ktks[0]),trkEta->at(ktks[0]),trkPhi->at(ktks[0]),PIMASS);
            track1par.SetPtEtaPhiM(trkPt->at(ktks[1]),trkEta->at(ktks[1]),trkPhi->at(ktks[1]),PIMASS);
            TKaon_DeltaR=track0par.DeltaR(track1par);

            TMatrixF covTrack0(3,3);
            float covTrack0Array[]={tppSQopQop->at(track0id),tppSQopLam->at(track0id),tppSQopPhi->at(track0id),tppSQopLam->at(track0id),tppSLamLam->at(track0id),tppSLamPhi->at(track0id),tppSQopPhi->at(track0id),tppSLamPhi->at(track0id),tppSPhiPhi->at(track0id)};
            covTrack0.SetMatrixArray(covTrack0Array);
            TMatrixF covTrack1(3,3);
            float covTrack1Array[]={tppSQopQop->at(track1id),tppSQopLam->at(track1id),tppSQopPhi->at(track1id),tppSQopLam->at(track1id),tppSLamLam->at(track1id),tppSLamPhi->at(track1id),tppSQopPhi->at(track1id),tppSLamPhi->at(track1id),tppSPhiPhi->at(track1id)};
            covTrack1.SetMatrixArray(covTrack1Array);
            // building jacobian matrix following this document:
            // http://cms.cern.ch/iCMS/jsp/openfile.jsp?type=NOTE&year=2006&files=NOTE2006_001.pdf" A. Strandlie, W. Wittek, "Propagation of Covariance Matrices...", CMS Note 2006/001
            TMatrixF jacTrack0(3,3);
            float jacTrack0Array[]={-trk0Q*trk0P*trk0P*cos(trk0Lam) * cos(trk0Phi),-trk0P*sin(trk0Lam)*cos(trk0Phi),-trk0P*cos(trk0Lam)*sin(trk0Phi),-trk0Q*trk0P*trk0P*cos(trk0Lam),-trk0P*sin(trk0Lam)*sin(trk0Phi),trk0P*cos(trk0Lam)*cos(trk0Phi),-trk0Q*trk0P*trk0P*sin(trk0Lam),trk0P*cos(trk0Lam),0};
            jacTrack0.SetMatrixArray(jacTrack0Array);
            TMatrixF jacTrack1(3,3);
            float jacTrack1Array[]={-trk1Q*trk1P*trk1P*cos(trk1Lam)*cos(trk1Phi),-trk1P*sin(trk1Lam)*cos(trk1Phi),-trk1P*cos(trk1Lam)*sin(trk1Phi),-trk1Q*trk1P*trk1P*cos(trk1Lam),-trk1P*sin(trk1Lam)*sin(trk1Phi),trk1P*cos(trk1Lam)*cos(trk1Phi),-trk1Q*trk1P*trk1P*sin(trk1Lam),trk1P*cos(trk1Lam),0};
            jacTrack1.SetMatrixArray(jacTrack1Array);

            covP_Track0=jacTrack0*covTrack0*jacTrack0.T();
            covP_Track1=jacTrack1*covTrack1*jacTrack1.T();

            hPhi->Fill(svtMass->at(sub[k]));
            TPhiMassFromSV=svtMass->at(sub[k]);
          
            VProbPhi=TMath::Prob(svtChi2->at(sub[k]),svtNDOF->at(sub[k]));
          
            if (debug){
              cout << "Phi mass (from vertex fit)" << svtMass->at(sub[k]) << " n tracks  in SV " << ktks.size()  << " n tracks saved " << k_idx.size() << " n tracks refitted " << k_idx_refit.size() << endl;
              cout << "Kaon indices norefit ";
              if (debug) cout << "Size of kaon collection " << ktks.size() << endl;
              for (uint t=0; t< ktks.size(); ++t) cout << ktks[t] << " ";
              cout << endl;
              cout << "Kaon indices refit ";
              for (uint t=0; t< ktks.size(); ++t) cout << k_idx_refit[t] << " ";
              cout << endl;
              cout << endl;
            }
          }
        }
      
        if ( (k_idx_refit[0] >= int(tvpPt->size()) ||  k_idx_refit[1] >= int(tvpPt->size()))  ||  (mu_idx_refit[0] >= int(tvpPt->size()) ||  mu_idx_refit[1] >= int(tvpPt->size())) ){
          cout << "PROBLEM IN THE vpIndex FUNCTION in run " << runNumber << " event " << eventNumber << " lumi " << lumiSection << endl;
          cout << "Size of the full track collection " << nTracks << " N refitted " << tvpPt->size() << endl;
          cout << " refitted indices k1 " <<  k_idx_refit[0] << " k2 " <<  k_idx_refit[1] << " mu 1 " << mu_idx_refit[0] << " m2 " <<  mu_idx_refit[1] << endl;
          return false;
        }

        //analyze the muons and tracks now
          
        if (k_idx.size()!=2 || mu_idx.size()!=2 || k_idx_refit.size()!=2 || mu_idx_refit.size()!=2 ) continue;

        if (trkCharge->at(mu_idx[0]) >0){
          mu_p_norefit.SetPtEtaPhiM(trkPt->at(mu_idx[0]),trkEta->at(mu_idx[0]),trkPhi->at(mu_idx[0]),MUMASS);
          mu_m_norefit.SetPtEtaPhiM(trkPt->at(mu_idx[1]),trkEta->at(mu_idx[1]),trkPhi->at(mu_idx[1]),MUMASS);
          mupH=trkHitPattern->at(mu_idx[0]); mumH=trkHitPattern->at(mu_idx[1]);
        } else if (trkCharge->at(mu_idx[0]) < 0){
          mu_m_norefit.SetPtEtaPhiM(trkPt->at(mu_idx[0]),trkEta->at(mu_idx[0]),trkPhi->at(mu_idx[0]),MUMASS);
          mu_p_norefit.SetPtEtaPhiM(trkPt->at(mu_idx[1]),trkEta->at(mu_idx[1]),trkPhi->at(mu_idx[1]),MUMASS);
          mupH=trkHitPattern->at(mu_idx[1]); mumH=trkHitPattern->at(mu_idx[0]);
        }
          
        if (trkCharge->at(k_idx[0]) >0){
          k_p_norefit.SetPtEtaPhiM(trkPt->at(k_idx[0]),trkEta->at(k_idx[0]),trkPhi->at(k_idx[0]),KMASS);
          k_m_norefit.SetPtEtaPhiM(trkPt->at(k_idx[1]),trkEta->at(k_idx[1]),trkPhi->at(k_idx[1]),KMASS);
          kpH=trkHitPattern->at(k_idx[0]); kmH=trkHitPattern->at(mu_idx[1]);
        } else if (trkCharge->at(k_idx[0]) < 0){
          k_m_norefit.SetPtEtaPhiM(trkPt->at(k_idx[0]),trkEta->at(k_idx[0]),trkPhi->at(k_idx[0]),KMASS);
          k_p_norefit.SetPtEtaPhiM(trkPt->at(k_idx[1]),trkEta->at(k_idx[1]),trkPhi->at(k_idx[1]),KMASS);
          kpH=trkHitPattern->at(k_idx[1]); kmH=trkHitPattern->at(mu_idx[0]);    
        }
      
        //adding refitted tracks
          
        if (trkCharge->at(mu_idx[0]) >0){
          mu_p.SetPtEtaPhiM(tvpPt->at(mu_idx_refit[0]),tvpEta->at(mu_idx_refit[0]),tvpPhi->at(mu_idx_refit[0]),MUMASS);
          mu_m.SetPtEtaPhiM(tvpPt->at(mu_idx_refit[1]),tvpEta->at(mu_idx_refit[1]),tvpPhi->at(mu_idx_refit[1]),MUMASS);
        } else if (trkCharge->at(mu_idx[0]) < 0){
          mu_m.SetPtEtaPhiM(tvpPt->at(mu_idx_refit[0]),tvpEta->at(mu_idx_refit[0]),tvpPhi->at(mu_idx_refit[0]),MUMASS);
          mu_p.SetPtEtaPhiM(tvpPt->at(mu_idx_refit[1]),tvpEta->at(mu_idx_refit[1]),tvpPhi->at(mu_idx_refit[1]),MUMASS);
        }
          
        if (trkCharge->at(k_idx[0]) >0){
          k_p.SetPtEtaPhiM(tvpPt->at(k_idx_refit[0]),tvpEta->at(k_idx_refit[0]),tvpPhi->at(k_idx_refit[0]),KMASS);
          k_m.SetPtEtaPhiM(tvpPt->at(k_idx_refit[1]),tvpEta->at(k_idx_refit[1]),tvpPhi->at(k_idx_refit[1]),KMASS);
        } else if (trkCharge->at(k_idx[0]) < 0){
          k_m.SetPtEtaPhiM(tvpPt->at(k_idx_refit[0]),tvpEta->at(k_idx_refit[0]),tvpPhi->at(k_idx_refit[0]),KMASS);
          k_p.SetPtEtaPhiM(tvpPt->at(k_idx_refit[1]),tvpEta->at(k_idx_refit[1]),tvpPhi->at(k_idx_refit[1]),KMASS);
        }

        if (debug) {
          cout << "+++++++++++++ Building Candidates ++++++++++++++" << endl;
          cout << "From refitted tracks: J/psi mass " << (mu_p+mu_m).Mag() << " Phi mass " << (k_p+k_m).Mag() << endl;
          cout << endl;
          cout << "From native tracks: J/psi mass " << (mu_p_norefit+mu_m_norefit).Mag() << " Phi mass " << (k_p_norefit+k_m_norefit).Mag() << endl;
          cout << endl;
        }

        Jpsi=mu_p + mu_m;
        
        Phi=k_p + k_m;
          
        Bs = Jpsi + Phi;
          
        JpsiNoRefit=mu_p_norefit+mu_m_norefit;
        
        PhiNoRefit=k_p_norefit + k_m_norefit;
          
        BsNoRefit = JpsiNoRefit + PhiNoRefit;
          
        int pvIndex=-999;
        //selecting the PV and computing Lxy and Error
          
        FindPV(SVpos, PVpos, Bs, pvIndex);

        if (pvIndex < 0){
          if (debug) cout << "NO PV FOUND FOR THIS EVENT" << endl;
          return false; //return if no PV is identified
        }

        int ipv = svtPVtx->at(SelectedSV);
        if(ipv<0) return false;

        int ipv_refitANBS = -1;
        int ipv_refitAT = -1;
        int ipv_refitCNBS = -1;
        int ipv_refitCO = -1;
        for(int isub=0; isub<nCompVts; ++isub){
          if(subPart->at(isub)!=SelectedSV) continue;
          int ivtx = subSVtx->at(isub);
          if(svtType->at(ivtx)==PDEnumString::refittedPVAllTrkNoBS) {ipv_refitANBS = ivtx; } 
          if(svtType->at(ivtx)==PDEnumString::refittedPVAllTrk) {ipv_refitAT = ivtx; } 
          if(svtType->at(ivtx)==PDEnumString::refittedPVRCOnlyNoBS) {ipv_refitCNBS = ivtx; }
          if(svtType->at(ivtx)==PDEnumString::refittedPVRCOnly) {ipv_refitCO = ivtx; }
        }
        //cout<<"PV "<<pvtX->at(ipv)<<" "<<pvtY->at(ipv)<<" "<<pvtZ->at(ipv)<<endl;
        //cout<<"PV refit all "<<svtX->at(ipv_refitANBS)<<" "<<svtY->at(ipv_refitANBS)<<" "<<svtZ->at(ipv_refitANBS)<<endl;
        //cout<<"PV refit all with BS "<<svtX->at(ipv_refitAT)<<" "<<svtY->at(ipv_refitAT)<<" "<<svtZ->at(ipv_refitAT)<<endl;
        //cout<<"PV refit reco tracks "<<svtX->at(ipv_refitCNBS)<<" "<<svtY->at(ipv_refitCNBS)<<" "<<svtZ->at(ipv_refitCNBS)<<endl;
        //cout<<"PV refit reco track with BS "<<svtX->at(ipv_refitCO)<<" "<<svtY->at(ipv_refitCO)<<" "<<svtZ->at(ipv_refitCO)<<endl;
        TVector3 PVposPD;
        PVposPD.SetXYZ(pvtX->at(ipv),pvtY->at(ipv),pvtZ->at(ipv));
        TVector3 PVSVdistPD=SVpos-PVposPD;
        TVector3 PVposANBS;
        if (ipv_refitANBS>=0) PVposANBS.SetXYZ(svtX->at(ipv_refitANBS),svtY->at(ipv_refitANBS),svtZ->at(ipv_refitANBS));
        TVector3 PVSVdistANBS=SVpos-PVposANBS;
        TVector3 PVposAT;
        if (ipv_refitAT>=0) PVposAT.SetXYZ(svtX->at(ipv_refitAT),svtY->at(ipv_refitAT),svtZ->at(ipv_refitAT));
        TVector3 PVSVdistAT=SVpos-PVposAT;
        TVector3 PVposCNBS;
        if (ipv_refitCNBS>=0) PVposCNBS.SetXYZ(svtX->at(ipv_refitCNBS),svtY->at(ipv_refitCNBS),svtZ->at(ipv_refitCNBS));
        TVector3 PVSVdistCNBS=SVpos-PVposCNBS;
        TVector3 PVposCO;
        if (ipv_refitCO>=0) PVposCO.SetXYZ(svtX->at(ipv_refitCO),svtY->at(ipv_refitCO),svtZ->at(ipv_refitCO));
        TVector3 PVSVdistCO=SVpos-PVposCO;


        inizializeOsMuonTagVars();
        setVtxOsMuonTag(iSV, pvIndex); //set PV and SVT for tagging class
        bool istagged = makeOsMuonTagging();

        Ttag = getOsMuonTag(); //get tag decision -1, +1 or 0 
   
        if(istagged){
          Tmtag = getOsMuonTagMistagProbRaw();
          TmtagCal = getOsMuonTagMistagProbCalProcess();
          TmtagCalBs = getOsMuonTagMistagProbCalProcessBuBs();
        }



        float Lxy, Lxy3D;
        float LxyPD,LxyANBS,LxyAT,LxyCNBS,LxyCO;      

        TVector3 BsMomentum=Bs.Vect();//3Momentum of Bs
        TVector3 PVSVdist=SVpos-PVpos;
          
        Lxy3D=BSMASS/Bs.P()*PVSVdist.Dot(BsMomentum)/BsMomentum.Mag();//(PVpos-SVpos).Mag();
        
        BsMomentum.SetZ(0.);
        PVSVdist.SetZ(0.);
      
        Lxy=BSMASS*PVSVdist.Dot(BsMomentum)/BsMomentum.Mag2(); //Transverse ctau after setting to 0 the Z components of the vectors above
        LxyPD=BSMASS*PVSVdistPD.Dot(BsMomentum)/BsMomentum.Mag2();      
        if (ipv_refitANBS>=0) { LxyANBS=BSMASS/Bs.Pt()*PVSVdistANBS.Dot(BsMomentum)/BsMomentum.Mag(); } else { LxyANBS=-99999;}
        if (ipv_refitAT>=0) {LxyAT=BSMASS/Bs.Pt()*PVSVdistAT.Dot(BsMomentum)/BsMomentum.Mag(); } else {LxyAT=-99999;}
        if (ipv_refitCNBS>=0) {LxyCNBS=BSMASS/Bs.Pt()*PVSVdistCNBS.Dot(BsMomentum)/BsMomentum.Mag(); } else {LxyCNBS=-99999;}      
        if (ipv_refitCO>=0) {LxyCO=BSMASS/Bs.Pt()*PVSVdistCO.Dot(BsMomentum)/BsMomentum.Mag(); } else {LxyCO=-99999;}      

        if (debug) cout << "Lxy " << Lxy << " Lxy 3D " << Lxy3D << endl;
          
        TMatrixF covSV(3,3);
        float covSVArray[]={svtSxx->at(iSV),svtSxy->at(iSV),svtSxz->at(iSV),svtSxy->at(iSV),svtSyy->at(iSV),svtSyz->at(iSV), svtSxz->at(iSV),svtSyz->at(iSV),svtSzz->at(iSV)};
        covSV.SetMatrixArray(covSVArray);
          
        TMatrixF covPV(3,3);
        float covPVArray[]={pvtSxx->at(pvIndex),pvtSxy->at(pvIndex),pvtSxz->at(pvIndex),pvtSxy->at(pvIndex),pvtSyy->at(pvIndex),pvtSyz->at(pvIndex), pvtSxz->at(pvIndex),pvtSyz->at(pvIndex),pvtSzz->at(pvIndex)};
        covPV.SetMatrixArray(covPVArray);
          
        TMatrixF covTot= covSV+covPV;
          
        float distArray[]={float((SVpos-PVpos).X()),float((SVpos-PVpos).Y()),float((SVpos-PVpos).Z())};
        TVectorF diff(3,distArray);

        float distArray2D[]={float((SVpos-PVpos).X()),float((SVpos-PVpos).Y()),0.};
        TVectorF diff2D(3,distArray2D);
          
        float LxyErr=-999,Lxy3DErr=-999;
          
        if (diff2D.Norm2Sqr()==0 || diff.Norm2Sqr()==0) continue; //if the secondary vertex is exactly the same as PV continue
          
        TVector3 diff2DV(diff2D[0],diff2D[1],diff2D[2]); 
        LxyErr=TMath::Abs(BSMASS/BsMomentum.Mag()*sqrt(covTot.Similarity(diff2D))/sqrt(diff2D.Norm2Sqr())*cos(diff2DV.Angle(BsMomentum))); 
        Lxy3DErr= BSMASS/Bs.P()*sqrt(covTot.Similarity(diff))/sqrt(diff.Norm2Sqr()); 
          
        float temparray[]={float(BsMomentum[0]),float(BsMomentum[1]),float(BsMomentum[2])}; 
        TVectorF BsMomentum2DF(3,temparray);
        
        float LxyErrBis= sqrt(pow(LxyErr,2)+pow(B0MASS/BsMomentum.Mag2()*diff2D.Norm1()*sqrt(TMath::Abs(covP_Track0.Similarity(BsMomentum2DF)))/sqrt(BsMomentum2DF.Norm2Sqr())*cos(diff2DV.Angle(BsMomentum)),2)); 
        cout<< Lxy << " " <<  LxyErr<< " "<<LxyErrBis<<" "<< (LxyErr-LxyErrBis)/LxyErr*100<<endl;

        if (debug) cout << "LxyErr " << LxyErr << " Lxy 3DErr " << Lxy3DErr << endl;
        
        //adding a pion to search for Bc->Bs pi      
        int PionIndex=-999;

        float dRtmp=1;
        float pttmp=1;

        for (int i=0; i<nTracks; ++i){
          TLorentzVector pi;
          pi.SetPtEtaPhiM(trkPt->at(i),trkEta->at(i),trkPhi->at(i),PIMASS);
          if (pi.Pt() < 1) continue;
          if (pi.DeltaR(mu_p) < 1.e-3 || pi.DeltaR(mu_m) < 1.e-3 || pi.DeltaR(k_p) < 1.e-3 || pi.DeltaR(k_m) < 1.e-3) continue;
          if ((Bs+pi).M() < 5.0 || (Bs+pi).M() > 6.6) continue;
          float dR=Bs.DeltaR(pi);
          if (dR < dRtmp && pi.Pt()> pttmp){
            PionIndex=i;
            pttmp=pi.Pt();
            dRtmp=dR;
          }
        }

        if (PionIndex>=0){
          TLorentzVector pi;
          pi.SetPtEtaPhiM(trkPt->at(PionIndex),trkEta->at(PionIndex),trkPhi->at(PionIndex),PIMASS);
          TPiPt=trkPt->at(PionIndex);
          TPiPhi=trkPhi->at(PionIndex);
          TPiEta=trkEta->at(PionIndex);
          TPiHits=trkHitPattern->at(PionIndex);
          TBcMass=(Bs+pi).M();
        }

        //computing angles
          
        // the betas for the boost
        TVector3 p3_JPsi;
        p3_JPsi = Jpsi.Vect();
        p3_JPsi *= -1./Jpsi.E();

        // the boost matrix
        TLorentzRotation boost_jpsi(p3_JPsi);
        TLorentzVector p_JPsi_JPsi = boost_jpsi.VectorMultiplication(Jpsi);
          
        // the different momenta in the new frame                                                                                                       
        TLorentzVector p_JPsi_muplus = boost_jpsi.VectorMultiplication(mu_p); 
        TLorentzVector p_JPsi_Kplus = boost_jpsi.VectorMultiplication(k_p); 
        TLorentzVector p_JPsi_phi = boost_jpsi.VectorMultiplication(Phi);                                                                       
          
        // the 3-momenta
        TVector3 p3_JPsi_muplus = p_JPsi_muplus.Vect();
        TVector3 p3_JPsi_Kplus = p_JPsi_Kplus.Vect();
        TVector3 p3_JPsi_phi = p_JPsi_phi.Vect();
      
        // coordinate system
        TVector3 x,y,z;
        x = p3_JPsi_phi.Unit();
        y = p3_JPsi_Kplus.Unit() - (x * (x * p3_JPsi_Kplus.Unit()));
        y = y.Unit();
        z = x.Cross(y);
          
        // Transversity Basis
        angle_costheta = p3_JPsi_muplus.Unit() * z;
          
        double cos_phi = p3_JPsi_muplus.Unit() * x / TMath::Sqrt(1 - angle_costheta*angle_costheta);
        double sin_phi = p3_JPsi_muplus.Unit() * y / TMath::Sqrt(1 - angle_costheta*angle_costheta);
        angle_phi = TMath::ACos(cos_phi);
        if (sin_phi < 0){
          angle_phi =  -angle_phi;
        }
      
        // boosting in phi restframe                                                                                                          
        // the betas for the boost
        TVector3 p3_phi;
        p3_phi = Phi.Vect();
        p3_phi *= -1./Phi.E();
          
        // the boost matrix
        TLorentzRotation boost_phi(p3_phi);
        TLorentzVector p_phi_phi;
        p_phi_phi = boost_phi.VectorMultiplication(Phi);
          
        // the different momenta in the new frame
        TLorentzVector p_phi_Kplus = boost_phi.VectorMultiplication(k_p);
        TLorentzVector p_phi_JPsi = boost_phi.VectorMultiplication(Jpsi);
        TLorentzVector p_phi_Bs  = boost_phi.VectorMultiplication(Bs);
          
        // the 3-momenta
        TVector3 p3_phi_Kplus = p_phi_Kplus.Vect();
        TVector3 p3_phi_JPsi = p_phi_JPsi.Vect();
        TVector3 p3_phi_Bs = p_phi_Bs.Vect();
          
        angle_cospsi = -1 * p3_phi_Kplus.Unit() * p3_phi_JPsi.Unit();
        
        hCosPsi->Fill(angle_cospsi); 
        hCosTheta->Fill(angle_costheta);
        hAnglePhi->Fill(angle_phi);
          
        if (debug) cout << "CosPsi " << angle_cospsi << " CosTheta " << angle_costheta << " Angle Phi " << angle_phi << endl;
        nCandidates++;
            
        //MC Truth here
        if (use_gen){
        
          if (debug) cout << "-------------------- GEN info now" << endl;
        
          int BsIndex=-999, JpsiIndex=-999, PhiIndex=-999, KpIndex=-999, KmIndex=-999, MupIndex=-999, MumIndex=-999;

          //Search Bs
          for (int p=0; p< nGenP; ++p){
            BsIndex   = -999;
            JpsiIndex = -999;
            PhiIndex  = -999;
            KpIndex   = -999;
            KmIndex   = -999;
            MupIndex  = -999;
            MumIndex  = -999;

            if(abs(genId->at(p))!=531) continue;
            int tagMix = GetMixStatus( p );
            if(tagMix == 2) continue; // 2 -> Bs that are going to oscillate, 1 -> Bs have oscillated, 0 -> no oscillation

            const vector <int>& aD = allDaughters(p);
            for(int it:aD){
              if(abs(genId->at(it))==443) JpsiIndex=it;
              if(abs(genId->at(it))==333) PhiIndex=it;
            }
            if( JpsiIndex==-999 || PhiIndex==-999 ) continue;

            const vector <int>& aDjpsi = allDaughters(JpsiIndex);
            const vector <int>& aDphi = allDaughters(PhiIndex);

            //Jpsi muons
            for(int it:aDjpsi){
              if(genId->at(it)==13) MumIndex=it;
              if(genId->at(it)==-13) MupIndex=it;
            }
            if( MumIndex==-999 || MupIndex==-999 ) continue;

            //Phi Kaons
            for(int it:aDphi){
              if(genId->at(it)==321) KpIndex=it;
              if(genId->at(it)==-321) KmIndex=it;
            }
            if( KpIndex==-999 || KmIndex==-999 ) continue;

            BsIndex=p;
            if(genId->at(p)>0) TgenFlavour=1;
            else  TgenFlavour=-1;

            if(tagMix == 1) TgenFlavour*=-1;

            break;
          }

          if (debug) cout << "Found MC event, tot particles " << nGenP << " indexes:" << BsIndex << " " <<  JpsiIndex << " " <<  PhiIndex << " " << MupIndex  << " " << MumIndex << " " << KpIndex << " "<< KmIndex << endl;      

          if (BsIndex==-999) return false;
    
          TVector3 GenSV, GenPV, GenDiff, BsMomentumGen;
          GenSV.SetXYZ(genVx->at(JpsiIndex),genVy->at(JpsiIndex),genVz->at(JpsiIndex));//checked that the Jpsi SV is exactly the same as the Phi SV
          if (GetMixStatus(BsIndex)==0) {
            GenPV.SetXYZ(genVx->at(BsIndex),genVy->at(BsIndex),genVz->at(BsIndex));
          } else {
            const vector <int>& aM = allMothers(BsIndex);
            if( aM.size()>1)  cout<<"Bs with more than one mother"<<endl;
            GenPV.SetXYZ(genVx->at(aM[0]),genVy->at(aM[0]),genVz->at(aM[0]));
          }

          TLorentzVector BsGen, JpsiGen, PhiGen, MupGen, MumGen, KpGen, KmGen;
          
          BsGen.SetPtEtaPhiM(genPt->at(BsIndex), genEta->at(BsIndex), genPhi->at(BsIndex), genMass->at(BsIndex));
        
          BsMomentumGen=BsGen.Vect();
        
          GenDiff=GenSV-GenPV;
        
          float GenLxy3D=BsGen.M()/BsGen.P()*GenDiff.Dot(BsMomentumGen)/BsMomentumGen.Mag();
        
          GenDiff.SetZ(0.);//Lxy in 2D now
          BsMomentumGen.SetZ(0.);
        
          float GenLxy=BsGen.M()/BsGen.Pt()*GenDiff.Dot(BsMomentumGen)/BsMomentumGen.Mag();
          
          if (debug) cout << "Gen Lxy3D " << GenLxy3D << " Gen Lxy " << GenLxy << " BsIndex " << BsIndex << " JpsiIndex " << JpsiIndex << endl;
        
          JpsiGen.SetPtEtaPhiM(genPt->at(JpsiIndex), genEta->at(JpsiIndex), genPhi->at(JpsiIndex), genMass->at(JpsiIndex));
          PhiGen.SetPtEtaPhiM(genPt->at(PhiIndex), genEta->at(PhiIndex), genPhi->at(PhiIndex), genMass->at(PhiIndex));
        
          MupGen.SetPtEtaPhiM(genPt->at(MupIndex), genEta->at(MupIndex), genPhi->at(MupIndex), genMass->at(MupIndex));
          MumGen.SetPtEtaPhiM(genPt->at(MumIndex), genEta->at(MumIndex), genPhi->at(MumIndex), genMass->at(MumIndex));
        
          KpGen.SetPtEtaPhiM(genPt->at(KpIndex), genEta->at(KpIndex), genPhi->at(KpIndex), genMass->at(KpIndex));
          KmGen.SetPtEtaPhiM(genPt->at(KmIndex), genEta->at(KmIndex), genPhi->at(KmIndex), genMass->at(KmIndex));
        
          if (debug) cout << " Gen Masses " << BsGen.Mag() << " " << JpsiGen.Mag() << " " << PhiGen.Mag() << " " << MupGen.Mag() << " " << MumGen.Mag() << " " << KpGen.Mag() << " " << KmGen.Mag() << endl;
        
          //mc matching
          if (KpGen.DeltaR(k_p) < 0.005 && KmGen.DeltaR(k_m)< 0.005 && MumGen.DeltaR(mu_m) < 0.005 &&  MupGen.DeltaR(mu_p) < 0.005) Tmatched=1;
        
          if (debug) cout << "Distance GEN reco: Kp " << KpGen.DeltaR(k_p) << " Km " <<  KmGen.DeltaR(k_m) << " Mum " << MumGen.DeltaR(mu_m) << " Mup " << MupGen.DeltaR(mu_p) << endl;
        
          //Computing Gen angles
        
          double angle_costhetaGen;
          double angle_phiGen;
          double angle_cospsiGen;
        
          TVector3 p3_JPsiGen;
          p3_JPsiGen = JpsiGen.Vect();
          p3_JPsiGen *= -1./JpsiGen.E();
        
          // the boost matrix
          TLorentzRotation boost_jpsiGen(p3_JPsiGen);
          TLorentzVector p_JPsi_JPsiGen = boost_jpsiGen.VectorMultiplication(JpsiGen);
        
          // the different momenta in the new frame                                                                                                       
          TLorentzVector p_JPsi_muplusGen = boost_jpsiGen.VectorMultiplication(MupGen); 
          TLorentzVector p_JPsi_KplusGen = boost_jpsiGen.VectorMultiplication(KpGen); 
          TLorentzVector p_JPsi_phiGen = boost_jpsiGen.VectorMultiplication(PhiGen);      
        
          // the 3-momenta
          TVector3 p3_JPsi_muplusGen = p_JPsi_muplusGen.Vect();
          TVector3 p3_JPsi_KplusGen = p_JPsi_KplusGen.Vect();
          TVector3 p3_JPsi_phiGen = p_JPsi_phiGen.Vect();
        
          // coordinate syste`m
          TVector3 xg,yg,zg;
          xg = p3_JPsi_phiGen.Unit();
          yg = p3_JPsi_KplusGen.Unit() - (xg * (xg * p3_JPsi_KplusGen.Unit()));
          yg = yg.Unit();
          zg = xg.Cross(yg);
        
          // Transversity Basis
          angle_costhetaGen = p3_JPsi_muplusGen.Unit() * zg;
        
          double cos_phiGen = p3_JPsi_muplusGen.Unit() * xg / TMath::Sqrt(1 - angle_costhetaGen*angle_costhetaGen);
          double sin_phiGen = p3_JPsi_muplusGen.Unit() * yg / TMath::Sqrt(1 - angle_costhetaGen*angle_costhetaGen);
          angle_phiGen = TMath::ACos(cos_phiGen);
          if (sin_phiGen < 0){
            angle_phiGen =  -angle_phiGen;
          }
        
          // boosting in phi restframe                                                                                                          
          // the betas for the boost
          TVector3 p3_phiGen;
          p3_phiGen = PhiGen.Vect();
          p3_phiGen *= -1./PhiGen.E();
        
          // the boost matrix
          TLorentzRotation boost_phiGen(p3_phiGen);
          TLorentzVector p_phi_phiGen;
          p_phi_phiGen = boost_phiGen.VectorMultiplication(PhiGen);
        
          // the different momenta in the new frame
          TLorentzVector p_phi_KplusGen = boost_phiGen.VectorMultiplication(KpGen);
          TLorentzVector p_phi_JPsiGen = boost_phiGen.VectorMultiplication(JpsiGen);
          TLorentzVector p_phi_BsGen  = boost_phiGen.VectorMultiplication(BsGen);
        
          // the 3-momenta
          TVector3 p3_phi_KplusGen = p_phi_KplusGen.Vect();
          TVector3 p3_phi_JPsiGen = p_phi_JPsiGen.Vect();
          TVector3 p3_phi_BsGen = p_phi_BsGen.Vect();
        
          angle_cospsiGen = -1 * p3_phi_KplusGen.Unit() * p3_phi_JPsiGen.Unit();
        
          if (debug) cout << "CosPsi GEN " << angle_cospsiGen << " CosTheta GEN " << angle_costhetaGen << " Angle Phi GEN " << angle_phiGen << endl;
        
          //GEN variables for the Tree
          Tangle_cospsiGen=angle_cospsiGen;
          Tangle_costhetaGen=angle_costhetaGen;
          Tangle_phiGen=angle_phiGen;
        
          TLxyGen=GenLxy;
          TLxyzGen=GenLxy3D;
        
          TBsMassGen=BsGen.M();
          TBsPtGen=BsGen.Pt();
          TBsEtaGen=BsGen.Eta();
          TBsPhiGen=BsGen.Phi();
        
          TJpsiMassGen=JpsiGen.M();
          TJpsiPtGen=JpsiGen.Pt();
          TJpsiEtaGen=JpsiGen.Eta();
          TJpsiPhiGen=JpsiGen.Phi();
        
          TPhiMassGen=PhiGen.M();
          TPhiPtGen=PhiGen.Pt();
          TPhiEtaGen=PhiGen.Eta();
          TPhiPhiGen=PhiGen.Phi();
        
          TMupPtGen=MupGen.Pt();
          TMupPhiGen=MupGen.Phi();
          TMupEtaGen=MupGen.Eta();
          
          TMumPtGen=MumGen.Pt();
          TMumPhiGen=MumGen.Phi();
          TMumEtaGen=MumGen.Eta();
          
          TKpPtGen=KpGen.Pt();
          TKpPhiGen=KpGen.Phi();
          TKpEtaGen=KpGen.Eta();
          
          TKmPtGen=KmGen.Pt();
          TKmPhiGen=KmGen.Phi();
          TKmEtaGen=KmGen.Eta();
        }
      
        //Fill TTree here
          
        if (debug) cout << "Filling TTree now" << endl;

        Tangle_cospsi=angle_cospsi;
        Tangle_costheta=angle_costheta;
        Tangle_phi=angle_phi;
        
        TLxy=Lxy;
        TLxyPD=LxyPD;
        TLxyANBS=LxyANBS;
        TLxyAT=LxyAT;
        TLxyCNBS=LxyCNBS;
        TLxyCO=LxyCO;
        TLxyz=Lxy3D;
        TLxyErr=LxyErr;
        TLxyzErr=Lxy3DErr;
          
        TVprob=VProb;
        TVprobJpsi=VProbJpsi;
        TVprobPhi=VProbPhi;
          
        TBsMass=Bs.M();
        TBsPt=Bs.Pt();
        TBsEta=Bs.Eta();
        TBsPhi=Bs.Phi();
          
        TBsMassNoRefit=BsNoRefit.M();
        TJpsiMassNoRefit=JpsiNoRefit.M();
        TPhiMassNoRefit=PhiNoRefit.M();
          
        TJpsiMass=Jpsi.M();
        TJpsiPt=Jpsi.Pt();
        TJpsiEta=Jpsi.Eta();
        TJpsiPhi=Jpsi.Phi();
          
        TPhiMass=Phi.M();
        TPhiPt=Phi.Pt();
        TPhiEta=Phi.Eta();
        TPhiPhi=Phi.Phi();
          
        TMupPt=mu_p.Pt();
        TMupPhi=mu_p.Phi();
        TMupEta=mu_p.Eta();
        TMupHits=mupH;
          
        TMumPt=mu_m.Pt();
        TMumPhi=mu_m.Phi();
        TMumEta=mu_m.Eta();
        TMumHits=mumH;
          
        TKpPt=k_p.Pt();
        TKpPhi=k_p.Phi();
        TKpEta=k_p.Eta();
        TKpHits=kpH;
          
        TKmPt=k_m.Pt();
        TKmPhi=k_m.Phi();
        TKmEta=k_m.Eta();
        TKmHits=kmH;

        outTree->Fill();

      }

      if (TisBs==0){

        std::vector<int> mu_idx, k_idx, mu_idx_refit, k_idx_refit;
        
        TLorentzVector mu_p, mu_m, k_p, k_m, mu_p_norefit, mu_m_norefit, k_p_norefit, k_m_norefit, Bs, Jpsi, Phi, BsNoRefit, JpsiNoRefit, PhiNoRefit;
        TVector3 SVpos, PVpos;
          
        float VProb=-999, VProbJpsi=-999, VProbPhi=-999;
        float mupH=-999,mumH=-999,kpH=-999, kmH=-999;
          
        double angle_costheta;
        double angle_phi;
        double angle_cospsi;
        //    double AngleBsDecayLength;
          
        const vector<int>& tks = tracksFromSV( iSV );
        int n = tks.size();
          
        if (!n) continue;
          
        SVpos.SetXYZ(svtX->at( iSV ),svtY->at( iSV ),svtZ->at( iSV ));
          
        //filling some histos
        hBs->Fill( svtMass->at( iSV ) );
        VProb=TMath::Prob(svtChi2->at(iSV),svtNDOF->at(iSV));    
          
        TBsMassFromSV=svtMass->at( iSV );
      
        //Select the sub-secondary vertices (Jpsi and Phi)
        const vector<int>& sub = subVtxFromSV( iSV );
          
        if (sub.size() !=2) continue;   
      
        for (uint k=0; k<sub.size(); k++){
        
          if (svtType->at(sub[k]) == PDEnumString::svtJPsi){
            const vector<int>& mutks = tracksFromSV( sub[k] );
            for (uint t=0; t< mutks.size(); ++t) {mu_idx.push_back(mutks[t]); mu_idx_refit.push_back(vpIndex(mutks[t],sub[k]));}
          
            hJpsi->Fill(svtMass->at(sub[k]));
            TJpsiMassFromSV=svtMass->at(sub[k]);
          
            VProbJpsi=TMath::Prob(svtChi2->at(sub[k]),svtNDOF->at(sub[k]));
          
            if (debug) {
              cout << "J/psi mass (from vertex fit) " << svtMass->at(sub[k]) << " n tracks in SV " << mutks.size()  << " n tracks saved " << mu_idx.size() << " n tracks refitted " << mu_idx_refit.size() << endl;
              cout << "Muon indices norefit ";
              for (uint t=0; t< mutks.size(); ++t) cout << mutks[t] << " ";
              cout << endl;
              cout << "Muon indices refit ";
              for (uint t=0; t< mutks.size(); ++t) cout << mu_idx_refit[t] << " ";
              cout << endl;
              cout << endl;
            }   
          
          } else if (svtType->at(sub[k]) == PDEnumString::svtKx0){
            const vector<int>& ktks = tracksFromSV( sub[k] );
            for (uint t=0; t< ktks.size(); ++t) {
              k_idx.push_back(ktks[t]); 
              k_idx_refit.push_back(vpIndex(ktks[t],sub[k]));
            } 
          
            hPhi->Fill(svtMass->at(sub[k]));
            TPhiMassFromSV=svtMass->at(sub[k]);
            
            VProbPhi=TMath::Prob(svtChi2->at(sub[k]),svtNDOF->at(sub[k]));
          
            if (debug){
              cout << "Kstar mass (from vertex fit)" << svtMass->at(sub[k]) << " n tracks  in SV " << ktks.size()  << " n tracks saved " << k_idx.size() << " n tracks refitted " << k_idx_refit.size() << endl;
              cout << "Kaon indices norefit ";
              for (uint t=0; t< ktks.size(); ++t) cout << ktks[t] << " ";
              cout << endl;
              cout << "Kaon indices refit ";
              for (uint t=0; t< ktks.size(); ++t) cout << k_idx_refit[t] << " ";
              cout << endl;
              cout << endl;
              cout << "Size of the full track collection " << nTracks << " N refitted " << tvpPt->size() << endl;
              if ( k_idx_refit[0] >= int(tvpPt->size()) ||  k_idx_refit[1] >= int(tvpPt->size())){
                cout << "PROBLEM IN THE vpIndex FUNCTION" << endl;
                return false;
              }
            }       
          }
        }

        //analyze the muons and tracks now
          
        if (k_idx.size()!=2 || mu_idx.size()!=2 || k_idx_refit.size()!=2 || mu_idx_refit.size()!=2 ) continue;
          
        if (trkCharge->at(mu_idx[0]) >0){
          mu_p_norefit.SetPtEtaPhiM(trkPt->at(mu_idx[0]),trkEta->at(mu_idx[0]),trkPhi->at(mu_idx[0]),MUMASS);
          mu_m_norefit.SetPtEtaPhiM(trkPt->at(mu_idx[1]),trkEta->at(mu_idx[1]),trkPhi->at(mu_idx[1]),MUMASS);
          mupH=trkHitPattern->at(mu_idx[0]); mumH=trkHitPattern->at(mu_idx[1]);
        } else if (trkCharge->at(mu_idx[0]) < 0){
          mu_m_norefit.SetPtEtaPhiM(trkPt->at(mu_idx[0]),trkEta->at(mu_idx[0]),trkPhi->at(mu_idx[0]),MUMASS);
          mu_p_norefit.SetPtEtaPhiM(trkPt->at(mu_idx[1]),trkEta->at(mu_idx[1]),trkPhi->at(mu_idx[1]),MUMASS);
          mupH=trkHitPattern->at(mu_idx[1]); mumH=trkHitPattern->at(mu_idx[0]);
        }
      
        //now we have K* here, so the two tracks must have a different mass hypotesis
        TLorentzVector cand1, cand2;

        k_m_norefit.SetPtEtaPhiM(trkPt->at(k_idx[0]),trkEta->at(k_idx[0]),trkPhi->at(k_idx[0]),KMASS);
        k_p_norefit.SetPtEtaPhiM(trkPt->at(k_idx[1]),trkEta->at(k_idx[1]),trkPhi->at(k_idx[1]),PIMASS);
          
        cand1= k_p_norefit+ k_m_norefit;

          
        k_p_norefit.SetPtEtaPhiM(trkPt->at(k_idx[0]),trkEta->at(k_idx[0]),trkPhi->at(k_idx[0]),PIMASS);
        k_m_norefit.SetPtEtaPhiM(trkPt->at(k_idx[1]),trkEta->at(k_idx[1]),trkPhi->at(k_idx[1]),KMASS);
          
        cand2= k_p_norefit+ k_m_norefit;

        if (abs(cand1.M()-svtMass->at( iSV )) <= abs(cand2.M()-svtMass->at( iSV ))){
          k_m_norefit.SetPtEtaPhiM(trkPt->at(k_idx[0]),trkEta->at(k_idx[0]),trkPhi->at(k_idx[0]),KMASS);
          k_p_norefit.SetPtEtaPhiM(trkPt->at(k_idx[1]),trkEta->at(k_idx[1]),trkPhi->at(k_idx[1]),PIMASS);
          kpH=trkHitPattern->at(k_idx[0]); kmH=trkHitPattern->at(mu_idx[1]);  
        }
        else {
          k_p_norefit.SetPtEtaPhiM(trkPt->at(k_idx[0]),trkEta->at(k_idx[0]),trkPhi->at(k_idx[0]),PIMASS);
          k_m_norefit.SetPtEtaPhiM(trkPt->at(k_idx[1]),trkEta->at(k_idx[1]),trkPhi->at(k_idx[1]),KMASS);
          kpH=trkHitPattern->at(k_idx[1]); kmH=trkHitPattern->at(mu_idx[0]);    
        }
      
        //adding refitted tracks
          
        if (trkCharge->at(mu_idx[0]) >0){
          mu_p.SetPtEtaPhiM(tvpPt->at(mu_idx_refit[0]),tvpEta->at(mu_idx_refit[0]),tvpPhi->at(mu_idx_refit[0]),MUMASS);
          mu_m.SetPtEtaPhiM(tvpPt->at(mu_idx_refit[1]),tvpEta->at(mu_idx_refit[1]),tvpPhi->at(mu_idx_refit[1]),MUMASS);
        } else if (trkCharge->at(mu_idx[0]) < 0){
          mu_m.SetPtEtaPhiM(tvpPt->at(mu_idx_refit[0]),tvpEta->at(mu_idx_refit[0]),tvpPhi->at(mu_idx_refit[0]),MUMASS);
          mu_p.SetPtEtaPhiM(tvpPt->at(mu_idx_refit[1]),tvpEta->at(mu_idx_refit[1]),tvpPhi->at(mu_idx_refit[1]),MUMASS);
        }

        //select the mass hypotesis closer to PDG value

        k_m.SetPtEtaPhiM(tvpPt->at(k_idx_refit[0]),tvpEta->at(k_idx_refit[0]),tvpPhi->at(k_idx_refit[0]),KMASS);
        k_p.SetPtEtaPhiM(tvpPt->at(k_idx_refit[1]),tvpEta->at(k_idx_refit[1]),tvpPhi->at(k_idx_refit[1]),PIMASS);
        cand1= k_p + k_m;


        k_p.SetPtEtaPhiM(tvpPt->at(k_idx_refit[0]),tvpEta->at(k_idx_refit[0]),tvpPhi->at(k_idx_refit[0]),PIMASS);
        k_m.SetPtEtaPhiM(tvpPt->at(k_idx_refit[1]),tvpEta->at(k_idx_refit[1]),tvpPhi->at(k_idx_refit[1]),KMASS);
        cand2= k_p + k_m;

        if (abs(cand1.M()-svtMass->at( iSV )) <= abs(cand2.M()-svtMass->at( iSV ))){
          k_m.SetPtEtaPhiM(tvpPt->at(k_idx_refit[0]),tvpEta->at(k_idx_refit[0]),tvpPhi->at(k_idx_refit[0]),KMASS);
          k_p.SetPtEtaPhiM(tvpPt->at(k_idx_refit[1]),tvpEta->at(k_idx_refit[1]),tvpPhi->at(k_idx_refit[1]),PIMASS);
        }
        else {
          k_p.SetPtEtaPhiM(tvpPt->at(k_idx_refit[0]),tvpEta->at(k_idx_refit[0]),tvpPhi->at(k_idx_refit[0]),PIMASS);
          k_m.SetPtEtaPhiM(tvpPt->at(k_idx_refit[1]),tvpEta->at(k_idx_refit[1]),tvpPhi->at(k_idx_refit[1]),KMASS);
        }
      
        //float pt_mu=trkPt->at(mu_idx[0]);

        if (debug) {
          cout << "+++++++++++++ Building Candidates ++++++++++++++" << endl;
          cout << "From refitted tracks: J/psi mass " << (mu_p+mu_m).Mag() << " Kstar mass " << (k_p+k_m).Mag() << endl;
          cout << endl;
          cout << "From native tracks: J/psi mass " << (mu_p_norefit+mu_m_norefit).Mag() << " Star mass " << (k_p_norefit+k_m_norefit).Mag() << endl;
          cout << endl;
        }
      
        Jpsi=mu_p + mu_m;
        
        Phi=k_p + k_m;
          
        Bs = Jpsi + Phi;
          
        JpsiNoRefit=mu_p_norefit+mu_m_norefit;
          
        PhiNoRefit=k_p_norefit + k_m_norefit;
        
        BsNoRefit = JpsiNoRefit + PhiNoRefit;
          
        int pvIndex=-999;
        //selecting the PV and computing Lxy and Error
          
        FindPV(SVpos, PVpos, Bs, pvIndex);
        if (pvIndex < 0){
          if (debug) cout << "NO PV FOUND FOR THIS EVENT" << endl;
          return false; //return if no PV is identified
        }
            
        inizializeOsMuonTagVars();
        setVtxOsMuonTag(iSV, pvIndex); //set PV and SVT for tagging class
        bool istagged = makeOsMuonTagging();

        Ttag = int(getOsMuonTag()); //get tag decision -1, +1 or 0 

        if(istagged){
          Tmtag = getOsMuonTagMistagProbRaw();
          TmtagCal = getOsMuonTagMistagProbCalProcess();
          TmtagCalBs = getOsMuonTagMistagProbCalProcessBuBs();
        }

        float Lxy, Lxy3D;
          
        TVector3 BsMomentum=Bs.Vect();//3Momentum of Bs
        TVector3 PVSVdist=SVpos-PVpos;
          
        Lxy3D=B0MASS*PVSVdist.Dot(BsMomentum)/BsMomentum.Mag2();//(PVpos-SVpos).Mag();
        
        BsMomentum.SetZ(0.);
        PVSVdist.SetZ(0.);
          
        Lxy=B0MASS*PVSVdist.Dot(BsMomentum)/BsMomentum.Mag2(); //Transverse decay lenght after setting to 0 the Z components of the vectors above
          
        if (debug) cout << "Lxy " << Lxy << " Lxy 3D " << Lxy3D << endl;

        TMatrixF covSV(3,3);
        float covSVArray[]={svtSxx->at(iSV),svtSxy->at(iSV),svtSxz->at(iSV),svtSxy->at(iSV),svtSyy->at(iSV),svtSyz->at(iSV), svtSxz->at(iSV),svtSyz->at(iSV),svtSzz->at(iSV)};
        covSV.SetMatrixArray(covSVArray);
          
        TMatrixF covPV(3,3);
        float covPVArray[]={pvtSxx->at(pvIndex),pvtSxy->at(pvIndex),pvtSxz->at(pvIndex),pvtSxy->at(pvIndex),pvtSyy->at(pvIndex),pvtSyz->at(pvIndex), pvtSxz->at(pvIndex),pvtSyz->at(pvIndex),pvtSzz->at(pvIndex)};
        covPV.SetMatrixArray(covPVArray);
          
        TMatrixF covTot= covSV+covPV;
          
        float distArray[]={float((SVpos-PVpos).X()),float((SVpos-PVpos).Y()),float((SVpos-PVpos).Z())};
        TVectorF diff(3,distArray);

        float distArray2D[]={float((SVpos-PVpos).X()),float((SVpos-PVpos).Y()),0.};
        TVectorF diff2D(3,distArray2D);
        
            
        float LxyErr=-999,Lxy3DErr=-999;
             
        if (diff2D.Norm2Sqr()==0 || diff.Norm2Sqr()==0) continue; //if the secondary vertex is exactly the same as PV continue
          
        LxyErr=B0MASS/Bs.Pt()*sqrt(covTot.Similarity(diff2D))/sqrt(diff2D.Norm2Sqr()); 
        Lxy3DErr=B0MASS/Bs.P()*sqrt(covTot.Similarity(diff))/sqrt(diff.Norm2Sqr()); 
        
        
        // the betas for the boost
        TVector3 p3_JPsi;
        p3_JPsi = Jpsi.Vect();
        p3_JPsi *= -1./Jpsi.E();

        // the boost matrix
        TLorentzRotation boost_jpsi(p3_JPsi);
        TLorentzVector p_JPsi_JPsi = boost_jpsi.VectorMultiplication(Jpsi);
          
        // the different momenta in the new frame                                                                                                       
        TLorentzVector p_JPsi_muplus = boost_jpsi.VectorMultiplication(mu_p); 
        TLorentzVector p_JPsi_Kplus = boost_jpsi.VectorMultiplication(k_p); 
        TLorentzVector p_JPsi_phi = boost_jpsi.VectorMultiplication(Phi);                                                                       
          
        // the 3-momenta
        TVector3 p3_JPsi_muplus = p_JPsi_muplus.Vect();
        TVector3 p3_JPsi_Kplus = p_JPsi_Kplus.Vect();
        TVector3 p3_JPsi_phi = p_JPsi_phi.Vect();
          
        // coordinate system
        TVector3 x,y,z;
        x = p3_JPsi_phi.Unit();
        y = p3_JPsi_Kplus.Unit() - (x * (x * p3_JPsi_Kplus.Unit()));
        y = y.Unit();
        z = x.Cross(y);
          
        // Transversity Basis
        angle_costheta = p3_JPsi_muplus.Unit() * z;
          
        double cos_phi = p3_JPsi_muplus.Unit() * x / TMath::Sqrt(1 - angle_costheta*angle_costheta);
        double sin_phi = p3_JPsi_muplus.Unit() * y / TMath::Sqrt(1 - angle_costheta*angle_costheta);
        angle_phi = TMath::ACos(cos_phi);
        if (sin_phi < 0){
          angle_phi =  -angle_phi;
        }
          
        // boosting in phi restframe                                                                                                          
        // the betas for the boost
        TVector3 p3_phi;
        p3_phi = Phi.Vect();
        p3_phi *= -1./Phi.E();
          
        // the boost matrix
        TLorentzRotation boost_phi(p3_phi);
        TLorentzVector p_phi_phi;
        p_phi_phi = boost_phi.VectorMultiplication(Phi);
          
        // the different momenta in the new frame
        TLorentzVector p_phi_Kplus = boost_phi.VectorMultiplication(k_p);
        TLorentzVector p_phi_JPsi = boost_phi.VectorMultiplication(Jpsi);
        TLorentzVector p_phi_Bs  = boost_phi.VectorMultiplication(Bs);
          
        // the 3-momenta
        TVector3 p3_phi_Kplus = p_phi_Kplus.Vect();
        TVector3 p3_phi_JPsi = p_phi_JPsi.Vect();
        TVector3 p3_phi_Bs = p_phi_Bs.Vect();
          
        angle_cospsi = -1 * p3_phi_Kplus.Unit() * p3_phi_JPsi.Unit();
        
        hCosPsi->Fill(angle_cospsi); 
        hCosTheta->Fill(angle_costheta);
        hAnglePhi->Fill(angle_phi);
        
        if (debug) cout << "CosPsi " << angle_cospsi << " CosTheta " << angle_costheta << " Angle Phi " << angle_phi << endl;
        nCandidates++;
    
        //MC Truth here
        if (use_gen){

          if (debug) cout << "-------------------- GEN info now" << endl;

          int BsIndex=-999, JpsiIndex=-999, PhiIndex=-999, KpIndex=-999, KmIndex=-999, MupIndex=-999, MumIndex=-999;
      
          //Search for B0
          for (int p=0; p< nGenP; ++p){
            BsIndex=-999;
            JpsiIndex=-999;
            PhiIndex=-999;
            KpIndex=-999;
            KmIndex=-999;
            MupIndex=-999;
            MumIndex=-999;

            if(abs(genId->at(p))!=511) continue;
            int tagMix = GetMixStatus( p );
            if(tagMix == 2) continue; // 2 -> B that are going to oscillate, 1 -> B have oscillated, 0 -> no oscillation

            const vector <int>& aD = allDaughters(p);
            for(int it:aD){
              if(abs(genId->at(it))==443) JpsiIndex=it;
              if(abs(genId->at(it))==313) PhiIndex=it;
            }
            if( JpsiIndex==-999 || PhiIndex==-999 ) continue;

            const vector <int>& aDjpsi = allDaughters(JpsiIndex);
            // const vector <int>& aDK0 = allDaughters(PhiIndex);

            //Jpsi muons
            for(int it:aDjpsi){
              if(genId->at(it)==13) MumIndex=it;
              if(genId->at(it)==-13) MupIndex=it;
            }
            if( MumIndex==-999 || MupIndex==-999 ) continue;

            // //K*0
            // for(int it:aDK0){
            //   if(abs(genId->at(it))==321) KpIndex=it;
            //   if(abs(genId->at(it))==211) KmIndex=it;
            // }
            // if( KpIndex==-999 || KmIndex==-999 ) continue;

            BsIndex=p;
            if(genId->at(p)>0) TgenFlavour=1;
            else  TgenFlavour=-1;

            if(tagMix == 1) TgenFlavour*=-1;

            break;
          }

          if (BsIndex==-999) return false;

          if (debug) cout << "Found MC event, tot particles " << nGenP << " indexes:" << BsIndex << " " <<  JpsiIndex << " " <<  PhiIndex << " " << MupIndex  << " " << MumIndex << " " << KpIndex << " "<< KmIndex << endl;   
            
          TVector3 GenSV, GenPV, GenDiff, BsMomentumGen;
          GenSV.SetXYZ(genVx->at(JpsiIndex),genVy->at(JpsiIndex),genVz->at(JpsiIndex));//checked that the Jpsi SV is exactly the same as the Phi SV
          if (GetMixStatus(BsIndex)==0) {
            GenPV.SetXYZ(genVx->at(BsIndex),genVy->at(BsIndex),genVz->at(BsIndex));
          } else {
            const vector <int>& aM = allMothers(BsIndex);
            if( aM.size()>1)  cout<<"Bs with more than one mother"<<endl;
            GenPV.SetXYZ(genVx->at(aM[0]),genVy->at(aM[0]),genVz->at(aM[0]));
          }
        
          TLorentzVector BsGen, JpsiGen, PhiGen, MupGen, MumGen, KpGen, KmGen;
          
          BsGen.SetPtEtaPhiM(genPt->at(BsIndex), genEta->at(BsIndex), genPhi->at(BsIndex), genMass->at(BsIndex));
        
          BsMomentumGen=BsGen.Vect();
        
          GenDiff=GenSV-GenPV;
        
          float GenLxy3D=BsGen.M()/BsGen.P()*GenDiff.Dot(BsMomentumGen)/BsMomentumGen.Mag();
        
          GenDiff.SetZ(0.);//Lxy in 2D now
          BsMomentumGen.SetZ(0.);
        
          float GenLxy=BsGen.M()/BsGen.Pt()*GenDiff.Dot(BsMomentumGen)/BsMomentumGen.Mag();
        
          if (debug) cout << "Gen Lxy3D " << GenLxy3D << " Gen Lxy " << GenLxy << " BsIndex " << BsIndex << " JpsiIndex " << JpsiIndex << endl;
        
          JpsiGen.SetPtEtaPhiM(genPt->at(JpsiIndex), genEta->at(JpsiIndex), genPhi->at(JpsiIndex), genMass->at(JpsiIndex));
          PhiGen.SetPtEtaPhiM(genPt->at(PhiIndex), genEta->at(PhiIndex), genPhi->at(PhiIndex), genMass->at(PhiIndex));
        
          MupGen.SetPtEtaPhiM(genPt->at(MupIndex), genEta->at(MupIndex), genPhi->at(MupIndex), genMass->at(MupIndex));
          MumGen.SetPtEtaPhiM(genPt->at(MumIndex), genEta->at(MumIndex), genPhi->at(MumIndex), genMass->at(MumIndex));
        
          KpGen.SetPtEtaPhiM(k_p.Pt(), k_p.Eta(), k_p.Phi(), k_p.M());//GEN tracks are set to the reco ones since they are not included in the collection
          KmGen.SetPtEtaPhiM(k_m.Pt(), k_m.Eta(), k_m.Phi(), k_m.M());
        
          if (debug) cout << " Gen Masses " << BsGen.Mag() << " " << JpsiGen.Mag() << " " << PhiGen.Mag() << " " << MupGen.Mag() << " " << MumGen.Mag() << " " << KpGen.Mag() << " " << KmGen.Mag() << endl;
        
          //mc matching
          if (PhiGen.DeltaR(k_p+k_m) < 0.005 && MumGen.DeltaR(mu_m) < 0.005 &&  MupGen.DeltaR(mu_p) < 0.005) Tmatched=1;
        
          if (debug) cout << "Distance GEN reco: Kstar " << PhiGen.DeltaR(k_p+k_m)  << " Mum " << MumGen.DeltaR(mu_m) << " Mup " << MupGen.DeltaR(mu_p) << endl;
        
          //Computing Gen angles
        
          double angle_costhetaGen;
          double angle_phiGen;
          double angle_cospsiGen;
        
          TVector3 p3_JPsiGen;
          p3_JPsiGen = JpsiGen.Vect();
          p3_JPsiGen *= -1./JpsiGen.E();
        
          // the boost matrix
          TLorentzRotation boost_jpsiGen(p3_JPsiGen);
          TLorentzVector p_JPsi_JPsiGen = boost_jpsiGen.VectorMultiplication(JpsiGen);
        
          // the different momenta in the new frame                                                                                                       
          TLorentzVector p_JPsi_muplusGen = boost_jpsiGen.VectorMultiplication(MupGen); 
          TLorentzVector p_JPsi_KplusGen = boost_jpsiGen.VectorMultiplication(KpGen); 
          TLorentzVector p_JPsi_phiGen = boost_jpsiGen.VectorMultiplication(PhiGen);      
        
          // the 3-momenta
          TVector3 p3_JPsi_muplusGen = p_JPsi_muplusGen.Vect();
          TVector3 p3_JPsi_KplusGen = p_JPsi_KplusGen.Vect();
          TVector3 p3_JPsi_phiGen = p_JPsi_phiGen.Vect();
        
          // coordinate syste`m
          TVector3 xg,yg,zg;
          xg = p3_JPsi_phiGen.Unit();
          yg = p3_JPsi_KplusGen.Unit() - (xg * (xg * p3_JPsi_KplusGen.Unit()));
          yg = yg.Unit();
          zg = xg.Cross(yg);
        
          // Transversity Basis
          angle_costhetaGen = p3_JPsi_muplusGen.Unit() * zg;
        
          double cos_phiGen = p3_JPsi_muplusGen.Unit() * xg / TMath::Sqrt(1 - angle_costhetaGen*angle_costhetaGen);
          double sin_phiGen = p3_JPsi_muplusGen.Unit() * yg / TMath::Sqrt(1 - angle_costhetaGen*angle_costhetaGen);
          angle_phiGen = TMath::ACos(cos_phiGen);
          if (sin_phiGen < 0){
            angle_phiGen =  -angle_phiGen;
          }
        
          // boosting in phi restframe                                                                                                          
          // the betas for the boost
          TVector3 p3_phiGen;
          p3_phiGen = PhiGen.Vect();
          p3_phiGen *= -1./PhiGen.E();
        
          // the boost matrix
          TLorentzRotation boost_phiGen(p3_phiGen);
          TLorentzVector p_phi_phiGen;
          p_phi_phiGen = boost_phiGen.VectorMultiplication(PhiGen);
        
          // the different momenta in the new frame
          TLorentzVector p_phi_KplusGen = boost_phiGen.VectorMultiplication(KpGen);
          TLorentzVector p_phi_JPsiGen = boost_phiGen.VectorMultiplication(JpsiGen);
          TLorentzVector p_phi_BsGen  = boost_phiGen.VectorMultiplication(BsGen);
        
          // the 3-momenta
          TVector3 p3_phi_KplusGen = p_phi_KplusGen.Vect();
          TVector3 p3_phi_JPsiGen = p_phi_JPsiGen.Vect();
          TVector3 p3_phi_BsGen = p_phi_BsGen.Vect();
        
          angle_cospsiGen = -1 * p3_phi_KplusGen.Unit() * p3_phi_JPsiGen.Unit();
        
          if (debug) cout << "CosPsi GEN " << angle_cospsiGen << " CosTheta GEN " << angle_costhetaGen << " Angle Phi GEN " << angle_phiGen << endl;
        
          //GEN variables for the Tree
          Tangle_cospsiGen=angle_cospsiGen;
          Tangle_costhetaGen=angle_costhetaGen;
          Tangle_phiGen=angle_phiGen;
        
          TLxyGen=GenLxy;
          TLxyzGen=GenLxy3D;
        
          TBsMassGen=BsGen.M();
          TBsPtGen=BsGen.Pt();
          TBsEtaGen=BsGen.Eta();
          TBsPhiGen=BsGen.Phi();
        
          TJpsiMassGen=JpsiGen.M();
          TJpsiPtGen=JpsiGen.Pt();
          TJpsiEtaGen=JpsiGen.Eta();
          TJpsiPhiGen=JpsiGen.Phi();
        
          TPhiMassGen=PhiGen.M();
          TPhiPtGen=PhiGen.Pt();
          TPhiEtaGen=PhiGen.Eta();
          TPhiPhiGen=PhiGen.Phi();
        
          TMupPtGen=MupGen.Pt();
          TMupPhiGen=MupGen.Phi();
          TMupEtaGen=MupGen.Eta();
          
          TMumPtGen=MumGen.Pt();
          TMumPhiGen=MumGen.Phi();
          TMumEtaGen=MumGen.Eta();
        
          TKpPtGen=KpGen.Pt();
          TKpPhiGen=KpGen.Phi();
          TKpEtaGen=KpGen.Eta();
        
          TKmPtGen=KmGen.Pt();
          TKmPhiGen=KmGen.Phi();
          TKmEtaGen=KmGen.Eta();
        }
      
        //Fill TTree here
          
        Tangle_cospsi=angle_cospsi;
        Tangle_costheta=angle_costheta;
        Tangle_phi=angle_phi;
          
        TLxy=Lxy;
        TLxyz=Lxy3D;
        TLxyErr=LxyErr;
        TLxyzErr=Lxy3DErr;
          
        TVprob=VProb;
        TVprobJpsi=VProbJpsi;
        TVprobPhi=VProbPhi;
          
        TBsMass=Bs.M();
        TBsPt=Bs.Pt();
        TBsEta=Bs.Eta();
        TBsPhi=Bs.Phi();
          
        TBsMassNoRefit=BsNoRefit.M();
        TJpsiMassNoRefit=JpsiNoRefit.M();
        TPhiMassNoRefit=PhiNoRefit.M();
          
        TJpsiMass=Jpsi.M();
        TJpsiPt=Jpsi.Pt();
        TJpsiEta=Jpsi.Eta();
        TJpsiPhi=Jpsi.Phi();
          
        TPhiMass=Phi.M();
        TPhiPt=Phi.Pt();
        TPhiEta=Phi.Eta();
        TPhiPhi=Phi.Phi();
          
        TMupPt=mu_p.Pt();
        TMupPhi=mu_p.Phi();
        TMupEta=mu_p.Eta();
        TMupHits=mupH;
          
        TMumPt=mu_m.Pt();
        TMumPhi=mu_m.Phi();
        TMumEta=mu_m.Eta();
        TMumHits=mumH;
          
        TKpPt=k_p.Pt();
        TKpPhi=k_p.Phi();
        TKpEta=k_p.Eta();
        TKpHits=kpH;
          
        TKmPt=k_m.Pt();
        TKmPhi=k_m.Phi();
        TKmEta=k_m.Eta();
        TKmHits=kmH;
          
        outTree->Fill();

      }
    }
  }


  if (genOnly){//simply copyed and pasted the GEN level analysis, it works only with Bs!!!

    if (debug) cout << "-------------------- ONLY GEN info will be analyzed, THIS MESSAGE SHOULD NEVER APPEAR RUNNING ON DATA" << endl;

    int BsIndex=-999, JpsiIndex=-999, PhiIndex=-999, KpIndex=-999, KmIndex=-999, MupIndex=-999, MumIndex=-999;

    /*int nBs=0;
    for (int p=0; p< nGenP; ++p){
      if (abs(genId->at(p))==531 && (abs(genId->at(genMother->at(p)))==533 || abs(genId->at(genMother->at(p)))==21 || abs(genId->at(genMother->at(p)))<10)) nBs++;
    }  

    if (nBs>1){if (debug) cout << nBs << " GenBs, returning " << endl; return false;} // consider only MC events with a single Bs
    */  //Why do remove events with two Bs?
    //Search Bs
    for (int p=0; p< nGenP; ++p){
      BsIndex   = -999;
      JpsiIndex = -999;
      PhiIndex  = -999;
      KpIndex   = -999;
      KmIndex   = -999;
      MupIndex  = -999;
      MumIndex  = -999;

      if(abs(genId->at(p))!=531) continue;
      int tagMix = GetMixStatus( p );
      if(tagMix == 2) continue; // 2 -> Bs that are going to oscillate, 1 -> Bs have oscillated, 0 -> no oscillation

      const vector <int>& aD = allDaughters(p);
      for(int it:aD){
        if(abs(genId->at(it))==443) JpsiIndex=it;
        if(abs(genId->at(it))==333) PhiIndex=it;
      }
      if( JpsiIndex==-999 || PhiIndex==-999 ) continue;

      const vector <int>& aDjpsi = allDaughters(JpsiIndex);
      const vector <int>& aDphi = allDaughters(PhiIndex);

      //Jpsi muons
      for(int it:aDjpsi){
        if(genId->at(it)==13) MumIndex=it;
        if(genId->at(it)==-13) MupIndex=it;
      }
      if( MumIndex==-999 || MupIndex==-999 ) continue;

      //Phi Kaons
      for(int it:aDphi){
        if(genId->at(it)==321) KpIndex=it;
        if(genId->at(it)==-321) KmIndex=it;
      }
      if( KpIndex==-999 || KmIndex==-999 ) continue;

      BsIndex=p;
      if(genId->at(p)>0) TgenFlavour=1;
      else  TgenFlavour=-1;

      if(tagMix == 1) TgenFlavour*=-1;

      break;
    }

    if (BsIndex==-999) return false;

    if (debug) cout << "Found MC event, tot particles " << nGenP << " indexes:" << BsIndex << " " <<  JpsiIndex << " " <<  PhiIndex << " " << MupIndex  << " " << MumIndex << " " << KpIndex << " "<< KmIndex << endl;   

    TVector3 GenSV, GenPV, GenDiff, BsMomentumGen;
    GenSV.SetXYZ(genVx->at(JpsiIndex),genVy->at(JpsiIndex),genVz->at(JpsiIndex));//checked that the Jpsi SV is exactly the same as the Phi SV
    if (GetMixStatus(BsIndex)==0) {
      GenPV.SetXYZ(genVx->at(BsIndex),genVy->at(BsIndex),genVz->at(BsIndex));
    } else {
      const vector <int>& aM = allMothers(BsIndex);
      if( aM.size()>1)  cout<<"Bs with more than one mother"<<endl;
      GenPV.SetXYZ(genVx->at(aM[0]),genVy->at(aM[0]),genVz->at(aM[0]));
    }
 
    TLorentzVector BsGen, JpsiGen, PhiGen, MupGen, MumGen, KpGen, KmGen;

    BsGen.SetPtEtaPhiM(genPt->at(BsIndex), genEta->at(BsIndex), genPhi->at(BsIndex), genMass->at(BsIndex));

    BsMomentumGen=BsGen.Vect();

    GenDiff=GenSV-GenPV;

    float GenLxy3D=BsGen.M()/BsGen.P()*GenDiff.Mag();//Dot(BsMomentumGen)/BsMomentumGen.Mag();

    GenDiff.SetZ(0.);//Lxy in 2D now
    BsMomentumGen.SetZ(0.);

    float GenLxy=BsGen.M()/BsGen.Pt()*GenDiff.Mag();//Dot(BsMomentumGen)/BsMomentumGen.Mag();

    //if (GenLxy<0) {cout << " Lxy " << GenLxy << " Bs eta " << BsGen.Eta() << " Bs pt " << BsGen.Pt() << " Bs Mass " << BsGen.M() << " "  << BsMomentumGen.X() << " " << BsMomentumGen.Y() << " " << BsMomentumGen.Z() << " SV " << GenSV.X() << " " << GenSV.Y() << " " << GenSV.Z() << " " << " PV " << GenPV.X() << " " << GenPV.Y() << " " << GenPV.Z() << " " << " Mother " << genId->at(genMother->at(BsIndex))<< " Index " << BsIndex <<endl;}
    //if (GenLxy>0) return false;

    if (debug) cout << "Gen Lxy3D " << GenLxy3D << " Gen Lxy " << GenLxy << " BsIndex " << BsIndex << " JpsiIndex " << JpsiIndex << endl;
    
    JpsiGen.SetPtEtaPhiM(genPt->at(JpsiIndex), genEta->at(JpsiIndex), genPhi->at(JpsiIndex), genMass->at(JpsiIndex));
    PhiGen.SetPtEtaPhiM(genPt->at(PhiIndex), genEta->at(PhiIndex), genPhi->at(PhiIndex), genMass->at(PhiIndex));
    
    MupGen.SetPtEtaPhiM(genPt->at(MupIndex), genEta->at(MupIndex), genPhi->at(MupIndex), genMass->at(MupIndex));
    MumGen.SetPtEtaPhiM(genPt->at(MumIndex), genEta->at(MumIndex), genPhi->at(MumIndex), genMass->at(MumIndex));
    
    KpGen.SetPtEtaPhiM(genPt->at(KpIndex), genEta->at(KpIndex), genPhi->at(KpIndex), genMass->at(KpIndex));
    KmGen.SetPtEtaPhiM(genPt->at(KmIndex), genEta->at(KmIndex), genPhi->at(KmIndex), genMass->at(KmIndex));
  
    if (debug) cout << " Gen Masses " << BsGen.Mag() << " " << JpsiGen.Mag() << " " << PhiGen.Mag() << " " << MupGen.Mag() << " " << MumGen.Mag() << " " << KpGen.Mag() << " " << KmGen.Mag() << endl;

    //Computing Gen angles

    double angle_costhetaGen;
    double angle_phiGen;
    double angle_cospsiGen;

    TVector3 p3_JPsiGen;
    p3_JPsiGen = JpsiGen.Vect();
    p3_JPsiGen *= -1./JpsiGen.E();

    // the boost matrix
    TLorentzRotation boost_jpsiGen(p3_JPsiGen);
    TLorentzVector p_JPsi_JPsiGen = boost_jpsiGen.VectorMultiplication(JpsiGen);

    // the different momenta in the new frame                                                                                                       
    TLorentzVector p_JPsi_muplusGen = boost_jpsiGen.VectorMultiplication(MupGen); 
    TLorentzVector p_JPsi_KplusGen = boost_jpsiGen.VectorMultiplication(KpGen); 
    TLorentzVector p_JPsi_phiGen = boost_jpsiGen.VectorMultiplication(PhiGen);      

    // the 3-momenta
    TVector3 p3_JPsi_muplusGen = p_JPsi_muplusGen.Vect();
    TVector3 p3_JPsi_KplusGen = p_JPsi_KplusGen.Vect();
    TVector3 p3_JPsi_phiGen = p_JPsi_phiGen.Vect();

    // coordinate syste`m
    TVector3 xg,yg,zg;
    xg = p3_JPsi_phiGen.Unit();
    yg = p3_JPsi_KplusGen.Unit() - (xg * (xg * p3_JPsi_KplusGen.Unit()));
    yg = yg.Unit();
    zg = xg.Cross(yg);

    // Transversity Basis
    angle_costhetaGen = p3_JPsi_muplusGen.Unit() * zg;

    double cos_phiGen = p3_JPsi_muplusGen.Unit() * xg / TMath::Sqrt(1 - angle_costhetaGen*angle_costhetaGen);
    double sin_phiGen = p3_JPsi_muplusGen.Unit() * yg / TMath::Sqrt(1 - angle_costhetaGen*angle_costhetaGen);
    angle_phiGen = TMath::ACos(cos_phiGen);
    if (sin_phiGen < 0){
      angle_phiGen =  -angle_phiGen;
    }

    // boosting in phi restframe                                                                                                          
    // the betas for the boost
    TVector3 p3_phiGen;
    p3_phiGen = PhiGen.Vect();
    p3_phiGen *= -1./PhiGen.E();

    // the boost matrix
    TLorentzRotation boost_phiGen(p3_phiGen);
    TLorentzVector p_phi_phiGen;
    p_phi_phiGen = boost_phiGen.VectorMultiplication(PhiGen);

    // the different momenta in the new frame
    TLorentzVector p_phi_KplusGen = boost_phiGen.VectorMultiplication(KpGen);
    TLorentzVector p_phi_JPsiGen = boost_phiGen.VectorMultiplication(JpsiGen);
    TLorentzVector p_phi_BsGen  = boost_phiGen.VectorMultiplication(BsGen);

    // the 3-momenta
    TVector3 p3_phi_KplusGen = p_phi_KplusGen.Vect();
    TVector3 p3_phi_JPsiGen = p_phi_JPsiGen.Vect();
    TVector3 p3_phi_BsGen = p_phi_BsGen.Vect();

    angle_cospsiGen = -1 * p3_phi_KplusGen.Unit() * p3_phi_JPsiGen.Unit();

    if (debug) cout << "CosPsi GEN " << angle_cospsiGen << " CosTheta GEN " << angle_costhetaGen << " Angle Phi GEN " << angle_phiGen << endl;

    //GEN variables for the Tree
    Tangle_cospsiGen=angle_cospsiGen;
    Tangle_costhetaGen=angle_costhetaGen;
    Tangle_phiGen=angle_phiGen;

    TLxyGen=GenLxy;
    TLxyzGen=GenLxy3D;

    TBsMassGen=BsGen.M();
    TBsPtGen=BsGen.Pt();
    TBsEtaGen=BsGen.Eta();
    TBsPhiGen=BsGen.Phi();

    TJpsiMassGen=JpsiGen.M();
    TJpsiPtGen=JpsiGen.Pt();
    TJpsiEtaGen=JpsiGen.Eta();
    TJpsiPhiGen=JpsiGen.Phi();

    TPhiMassGen=PhiGen.M();
    TPhiPtGen=PhiGen.Pt();
    TPhiEtaGen=PhiGen.Eta();
    TPhiPhiGen=PhiGen.Phi();

    TMupPtGen=MupGen.Pt();
    TMupPhiGen=MupGen.Phi();
    TMupEtaGen=MupGen.Eta();

    TMumPtGen=MumGen.Pt();
    TMumPhiGen=MumGen.Phi();
    TMumEtaGen=MumGen.Eta();

    TKpPtGen=KpGen.Pt();
    TKpPhiGen=KpGen.Phi();
    TKpEtaGen=KpGen.Eta();

    TKmPtGen=KmGen.Pt();
    TKmPhiGen=KmGen.Phi();
    TKmEtaGen=KmGen.Eta();
    
    outTree->Fill();
  }

  flag = true;
  return flag;

}



int PDAnalyzer::HLTMatch(int trackIdx,int isMu){
  
  int hltIdx= -999;

  float mass;//to be set properly for the two objects
  if (isMu==1) mass=0.105658;
  else mass=0.493677;;

  TLorentzVector trk, hltObj;
  trk.SetPtEtaPhiM(trkPt->at(trackIdx),trkEta->at(trackIdx),trkPhi->at(trackIdx),mass);
  
  float minDist=0.1; 

  // type 5 = track, 2 = muon
  for (int k=0; k< nHLTObjects; k++){

    if (hltObjType->at(k) != 2 && isMu==1) continue; // check only matches between trigger muons and muons
    if (hltObjType->at(k) != 5 && isMu==0) continue; // check only matches between trigger tracks and kaons

    hltObj.SetPtEtaPhiM(hltPt->at(k),hltEta->at(k),hltPhi->at(k),mass);
    
    float dRtmp=trk.DeltaR(hltObj);
    
    if (dRtmp < minDist){
      minDist=dRtmp;
      hltIdx=k;
    }
  }
  return hltIdx;
}


int PDAnalyzer::GetMixStatus( unsigned int genindex )
{
  int Code = genId->at( genindex );

  const vector <int>& aD = allDaughters(genindex);
  if( aD.size()>0 && genId->at(aD[0]) == -Code ) return 2;

  const vector <int>& aM = allMothers(genindex); 
  if( aM.size()>0 && genId->at(aM[0]) == -Code ) return 1;

  return 0;
}



void PDAnalyzer::FindPV(TVector3 & sv, TVector3 & pv, TLorentzVector & BsP4, int & index) {
  
  TVector3 BsP3=BsP4.Vect();
  TVector3 tmpPV;
  float tmpCos=-999.;

  int iPV;
  for ( iPV = 0; iPV < nPVertices; ++iPV ) {
    tmpPV.SetXYZ(pvtX->at(iPV),pvtY->at(iPV),pvtZ->at(iPV));
    TVector3 diff= sv-tmpPV; //vector pointing from PV to SV
    if (abs(sv.Z()-tmpPV.Z())>0.5) continue; //if the PV has a z distance more than 0.5 cm just skip the vertex
    //diff.SetZ(0.);
    //BsP3.SetZ(0.);
    float cosPoint=cos(diff.Angle(BsP3));
    if (cosPoint > tmpCos){
      tmpCos=cosPoint;
      TcosPoint=cosPoint;
      pv.SetXYZ(pvtX->at(iPV),pvtY->at(iPV),pvtZ->at(iPV));
      index=iPV;
    }
  }  

  return;
}

void PDAnalyzer::resetTree() {

  if (!genOnly){

    Ttag=-999;
    Tmtag=-999;
    TmtagCal=-999;
    TmtagCalBs=-999;

    Tangle_cospsi=-999;
    Tangle_costheta=-999;
    Tangle_phi=-999;

    TLxy=-999;
    TLxyz=-999;
    TLxyErr=-999;
    TLxyzErr=-999;

    TVprob=-999;
    TVprobPhi=-999;
    TVprobJpsi=-999;

    TBsMass=-999; 
    TBsPt=-999;
    TBsPhi=-999;
    TBsEta=-999;

    TJpsiMass=-999; 
    TJpsiPt=-999;
    TJpsiPhi=-999;
    TJpsiEta=-999;

    TPhiMass=-999; 
    TPhiPt=-999;
    TPhiPhi=-999;
    TPhiEta=-999;

    TMupPt=-999;
    TMupPhi=-999;
    TMupEta=-999;
    TMupHits=-999;

    TMumPt=-999;
    TMumPhi=-999;
    TMumEta=-999;
    TMumHits=-999;

    TKpPt=-999;
    TKpPhi=-999;
    TKpEta=-999;
    TKpHits=-999;

    TPiPt=-999;
    TPiPhi=-999;
    TPiEta=-999;
    TPiHits=-999;
  
    TNpv=-999;
    TBcMass=-999;

    TKmPt=-999;
    TKmPhi=-999;
    TKmEta=-999;
    TKmHits=-999;

    TBsMassNoRefit=-999;
    TJpsiMassNoRefit=-999;
    TPhiMassNoRefit=-999;

    TBsMassFromSV=-999;
    TJpsiMassFromSV=-999;
    TPhiMassFromSV=-999;
    
    TcosPoint=-999;
    
    Trun=-999;
    Tevt=-999;
    Tlumi=-999;
    
    THLT_Jtktk=-999; 
    THLT_Jmu=-999; 
    THLT_Jtk=-999;
  }

  if (use_gen || genOnly){

    Tmatched=-999;//mc matching

    TgenFlavour=-999;

    Tangle_cospsiGen=-999;
    Tangle_costhetaGen=-999;
    Tangle_phiGen=-999;

    TLxyGen=-999;
    TLxyzGen=-999;

    TBsMassGen=-999; 
    TBsPtGen=-999;
    TBsPhiGen=-999;
    TBsEtaGen=-999;

    TJpsiMassGen=-999; 
    TJpsiPtGen=-999;
    TJpsiPhiGen=-999;
    TJpsiEtaGen=-999;

    TPhiMassGen=-999; 
    TPhiPtGen=-999;
    TPhiPhiGen=-999;
    TPhiEtaGen=-999;

    TMupPtGen=-999;
    TMupPhiGen=-999;
    TMupEtaGen=-999;
       
    TMumPtGen=-999;
    TMumPhiGen=-999;
    TMumEtaGen=-999;
  
    TKpPtGen=-999;
    TKpPhiGen=-999;
    TKpEtaGen=-999;

    TKmPtGen=-999;
    TKmPhiGen=-999;
    TKmEtaGen=-999;

  }

  return;
}

void PDAnalyzer::endJob() {
// to skim the N-tuple "uncomment" the following line
//  closeSkim();

// additional features
//  DataSetFilter::endJob();                       // dataset filter
//  tWriter->close();                              // second ntuple
    return;
}


void PDAnalyzer::save() {
#  if UTIL_USE == FULL
  // explicit saving not necessary for "autoSavedObjects"
  autoSave();
#elif UTIL_USE == BARE
  // explicit save histos when not using the full utility
  hptmumax->Write();
  hptmu2nd->Write();
  hptmu   ->Write();
#endif

  outTree->Write();

  return;
}


