// -*- C++ -*-
//
// Package:    Eff_Directory/effAnalyzer
// Class:      effAnalyzer
//
/**\class effAnalyzer effAnalyzer.cc Eff_Directory/effAnalyzer/plugins/effAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Muhammad Muhammad Alibordi
//         Created:  Mon, 20 Apr 2020 10:46:25 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "FWCore/ParameterSet/interface/ParameterSetfwd.h"
#include <iostream>
#include <memory>
#include <vector>
#include <iomanip>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TTree.h>
#include <TMath.h>
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "PhysicsTools/PatAlgos/plugins/PATJetProducer.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include <iostream>
#include <algorithm>
#include <bitset>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupMixingContent.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupVertexContent.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "MuonAnalysis/MuonAssociators/interface/L1MuonMatcherAlgo.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/Association.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/TriggerReport.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"   
#include "DataFormats/MuonReco/interface/MuonQuality.h" 
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"
#include "DataFormats/PatCandidates/interface/TriggerCondition.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"   
using namespace edm;
    using namespace std;
    using namespace reco;
    using namespace pat;






class effAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit effAnalyzer(const edm::ParameterSet&);
      ~effAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    typedef math::XYZTLorentzVectorF LorentzVector;
   

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
    
   void bookHists(edm::Service<TFileService>& fs, const std::string& suffix);
   void fillHists(const LorentzVector& lv, const std::string& suffix); 

      // ----------member data ---------------------------
    
    double nominalMuonMass = 0.1056583;
    double nominalJpsiMass = 3.096916;
    bool isMCstudy_;
    double JpsiMassWindow; 
    edm::InputTag muon;
    edm::EDGetTokenT<edm::View<reco::Muon>> muonTok;
    edm::InputTag primaryvertex;
    edm::EDGetTokenT<edm::View<reco::Vertex>> primaryvertexTok;
    edm::InputTag track;
    edm::EDGetTokenT<edm::View<reco::Track>> trackTok;
    edm::InputTag geninfo;
    edm::EDGetTokenT<edm::View<reco::GenParticle>> geninfoTok;
    edm::InputTag triggerresults;
    edm::EDGetTokenT<edm::TriggerResults>triggerresultsTok;
    edm::InputTag trigsummary;
    edm::EDGetTokenT<trigger::TriggerEvent>trigsummaryTok;
    edm::InputTag filtertag;
    HLTConfigProvider hltConfiguration;
    edm::InputTag offlineBS;
    edm::EDGetTokenT<reco::BeamSpot> offlineBSTok;
    
     std::map<std::string,TH1F*> hists_1d_;  
    //-------------Tree and variables -------------------
     TTree * eff_muon;
    Double_t trkmuon_pt_all = -999999;
             Double_t trkmuon_eta_all = -999999;
             Double_t trkmuon_phi_all = -999999;
             Double_t trkmuon_pt_pass = -999999;
             Double_t trkmuon_eta_pass = -999999;
             Double_t trkmuon_phi_pass = -999999;
             Double_t trkmuon_pt_fail = -999999;
             Double_t trkmuon_eta_fail = -999999;
             Double_t trkmuon_phi_fail = -999999;
             Double_t trigmuon_pt_all = -999999;
             Double_t trigmuon_eta_all = -999999;
             Double_t trigmuon_phi_all = -999999;
             Double_t trigmuon_pt_pass = -999999;
             Double_t trigmuon_eta_pass = -999999;
             Double_t trigmuon_phi_pass = -999999;
             Double_t trigmuon_pt_fail = -999999;
             Double_t trigmuon_eta_fail = -999999;
             Double_t trigmuon_phi_fail = -999999;
             Double_t trkDiMuon_all = -999999;
             Double_t trkDiMuon_pass = -999999;
             Double_t trkDiMuon_fail = -999999;
             Double_t trigDiMuon_all = -999999;
             Double_t trigDiMuon_pass = -999999;
             Double_t trigDiMuon_fail = -999999;
    Double_t muon_gen_pt, muon_gen_eta, muon_gen_phi; //genmatched tracking soft muons
   bool triggerbit_HLT_UPC_SingleMu0 ;
   const MagneticField* magneticField_; 
  double BSx         = -9999999.;
  double BSy         = -9999999.;
  double BSz         = -9999999.;
  double BSdx        = -9999999.;
  double BSdy        = -9999999.;
  double BSdz        = -9999999.;
  double BSdxdz      = -9999999.;
  double BSdydz      = -9999999.;
  double BSsigmaZ    = -9999999.;
  double BSdsigmaZ   = -9999999.;
  Double_t MuonsDCA=-999999;   
  Double_t vtxProb_dimuon = -999999;
  Double_t CosAlpha = -999999; 
  Double_t deltaR_Trigobj_Tag = -999999;
  Double_t deltaR_Trk_Tag = -999999;
  Double_t deltaR_Trk_Probe = -999999;


};


effAnalyzer::effAnalyzer(const edm::ParameterSet& iConfig):
  isMCstudy_(iConfig.getParameter<bool>("isMCstudy")),
  JpsiMassWindow(iConfig.getParameter<double>("JpsiWindow")),
  muon(iConfig.getUntrackedParameter<edm::InputTag>("muon")),
  muonTok(consumes<edm::View<reco::Muon>>(muon)),
  primaryvertex(iConfig.getUntrackedParameter<edm::InputTag>("primaryvertex")),
  primaryvertexTok(consumes<edm::View<reco::Vertex>>(primaryvertex)),
  track(iConfig.getUntrackedParameter<edm::InputTag>("track")),
  trackTok(consumes<edm::View<reco::Track>>(track)),
  geninfo(iConfig.getUntrackedParameter<edm::InputTag>("geninfo")),
  geninfoTok(consumes<edm::View<reco::GenParticle>>(geninfo)),
  triggerresults(iConfig.getUntrackedParameter<edm::InputTag>("triggerresults")),
  triggerresultsTok(consumes<edm::TriggerResults>(triggerresults)),
  trigsummary(iConfig.getUntrackedParameter<edm::InputTag>("trigsummary")),
  trigsummaryTok(consumes<trigger::TriggerEvent>(trigsummary)),
  filtertag(iConfig.getUntrackedParameter<edm::InputTag>("filtertag")),
  offlineBS(iConfig.getUntrackedParameter<edm::InputTag>("offlineBS")),
  offlineBSTok(consumes<reco::BeamSpot>(offlineBS))

{
                               usesResource("TFileService");
                               edm::Service<TFileService> fs;
                               eff_muon = fs->make<TTree>("eff_muon","eff_muon");
    
             eff_muon->Branch("muon_gen_pt", &muon_gen_pt); eff_muon->Branch("muon_gen_eta", &muon_gen_eta);
             eff_muon->Branch("muon_gen_phi", &muon_gen_phi);
             eff_muon->Branch("trkmuon_pt_all", &trkmuon_pt_all); eff_muon->Branch("trkmuon_eta_all", &trkmuon_eta_all);
             eff_muon->Branch("trkmuon_phi_all", &trkmuon_phi_all);
             eff_muon->Branch("trkmuon_pt_pass", &trkmuon_pt_pass); eff_muon->Branch("trkmuon_eta_pass", &trkmuon_eta_pass);
             eff_muon->Branch("trkmuon_phi_pass", &trkmuon_phi_pass);
             eff_muon->Branch("trkmuon_pt_fail", &trkmuon_pt_fail); eff_muon->Branch("trkmuon_eta_fail", &trkmuon_eta_fail);
             eff_muon->Branch("trkmuon_phi_fail", &trkmuon_phi_fail);
             eff_muon->Branch("trkDiMuon_all",&trkDiMuon_all);
             eff_muon->Branch("trkDiMuon_pass",&trkDiMuon_pass);
             eff_muon->Branch("trkDiMuon_fail",&trkDiMuon_fail);
             eff_muon->Branch("trigmuon_pt_all", &trigmuon_pt_all); eff_muon->Branch("trigmuon_eta_all", &trigmuon_eta_all);
             eff_muon->Branch("trigmuon_phi_all", &trigmuon_phi_all);
             eff_muon->Branch("trigmuon_pt_pass", &trigmuon_pt_pass); eff_muon->Branch("trigmuon_eta_pass", &trigmuon_eta_pass);
             eff_muon->Branch("trigmuon_phi_pass", &trigmuon_phi_pass);
             eff_muon->Branch("trigmuon_pt_fail", &trigmuon_pt_fail); eff_muon->Branch("trigmuon_eta_fail", &trigmuon_eta_fail);
             eff_muon->Branch("trigmuon_phi_fail", &trigmuon_phi_fail);
             eff_muon->Branch("trigDiMuon_all",&trigDiMuon_all);
             eff_muon->Branch("trigDiMuon_pass",&trigDiMuon_pass);
             eff_muon->Branch("trigDiMuon_fail",&trigDiMuon_fail);
             eff_muon->Branch("MuonsDCA",&MuonsDCA);
             eff_muon->Branch("vtxProb_dimuon",&vtxProb_dimuon);
             eff_muon->Branch("CosAlpha",&CosAlpha);
             eff_muon->Branch("deltaR_Trigobj_Tag",&deltaR_Trigobj_Tag); 
             eff_muon->Branch("deltaR_Trk_Tag",&deltaR_Trk_Tag);
             eff_muon->Branch("deltaR_Trk_Probe",&deltaR_Trk_Probe);  
             hists_1d_["h_passtrig"] = fs->make<TH1F>("h_passtrig" , "; passed trigger" , 2 , 0. , 2. );
             hists_1d_["h_trigmll_all"] = fs->make<TH1F>("h_trigmll_all" , "; m_{ll} (GeV)" , 100 , 2.7 , 3.5 );
             hists_1d_["h_trigmll_pass"] = fs->make<TH1F>("h_trigmll_pass" , "; m_{ll} (GeV)" , 100 , 2.7 , 3.5 );
             hists_1d_["h_trigmll_fail"] = fs->make<TH1F>("h_trigmll_fail" , "; m_{ll} (GeV)" , 100 , 2.7 , 3.5 );
             hists_1d_["h_trkmll_all"] = fs->make<TH1F>("h_trkmll_all", "; M_{#mu^{+}#mu^{-}} (GeV)", 100, 2.7,3.5);
             hists_1d_["h_trkmll_pass"] = fs->make<TH1F>("h_trkmll_pass", "; M_{#mu^{+}#mu^{-}} (GeV)", 100, 2.7, 3.5);
             hists_1d_["h_trkmll_fail"] = fs->make<TH1F>("h_trkmll_fail", "; M_{#mu^{+}#mu^{-}} (GeV)", 100, 2.7,3.5);
             bookHists(fs,"probe_all");
             bookHists(fs,"probe_pass");
             bookHists(fs,"probe_fail");
 
      
    

}


effAnalyzer::~effAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}




// ------------ method called for each event  ------------
void
effAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace std;
    using namespace edm;
    using namespace reco;
    using namespace trigger;

 
    edm::Handle< View<reco::Muon>> muons;
    iEvent.getByToken(muonTok, muons);
    
    edm::Handle<View<reco::Vertex>> pv;
    iEvent.getByToken(primaryvertexTok, pv);
    
    edm::Handle<edm::View<reco::Track>> trackinfo;
    iEvent.getByToken(trackTok, trackinfo);
      
                 edm::Handle<edm::View<reco::GenParticle>> genpart;
                      if(isMCstudy_)
                           {
                                 iEvent.getByToken(geninfoTok, genpart);
                                 for(unsigned int n = 0; n< genpart->size(); ++n){
                                 const reco::GenParticle & gen = (*genpart)[n];
                                 if (std::abs(gen.pdgId()) != 13 || gen.status() != 1) continue; 
                                 muon_gen_pt = gen.pt(); muon_gen_eta = gen.eta(); muon_gen_phi = gen.phi();
                              }
                          }
    
    edm::Handle<edm::TriggerResults> hltresults;
    iEvent.getByToken(triggerresultsTok, hltresults);
    
    edm::Handle<trigger::TriggerEvent> trigsummaryinfo;
    iEvent.getByToken(trigsummaryTok, trigsummaryinfo);

                   edm::ESHandle<MagneticField> magneticField;
                   iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
                   magneticField_ = &*magneticField;


   edm::Handle<reco::BeamSpot> vertexBeamSpot ;
   iEvent.getByToken(offlineBSTok,vertexBeamSpot);

         BSx = vertexBeamSpot->x0(); 
         BSy = vertexBeamSpot->y0();
         BSz = vertexBeamSpot->z0();
         BSdxdz = vertexBeamSpot->dxdz();
         BSdydz = vertexBeamSpot->dydz();

     std::vector<TLorentzVector> trackingMuons_pos;
     std::vector<TLorentzVector> trackingMuons_neg;
    for(unsigned int k =0; k < trackinfo->size() ; ++k)
       {  
          const reco::Track & track1 = (*trackinfo)[k];
          if(track1.charge() < 0)continue;
          if(track1.pt() < 1.0) continue;
          if(fabs(track1.eta()) > 2.4) continue;
          if(!track1.quality(reco::TrackBase::highPurity))continue;
          if(track1.numberOfValidHits() < 1)continue;
          TLorentzVector track1_lv;
          track1_lv.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi(), nominalMuonMass);
          trackingMuons_pos.push_back(TLorentzVector(track1_lv));     
         

         for(unsigned int l =0; l < trackinfo->size() ; ++l)
             {
                const reco::Track & track2 = (*trackinfo)[l];
          if(track2.charge() > 0)continue;
          if(track2.pt() < 1.0) continue;
          if(fabs(track2.eta()) > 2.4) continue;
          if(!track2.quality(reco::TrackBase::highPurity))continue;
          if(track2.numberOfValidHits() < 1)continue;
          TLorentzVector track2_lv;
          track2_lv.SetPtEtaPhiM(track2.pt(), track2.eta(), track2.phi(), nominalMuonMass);
          trackingMuons_neg.push_back(TLorentzVector(track2_lv));          

      
}             
      }


    

        for (unsigned int j =0 ; j< muons->size(); ++j)
            {
                 const reco::Muon & muon1_trk = (*muons)[j];
                 if(muon1_trk.charge()< 0)continue;
                 if(muon1_trk.pt() < 1.0) continue;
                 if(fabs(muon1_trk.eta()) > 2.4) continue;
                 TrackRef intracktag = muon1_trk.innerTrack();
                 if(intracktag.isNull())continue;
                 if(!muon1_trk.isGlobalMuon())continue;
                 if(!(muon1_trk.passed(reco::Muon::SoftCutBasedId)))continue;
                 if(!(muon1_trk.passed(reco::Muon::TkIsoLoose)))continue;
                 bool trackingMatch_tag = false;
             for (unsigned int trackindex=0; trackindex< trackingMuons_pos.size(); ++trackindex)
                {
                        deltaR_Trk_Tag = ROOT::Math::VectorUtil::DeltaR(muon1_trk.p4(),trackingMuons_pos.at(trackindex));
                        if (ROOT::Math::VectorUtil::DeltaR(muon1_trk.p4(),trackingMuons_pos.at(trackindex)) < 0.7) trackingMatch_tag = true;
                 }
                 if (!trackingMatch_tag) continue;



       for (unsigned int l =0 ; l< muons->size(); ++l)
         {
              const reco::Muon & muon2_trk = (*muons)[l];
              if(muon2_trk.charge() > 0)continue;             
              if(muon2_trk.pt() < 1.0) continue;
              if(fabs(muon2_trk.eta()) > 2.4) continue;
              TrackRef intracktag2 = muon2_trk.innerTrack();
              if(intracktag2.isNull())continue;
              //if(!muon2_trk.isGlobalMuon())continue;
              if(!muon2_trk.isStandAloneMuon())continue;
              if(!(muon1_trk.passed(reco::Muon::SoftCutBasedId)))continue;
              trkmuon_pt_all = muon2_trk.pt(); trkmuon_eta_all = muon2_trk.eta(); trkmuon_phi_all = muon2_trk.phi();
               //std::cout<<"pT :: "<<trkmuon_pt_all<<"\t"<<"eta :: "<<trkmuon_eta_all<<"\t"<<"phi :: "<<trkmuon_phi_all<<"\n";
              TLorentzVector tracktagleg, trackprobeleg;
              tracktagleg.SetPtEtaPhiM(muon1_trk.pt(), muon1_trk.eta(), muon1_trk.phi(), nominalMuonMass);
              trackprobeleg.SetPtEtaPhiM(muon2_trk.pt(), muon2_trk.eta(), muon2_trk.phi(), nominalMuonMass);
              TLorentzVector trk_invdimuon = tracktagleg + trackprobeleg;
              if(abs(trk_invdimuon.M()-nominalJpsiMass) > JpsiMassWindow)continue;
              trkDiMuon_all = trk_invdimuon.M();
              hists_1d_["h_trkmll_all"]->Fill(trk_invdimuon.M());

             bool trackingMatch_probe = false;
             for (unsigned int trackindex2=0; trackindex2< trackingMuons_neg.size(); ++trackindex2)
               {
                  deltaR_Trk_Probe = ROOT::Math::VectorUtil::DeltaR(muon2_trk.p4(),trackingMuons_neg.at(trackindex2));
                  if (fabs(muon2_trk.eta() - trackingMuons_neg.at(trackindex2).Eta())> 0.2)continue;
                  if (ROOT::Math::VectorUtil::DeltaR(muon2_trk.p4(),trackingMuons_neg.at(trackindex2)) < 0.7) trackingMatch_probe = true;
             }
       
           if(trackingMatch_probe)
                 { 
                   trkmuon_pt_pass = muon2_trk.pt(); 
                   trkmuon_eta_pass = muon2_trk.eta(); 
                   trkmuon_phi_pass = muon2_trk.phi();
                   trkDiMuon_pass = trk_invdimuon.M();
                   hists_1d_["h_trkmll_pass"]->Fill(trk_invdimuon.M());
                  }
            else
               {
                 trkmuon_pt_fail = muon2_trk.pt(); 
                 trkmuon_eta_fail = muon2_trk.eta(); 
                 trkmuon_phi_fail = muon2_trk.phi();
                 trkDiMuon_fail = trk_invdimuon.M();
                 hists_1d_["h_trkmll_fail"]->Fill(trk_invdimuon.M());
               }

   
       }
    }


        trigger::TriggerObjectCollection allTriggerObjects;
        trigger::TriggerObjectCollection selectedObjects;
    
    
    TLorentzVector mtrk;
    std::vector<TLorentzVector> trigMuons;
    const  edm::TriggerNames & triggernames = iEvent.triggerNames(*hltresults);
    int ntrigs = hltresults->size();
   
        size_t filterIndex = (*trigsummaryinfo).filterIndex(filtertag);
        allTriggerObjects = trigsummaryinfo->getObjects();
        if (filterIndex < (*trigsummaryinfo).sizeFilters())
        { 
          const trigger::Keys &keys = (*trigsummaryinfo).filterKeys(filterIndex);
          for (size_t j = 0; j < keys.size(); j++)
           {
            trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
            selectedObjects.push_back(foundObject);
            mtrk.SetPtEtaPhiM(selectedObjects[j].pt(), selectedObjects[j].eta(), selectedObjects[j].phi(), nominalMuonMass);
            trigMuons.push_back(TLorentzVector(mtrk));
            //std::cout<<"  -pt "<< selectedObjects[j].pt()<<" -eta "<<selectedObjects[j].eta()<<" -phi "<<selectedObjects[j].phi()<<"\n";
           }
        }



    for (int itrig = 0; itrig != ntrigs; ++itrig)
    {
                                
        TString trigName = triggernames.triggerName(itrig);
        if (trigName=="HLT_HIUPC_SingleMu0_NotMBHF2AND_v1") triggerbit_HLT_UPC_SingleMu0 = hltresults->accept(itrig);
            
    }
    
    edm::View<reco::Vertex>::const_iterator vertex ;
    for ( edm::View<reco::Vertex>::const_iterator vtx = pv->begin(); vtx != pv->end(); ++vtx ) {
      if (  !vtx->isFake() && vtx->ndof()>=4. && vtx->position().Rho()<=6.0 && fabs(vtx->position().Z())<=24.0 ) {
        if (vertex == pv->end()) {
      vertex = vtx;
      
      break;
        }
      }
    }
    //const float dr_trigmatch = 0.4;
    const float dr_trigmatch_wider = 0.7;
    
     for (unsigned int k =0 ; k< muons->size(); ++k)
{
    const reco::Muon & muon_tag = (*muons)[k];
    
    //std::cout<<" - tag cand: pt: "<<muon_tag.pt()<<", eta: "<< muon_tag.eta()<<", phi: "<< muon_tag.phi()<<"\n";
    //double pfiso04 = (muon_tag->pfIsolationR04().sumChargedHadronPt + max(0.,muon_tag->pfIsolationR04().sumNeutralHadronEt + muon_tag->pfIsolationR04().sumPhotonEt - 0.5*muon_tag->pfIsolationR04().sumPUPt))/muon_tag->pt();

          
         //---------------------------------------tag muon
                         if (muon_tag.pt() < 1.0) continue;
                         if (fabs(muon_tag.eta()) > 2.4) continue;  
                         TrackRef intracktag = muon_tag.innerTrack();
                         if(intracktag.isNull())continue;
                         if(!muon_tag.isGlobalMuon())continue;
                         if(!(muon_tag.passed(reco::Muon::SoftCutBasedId)))continue;
                         if(!(muon_tag.passed(reco::Muon::TkIsoLoose)))continue;
                         


                        /* if(!(muon_tag.passed(reco::Muon::CutBasedIdLoose)))continue;
                         if(muon_tag.passed(reco::Muon::CutBasedIdMedium ))continue;
                         if(muon_tag.passed(reco::Muon::CutBasedIdMediumPrompt))continue;
                         if(muon_tag.passed(reco::Muon::CutBasedIdTight))continue;
                         if(muon_tag.passed(reco::Muon::CutBasedIdGlobalHighPt))continue;
                         if(muon_tag.passed(reco::Muon::CutBasedIdTrkHighPt))continue;
                         if(muon_tag.passed(reco::Muon::PFIsoVeryLoose))continue;
                         if(muon_tag.passed(reco::Muon::PFIsoLoose))continue;
                         if(muon_tag.passed(reco::Muon::PFIsoMedium))continue;
                         if(muon_tag.passed(reco::Muon::PFIsoTight))continue;
                         if(muon_tag.passed(reco::Muon::PFIsoVeryTight))continue;
                         if(muon_tag.passed(reco::Muon::TkIsoLoose))continue;
                         if(muon_tag.passed(reco::Muon::TkIsoTight))continue;
                         if(!(muon_tag.passed(reco::Muon::SoftCutBasedId)))continue;
                         if(!(muon_tag.passed(reco::Muon::SoftMvaId )))continue;
                         if(muon_tag.passed(reco::Muon::MvaLoose))continue;
                         if(muon_tag.passed(reco::Muon::MvaMedium))continue;
                         if(muon_tag.passed(reco::Muon::MvaTight))continue;
                         if(muon_tag.passed(reco::Muon::MiniIsoLoose))continue;
                         if(muon_tag.passed(reco::Muon:: MiniIsoMedium))continue;
                         if(muon_tag.passed(reco::Muon::MiniIsoTight))continue;
                         if(muon_tag.passed(reco::Muon::MiniIsoVeryTight))continue;*/
                         
                            


           bool trigmatch_tag = false;
          for (unsigned int trigindex=0; trigindex< trigMuons.size(); ++trigindex) {
           deltaR_Trigobj_Tag = ROOT::Math::VectorUtil::DeltaR(muon_tag.p4(),trigMuons.at(trigindex));
           if (ROOT::Math::VectorUtil::DeltaR(muon_tag.p4(),trigMuons.at(trigindex)) < dr_trigmatch_wider) trigmatch_tag = true;
         }
          if (!trigmatch_tag) continue;
   

//---------------------------------------------probemuon
     for ( unsigned int l =0 ; l< muons->size(); ++l ) {
         
                         const reco::Muon & muon_probe = (*muons)[l];
                         if( muon_tag.charge() == muon_probe.charge());
                         if (ROOT::Math::VectorUtil::DeltaR(muon_tag.p4(),muon_probe.p4()) < 0.01) continue;//tag and probe should be differenrt from each other 
                         if (muon_probe.pt() < 1.0) continue;
                         if (fabs(muon_probe.eta()) > 2.4) continue; 
                         TrackRef intrackprobe = muon_probe.innerTrack();
                         if(intrackprobe.isNull())continue;
                         if(!muon_probe.isStandAloneMuon())continue; 
                         //if(!muon_probe.isGlobalMuon())continue;       
                         if(!(muon_probe.passed(reco::Muon::SoftCutBasedId)))continue;

                               

                              //std::cout<<"hout not the fucking faded gustavo"<<"\n";
                              TrackRef tagTrkRef = muon_tag.get<TrackRef>();
                              TrackRef probeTrkRef = muon_probe.get<TrackRef>();
                              TransientTrack tagTrkTT(tagTrkRef, magneticField_);
                              TransientTrack probeTrkTT(probeTrkRef, magneticField_);
                              if(tagTrkTT.numberOfValidHits() <= 0) continue;
                              if(probeTrkTT.numberOfValidHits() <= 0) continue;
                              std::vector<TransientTrack> allTrk;
                              allTrk.push_back(tagTrkTT);
                              allTrk.push_back(probeTrkTT);
                              if(allTrk.size() <  2 )continue;
                              KalmanVertexFitter kvf(true);
                              TransientVertex tv = kvf.vertex(allTrk);

                              Vertex vertex = tv;
                              vtxProb_dimuon = TMath::Prob(vertex.chi2(),(int)vertex.ndof());
                              math::XYZVector  pperp(muon_tag.px() + muon_probe.px(), muon_tag.py() + muon_probe.py(), 0.);
                              reco::Vertex::Point  vpoint = vertex.position(); 
                              GlobalPoint secondaryVertex (vpoint.x(), vpoint.y(), vpoint.z());
                            // std::cout<<"Meku tastar guitara knhor aii aiia ii ai mia more"<<vpoint.x()<<"\n";
              GlobalPoint displacementFromBeamspot( -1*((BSx -  secondaryVertex.x()) +  (secondaryVertex.z() - BSz) * BSdxdz),-1*((BSy - secondaryVertex.y())+  (secondaryVertex.z() - BSz) * BSdydz), 0);
                              reco::Vertex::Point vperp(displacementFromBeamspot.x(),displacementFromBeamspot.y(),0.);
                              CosAlpha = vperp.Dot(pperp)/(vperp.R()*pperp.R());        
                              TrajectoryStateClosestToPoint tagTS = tagTrkTT.impactPointTSCP();
                              TrajectoryStateClosestToPoint probeTS = probeTrkTT.impactPointTSCP(); 
                             if (tagTS.isValid() && probeTS.isValid()) {
       
                                                                                ClosestApproachInRPhi cApp;
                                                                                cApp.calculate(tagTS.theState(), probeTS.theState());
                                                                                MuonsDCA = cApp.distance();
                                                                                //std::cout<<"Bailando, bailando.. "<<MuonsDCA<<"\n";
                                                                        }

         

                    //reco::CompositeCandidate dimuon;
                    //dimuon.addDaughter(muon_tag);
                    //dimuon.addDaughter(muon_probe);
                    //AddFourMomenta fourvec;
                    //fourvec.set(dimuon); 
         TLorentzVector tagleg, probeleg;
         tagleg.SetPtEtaPhiM(muon_tag.pt(), muon_tag.eta(), muon_tag.phi(), nominalMuonMass);
         probeleg.SetPtEtaPhiM(muon_probe.pt(), muon_probe.eta(), muon_probe.phi(), nominalMuonMass);
         TLorentzVector dimuon = tagleg+probeleg;
         if(abs(dimuon.M()-nominalJpsiMass) > JpsiMassWindow)continue;
         trigDiMuon_all = dimuon.M();
         hists_1d_["h_trigmll_all"]->Fill(dimuon.M());
         fillHists(LorentzVector(muon_probe.p4()),"probe_all"); 
         trigmuon_pt_all = muon_probe.pt(); trigmuon_eta_all = muon_probe.eta(); trigmuon_phi_all = muon_probe.phi();
           

         bool trigmatch_probe = false;
          
         for (unsigned int trigindex=0; trigindex < trigMuons.size(); ++trigindex) {
           if (ROOT::Math::VectorUtil::DeltaR(muon_probe.p4(),trigMuons.at(trigindex)) < dr_trigmatch_wider) trigmatch_probe = true;
           }
         
         
      if (trigmatch_probe)
               {
                           trigmuon_pt_pass = muon_probe.pt(); 
                           trigmuon_eta_pass = muon_probe.eta(); 
                           trigmuon_phi_pass = muon_probe.phi();
                           trigDiMuon_pass = dimuon.M();
                           fillHists(LorentzVector(muon_probe.p4()),"probe_pass");
                           hists_1d_["h_trigmll_pass"]->Fill(dimuon.M());
               }
        else 
              { 
                       trigmuon_pt_fail = muon_probe.pt();
                       trigmuon_eta_fail = muon_probe.eta(); 
                       trigmuon_phi_fail = muon_probe.phi();
                       trigDiMuon_fail = dimuon.M();
                       fillHists(LorentzVector(muon_probe.p4()),"probe_fail");
                       hists_1d_["h_trigmll_fail"]->Fill(dimuon.M());
              }
     }
     }
    
    eff_muon->Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
    
    
}
void effAnalyzer::bookHists(edm::Service<TFileService>& fs, const std::string& suffix) {

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_pt"+suf] = fs->make<TH1F>(Form("h_pt%s",suf.c_str()) , "; p_{T} [GeV]" , 100 , 0. , 20. );
  hists_1d_["h_eta"+suf] = fs->make<TH1F>(Form("h_eta%s",suf.c_str()) , "; #eta" , 100 , -2.4 , 2.4 );
  hists_1d_["h_phi"+suf] = fs->make<TH1F>(Form("h_phi%s",suf.c_str()) , "; #phi" , 100 , -TMath::Pi(),TMath::Pi() );

  return;
}

void effAnalyzer::fillHists(const LorentzVector& lv, const std::string& suffix) {

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_pt"+suf]->Fill(lv.pt());
  hists_1d_["h_eta"+suf]->Fill(lv.eta());
  hists_1d_["h_phi"+suf]->Fill(lv.phi());

  return;
}

// ------------ method called once each job just before starting event loop  ------------
void
effAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
effAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
effAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(effAnalyzer);
