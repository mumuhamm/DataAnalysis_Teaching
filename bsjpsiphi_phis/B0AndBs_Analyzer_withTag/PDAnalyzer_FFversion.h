#ifndef PDAnalyzer_H
#define PDAnalyzer_H

#include "TH1.h"
#include "PDAnalyzerUtil.h"
#include "PDAnalysis/Ntu/interface/PDGenHandler.h"
#include "PDMuonVar.h"
#include "PDSoftMuonMvaEstimator.h"
#include "AlbertoUtil.h"
#include "OSMuonMvaTag.h"
// to skim the N-tuple "uncomment" the following line
//#include "NtuTool/Common/interface/TreeFilter.h"
#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "Math/LorentzVector.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TVector3.h"
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>


// additional features
//#include "DataSetFilter.h"                       // dataset filter
// class PDSecondNtupleWriter;                      // second ntuple

// to skim the N-tuple replace the the following line
// with the "commented" ones
class PDAnalyzer: public virtual PDAnalyzerUtil
,                 public virtual PDGenHandler
,                 public virtual PDMuonVar
,                 public virtual PDSoftMuonMvaEstimator
,                 public virtual AlbertoUtil
,                 public virtual OSMuonMvaTag
// to skim the N-tuple "uncomment" the following line
//,                 public virtual TreeFilter
 {

 public:

  PDAnalyzer();
  virtual ~PDAnalyzer();

  // function called before starting the analysis
  virtual void beginJob();

  // functions to book the histograms
  void book();

  // functions called for each event
  // function to reset class content before reading from file
  virtual void reset();
  // function to do event-by-event analysis,
  // return value "true" for accepted events
  virtual bool analyze( int entry, int event_file, int event_tot );

  // function called at the end of the analysis
  virtual void endJob();

  // functions called at the end of the event loop
// to plot some histogram immediately after the ntuple loop
// "uncomment" the following line
//  virtual void plot();   // plot the histograms on the screen
  virtual void save();   // save the histograms on a ROOT file

  virtual void resetTree();

  virtual void FindPV(TVector3 &, TVector3 &, TLorentzVector &, int &);
  virtual int  HLTMatch(int, int);
  int GetMixStatus( unsigned int genindex );


  bool verbose;

 protected:

  double ptCut;
  double maxAngle;

  //PARAMETERS
  TString outputFile;
  TString process;

// additional features: second ntuple
  // PDSecondNtupleWriter* tWriter;                 // second ntuple

 private:

  TH1F* hAnglePhi;
  TH1F* hCosTheta;
  TH1F* hCosPsi;

  TH1F* hPhi;
  TH1F* hJpsi;
  TH1F* hBs;

  TH1F* hPhiHLT;
  TH1F* hJpsiHLT;

  bool genOnly;
  bool debug;

  TTree* outTree;
  float  TVprob,TVprobJpsi, TVprobPhi, Tmatched, Tangle_cospsiGen,  Tangle_costhetaGen, Tangle_phiGen, Tangle_cospsi,  Tangle_costheta, Tangle_phi, TLxy, TLxyPD, TLxyANBS, TLxyAT, TLxyCNBS, TLxyCO,  TLxyz, TLxyErr, TLxyzErr, TLxyGen, TLxyzGen;
  int TnBs;

  //the Bs
  float TBsMass, TBsPt, TBsPhi, TBsEta, TBsMassGen, TBsPtGen, TBsPhiGen, TBsEtaGen;
  float TBsMassNoRefit, TJpsiMassNoRefit, TPhiMassNoRefit, TBsMassFromSV, TJpsiMassFromSV, TPhiMassFromSV;

  //Jpsi and Phi
  float TJpsiMass, TJpsiPt, TJpsiPhi, TJpsiEta, TJpsiMassGen, TJpsiPtGen, TJpsiPhiGen, TJpsiEtaGen;
  float TPhiMass, TPhiPt, TPhiPhi, TPhiEta, TPhiMassGen, TPhiPtGen, TPhiPhiGen, TPhiEtaGen;

  //Muons and Kaons
  float TMupPt, TMupPhi, TMupEta, TMupHits, TMupPtGen, TMupPhiGen, TMupEtaGen;
  float TMumPt, TMumPhi, TMumEta, TMumHits, TMumPtGen, TMumPhiGen, TMumEtaGen;

  float  TKpPt, TKpPhi, TKpEta,TKpHits, TKpPtGen, TKpPhiGen, TKpEtaGen;
  float  TKmPt, TKmPhi, TKmEta, TKmHits, TKmPtGen, TKmPhiGen, TKmEtaGen;

  float  TPiPt, TPiPhi, TPiEta, TPiHits, TNpv, TBcMass;

  float TcosPoint;

  //event
  float Tmtag, TmtagCal, TmtagCalBs, TmtagErr;
  float Trun, Tlumi, THLT_Jtktk, THLT_Jmu, THLT_Jtk, TisBs;  
  int   TgenFlavour,Ttag;
  float TmatchHLTmu, TmatchHLTk;
  ULong64_t Tevt;

  float TKaon_Deltadz,TKaon_DeltaR;
  // dummy copy constructor and assignment
  PDAnalyzer           ( const PDAnalyzer& );
  PDAnalyzer& operator=( const PDAnalyzer& );

};


#endif

