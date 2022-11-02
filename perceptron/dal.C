// Author :  Muhammad Alibordi

//
#include "TGenPhaseSpace.h"
#include "TH2.h"
#include "TH1.h"
#include "TRandom.h"
#include "math.h"
#include "TF1.h"
#include "TComplex.h"

using namespace RooFit; 
using namespace std;

void dal(){
  

 


 
   TLorentzVector P(0,0,0,5.279);
   TLorentzVector p1,p2,p3,p12,p13,p23;
   
   Int_t n=3;
   Double_t m[3] = {0.493677,0.497672,0.497672};
   Double_t m12,m13,m23;
   
   TGenPhaseSpace S;
   
   S.SetDecay(P,n,m,"default");
   
   TH2F* dalitz = new TH2F("dalitz","Dalitz plot; M^{2}(K^{+}K_{s}) (GeV)^{2} ; M^{2}(K^{+}K_{s}) (GeV)^{2}",100,0,25,100,0,25);
  Double_t w;
  
  for(int i=0;i<500000;i++) {
    w = S.Generate();
    TLorentzVector *p1 = S.GetDecay(0);
    TLorentzVector *p2 = S.GetDecay(1);
    TLorentzVector *p3 = S.GetDecay(2);
    TLorentzVector p12= *p1+ *p2;
    TLorentzVector p13= *p1+ *p3;
    TLorentzVector p23= *p2+ *p3;
      
     
    m12 = p12.M2();
    m13 = p13.M2();
    m23 = p23.M2();
      
    dalitz->Fill(m12,m13,w);
     
    // cout<<w<<" ; "<<S.GetWtMax() <<endl;
  }
  TCanvas *c = new TCanvas("c", "c", 0, 0, 800,600);
  dalitz->Draw("COLZ");

}
