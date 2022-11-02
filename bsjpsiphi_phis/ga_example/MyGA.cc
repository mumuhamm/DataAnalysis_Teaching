
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"

using namespace std;
using namespace TMVA;

class MyEvent {
public:
    float ptMup, ptMum, etaMup, etaMum, ptKp, ptKm, etaKp, etaKm, Lxy, ptBs, etaBs, vprobBs;
};

std::vector<MyEvent> vec_signal, vec_data;

void importEvents(std::vector<MyEvent>& vec, TString filename)
{
    std::cout << "Loading " << filename << std::endl;
    
    TFile *fin = new TFile(filename);
    TTree *tree = (TTree*)fin->Get("OutTree_PlusCuts");
    
    float _ptMup,_ptMum, _etaMup, _etaMum, _ptKp, _ptKm, _etaKp, _etaKm, _Lxy, _ptBs, _etaBs, _vprobBs, _LxyErr;

    tree->SetBranchAddress("Lxy",&_Lxy);
    tree->SetBranchAddress("LxyErr",&_LxyErr);
    tree->SetBranchAddress("Bs_Pt",&_ptBs);
    tree->SetBranchAddress("Bs_Eta",&_etaBs);
    tree->SetBranchAddress("Bs_VProb",&_vprobBs);

    tree->SetBranchAddress("Mum_Pt",&_ptMum);
    tree->SetBranchAddress("Mum_Eta",&_etaMum);

    tree->SetBranchAddress("Mup_Pt",&_ptMup);
    tree->SetBranchAddress("Mup_Eta",&_etaMup);

    tree->SetBranchAddress("Kp_Pt",&_ptKp);
    tree->SetBranchAddress("Kp_Eta",&_etaKp);

    tree->SetBranchAddress("Km_Pt",&_ptKm);
    tree->SetBranchAddress("Km_Eta",&_etaKm);

    
    for (int evt=0; evt<tree->GetEntries(); evt++) {
        tree->GetEntry(evt);
        
        if (fabs(_etaMum)<2.5 && fabs(_etaMup)<2.5 && fabs(_etaKm)<2.5 && fabs(_etaKp)<2.5 && _ptMum>3.5 && _ptMup>3.5 && _ptKm>0.6 && _ptKp >0.6 && (_Lxy/_LxyErr)>3) {
            MyEvent event;

	    event.Lxy=_Lxy;
	    event.ptBs=_ptBs;
	    event.etaBs=fabs(_etaBs);
	    event.vprobBs=_vprobBs;

	    event.ptKm = _ptKm;
            event.ptKp   = _ptKp;
            event.etaKm  = fabs(_etaKm);
	    event.etaKp  = fabs(_etaKp);

            event.ptMum = _ptMum;
            event.ptMup   = _ptMup;
            event.etaMum  = fabs(_etaMum);
	    event.etaMup  = fabs(_etaMup);

            
            vec.push_back(event);
        }
    }
    
    std::cout << vec.size() << " events stored." << std::endl;
}

class MyFitness : public IFitterTarget {
public:
    MyFitness() : IFitterTarget() {
    }
    
    Double_t EstimatorFunction( std::vector<Double_t> & factors ){
        
        double CUT_LXY  = factors.at(0);
        double CUT_PTBS = factors.at(1);
	double CUT_ETABS = factors.at(2);
	double CUT_VPROBBS = factors.at(3);

	double CUT_PTMU  = factors.at(4);
	double CUT_ETAMU = factors.at(5);

	double CUT_PTK  = factors.at(6);
        double CUT_ETAK = factors.at(7);	
        
        double S = 0., N = 0.;
        
        // loop over signal sample
        for( std::vector<MyEvent>::iterator it=vec_signal.begin(); it!=vec_signal.end(); it++) {
            
            if (it->Lxy<=CUT_LXY) continue;
            if (it->ptBs<=CUT_PTBS) continue;
            if (it->etaBs>=CUT_ETABS) continue;
	    if (it->vprobBs<=CUT_VPROBBS) continue;

	    if (it->ptMup<=CUT_PTMU) continue;
	    if (it->ptMum<=CUT_PTMU) continue;

	    if (it->etaMum>=CUT_ETAMU) continue;
	    if (it->etaMup>=CUT_ETAMU) continue;

	    if (it->ptKm<=CUT_PTK) continue;
	    if (it->ptKp<=CUT_PTK) continue;

	    if (it->etaKm>=CUT_ETAK) continue;
	    if (it->etaKp>=CUT_ETAK) continue;
            
            S += 1.;
        }
        S =S*0.1; // scale factor between the data fit and the signal MC
        
        // loop over data
        for( std::vector<MyEvent>::iterator it=vec_data.begin(); it!=vec_data.end(); it++) {
            
	  if (it->Lxy<=CUT_LXY) continue;
	  if (it->ptBs<=CUT_PTBS) continue;
	  if (it->etaBs>=CUT_ETABS) continue;
	  if (it->vprobBs<=CUT_VPROBBS) continue;

	  if (it->ptMup<=CUT_PTMU) continue;
	  if (it->ptMum<=CUT_PTMU) continue;

	  if (it->ptKp<=CUT_PTK) continue;
	  if (it->ptKm<=CUT_PTK) continue;

	  if (it->etaMum>=CUT_ETAMU) continue;
	  if (it->etaMup>=CUT_ETAMU) continue;

	  if (it->etaKm>=CUT_ETAK) continue;
	  if (it->etaKp>=CUT_ETAK) continue;
            
	  N += 1.;
        }

	N=N*1.64;//factor to have the full side bands considered

        double significance = S/sqrt(S+N);
	
        return -significance;
    }
};

int main()
{
    importEvents(vec_signal,"/lustre/cmswork/ffiori/AnalysisTrees/TreesForCutOptimization/Plots_Tot_Bs_MC2017_JpsiTkTk_Only.root");
    importEvents(vec_data,"/lustre/cmswork/ffiori/AnalysisTrees/TreesForCutOptimization/Plots_Tot_Bs_2017-2018_JpsiTkTk_Only_SB.root");
    
    std::vector<Interval*> ranges;
    ranges.push_back( new Interval(0.0001,0.1,101) );//Lxy
    ranges.push_back( new Interval(7.,15.,45) );//Bs pt
    ranges.push_back( new Interval(0.01,2.5,100) );//eta range for everything
    ranges.push_back( new Interval(0.001,1,101) );//Vprob Range

    ranges.push_back( new Interval(3.5,5,151) );//ptMup
    //    ranges.push_back( new Interval(3.5,5,151) );//ptMum

    ranges.push_back( new Interval(0.01,2.5,100) );//etaMum
    // ranges.push_back( new Interval(0.0,2.5,52) );//etaMupd

    ranges.push_back( new Interval(0.7,1.5,91) );//ptKm
    //ranges.push_back( new Interval(0.6,1.5,91) );//ptKp

    ranges.push_back( new Interval(0.01,2.5,100) );//etaKm
    //ranges.push_back( new Interval(0.0,2.5,52) );//etaKp

    std::cout << "\nInitial ranges:" << std::endl;
    for( std::vector<Interval*>::iterator it = ranges.begin(); it != ranges.end(); it++ ){
        std::cout << " range: " << (*it)->GetMin() << "   " << (*it)->GetMax() << std::endl;
    }
    
    IFitterTarget* myFitness = new MyFitness();
    
    GeneticAlgorithm mg( *myFitness, 100, ranges );
    
#define CONVSTEPS 20
#define CONVCRIT 0.0001
#define SCSTEPS 10
#define SCRATE 5
#define SCFACTOR 0.95
    
    do {
        mg.Init();
        
        mg.CalculateFitness();
        
        mg.GetGeneticPopulation().Print(0);

        mg.GetGeneticPopulation().TrimPopulation();
        
        mg.SpreadControl( SCSTEPS, SCRATE, SCFACTOR );
        
    } while (!mg.HasConverged( CONVSTEPS, CONVCRIT ));
    
    std::cout << "\nBest factors:" << std::endl;
    
    GeneticGenes* genes = mg.GetGeneticPopulation().GetGenes( 0 );
    std::vector<Double_t> gvec;
    gvec = genes->GetFactors();
    int n = 0;
    for( std::vector<Double_t>::iterator it = gvec.begin(); it<gvec.end(); it++ ){
        std::cout << "FACTOR " << n << " : " << (*it) << std::endl;
        n++;
    }
}
