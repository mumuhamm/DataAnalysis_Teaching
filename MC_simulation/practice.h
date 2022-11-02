#ifndef PRACTICE_H_HAS_BEEN_INCLUDED
#define PRACTICE_H_HAS_BEEN_INCLUDED
#pragma once
//--------------------------------
#include "TChain.h"
//#include "TFile.h"
#include <sstream>
#include <iostream>
#include <string>
#include <fstream> 
#include <vector>
#include <map>
#include <cmath>
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TBuffer.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"
#include "TRandom.h"


class analysis
{
   

 public:
   float px, py;
   float PT(float px, float py){ return(TMath::Sqrt(px*px + py*py));}
      
   };

#endif

