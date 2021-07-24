//Continumm Suppression TMVA classification :: Author Md. Alibordi
// Adopted from main TMVA classification code
//root -l -e 'TMVA::TMVAGui("BToD0_TMVA_Output.root")'
// When you get the output file
//===================================================================



#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;
int CS_TMVAClassification_BToD0( TString myMethodList = "" )
{

   
    TMVA::Tools::Instance();
    std::map<std::string,int> Use;
   
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   
   //naive Bayes estimator
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0;
   Use["LikelihoodPCA"]   = 0;
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   
   // Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0;
   Use["KNN"]             = 0;
                              
   //Discriminant Analysis & Function Discriminant analysis
   Use["LD"]              = 0;
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0;
   Use["HMatrix"]         = 0;
   Use["FDA_GA"]          = 0;
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   
   Use["MLP"]             = 0;
   Use["MLPBFGS"]         = 0;
   Use["MLPBNN"]          = 1;
   Use["CFMlpANN"]        = 0;
   Use["TMlpANN"]         = 0;
   Use["DNN_CPU"]         = 1;
   Use["DNN_GPU"]         = 0;
                              
   Use["SVM"]             = 0;
   Use["BDT"]             = 0;
   Use["BDTG"]            = 1;
   Use["BDTB"]            = 0;
   Use["BDTD"]            = 0;
   Use["BDTF"]            = 0;
   
      
   Use["RuleFit"]         = 0;
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;
                               
   std::cout << std::endl;
   std::cout << "==> Start continumsuppression_TMVA_BToD0" << std::endl;
   
     
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
      
      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);
         
         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return 1;
         }
         Use[regMethod] = 1;
      }
   }
  
   TString outfileName( "BToD0_TMVA_Output.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification_BToD0", outputFile,
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
   
   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
   
   dataloader->AddVariable("deltaE",'D');
   dataloader->AddVariable("Mbc",'D');
   dataloader->AddVariable("cosTBTO",'D');
   dataloader->AddVariable("DeltaZ",'D');
   dataloader->AddVariable("cosTBz",'D');
   dataloader->AddVariable("thrustBm",'D');
   dataloader->AddVariable("thrustOm",'D');
   dataloader->AddVariable("D0_pt",'D');
   dataloader->AddVariable("D0_E",'D');
   dataloader->AddVariable("KS0_pt",'D');
   dataloader->AddVariable("KS0_E",'D');
   dataloader->AddVariable("Pi0_pt",'D');
   dataloader->AddVariable("Pi0_E",'D');
  /* dataloader->AddVariable("D0_K_S0_pt",'D');
   dataloader->AddVariable("D0_K_S0_E",'D');
   dataloader->AddVariable("D0_K_S0_cosTheta",'D');
   dataloader->AddVariable("D0_K_S0_pi_0_pt",'D');
   dataloader->AddVariable("D0_K_S0_pi_0_E",'D');
   dataloader->AddVariable("D0_K_S0_pi_0_cosTheta",'D');
   */
   TString fnames1 = "/Users/md/Documents/Bs_Review/con_suppression/signal-DPi-10k.root";
   TString fnameb1 = "/Users/md/Documents/Bs_Review/con_suppression/continuum-dpi.root";
   TString fnameb2 = "/Users/md/Documents/Bs_Review/con_suppression/mixed-dpi-00.root";
   TString fnameb3 = "/Users/md/Documents/Bs_Review/con_suppression/charged-dpi-00.root";
   
   
   TFile *inputs1 = TFile::Open(fnames1);
   TFile *inputb1 = TFile::Open(fnameb1);
   TFile *inputb2 = TFile::Open(fnameb2);
   TFile *inputb3 = TFile::Open(fnameb3);
   
   std::cout << "--- TMVAClassification       : Using input file: "<<inputs1->GetName()<<"\n";
   std::cout << "--- TMVAClassification       : Using input file: "<<inputb1->GetName()<<"\n";
   std::cout << "--- TMVAClassification       : Using input file: "<<inputb2->GetName()<<"\n";
   std::cout << "--- TMVAClassification       : Using input file: "<<inputb3->GetName()<<"\n";
   
   
   TTree *signal1     = (TTree*)inputs1->Get("bp3");
   TTree *background1 = (TTree*)inputb1->Get("bp3");
   TTree *background2 = (TTree*)inputb2->Get("bp3");
   TTree *background3 = (TTree*)inputb3->Get("bp3");
   
   Double_t signalWeight;//     = 1.0;
   Double_t backgroundWeight1,  backgroundWeight2, backgroundWeight3;//, = 1.0;
   
   TFile *sig = new TFile("/Users/md/Documents/Bs_Review/con_suppression/signal-DPi-10k.root");
   TTree* sigtree = (TTree*)sig->Get("bp3");
   Double_t sigw;
   sigtree->SetBranchAddress("__weight__",&sigw);for (Int_t i=0;i<sigtree->GetEntries();i++) {sigtree->GetEntry(i);signalWeight=sigw;}
   
   TFile *bkg1 = new TFile("/Users/md/Documents/Bs_Review/con_suppression/continuum-dpi.root");
   TTree* bkgtree1 = (TTree*)bkg1->Get("bp3");
   Double_t bkgw1;
   bkgtree1->SetBranchAddress("__weight__",&bkgw1);for (Int_t i=0;i<bkgtree1->GetEntries();i++) {bkgtree1->GetEntry(i);backgroundWeight1=bkgw1;}
   
   TFile *bkg2 = new TFile("/Users/md/Documents/Bs_Review/con_suppression/mixed-dpi-00.root");
   TTree* bkgtree2 = (TTree*)bkg2->Get("bp3");
   Double_t bkgw2;
   bkgtree2->SetBranchAddress("__weight__",&bkgw2);for (Int_t i=0;i<bkgtree2->GetEntries();i++) {bkgtree2->GetEntry(i);backgroundWeight2=bkgw2;}
   
   TFile *bkg3 = new TFile("/Users/md/Documents/Bs_Review/con_suppression/charged-dpi-00.root");
   TTree* bkgtree3 = (TTree*)bkg3->Get("bp3");
   Double_t bkgw3;
   bkgtree3->SetBranchAddress("__weight__",&bkgw3);for (Int_t i=0;i<bkgtree3->GetEntries();i++) {bkgtree3->GetEntry(i);backgroundWeight3=bkgw3;}

      
      
      
   
   dataloader->AddSignalTree(signal1,signalWeight);
   dataloader->AddBackgroundTree(background1,backgroundWeight1);
   dataloader->AddBackgroundTree(background2,backgroundWeight2);
   dataloader->AddBackgroundTree(background3,backgroundWeight3);
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";
   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                       "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   
   //================Classifier train and test
   
   
   if (Use["Cuts"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "Cuts",
                          "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
   
   if (Use["CutsD"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsD",
                          "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );
   
   if (Use["CutsPCA"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsPCA",
                          "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );
   
   if (Use["CutsGA"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsGA",
                          "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );
   
   if (Use["CutsSA"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsSA",
                          "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
   
   if (Use["Likelihood"])
      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "Likelihood",
                          "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );
   
   if (Use["LikelihoodD"])
      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodD",
                          "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );
   
   if (Use["LikelihoodPCA"])
      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodPCA",
                          "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" );

   if (Use["LikelihoodKDE"])
      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodKDE",
                          "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" );
 
   if (Use["LikelihoodMIX"])
      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodMIX",
                          "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" );
 
   if (Use["PDERS"])
      factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERS",
                          "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );
   
   if (Use["PDERSD"])
      factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERSD",
                          "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );
   
   if (Use["PDERSPCA"])
      factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERSPCA",
                          "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );
  
   if (Use["PDEFoam"])
      factory->BookMethod( dataloader, TMVA::Types::kPDEFoam, "PDEFoam",
                          "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );
   
   if (Use["PDEFoamBoost"])
      factory->BookMethod( dataloader, TMVA::Types::kPDEFoam, "PDEFoamBoost",
                          "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );
  
   if (Use["KNN"])
      factory->BookMethod( dataloader, TMVA::Types::kKNN, "KNN",
                          "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );
   
   if (Use["HMatrix"])
      factory->BookMethod( dataloader, TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );
   
   if (Use["LD"])
      factory->BookMethod( dataloader, TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
   
   if (Use["Fisher"])
      factory->BookMethod( dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
   
   if (Use["FisherG"])
      factory->BookMethod( dataloader, TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );
   
   if (Use["BoostedFisher"])
      factory->BookMethod( dataloader, TMVA::Types::kFisher, "BoostedFisher",
                          "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );
   
   if (Use["FDA_MC"])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MC",
                          "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );
   
   if (Use["FDA_GA"])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_GA",
                          "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=100:Cycles=2:Steps=5:Trim=True:SaveBestGen=1" );
   
   if (Use["FDA_SA"])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_SA",
                          "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
   
   if (Use["FDA_MT"])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MT",
                          "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );
   
   if (Use["FDA_GAMT"])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_GAMT",
                          "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );
   
   if (Use["FDA_MCMT"])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MCMT",
                          "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );
   
  
   if (Use["MLP"])
      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
   
   if (Use["MLPBFGS"])
      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );
   
   if (Use["MLPBNN"])
      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=60:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" );
   if (Use["DNN_CPU"] or Use["DNN_GPU"]) {
        
      TString layoutString ("Layout=TANH|128,TANH|128,TANH|128,LINEAR");
      
       
      TString training0("LearningRate=1e-1,Momentum=0.9,Repetitions=1,"
                        "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.5+0.5+0.5, Multithreading=True");
      TString training1("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                        "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
      TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
                        "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
      TString trainingStrategyString ("TrainingStrategy=");
      trainingStrategyString += training0 + "|" + training1 + "|" + training2;
      
    
      TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                          "WeightInitialization=XAVIERUNIFORM");
      dnnOptions.Append (":"); dnnOptions.Append (layoutString);
      dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);
      
     
      if (Use["DNN_GPU"]) {
         TString gpuOptions = dnnOptions + ":Architecture=GPU";
         factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN_GPU", gpuOptions);
      }
        
      if (Use["DNN_CPU"]) {
         TString cpuOptions = dnnOptions + ":Architecture=CPU";
         factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN_CPU", cpuOptions);
      }
   }
  
   if (Use["CFMlpANN"])
      factory->BookMethod( dataloader, TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N"  );
   if (Use["TMlpANN"])
      factory->BookMethod( dataloader, TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  );
   if (Use["SVM"])
      factory->BookMethod( dataloader, TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );
   
   if (Use["BDTG"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
                          "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
   
   if (Use["BDT"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                          "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
   
   if (Use["BDTB"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTB",
                          "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );
   
   if (Use["BDTD"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTD",
                          "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );
   
   if (Use["BDTF"])
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTF",
                          "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );
   
   if (Use["RuleFit"])
      factory->BookMethod( dataloader, TMVA::Types::kRuleFit, "RuleFit",
                          "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );
   
   //=========================================
   
   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();
   outputFile->Close();
   
   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;
   
   delete factory;
   delete dataloader;
   if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );
   
   return 0;
}

int main( int argc, char** argv )
{
      
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   return CS_TMVAClassification_BToD0(methodList);
}
   
   
   
   
   
   
   
   

