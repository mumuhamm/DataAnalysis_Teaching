# BsJpsiPhi analysis toolkit

better performance on lxplus for python scripts:

ssh -Y name@lxplus7.cern.ch

. /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.14.06/x86_64-centos7-gcc48-opt/bin/thisroot.sh

# B0AndBs_Analyzer

this is the analyzer to be used, it searches for Bs and B0 secondary vertices and dump relvant quantity in a TTree.
Priority is given to the Bs candidates (i.e if a B0 and a Bs are present in the same event, the B0 will be not saved), 
it select the one with the highest vertex probability with the basic requirement of HLT matching for the the J/psi (i.e the two muons)

The analyzer is build using the Padova code described here https://twiki.cern.ch/twiki/bin/view/CMS/PDGeneralNtuple , 
to use it please follow instruction here https://twiki.cern.ch/twiki/bin/view/CMS/PDGeneralNtuple#step_3_b to install the code, 
then copy the PDAnalyzer.cc/h files in the relevant directory (<YourArea>/src/PDAnalysis/Ntu/bin/.) and issue a scram b.

To run on some data you need a list of root files as those provided here https://cernbox.cern.ch/index.php/s/zvxnXqszJhNtHp2

A possible command can be:
pdTreeAnalyze <file list>.list YourHisto.root -v histoMode RECREATE -v use_gen f -v use_tracks t -v use_hlts t -v use_pvts t -v use_hlto t -n 1000

In the analyzer there are a cpl of relevant flags, at the moment to be set manually:

- debug if set to true will produce quite a lot of printouts to help in case of bugs
- genOny if set to true, only the GEN level information will be processed (working only for the Bs case atm)

