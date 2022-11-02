# BsJpsiPhiAnalyzer
This installation guide applies for PD/Legnaro network

## Installation Guide

Create a PD framework area
```
source /lustre/cmswork/share/ntuAna/NtuProdV13_02/build_read nameDirectory SCRAM_ARCH=slc6_amd64_gcc700 CMSSW_VERSION=CMSSW_10_3_0
```

Install tagging code in your nameDirectory/src
```
cp /lustre/cmswork/abragagn/BPH/albertoCode.tar.gz ./
tar -xvzf albertoCode.tar.gz
```
Copy from gitlab PDAnalyzer.cc/.h  and PDSecondNtupleData.h in yours 
`.../nameDirectory/src/PDAnalysis/Ntu/bin/`

Compile 
```
cd .../nameDirectory/src
scram b clean
scram b -j12
```

## How to Run
Command structure:
```
pdTreeAnalyze list histFile.root -v parameter value -n nEvt
```
Parameters string you should not touch:
```
-v histoMode RECREATE -v use_tracks t -v use_hlts t -v use_pvts t -v use_hlto t
```
Parameters you probably want to specify:
 - `outputFile` (string)
 -- output file name
 - `process` (string)
 - `tagCalChannel` (string)
 -- e.g. `BuJPsiKData2018`
 -- atm only `Bs` or `Bd`
 - `use_gen` (bool)
 -- `t` or `f` to switch on and off gen infos
 - `genOnly` (bool)
 - `debug` (bool)

File Lists locations:
`/lustre/cmswork/abragagn/ntuList/`

Example of a complete running command:
```
pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2018Lists/BsToJpsiPhi_2018_DCAP.list hist.root -v outputFile ntu.root -v process Bs -v histoMode RECREATE -v use_gen t -v genOnly f -v debug f -v use_tracks t -v use_hlts t -v use_pvts t -v use_hlto t -n 1000
```
