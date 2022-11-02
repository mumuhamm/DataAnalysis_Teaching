import ROOT as R
from array import array
import math
import json
import pickle

#
# Giacomo Fedi
#

LUMYXBIN=1000 #1/pb

json_data16=open("../Data/jsons/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_lumidict.json").read()
lumidict16 = json.loads(json_data16)
json_data17=open("../Data/jsons/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_lumidict.json").read()
lumidict17 = json.loads(json_data17)
json_data18=open("../Data/jsons/Cert_314472-321461_13TeV_PromptReco_Collisions18_JSON_lumidict.json").read()
lumidict18 = json.loads(json_data18)

lumitotal16=0.
lumitotal17=0.
lumitotal18=0.
for i in lumidict16:
  lumitotal16=lumitotal16+lumidict16[i]
for i in lumidict17:
  lumitotal17=lumitotal17+lumidict17[i]
for i in lumidict18:
  lumitotal18=lumitotal18+lumidict18[i]

print("Total lumy 2016=",lumitotal16," 1/pb")
print("Total lumy 2017=",lumitotal17," 1/pb")
print("Total lumy 2018=",lumitotal18," 1/pb")


json_crab16B=open("../Data/jsons/crab_Charmonium_Run2016B-07Aug17_ver2-v1_AOD_lumidict.json").read()
lumidict16B = json.loads(json_crab16B)
json_crab16C=open("../Data/jsons/crab_Charmonium_Run2016C-07Aug17-v1_AOD_lumidict.json").read()
lumidict16C = json.loads(json_crab16C)
json_crab16D=open("../Data/jsons/crab_Charmonium_Run2016D-07Aug17-v1_AOD_lumidict.json").read()
lumidict16D = json.loads(json_crab16D)
json_crab16E=open("../Data/jsons/crab_Charmonium_Run2016E-07Aug17-v1_AOD_lumidict.json").read()
lumidict16E = json.loads(json_crab16E)
json_crab16F=open("../Data/jsons/crab_Charmonium_Run2016F-07Aug17-v1_AOD_lumidict.json").read()
lumidict16F = json.loads(json_crab16F)
json_crab16G=open("../Data/jsons/crab_Charmonium_Run2016G-07Aug17-v1_AOD_lumidict.json").read()
lumidict16G = json.loads(json_crab16G)
json_crab16H=open("../Data/jsons/crab_Charmonium_Run2016H-07Aug17-v1_AOD_lumidict.json").read()
lumidict16H = json.loads(json_crab16H)

lumidict_crab16=dict(lumidict16B,**lumidict16C)
lumidict_crab16=dict(lumidict_crab16,**lumidict16D)
lumidict_crab16=dict(lumidict_crab16,**lumidict16E)
lumidict_crab16=dict(lumidict_crab16,**lumidict16F)
lumidict_crab16=dict(lumidict_crab16,**lumidict16G)
lumidict_crab16=dict(lumidict_crab16,**lumidict16H)


lumitotal_crab16=0.
for i in lumidict_crab16:
  lumitotal_crab16=lumitotal_crab16+lumidict_crab16[i]


print("Total lumy crab 2016=",lumitotal_crab16," 1/pb")



json_crab17B=open("../Data/jsons/crab_Charmonium_Run2017B-31Mar2018-v1_MINIAOD_lumidict.json").read()
lumidict17B = json.loads(json_crab17B)
json_crab17C=open("../Data/jsons/crab_Charmonium_Run2017C-31Mar2018-v1_MINIAOD_lumidict.json").read()
lumidict17C = json.loads(json_crab17C)
json_crab17D=open("../Data/jsons/crab_Charmonium_Run2017D-31Mar2018-v1_MINIAOD_lumidict.json").read()
lumidict17D = json.loads(json_crab17D)
json_crab17E=open("../Data/jsons/crab_Charmonium_Run2017E-31Mar2018-v1_MINIAOD_lumidict.json").read()
lumidict17E = json.loads(json_crab17E)
json_crab17F=open("../Data/jsons/crab_Charmonium_Run2017F-31Mar2018-v1_MINIAOD_lumidict.json").read()
lumidict17F = json.loads(json_crab17F)


lumidict_crab17=dict(lumidict17B,**lumidict17C) 
lumidict_crab17=dict(lumidict_crab17,**lumidict17D) 
lumidict_crab17=dict(lumidict_crab17,**lumidict17E) 
lumidict_crab17=dict(lumidict_crab17,**lumidict17F) 


lumitotal_crab17=0.
for i in lumidict_crab17:
  lumitotal_crab17=lumitotal_crab17+lumidict_crab17[i]


print("Total lumy crab 2017=",lumitotal_crab17," 1/pb")



json_crab18A1=open("../Data/jsons/crab_Charmonium_Run2018A-PromptReco-v1_MINIAOD_lumidict.json").read()
lumidict18A1 = json.loads(json_crab18A1)
json_crab18A2=open("../Data/jsons/crab_Charmonium_Run2018A-PromptReco-v2_MINIAOD_lumidict.json").read()
lumidict18A2 = json.loads(json_crab18A2)
json_crab18A3=open("../Data/jsons/crab_Charmonium_Run2018A-PromptReco-v3_MINIAOD_lumidict.json").read()
lumidict18A3 = json.loads(json_crab18A3)
json_crab18B1=open("../Data/jsons/crab_Charmonium_Run2018B-PromptReco-v1_MINIAOD_lumidict.json").read()
lumidict18B1 = json.loads(json_crab18B1)
json_crab18B2=open("../Data/jsons/crab_Charmonium_Run2018B-PromptReco-v2_MINIAOD_lumidict.json").read()
lumidict18B2 = json.loads(json_crab18B2)
json_crab18C1=open("../Data/jsons/crab_Charmonium_Run2018C-PromptReco-v1_MINIAOD_lumidict.json").read()
lumidict18C1 = json.loads(json_crab18C1)
json_crab18C2=open("../Data/jsons/crab_Charmonium_Run2018C-PromptReco-v2_MINIAOD_lumidict.json").read()
lumidict18C2 = json.loads(json_crab18C2)
json_crab18C3=open("../Data/jsons/crab_Charmonium_Run2018C-PromptReco-v3_MINIAOD_lumidict.json").read()
lumidict18C3 = json.loads(json_crab18C3)
json_crab18D=open("../Data/jsons/crab_Charmonium_Run2018D-PromptReco-v2_MINIAOD_lumidict.json").read()
lumidict18D = json.loads(json_crab18D)

lumidict_crab18=dict(lumidict18A1,**lumidict18A2) 
lumidict_crab18=dict(lumidict_crab18,**lumidict18A3)
lumidict_crab18=dict(lumidict_crab18,**lumidict18B1)
lumidict_crab18=dict(lumidict_crab18,**lumidict18B2)
lumidict_crab18=dict(lumidict_crab18,**lumidict18C1)
lumidict_crab18=dict(lumidict_crab18,**lumidict18C2)
lumidict_crab18=dict(lumidict_crab18,**lumidict18C3)
lumidict_crab18=dict(lumidict_crab18,**lumidict18D)

lumitotal_crab18=0.
for i in lumidict_crab18:
  lumitotal_crab18=lumitotal_crab18+lumidict_crab18[i]


print("Total lumy crab 2018=",lumitotal_crab18," 1/pb")


lumidict_run2=dict(lumidict_crab16,**lumidict_crab17)
lumidict_run2=dict(lumidict_run2,**lumidict_crab18)

data={}
idx=1
lumytemp=0
runstart=0
for run in range(270000,330000):
  srun=str(run)
  if srun not in lumidict_run2: continue
  if lumytemp==0 and lumidict_run2[srun]>0: 
    runstart=srun
    data[idx]={}
  lumytemp=lumytemp+lumidict_run2[srun]
  if lumytemp>LUMYXBIN:
    data[idx]["start"]=int(runstart)
    data[idx]["stop"]=int(srun)
    data[idx]["lumy"]=lumytemp
    lumytemp=0
    idx=idx+1
    
del data[idx]

print(data)
     
t=R.TChain("OutTreeB0")
t.Add("/afs/cern.ch/work/g/gfedi/B0samples/Tot_Run2016B.root")
t.Add("/afs/cern.ch/work/g/gfedi/B0samples/Tot_Run2016C.root")
t.Add("/afs/cern.ch/work/g/gfedi/B0samples/Tot_Run2016D.root")
t.Add("/afs/cern.ch/work/g/gfedi/B0samples/Tot_Run2016E.root")
t.Add("/afs/cern.ch/work/g/gfedi/B0samples/Tot_Run2016F.root")
t.Add("/afs/cern.ch/work/g/gfedi/B0samples/Tot_Run2016G.root")
t.Add("/afs/cern.ch/work/g/gfedi/B0samples/Tot_Run2016H.root")
t.Add("/afs/cern.ch/work/g/gfedi/B0samples/Tot_B0_2017-2018.root")

print("Total num events=",t.GetEntries())

selection=R.TEventList("selection","selection")
t.Draw(">>selection","B0_MassFromSV>5.15 && B0_MassFromSV<5.4 && Mum_Pt>4 && Mup_Pt>4 && K_Pt>0.8 && Pi_Pt>0.8 && Lxy>0.02 && abs(Jpsi_Mass-3.0969)<0.150 && B0_VProb>0.02 && (HLT_JpsiTkTk==1 || HLT_JpsiTk==1 || HLT_JpsiMu==1)")
t.SetEventList(selection)

tselection=t.CopyTree("")

print("Selected num events=",tselection.GetEntries())

for idx in data:
  data[idx]["data"]=[]


#massplot=R.TH1F("massplot","massplot",50,5.15,5.4)
#tselection.Draw("B0_MassFromSV>>massplot","")
#massplot.Draw()

for idx in data:
    if data[idx]['start']<294927:
      data[idx]["dataTk"]=tselection.CopyTree("run>="+str(data[idx]["start"])+" && run<="+str(data[idx]["stop"])+" && HLT_JpsiTk==1")
    else:
      data[idx]["dataTk"]=tselection.CopyTree("run>="+str(data[idx]["start"])+" && run<="+str(data[idx]["stop"])+" && HLT_JpsiTkTk==1")
    print("Bin=",idx," Num Tk=",data[idx]["dataTk"].GetEntries())
    data[idx]["dataMu"]=tselection.CopyTree("run>="+str(data[idx]["start"])+" && run<="+str(data[idx]["stop"])+" && HLT_JpsiMu==1")
    print("Bin=",idx," Num Mu=",data[idx]["dataMu"].GetEntries())

with open('/eos/user/g/gfedi/B0split.pik', 'w') as outfile:
  pickle.dump(data, outfile)
