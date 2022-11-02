import json
import sys
import subprocess
import csv


print "export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH (bash)or " 
print "export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH (bash, CERN only)"

#export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH (bash)
jsonfile=sys.argv[1]
#json_data=open("Cert_294927-301141_13TeV_PromptReco_Collisions17_JSON_MuonPhys.txt").read()
#json_data=open("Cert_294927-302343_13TeV_PromptReco_Collisions17_JSON_MuonPhys.txt").read()
#json_data=open("Cert_294927-305364_13TeV_PromptReco_Collisions17_JSON_MuonPhys.txt").read()
json_data=open(jsonfile).read()
data = json.loads(json_data)

lumidict={}

for run in data:
  bashCommand = "brilcalc lumi --output-style csv -i \"{ "+str(run)+" : "+str(data[run])+" }\""
  print bashCommand
  output = subprocess.check_output(['bash','-c', bashCommand])
  print output
  print float(output.splitlines()[2].split(",")[5])/1e6
  lumidict[run]=float(output.splitlines()[2].split(",")[5])/1e6

with open(".."+jsonfile.strip(".json")+'_lumidict.json', 'w') as outfile:
    json.dump(lumidict, outfile)
