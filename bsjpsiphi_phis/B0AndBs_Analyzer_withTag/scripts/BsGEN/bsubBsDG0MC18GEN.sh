#!/bin/bash
i=0;
max=16;
n=500000;
while [ "$i" -le "$max" ]; do
  mkdir "sGEN$i";
  cd "sGEN$i";
  skip=$(python -c "print 0+$n*$i");
  echo $'#!/bin/sh' > script.sh
  echo $'#BSUB -o test.log' >> script.sh
  echo $'eval `scram runtime -sh`' >> script.sh
  echo "pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2018Lists/BsToJpsiPhi_genlevel1B_DCAP.list hist.root -v outputFile ntu.root -v process Bs -v tagCalChannel BsJPsiPhiMC2018 -v histoMode RECREATE -v use_gen t -v genOnly t -v debug f -v use_tracks t -v use_hlts t -v use_pvts t -v use_hlto t -n $n -s $skip" >> script.sh
  echo "" >> script.sh
  bsub < script.sh;
  cd ..;
  i=`expr "$i" + 1`;
done
