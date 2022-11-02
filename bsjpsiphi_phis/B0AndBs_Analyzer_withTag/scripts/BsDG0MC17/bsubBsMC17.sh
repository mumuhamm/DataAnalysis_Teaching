#!/bin/bash
i=0;
max=66;
n=1000000;
while [ "$i" -le "$max" ]; do
  mkdir "s$i";
  cd "s$i";
  skip=$(python -c "print 0+$n*$i");
  echo $'#!/bin/sh' > script.sh
  echo $'#BSUB -o test.log' >> script.sh
  echo $'eval `scram runtime -sh`' >> script.sh
  echo "pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2017Lists/BsToJpsiPhi_2017_DCAP.list hist.root -v outputFile ntu.root -v process Bs -v tagCalChannel BsJPsiPhiMC2017 -v histoMode RECREATE -v use_gen t -v genOnly f -v debug f -v use_tracks t -v use_hlts t -v use_pvts t -v use_hlto t -n $n -s $skip" >> script.sh
  echo "" >> script.sh
  bsub < script.sh;
  cd ..;
  i=`expr "$i" + 1`;
done
