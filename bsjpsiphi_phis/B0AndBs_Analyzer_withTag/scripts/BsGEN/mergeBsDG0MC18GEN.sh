#!/bin/bash
i=0;
max=16;
echo $'#!/bin/sh' > park.sh
echo -n "hadd ntuBsGEN.root" >> park.sh
while [ "$i" -le "$max" ]; do
  echo -n " sGEN$i/ntu" >> park.sh
  echo -n ".root" >> park.sh
  i=`expr "$i" + 1`;
done
echo " " >> park.sh
bash park.sh;
rm park.sh;
mv ntuBsGEN.root ../
