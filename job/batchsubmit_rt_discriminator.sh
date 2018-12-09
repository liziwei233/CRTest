#!/bin/bash

#cat backup/output.C
echo "./batchsubmit_rt_discrimanator [your sub filename] [factype]"
source /data2/R710liziwei/.bashrc
a1=$1
a2=USELG_Alfoilground_lobe30_spike70_withAbs
b1=$2
b2=CFD
subname=${a1:-$a2}
factype=${b1:-$b2}
filename="rtbatch.condo"
echo $filename
echo $subname
echo $factype
#return 1
for fac in 0.1 0.15 0.2 0.25 0.3
do
cp rtbatch_template.condo $filename
#sed -n "/Arguments/p" test.condo 

sed -i "/Arguments/s/filename/$subname/" $filename
sed -i "/Arguments/s/factype/$factype/" $filename
sed -i "/Arguments/s/fac/$fac/" $filename
condor_submit $filename
#sed -i 's//../g' test.condo > test$.condo
done



