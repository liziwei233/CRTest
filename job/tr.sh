#!/bin/bash
#NAME=/baby/one/more/time
#echo $NAME>name.log
FULLNAME=$(cat name.log)
echo $FULLNAME
POS=$(dirname $FULLNAME)
echo $POS
NAME=$(basename $POS)
echo $NAME
for thrd in 1 3 5 10 20 30
do
echo $thrd
hadd $POS/thrd$thrd.root $POS/{0..39}data_thrd$thrd.root
root -b -q "MultiTiers_TACor_sep.C(\"$POS/thrd$thrd\")"
done
