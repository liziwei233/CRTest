#!/bin/bash
#NAME=/baby/one/more/time
#echo $NAME>name.log
FULLNAME=$(cat name.log)
echo $FULLNAME
POS=$(dirname $FULLNAME)
echo $POS
NAME=$(basename $POS)
echo $NAME
hadd $POS/$NAME.root $POS/{0..39}data.root
root -b -q "MultiTiers_TACor_sep.C(\"$POS/$NAME\")"
