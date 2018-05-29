#!/bin/bash
echo "starting configure condor environment"
echo " CWD: " $PWD
echo " ENV: "
which root
which geant4-config
mkdir -p `dirname $1`

NAME=$1
thrd=$2
echo $NAME > ../job/name.log
root -b -q "../job/MultiTiersOutputfun_SiPM.C(\"$NAME\",$thrd)"
