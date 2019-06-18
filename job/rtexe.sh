#!/bin/bash
echo "starting configure condor environment"
echo " CWD: " $PWD
echo " ENV: "
which root
which geant4-config
mkdir -p `dirname $1`

NAME=$1
cluster=$2
echo $NAME > ../job/name.log
echo `date` >> $(dirname $3)/cluster.log
echo -e "$cluster\n\n" >> $(dirname $3)/cluster.log
root -b -q "../job/MultiTiersOutputfun_SiPM.C(\"$NAME\")"
