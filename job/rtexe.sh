#!/bin/bash
echo "starting configure condor environment"
echo " CWD: " $PWD
echo " ENV: "
which root
which geant4-config
mkdir -p `dirname $1`

NAME=$1
id=$2
cluster=$3
echo $NAME > ../job/name.log
flag=0
if [ $id -eq 0 ]
then
echo `date` >> $(dirname $1)/cluster.log
echo -e "$cluster\n\n" >> $(dirname $1)/cluster.log
fi

root -b -q "../job/MultiTiersOutputfun_sample.C(\"$NAME\",9e4)"
#root -b -q "../job/MultiTiersOutputfun_SiPM.C(\"$NAME\")"
