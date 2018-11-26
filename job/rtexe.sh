#!/bin/bash
echo "starting configure condor environment"
echo " CWD: " $PWD
echo " ENV: "
which root
which geant4-config
mkdir -p `dirname $1`

NAME=$1
thrd=$2
cluster=$3
echo $NAME > ../job/name.log
echo `date` >> $(dirname $1)/cluster.log
echo -e "$cluster\n\n" >> $(dirname $1)/cluster.log
root -b -q "../job/Outputfun_MCP.C(\"$NAME\")"
