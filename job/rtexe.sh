#!/bin/bash
echo "starting configure condor environment"
echo " CWD: " $PWD
echo " ENV: "
which root
which geant4-config
mkdir -p `dirname $1`

NAME=$1
fac=$2
factype=$3
cluster=$4
echo $NAME > ../job/name.log
echo `date` >> $(dirname $1)/cluster.log
echo -e "$cluster\n\n" >> $(dirname $1)/cluster.log
root -b -q "../job/Outputfun_MCP_multiDsk.C(\"$NAME\",$fac,\"$factype\")"
#root -b -q "../job/Outputfun_MCP.C(\"$NAME\",$fac,\"$factype\")"
#root -b -q "../job/Outputfun_MCP_sample.C(\"$NAME\",$fac,\"$factype\")"
