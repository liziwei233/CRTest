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
echo `date` >> $(dirname $1)/cluster.log
echo -e "$cluster\n\n" >> $(dirname $1)/cluster.log
#root -b -q "../job/Outputfun_MCP_MultiCFD.C(\"$NAME\")"
#root -b -q "../job/Outputfun_MCP_MultiFIX.C(\"$NAME\")"
root -b -q "../job/Outputfun_MCP.C(\"$NAME\")"
